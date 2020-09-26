#include "stagbl/private/stagblsystemimpl.h"
#include "stagblsystemsimpleimpl.h"
#include "petscblaslapack.h"
#include <stdlib.h>

/* For DMStagStencilToIndexLocal(), which should be made public */
#include <petsc/private/dmstagimpl.h>

PetscErrorCode StagBLSystemDestroy_Simple(StagBLSystem system)
{
  PetscErrorCode      ierr;
  StagBLSystem_Simple *data = (StagBLSystem_Simple*) system->data;

  PetscFunctionBegin;
  if (data->mat) free(data->mat);
  if (data->rhs) {
    ierr = StagBLArrayDestroy(&data->rhs);CHKERRQ(ierr);
  }
  free(system->data);
  system->data = NULL;
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLSystemOperatorSetValuesStencil_Simple(StagBLSystem system,PetscInt nrows,const DMStagStencil *rows,PetscInt ncols,const DMStagStencil *cols, const PetscScalar *values)
{
  PetscErrorCode         ierr;
  StagBLSystem_Simple    *data = (StagBLSystem_Simple*) system->data;
  PetscInt               *row_indices_local,*col_indices_local;
  const PetscInt         *l2g_indices;
  DM                     dm;
  PetscInt               dim,global_offset;
  ISLocalToGlobalMapping l2g;
  PetscMPIInt            size;

  PetscFunctionBegin;

  ierr = StagBLGridPETScGetDM(system->grid,&dm);CHKERRQ(ierr);
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  row_indices_local = (PetscInt*) malloc(nrows * sizeof(PetscInt));
  col_indices_local = (PetscInt*) malloc(ncols * sizeof(PetscInt));
  ierr = DMStagStencilToIndexLocal(dm,dim,nrows,rows,row_indices_local);CHKERRQ(ierr);
  ierr = DMStagStencilToIndexLocal(dm,dim,ncols,cols,col_indices_local);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm),&size);CHKERRQ(ierr);
  if (size != 1) {
    StagBLError1(PetscObjectComm((PetscObject)dm),"%s only implemented for uniprocessor case",__func__);
  } else {
    global_offset = 0;
  }
  ierr = DMGetLocalToGlobalMapping(dm,&l2g);CHKERRQ(ierr);
  ierr = ISLocalToGlobalMappingGetIndices(l2g,&l2g_indices);CHKERRQ(ierr);
  for (PetscInt i=0; i<nrows; ++i) {
    const PetscInt row = l2g_indices[row_indices_local[i]] - global_offset;

    for (PetscInt j=0; j<ncols; ++j) {
      const PetscInt col = l2g_indices[col_indices_local[j]] - global_offset;

      data->mat[row + data->system_size * col] = values[j + nrows * i];
    }
  }
  ierr = ISLocalToGlobalMappingRestoreIndices(l2g,&l2g_indices);CHKERRQ(ierr);
  free(row_indices_local);
  free(col_indices_local);
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLSystemRHSSetConstant_Simple(StagBLSystem system,PetscScalar value)
{
  PetscErrorCode      ierr;
  StagBLSystem_Simple *data = (StagBLSystem_Simple*) system->data;

  PetscFunctionBegin;
  ierr = StagBLArraySetLocalConstant(data->rhs,value);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLSystemRHSSetValuesStencil_Simple(StagBLSystem system,PetscInt nrows,const DMStagStencil *rows, const PetscScalar *values)
{
  PetscErrorCode      ierr;
  StagBLSystem_Simple *data = (StagBLSystem_Simple*) system->data;

  PetscFunctionBegin;
  ierr = StagBLArraySetLocalValuesStencil(data->rhs,nrows,rows,values);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// FIXME make this static and a class function (and get rid of StagBLSolver)
PetscErrorCode StagBLSystemSolve_Simple(StagBLSystem system,StagBLArray sol)
{
  PetscErrorCode      ierr;
  StagBLSystem_Simple *data = (StagBLSystem_Simple*) system->data;
  PetscScalar         *rhs_raw,*sol_raw,*mat_copy;

  PetscFunctionBegin;
  // FIXME only do l2g if needed:
  ierr = StagBLArrayLocalToGlobal(data->rhs);CHKERRQ(ierr);
  ierr = StagBLArraySimpleGetGlobalRaw(data->rhs,&rhs_raw);CHKERRQ(ierr);
  ierr = StagBLArraySimpleGetGlobalRaw(sol,&sol_raw);CHKERRQ(ierr); // FIXME this definitely should do a type check

  /* Use LAPACK to do the solve */
  mat_copy = (PetscScalar*) malloc(data->system_size * data->system_size * sizeof(PetscScalar));
  for (PetscInt i=0; i<data->system_size * data->system_size; ++i) mat_copy[i] = data->mat[i];
  for (PetscInt i=0; i<data->system_size; ++i) sol_raw[i] = rhs_raw[i];
  {
    PetscBLASInt info, nrhs=1;
    PetscBLASInt *pivots = (PetscBLASInt*) malloc(data->system_size * sizeof(PetscBLASInt));
    PetscStackCallBLAS("LAPACKgesv",LAPACKgesv_(
          &data->system_size, // N
          &nrhs,              // NRHS
          data->mat,          // A
          &data->system_size, // LDA
          pivots,             // IPIV
          sol_raw,            // B
          &data->system_size, // LDB
          &info               // INFO
    ));

    if (info != 0) StagBLError1(PETSC_COMM_SELF,"dgsev failed with code %d",info);
    free(pivots);
  }
  free(mat_copy);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemCreate_Simple(StagBLSystem system)
{
  PetscErrorCode       ierr;
  StagBLSystem_Simple  *data;
  PetscMPIInt          size;
  DM                   dm;

  PetscFunctionBegin;
  system->ops->destroy = StagBLSystemDestroy_Simple;
  system->ops->operatorsetvaluesstencil = StagBLSystemOperatorSetValuesStencil_Simple;
  system->ops->rhssetconstant = StagBLSystemRHSSetConstant_Simple;
  system->ops->rhssetvaluesstencil = StagBLSystemRHSSetValuesStencil_Simple;

  system->solver_type = STAGBLSOLVERSIMPLE;

  ierr = StagBLGridPETScGetDM(system->grid,&dm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm),&size);CHKERRQ(ierr);
  if (size > 1) StagBLError(PetscObjectComm((PetscObject)dm),"Simple solver is only for uniprocessor/sequential use");
  system->data = (void*) malloc(sizeof(StagBLSystem_Simple));
  data = (StagBLSystem_Simple*) system->data;

  ierr = DMStagGetEntries(dm,&data->system_size);CHKERRQ(ierr);
  ierr = StagBLArrayCreate(system->grid,&data->rhs,STAGBLARRAYSIMPLE);CHKERRQ(ierr);
  data->mat = (PetscScalar*) calloc(data->system_size * data->system_size, sizeof(PetscScalar));
  PetscFunctionReturn(0);
}

