#include "stagbl/private/stagblsystemimpl.h"
#include "stagblsystempetscimpl.h"
#include <stdlib.h>

/* Note: this residual is Ax - b, not the usual b - Ax, so that the Jacobian is A, not -A */
static PetscErrorCode StagBLSystemPetscResidual_Default(SNES snes,Vec x, Vec f, void *ctx)
{
  PetscErrorCode             ierr;
  StagBLSystem               stagblsystem = (StagBLSystem) ctx;
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBegin;
  STAGBL_UNUSED(snes);
  ierr = MatAssemblyBegin(data->mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(data->rhs);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(data->mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatMult(data->mat,x,f);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(data->rhs);CHKERRQ(ierr);
  ierr = VecAXPY(f,-1.0,data->rhs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLSystemPetscJacobian_Default(SNES snes,Vec x, Mat Amat, Mat Pmat,void *ctx)
{
  PetscErrorCode             ierr;
  StagBLSystem               stagblsystem = (StagBLSystem) ctx;
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBegin;
  STAGBL_UNUSED(snes);
  STAGBL_UNUSED(x);
  ierr = MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCopy(data->mat,Amat,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  if (Pmat != Amat) {
    ierr = MatAssemblyBegin(Pmat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Pmat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatCopy(data->mat,Pmat,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemDestroy_PETSc(StagBLSystem stagblsystem)
{
  PetscErrorCode     ierr;
  StagBLSystem_PETSc *data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBegin;
  if (data->rhs) {
    ierr = VecDestroy(&data->rhs);CHKERRQ(ierr);
  }
  if (data->mat) {
    ierr = MatDestroy(&data->mat);CHKERRQ(ierr);
  }
  if (data->snes) {
    ierr = SNESDestroy(&data->snes);CHKERRQ(ierr);
  }
  free(stagblsystem->data);
  stagblsystem->data = NULL;
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLSystemOperatorSetValuesStencil_PETSc(StagBLSystem system,PetscInt nrows,const DMStagStencil *rows,PetscInt ncols,const DMStagStencil *cols, const PetscScalar *values)
{
  PetscErrorCode     ierr;
  StagBLSystem_PETSc *data = (StagBLSystem_PETSc*) system->data;
  DM                 dm;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(system->grid,&dm);CHKERRQ(ierr);
  ierr = DMStagMatSetValuesStencil(dm,data->mat,nrows,rows,ncols,cols,values,INSERT_VALUES);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLSystemRHSSetConstant_PETSc(StagBLSystem system,PetscScalar value)
{
  PetscErrorCode      ierr;
  StagBLSystem_PETSc *data = (StagBLSystem_PETSc*) system->data;

  PetscFunctionBegin;
  ierr = VecSet(data->rhs,value);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLSystemRHSSetValuesStencil_PETSc(StagBLSystem system,PetscInt nrows,const DMStagStencil *rows, const PetscScalar *values)
{
  PetscErrorCode     ierr;
  StagBLSystem_PETSc *data = (StagBLSystem_PETSc*) system->data;
  DM                 dm;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(system->grid,&dm);CHKERRQ(ierr);
  ierr = DMStagVecSetValuesStencil(dm,data->rhs,nrows,rows,values,INSERT_VALUES);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLSystemSolve_PETSc(StagBLSystem system,StagBLArray sol)
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) system->data;
  StagBLGrid                 grid;
  PetscErrorCode             ierr;
  Vec                        vec_sol;
  DM                         dm;

  PetscFunctionBegin;
  ierr = StagBLArrayGetStagBLGrid(sol,&grid);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(grid,&dm);CHKERRQ(ierr);

  ierr = StagBLArrayPETScGetGlobalVec(sol,&vec_sol);CHKERRQ(ierr);

  /* Create the solution Vec, if needbe */
  if (!vec_sol)
  {
    Vec *p_vec_sol;

    ierr = StagBLArrayPETScGetGlobalVecPointer(sol,&p_vec_sol);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm,p_vec_sol);CHKERRQ(ierr);
    vec_sol = *p_vec_sol;
  }

  /* Create the SNES object from the system, if needbe */
  if (!data->snes) {

    ierr = SNESCreate(PetscObjectComm((PetscObject)dm),&data->snes);CHKERRQ(ierr);
    ierr = SNESSetFunction(data->snes,NULL,data->residual_function,(void*)system);CHKERRQ(ierr);
    ierr = SNESSetJacobian(data->snes,NULL,NULL,data->jacobian_function,(void*)system);CHKERRQ(ierr);
    {
      PetscMPIInt size;
      KSP         ksp;
      PC          pc;

      ierr = SNESSetType(data->snes,SNESKSPONLY);CHKERRQ(ierr);
      ierr = SNESGetKSP(data->snes,&ksp);CHKERRQ(ierr);
      ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);
      ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

      ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm),&size);CHKERRQ(ierr);
      if (size == 1) {
        if (PetscDefined(HAVE_SUITESPARSE)) {
          ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
          ierr = PCFactorSetMatSolverType(pc,MATSOLVERUMFPACK);CHKERRQ(ierr);
        }
      } else {
        if (PetscDefined(HAVE_SUPERLU_DIST)) {
          ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
          ierr = PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU_DIST);CHKERRQ(ierr);
        }
      }
    }
    ierr = SNESSetFromOptions(data->snes);CHKERRQ(ierr);
  }

  ierr = SNESSolve(data->snes,NULL,vec_sol);CHKERRQ(ierr);
  {
    SNESConvergedReason reason;

    ierr = SNESGetConvergedReason(data->snes,&reason);CHKERRQ(ierr);
    if (reason < 0) SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_CONV_FAILED,"Solve failed: %s",SNESConvergedReasons[reason]);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemCreate_PETSc(StagBLSystem system)
{
  PetscErrorCode     ierr;
  StagBLSystem_PETSc *data;
  DM                 dm;

  PetscFunctionBegin;
  system->data = (void*) malloc(sizeof(StagBLSystem_PETSc));
  system->ops->destroy = StagBLSystemDestroy_PETSc;
  system->ops->operatorsetvaluesstencil = StagBLSystemOperatorSetValuesStencil_PETSc;
  system->ops->rhssetconstant = StagBLSystemRHSSetConstant_PETSc;
  system->ops->rhssetvaluesstencil = StagBLSystemRHSSetValuesStencil_PETSc;
  system->ops->solve = StagBLSystemSolve_PETSc;

  data = (StagBLSystem_PETSc*) system->data;
  data->residual_function = StagBLSystemPetscResidual_Default;
  data->jacobian_function = StagBLSystemPetscJacobian_Default;
  data->snes = NULL;

  ierr = StagBLGridPETScGetDM(system->grid,&dm);CHKERRQ(ierr);
  ierr = DMCreateMatrix(dm,&data->mat);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&data->rhs);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemPETScGetJacobianFunction(StagBLSystem stagblsystem,PetscErrorCode (**jacobian_function)(SNES,Vec,Mat,Mat,void*))
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBegin;
  *jacobian_function = data->jacobian_function;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemPETScGetMat(StagBLSystem stagblsystem,Mat *mat)
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBegin;
  *mat = data->mat;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemPETScGetMatPointer(StagBLSystem stagblsystem,Mat **mat)
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBegin;
  *mat = &(data->mat);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemPETScGetResidualFunction(StagBLSystem stagblsystem,PetscErrorCode (**residual_function)(SNES,Vec,Vec,void*))
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBegin;
  *residual_function = data->residual_function;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemPETScGetVec(StagBLSystem stagblsystem,Vec *vec)
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBegin;
  *vec = data->rhs;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemPETScGetVecPointer(StagBLSystem stagblsystem,Vec **vec)
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBegin;
  *vec = &(data->rhs);
  PetscFunctionReturn(0);
}
