#include "stagbl/private/stagblarrayimpl.h"
#include "stagblarraypetscimpl.h"
#include <stdlib.h>

/* For DMStagStencilToIndexLocal(), which should be made public */
#include <petsc/private/dmstagimpl.h>

static PetscErrorCode StagBLArrayPETScCreateGlobalVector_Private(StagBLArray);
static PetscErrorCode StagBLArrayPETScCreateLocalVector_Private(StagBLArray);

PetscErrorCode StagBLArrayDestroy_PETSc(StagBLArray stagblarray)
{
  StagBLArray_PETSc *data= (StagBLArray_PETSc*) stagblarray->data;
  if (data->local) {
    VecDestroy(&data->local);
  }
  if (data->global) {
    VecDestroy(&data->global);
  }
  free(stagblarray->data);
  stagblarray->data = NULL;
  return 0;
}

static PetscErrorCode StagBLArrayGetLocalValuesStencil_PETSc(StagBLArray array,PetscInt n,const DMStagStencil *pos,PetscScalar *values)
{
  PetscErrorCode    ierr;
  StagBLArray_PETSc *data = (StagBLArray_PETSc*) array->data;
  DM                dm;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  if (!data->local) StagBLError(PetscObjectComm((PetscObject)dm),"Local array data not defined");
  if (!array->current_local) StagBLError(PetscObjectComm((PetscObject)dm),"Local array data not current");
  ierr = DMStagVecGetValuesStencil(dm,data->local,n,pos,values);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLArrayGlobalToLocal_PETSc(StagBLArray array)
{
  PetscErrorCode    ierr;
  StagBLArray_PETSc *data = (StagBLArray_PETSc*) array->data;
  DM                dm;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  if (!data->global) StagBLError(PetscObjectComm((PetscObject)dm),"Cannot perform global-to-local before global data is defined");
  if (!array->current_global) StagBLError(PetscObjectComm((PetscObject)dm),"Refusing to perform global-to-local with stale global data");
  if (!data->local) {
    ierr = StagBLArrayPETScCreateLocalVector_Private(array);CHKERRQ(ierr);
  }
  ierr = DMGlobalToLocal(dm,data->global,INSERT_VALUES,data->local);CHKERRQ(ierr);
  array->current_local = PETSC_TRUE;
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLArrayLocalToGlobal_PETSc(StagBLArray array)
{
  PetscErrorCode    ierr;
  StagBLArray_PETSc *data = (StagBLArray_PETSc*) array->data;
  DM                dm;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  if (!data->local) StagBLError(PetscObjectComm((PetscObject)dm),"Cannot perform local-to-global before local data is defined");
  if (!array->current_local) StagBLError(PetscObjectComm((PetscObject)dm),"Refusing to perform global-to-local with stale global data");
  if (!data->global) {
    ierr = StagBLArrayPETScCreateGlobalVector_Private(array);CHKERRQ(ierr);
  }
  ierr = DMLocalToGlobal(dm,data->local,INSERT_VALUES,data->global);CHKERRQ(ierr);
  array->current_global = PETSC_TRUE;
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLArrayPrint_PETSc(StagBLArray array)
{
  PetscErrorCode    ierr;
  StagBLArray_PETSc *data = (StagBLArray_PETSc*) array->data;

  PetscFunctionBegin;
  if (data->local) {
    ierr = PetscPrintf(PetscObjectComm((PetscObject)data->local),"local Vec:\n");CHKERRQ(ierr);
    ierr = VecView(data->local,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  if (data->global) {
    ierr = PetscPrintf(PetscObjectComm((PetscObject)data->global),"global Vec:\n");CHKERRQ(ierr);
    ierr = VecView(data->global,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLArrayPETScCreateGlobalVector_Private(StagBLArray array)
{
  PetscErrorCode    ierr;
  StagBLArray_PETSc *data = (StagBLArray_PETSc*) array->data;

  PetscFunctionBegin;
  if (!data->global) {
    DM dm;

    ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm,&data->global);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLArrayPETScCreateLocalVector_Private(StagBLArray array)
{
  PetscErrorCode    ierr;
  StagBLArray_PETSc *data = (StagBLArray_PETSc*) array->data;

  PetscFunctionBegin;
  if (!data->local) {
    DM dm;

    ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(dm,&data->local);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLArraySetLocalValuesStencil_PETSc(StagBLArray array,PetscInt n,const DMStagStencil *pos,const PetscScalar *values)
{
  PetscErrorCode    ierr;
  StagBLArray_PETSc *data = (StagBLArray_PETSc*) array->data;
  DM                dm;
  PetscInt          *indices_local;
  PetscInt          dim;
  PetscScalar       *local_raw;

  PetscFunctionBegin;
  if (!data->local) {
    ierr = StagBLArrayPETScCreateLocalVector_Private(array);CHKERRQ(ierr);
    array->current_local = PETSC_TRUE;
  }
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  indices_local = (PetscInt*) malloc(n * sizeof(PetscInt));
  ierr = DMStagStencilToIndexLocal(dm,dim,n,pos,indices_local);CHKERRQ(ierr);
  ierr = VecGetArray(data->local,&local_raw);CHKERRQ(ierr);
  for (PetscInt i=0; i<n; ++i) local_raw[indices_local[i]] = values[i];
  ierr = VecRestoreArray(data->local,&local_raw);CHKERRQ(ierr);
  free(indices_local);
  PetscFunctionReturn(0);
}


static PetscErrorCode StagBLArraySetLocalConstant_PETSc(StagBLArray array, PetscScalar value)
{
  PetscErrorCode    ierr;
  StagBLArray_PETSc *data = (StagBLArray_PETSc*) array->data;

  PetscFunctionBegin;
  if (!data->local) {
    ierr = StagBLArrayPETScCreateLocalVector_Private(array);CHKERRQ(ierr);
  }
  ierr = VecSet(data->local,value);CHKERRQ(ierr);
  array->current_local = PETSC_TRUE;
  array->current_global = PETSC_FALSE;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayCreate_PETSc(StagBLArray stagblarray)
{
  StagBLArray_PETSc *data;
  stagblarray->data = (void*) malloc(sizeof(StagBLArray_PETSc));
  data = (StagBLArray_PETSc*) stagblarray->data;
  data->local  = NULL;
  data->global = NULL;
  stagblarray->ops->destroy = StagBLArrayDestroy_PETSc;
  stagblarray->ops->getlocalvaluesstencil = StagBLArrayGetLocalValuesStencil_PETSc;
  stagblarray->ops->globaltolocal = StagBLArrayGlobalToLocal_PETSc;
  stagblarray->ops->localtoglobal = StagBLArrayLocalToGlobal_PETSc;
  stagblarray->ops->print = StagBLArrayPrint_PETSc;
  stagblarray->ops->setlocalvaluesstencil = StagBLArraySetLocalValuesStencil_PETSc;
  stagblarray->ops->setlocalconstant = StagBLArraySetLocalConstant_PETSc;
  return 0;
}

PetscErrorCode StagBLArrayPETScGetGlobalVec(StagBLArray array,Vec *vec)
{
  PetscErrorCode            ierr;
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) array->data;

  PetscFunctionBeginHot;
  ierr = StagBLEnforceType(array->type,STAGBLARRAYPETSC);CHKERRQ(ierr);
  *vec = data->global;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayPETScGetLocalVec(StagBLArray array,Vec *vec)
{
  PetscErrorCode            ierr;
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) array->data;

  PetscFunctionBeginHot;
  ierr = StagBLEnforceType(array->type,STAGBLARRAYPETSC);CHKERRQ(ierr);
  *vec = data->local;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayPETScGetLocalVecPointer(StagBLArray array,Vec **vec)
{
  PetscErrorCode            ierr;
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) array->data;

  PetscFunctionBegin;
  ierr = StagBLEnforceType(array->type,STAGBLARRAYPETSC);CHKERRQ(ierr);
  *vec = &(data->local);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayPETScGetGlobalVecPointer(StagBLArray array,Vec **vec)
{
  PetscErrorCode            ierr;
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) array->data;

  PetscFunctionBegin;
  ierr = StagBLEnforceType(array->type,STAGBLARRAYPETSC);CHKERRQ(ierr);
  *vec = &(data->global);
  PetscFunctionReturn(0);
}
