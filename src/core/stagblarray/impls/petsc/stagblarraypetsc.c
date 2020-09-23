#include "stagbl/private/stagblarrayimpl.h"
#include "stagblarraypetscimpl.h"
#include <stdlib.h>

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


PetscErrorCode StagBLArraySetLocalConstant_PETSc(StagBLArray array, PetscScalar value)
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
  stagblarray->ops->print = StagBLArrayPrint_PETSc;
  stagblarray->ops->setlocalconstant = StagBLArraySetLocalConstant_PETSc;
  return 0;
}

PetscErrorCode StagBLArrayPETScGetGlobalVec(StagBLArray stagblarray,Vec *vec)
{
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) stagblarray->data;
  *vec = data->global;
  return 0;
}

PetscErrorCode StagBLArrayPETScGetLocalVec(StagBLArray stagblarray,Vec *vec)
{
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) stagblarray->data;
  *vec = data->local;
  return 0;
}

PetscErrorCode StagBLArrayPETScGetLocalVecPointer(StagBLArray stagblarray,Vec **vec)
{
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) stagblarray->data;
  *vec = &(data->local);
  return 0;
}

PetscErrorCode StagBLArrayPETScGetGlobalVecPointer(StagBLArray stagblarray,Vec **vec)
{
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) stagblarray->data;
  *vec = &(data->global);
  return 0;
}
