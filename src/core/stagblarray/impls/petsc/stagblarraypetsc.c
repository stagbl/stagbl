#include "stagbl/private/stagblarrayimpl.h"
#include "stagblarraypetscimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLArrayDestroy_PETSc(StagBLArray stagblarray)
{
  StagBLArray_PETSc *data= (StagBLArray_PETSc*) stagblarray->data;
  if (data->vecLocal) {
    VecDestroy(&data->vecLocal);
  }
  if (data->vecGlobal) {
    VecDestroy(&data->vecGlobal);
  }
  free(stagblarray->data);
  stagblarray->data = NULL;
  return 0;
}

PetscErrorCode StagBLArrayCreate_PETSc(StagBLArray stagblarray)
{
  StagBLArray_PETSc *data;
  stagblarray->data = (void*) malloc(sizeof(StagBLArray_PETSc));
  data = (StagBLArray_PETSc*) stagblarray->data;
  data->vecLocal  = NULL;
  data->vecGlobal = NULL;
  stagblarray->ops->destroy = StagBLArrayDestroy_PETSc;
  return 0;
}

PetscErrorCode StagBLArrayPETScGetGlobalVec(StagBLArray stagblarray,Vec *vec)
{
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) stagblarray->data;
  *vec = data->vecGlobal;
  return 0;
}

PetscErrorCode StagBLArrayPETScGetLocalVec(StagBLArray stagblarray,Vec *vec)
{
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) stagblarray->data;
  *vec = data->vecLocal;
  return 0;
}

PetscErrorCode StagBLArrayPETScGetLocalVecPointer(StagBLArray stagblarray,Vec **vec)
{
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) stagblarray->data;
  *vec = &(data->vecLocal);
  return 0;
}

PetscErrorCode StagBLArrayPETScGetGlobalVecPointer(StagBLArray stagblarray,Vec **vec)
{
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) stagblarray->data;
  *vec = &(data->vecGlobal);
  return 0;
}
