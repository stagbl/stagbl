#include "stagbl/private/stagblsystemimpl.h"
#include "stagblsystempetscimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLSystemDestroy_PETSc(StagBLSystem stagblsystem)
{
  PetscErrorCode     ierr;
  StagBLSystem_PETSc *data = (StagBLSystem_PETSc*) stagblsystem->data;

  PetscFunctionBeginUser;
  if (data->rhs) {
    ierr = VecDestroy(&data->rhs);CHKERRQ(ierr);
  }
  if (data->mat) {
    ierr = MatDestroy(&data->mat);CHKERRQ(ierr);
  }
  free(stagblsystem->data);
  stagblsystem->data = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemCreate_PETSc(StagBLSystem stagblsystem)
{
  StagBLSystem_PETSc *data;

  stagblsystem->data = (void*) malloc(sizeof(StagBLSystem_PETSc));
  data = (StagBLSystem_PETSc*) stagblsystem->data;
  data->mat = NULL;
  data->rhs = NULL;
  stagblsystem->ops->destroy = StagBLSystemDestroy_PETSc;
  return 0;
}

PetscErrorCode StagBLSystemPETScGetMat(StagBLSystem stagblsystem,Mat *mat)
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;
  *mat = data->mat;
  return 0;
}

PetscErrorCode StagBLSystemPETScGetMatPointer(StagBLSystem stagblsystem,Mat **mat)
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;
  *mat = &(data->mat);
  return 0;
}

PetscErrorCode StagBLSystemPETScGetVec(StagBLSystem stagblsystem,Vec *vec)
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;
  *vec = data->rhs;
  return 0;
}

PetscErrorCode StagBLSystemPETScGetVecPointer(StagBLSystem stagblsystem,Vec **vec)
{
  StagBLSystem_PETSc * const data = (StagBLSystem_PETSc*) stagblsystem->data;
  *vec = &(data->rhs);
  return 0;
}
