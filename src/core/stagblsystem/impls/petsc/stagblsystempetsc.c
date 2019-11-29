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
  ierr = MatMult(data->mat,x,f);CHKERRQ(ierr);
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
  // TODO this seems bogus, but we'll soon enough replace this with routines that actually compute the function/jacobian
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
  free(stagblsystem->data);
  stagblsystem->data = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemCreate_PETSc(StagBLSystem stagblsystem)
{
  StagBLSystem_PETSc *data;

  PetscFunctionBegin;
  stagblsystem->data = (void*) malloc(sizeof(StagBLSystem_PETSc));
  data = (StagBLSystem_PETSc*) stagblsystem->data;
  data->mat = NULL;
  data->rhs = NULL;
  data->residual_function = StagBLSystemPetscResidual_Default;
  data->jacobian_function = StagBLSystemPetscJacobian_Default;
  stagblsystem->ops->destroy = StagBLSystemDestroy_PETSc;
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
