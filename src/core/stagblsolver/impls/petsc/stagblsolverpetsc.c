#include "../../stagbllinearsolverimpl.h"
#include "stagbllinearsolverpetscimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLLinearSolverDestroy_PETSc(StagBLLinearSolver stagbllinearsolver)
{
  StagBLLinearSolver_PETSc *data = (StagBLLinearSolver_PETSc*) stagbllinearsolver->data;
  if (data->ksp) {
    KSPDestroy(&data->ksp);
  }
  free(stagbllinearsolver->data);
  stagbllinearsolver->data = NULL;

  return 0;
}

StagBLErrorCode StagBLLinearSolverCreate_PETSc(StagBLLinearSolver stagbllinearsolver)
{
  StagBLLinearSolver_PETSc *data;
  stagbllinearsolver->data = (void*) malloc(sizeof(StagBLLinearSolver_PETSc));
  data = (StagBLLinearSolver_PETSc*) stagbllinearsolver->data;
  data->ksp = NULL;
  stagbllinearsolver->ops->destroy = StagBLLinearSolverDestroy_PETSc;

  return 0;
}

StagBLErrorCode StagBLLinearSolverPETScGetKSPPointer(StagBLLinearSolver stagbllinearsolver,KSP **ksp)
{
  StagBLLinearSolver_PETSc * const data = (StagBLLinearSolver_PETSc*) stagbllinearsolver->data;
  *ksp = &(data->ksp);

  return 0;
}
