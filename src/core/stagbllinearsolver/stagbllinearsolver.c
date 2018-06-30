#include "stagbllinearsolverimpl.h"
#include <stdlib.h>

void StagBLLinearSolverCreate(StagBLLinearSolver *stagbllinearsolver)
{
  *stagbllinearsolver = malloc(sizeof(struct _p_StagBLLinearSolver));
  (*stagbllinearsolver)->ops = calloc(1,sizeof(struct _p_StagBLLinearSolverOps));

  // Setting Type and calling creation routine hard-coded for now
  (*stagbllinearsolver)->type = STAGBLLINEARSOLVERPETSC;
  (*stagbllinearsolver)->ops->create = StagBLLinearSolverCreate_PETSc; // Sets other ops
  ((*stagbllinearsolver)->ops->create)(*stagbllinearsolver);
}

void StagBLLinearSolverDestroy(StagBLLinearSolver *stagbllinearsolver)
{
  if ((*stagbllinearsolver)->ops->destroy) {
    ((*stagbllinearsolver)->ops->destroy)(*stagbllinearsolver);
  }
  free((*stagbllinearsolver)->ops);
  free(*stagbllinearsolver);
  *stagbllinearsolver = NULL;
}
