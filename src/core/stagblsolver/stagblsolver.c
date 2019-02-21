#include "stagbl/private/stagblsolverimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLSolverCreate(StagBLSystem system,StagBLSolver *stagblsolver)
{
  *stagblsolver = malloc(sizeof(struct _p_StagBLSolver));
  (*stagblsolver)->ops = calloc(1,sizeof(struct _p_StagBLSolverOps));

  (*stagblsolver)->system = system;

  // Setting Type and calling creation routine hard-coded for now
  (*stagblsolver)->type = STAGBLSOLVERPETSC;
  (*stagblsolver)->ops->create = StagBLSolverCreate_PETSc; // Sets other ops
  ((*stagblsolver)->ops->create)(*stagblsolver);
  return 0;
}

StagBLErrorCode StagBLSolverDestroy(StagBLSolver *solver)
{
  if (!*solver) return 0;
  if ((*solver)->ops->destroy) {
    ((*solver)->ops->destroy)(*solver);
  }
  free((*solver)->ops);
  free(*solver);
  *solver = NULL;
  return 0;
}

StagBLErrorCode StagBLSolverSolve(StagBLSolver stagblsolver, StagBLArray sol)
{
  // TODO check that solver's stored grid and sol's stored grid are compatible and have the same number of dof (or could even be stronger and assert that they must be the same object..)
  
  if (stagblsolver->ops->solve) {
    (stagblsolver->ops->solve)(stagblsolver,sol);
  } else {
    StagBLError(MPI_COMM_SELF,"StagBLSolverSolver not implemented for this type");
  }
  return 0;
}
