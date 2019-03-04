#include "stagbl/private/stagblsolverimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLSolverCreate(StagBLSystem system,StagBLSolver *stagblsolver)
{
  StagBLErrorCode ierr;

  *stagblsolver = malloc(sizeof(struct _p_StagBLSolver));
  (*stagblsolver)->ops = calloc(1,sizeof(struct _p_StagBLSolverOps));

  (*stagblsolver)->system = system;

  // Setting Type and calling creation routine hard-coded for now
  (*stagblsolver)->type = STAGBLSOLVERPETSC;
  (*stagblsolver)->ops->create = StagBLSolverCreate_PETSc; // Sets other ops
  ierr = ((*stagblsolver)->ops->create)(*stagblsolver);CHKERRQ(ierr);
  return 0;
}

StagBLErrorCode StagBLSolverGetSystem(StagBLSolver solver,StagBLSystem *system)
{
  *system = solver->system;
  return 0;
}

StagBLErrorCode StagBLSolverDestroy(StagBLSolver *solver)
{
  StagBLErrorCode ierr;

  if (!*solver) return 0;
  if ((*solver)->ops->destroy) {
    ierr = ((*solver)->ops->destroy)(*solver);CHKERRQ(ierr);
  }
  free((*solver)->ops);
  free(*solver);
  *solver = NULL;
  return 0;
}

StagBLErrorCode StagBLSolverSolve(StagBLSolver stagblsolver, StagBLArray sol)
{
  StagBLErrorCode ierr;
  // TODO check that solver's stored grid and sol's stored grid are compatible and have the same number of dof (or could even be stronger and assert that they must be the same object..)
  
  if (stagblsolver->ops->solve) {
    ierr = (stagblsolver->ops->solve)(stagblsolver,sol);CHKERRQ(ierr);
  } else {
    StagBLError(MPI_COMM_SELF,"StagBLSolverSolver not implemented for this type");
  }
  return 0;
}
