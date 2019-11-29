#include "stagbl/private/stagblsolverimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLSolverCreate(StagBLSystem system,StagBLSolver *stagblsolver)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,stagblsolver);CHKERRQ(ierr);
  ierr = PetscCalloc1(1,&(*stagblsolver)->ops);CHKERRQ(ierr);

  (*stagblsolver)->system = system;

  // Setting Type and calling creation routine hard-coded for now
  (*stagblsolver)->type = STAGBLSOLVERPETSC;
  (*stagblsolver)->ops->create = StagBLSolverCreate_PETSc; // Sets other ops
  ierr = ((*stagblsolver)->ops->create)(*stagblsolver);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSolverGetSystem(StagBLSolver solver,StagBLSystem *system)
{
  PetscFunctionBegin;
  *system = solver->system;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSolverDestroy(StagBLSolver *solver)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*solver) return 0;
  if ((*solver)->ops->destroy) {
    ierr = ((*solver)->ops->destroy)(*solver);CHKERRQ(ierr);
  }
  ierr = PetscFree((*solver)->ops);CHKERRQ(ierr);
  ierr = PetscFree(*solver);CHKERRQ(ierr);
  *solver = NULL;
  PetscFunctionReturn(0);
}

/**
  * Note: it doesn't make much sense to call this function if the array
  * and the solver's system don't correspond to compatible grids.
  */
PetscErrorCode StagBLSolverSolve(StagBLSolver stagblsolver, StagBLArray sol)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (stagblsolver->ops->solve) {
    ierr = (stagblsolver->ops->solve)(stagblsolver,sol);CHKERRQ(ierr);
  } else {
    StagBLError(MPI_COMM_SELF,"StagBLSolverSolver not implemented for this type");
  }
  PetscFunctionReturn(0);
}
