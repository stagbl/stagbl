#include "stagbl/private/stagblsolverimpl.h"
#include "stagblsolversimpleimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLSolverDestroy_Simple(StagBLSolver solver)
{
  PetscFunctionBegin;
  free(solver->data);
  solver->data = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSolverSolve_Simple(StagBLSolver solver,StagBLArray sol)
{
  PetscErrorCode              ierr;

  PetscFunctionBegin;
  // FIXME this is a kludge - I StagBLSolver should just be rolled into StagBLSystem for now.
  ierr = StagBLSystemSolve_Simple(solver->system,sol);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSolverCreate_Simple(StagBLSolver solver)
{
  //StagBLSolver_Simple *data;

  PetscFunctionBegin;
  solver->data = (void*) malloc(sizeof(StagBLSolver_Simple));
  //data = (StagBLSolver_Simple*) solver->data;
  solver->ops->destroy = StagBLSolverDestroy_Simple;
  solver->ops->solve = StagBLSolverSolve_Simple;
  PetscFunctionReturn(0);
}
