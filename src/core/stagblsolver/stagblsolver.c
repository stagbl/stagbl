#include "stagbl/private/stagblsolverimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLSolverCreate(StagBLSystem system,StagBLSolver *solver,StagBLSolverType type)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,solver);CHKERRQ(ierr);
  ierr = PetscCalloc1(1,&(*solver)->ops);CHKERRQ(ierr);

  (*solver)->type = type;
  (*solver)->system = system;

  /* Set the creation function and call it, which sets other ops */
  if (StagBLCheckType(type,STAGBLSOLVERPETSC)) {
      (*solver)->ops->create = StagBLSolverCreate_PETSc;
  } else if (StagBLCheckType(type,STAGBLSOLVERSIMPLE)) {
      (*solver)->ops->create = StagBLSolverCreate_Simple;
  } else StagBLError1(PETSC_COMM_WORLD,"System creation not implemented for type %s",type);
  ierr = ((*solver)->ops->create)(*solver);CHKERRQ(ierr);
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
  if (!*solver) PetscFunctionReturn(0);
  if ((*solver)->ops->destroy) {
    ierr = ((*solver)->ops->destroy)(*solver);CHKERRQ(ierr);
  }
  ierr = PetscFree((*solver)->ops);CHKERRQ(ierr);
  ierr = PetscFree(*solver);CHKERRQ(ierr);
  *solver = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSolverSolve(StagBLSolver solver, StagBLArray solution)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (solver->ops->solve) {
    ierr = (solver->ops->solve)(solver,solution);CHKERRQ(ierr);
  } else StagBLError(MPI_COMM_SELF,"StagBLSolverSolve not implemented for this type");
  PetscFunctionReturn(0);
}
