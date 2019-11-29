#include "stagbl/private/stagblsystemimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLSystemCreate(StagBLGrid grid,StagBLSystem *system)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,system);CHKERRQ(ierr);
  ierr = PetscCalloc1(1,&(*system)->ops);CHKERRQ(ierr);

  (*system)->grid = grid;

  // Setting Type and calling creation routine hard-coded for now
  (*system)->type = STAGBLSYSTEMPETSC;
  (*system)->ops->create = StagBLSystemCreate_PETSc; // Sets other ops
  ierr = ((*system)->ops->create)(*system);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemCreateStagBLSolver(StagBLSystem system,StagBLSolver *solver)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = StagBLSolverCreate(system,solver);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemDestroy(StagBLSystem *system)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*system) PetscFunctionReturn(0);
  if ((*system)->ops->destroy) {
    ierr = ((*system)->ops->destroy)(*system);CHKERRQ(ierr);
  }
  ierr = PetscFree((*system)->ops);CHKERRQ(ierr);
  ierr = PetscFree(*system);CHKERRQ(ierr);
  *system = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemGetGrid(StagBLSystem system,StagBLGrid *grid)
{
  PetscFunctionBegin;
  *grid = system->grid;
  PetscFunctionReturn(0);
}
