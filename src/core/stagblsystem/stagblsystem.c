#include "stagbl/private/stagblsystemimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLSystemCreate(StagBLGrid grid,StagBLSystem *system)
{
  PetscErrorCode ierr;

  *system = malloc(sizeof(struct _p_StagBLSystem));
  (*system)->ops = calloc(1,sizeof(struct _p_StagBLSystemOps));

  (*system)->grid = grid;

  // Setting Type and calling creation routine hard-coded for now
  (*system)->type = STAGBLSYSTEMPETSC;
  (*system)->ops->create = StagBLSystemCreate_PETSc; // Sets other ops
  ierr = ((*system)->ops->create)(*system);CHKERRQ(ierr);
  return 0;
}

PetscErrorCode StagBLSystemCreateStagBLSolver(StagBLSystem system,StagBLSolver *solver)
{
  PetscErrorCode ierr;
  ierr = StagBLSolverCreate(system,solver);CHKERRQ(ierr);
  return 0;
}

PetscErrorCode StagBLSystemDestroy(StagBLSystem *system)
{
  PetscErrorCode ierr;

  if (!*system) return 0;
  if ((*system)->ops->destroy) {
    ierr = ((*system)->ops->destroy)(*system);CHKERRQ(ierr);
  }
  free((*system)->ops);
  free(*system);
  *system = NULL;
  return 0;
}

PetscErrorCode StagBLSystemGetGrid(StagBLSystem system,StagBLGrid *grid)
{
  *grid = system->grid;
  return 0;
}
