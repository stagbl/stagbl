#include "stagbl/private/stagblsystemimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLSystemCreate(StagBLGrid grid,StagBLSystem *system)
{
  StagBLErrorCode ierr;

  *system = malloc(sizeof(struct _p_StagBLSystem));
  (*system)->ops = calloc(1,sizeof(struct _p_StagBLSystemOps));

  (*system)->grid = grid;

  // Setting Type and calling creation routine hard-coded for now
  (*system)->type = STAGBLSYSTEMPETSC;
  (*system)->ops->create = StagBLSystemCreate_PETSc; // Sets other ops
  ierr = ((*system)->ops->create)(*system);CHKERRQ(ierr);
  return 0;
}

StagBLErrorCode StagBLSystemCreateStagBLSolver(StagBLSystem system,StagBLSolver *solver)
{
  StagBLErrorCode ierr;
  ierr = StagBLSolverCreate(system,solver);CHKERRQ(ierr);
  return 0;
}

StagBLErrorCode StagBLSystemDestroy(StagBLSystem *system)
{
  StagBLErrorCode ierr;

  if (!*system) return 0;
  if ((*system)->ops->destroy) {
    ierr = ((*system)->ops->destroy)(*system);CHKERRQ(ierr);
  }
  free((*system)->ops);
  free(*system);
  *system = NULL;
  return 0;
}

StagBLErrorCode StagBLSystemGetGrid(StagBLSystem system,StagBLGrid *grid)
{
  *grid = system->grid;
  return 0;
}
