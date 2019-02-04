#include "stagbl/private/stagblsystemimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLSystemCreate(StagBLGrid grid,StagBLSystem *system)
{
  *system = malloc(sizeof(struct _p_StagBLSystem));
  (*system)->ops = calloc(1,sizeof(struct _p_StagBLSystemOps));

  (*system)->grid = grid;

  // Setting Type and calling creation routine hard-coded for now
  (*system)->type = STAGBLSYSTEMPETSC;
  (*system)->ops->create = StagBLSystemCreate_PETSc; // Sets other ops
  ((*system)->ops->create)(*system);
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
  if ((*system)->ops->destroy) {
    ((*system)->ops->destroy)(*system);
  }
  free((*system)->ops);
  free(*system);
  *system = NULL;
  return 0;
}
