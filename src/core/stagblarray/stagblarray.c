#include "stagbl/private/stagblarrayimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLArrayCreate(StagBLGrid grid, StagBLArray *stagblarray)
{
  PetscErrorCode ierr;
  *stagblarray = malloc(sizeof(struct _p_StagBLArray));
  (*stagblarray)->ops = calloc(1,sizeof(struct _p_StagBLArrayOps));

  (*stagblarray)->grid = grid;

  // Setting Type and calling creation routine hard-coded for now
  (*stagblarray)->type = STAGBLARRAYPETSC;
  (*stagblarray)->ops->create = StagBLArrayCreate_PETSc; // Sets other ops
  ierr = ((*stagblarray)->ops->create)(*stagblarray);CHKERRQ(ierr);
  return 0;
}

PetscErrorCode StagBLArrayDestroy(StagBLArray *stagblarray)
{
  PetscErrorCode ierr;
  if (!*stagblarray) return 0;
  if ((*stagblarray)->ops->destroy) {
    ierr = ((*stagblarray)->ops->destroy)(*stagblarray);CHKERRQ(ierr);
  }
  free((*stagblarray)->ops);
  free(*stagblarray);
  *stagblarray = NULL;
  return 0;
}

PetscErrorCode StagBLArrayGetStagBLGrid(StagBLArray stagblarray,StagBLGrid *grid)
{
  *grid = stagblarray->grid;
  return 0;
}
