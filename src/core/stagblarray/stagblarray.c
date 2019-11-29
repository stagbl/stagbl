#include "stagbl/private/stagblarrayimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLArrayCreate(StagBLGrid grid, StagBLArray *stagblarray)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,stagblarray);CHKERRQ(ierr);
  ierr = PetscCalloc1(1,&(*stagblarray)->ops);CHKERRQ(ierr);

  (*stagblarray)->grid = grid;

  // Setting Type and calling creation routine hard-coded for now
  (*stagblarray)->type = STAGBLARRAYPETSC;
  (*stagblarray)->ops->create = StagBLArrayCreate_PETSc; // Sets other ops
  ierr = ((*stagblarray)->ops->create)(*stagblarray);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayDestroy(StagBLArray *stagblarray)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*stagblarray) PetscFunctionReturn(0);
  if ((*stagblarray)->ops->destroy) {
    ierr = ((*stagblarray)->ops->destroy)(*stagblarray);CHKERRQ(ierr);
  }
  ierr = PetscFree((*stagblarray)->ops);CHKERRQ(ierr);
  ierr = PetscFree(*stagblarray);CHKERRQ(ierr);
  *stagblarray = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayGetStagBLGrid(StagBLArray stagblarray,StagBLGrid *grid)
{
  PetscFunctionBegin;
  *grid = stagblarray->grid;
  PetscFunctionReturn(0);
}
