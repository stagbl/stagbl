#include "stagbl/private/stagblarrayimpl.h"
#include <stdlib.h>

/**
  * Note: usually one would use StagBLGridCreateArray(), which calls this function.
  */
PetscErrorCode StagBLArrayCreate(StagBLGrid grid, StagBLArray *stagblarray, StagBLArrayType array_type)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,stagblarray);CHKERRQ(ierr);
  ierr = PetscCalloc1(1,&(*stagblarray)->ops);CHKERRQ(ierr);

  (*stagblarray)->grid = grid;
  (*stagblarray)->type = array_type;

  /* Set the creation function and call it, which sets other ops */
  if (StagBLCheckType(array_type,STAGBLARRAYPETSC)) {
      (*stagblarray)->ops->create = StagBLArrayCreate_PETSc;
  } else StagBLError1(PETSC_COMM_WORLD,"Array creation not implemented for type %s",array_type);
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
