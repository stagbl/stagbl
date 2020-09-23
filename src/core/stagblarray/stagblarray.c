#include "stagbl/private/stagblarrayimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLArrayCreate(StagBLGrid grid, StagBLArray *array, StagBLArrayType array_type)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,array);CHKERRQ(ierr);
  ierr = PetscCalloc1(1,&(*array)->ops);CHKERRQ(ierr);

  (*array)->type = array_type;
  (*array)->grid = grid;
  (*array)->current_local = PETSC_FALSE;
  (*array)->current_global = PETSC_FALSE;

  /* Set the creation function and call it, which sets other ops */
  if (StagBLCheckType(array_type,STAGBLARRAYPETSC)) {
      (*array)->ops->create = StagBLArrayCreate_PETSc;
  } else if (StagBLCheckType(array_type,STAGBLARRAYSIMPLE)) {
      (*array)->ops->create = StagBLArrayCreate_Simple;
  } else StagBLError1(PETSC_COMM_WORLD,"Array creation not implemented for type %s",array_type);
  ierr = ((*array)->ops->create)(*array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayDestroy(StagBLArray *array)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*array) PetscFunctionReturn(0);
  if ((*array)->ops->destroy) {
    ierr = ((*array)->ops->destroy)(*array);CHKERRQ(ierr);
  }
  ierr = PetscFree((*array)->ops);CHKERRQ(ierr);
  ierr = PetscFree(*array);CHKERRQ(ierr);
  *array = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayGetStagBLGrid(StagBLArray stagblarray,StagBLGrid *grid)
{
  PetscFunctionBegin;
  *grid = stagblarray->grid;
  PetscFunctionReturn(0);
}

