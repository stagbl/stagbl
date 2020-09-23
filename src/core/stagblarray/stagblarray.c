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

PetscErrorCode StagBLArrayGlobalToLocal(StagBLArray array)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (array->ops->globaltolocal) {
    ierr = (array->ops->globaltolocal)(array);CHKERRQ(ierr);
  } else StagBLError2(PETSC_COMM_WORLD,"%s not implemented for StagBLArray object of type %s",__func__,array->type);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayLocalToGlobal(StagBLArray array)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (array->ops->localtoglobal) {
    ierr = (array->ops->localtoglobal)(array);CHKERRQ(ierr);
  } else StagBLError2(PETSC_COMM_WORLD,"%s not implemented for StagBLArray object of type %s",__func__,array->type);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArraySetLocalConstant(StagBLArray array, PetscScalar value)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (array->ops->setlocalconstant) {
    ierr = (array->ops->setlocalconstant)(array,value);CHKERRQ(ierr);
  } else StagBLError2(PETSC_COMM_WORLD,"%s not implemented for StagBLArray object of type %s",__func__,array->type);
  PetscFunctionReturn(0);
}

/**
 Print the array contents to stdout, for debugging and diagnostic purposes only.
 */
PetscErrorCode StagBLArrayPrint(StagBLArray array)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (array->ops->print) {
    ierr = (array->ops->print)(array);CHKERRQ(ierr);
  } else  {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Printing not implemented\n");CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
