#include "stagbl/private/stagblgridimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLGridCreate(StagBLGrid *stagblgrid, StagBLGridType grid_type)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,stagblgrid);CHKERRQ(ierr);
  ierr = PetscCalloc1(1,&(*stagblgrid)->ops);CHKERRQ(ierr);

  (*stagblgrid)->type = grid_type;
  if (!grid_type  || StagBLCheckType(grid_type, STAGBLGRIDPETSC)) {
    (*stagblgrid)->ops->create = StagBLGridCreate_PETSc;
  } else StagBLError1(PETSC_COMM_WORLD,"Unsupport StagBLGridType %s",grid_type);
  ierr = ((*stagblgrid)->ops->create)(*stagblgrid);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLGridCreateCompatibleStagBLGrid(StagBLGrid grid,PetscInt dof0,PetscInt dof1,PetscInt dof2,PetscInt dof3,StagBLGrid *newgrid)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (grid->ops->createcompatiblestagblgrid) {
    ierr = (grid->ops->createcompatiblestagblgrid)(grid,dof0,dof1,dof2,dof3,newgrid);CHKERRQ(ierr);
  } else {
    StagBLError(MPI_COMM_SELF,"StagBLCreateCompatibleStagBLGrid not implemented for this type");
  }
  ierr = StagBLGridSetArrayType(*newgrid,grid->array_type);CHKERRQ(ierr);
  ierr = StagBLGridSetSystemType(*newgrid,grid->system_type);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLGridCreateStagBLArray(StagBLGrid grid,StagBLArray *array)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (grid->ops->createstagblarray) {
    ierr = (grid->ops->createstagblarray)(grid,array);CHKERRQ(ierr);
  } else {
    StagBLError(MPI_COMM_SELF,"StagBLGridCreateStagBLArray not implemented for this type");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLGridCreateStagBLSystem(StagBLGrid grid,StagBLSystem *system)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = StagBLSystemCreate(grid,system,grid->system_type);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLGridDestroy(StagBLGrid *stagblgrid)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*stagblgrid) PetscFunctionReturn(0);
  if ((*stagblgrid)->ops->destroy) {
    ierr = ((*stagblgrid)->ops->destroy)(*stagblgrid);CHKERRQ(ierr);
  }
  ierr = PetscFree((*stagblgrid)->ops);CHKERRQ(ierr);
  ierr = PetscFree(*stagblgrid);CHKERRQ(ierr);
  *stagblgrid = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLGridSetArrayType(StagBLGrid grid, StagBLArrayType array_type)
{
  PetscFunctionBegin;
  grid->array_type = array_type;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLGridSetSystemType(StagBLGrid grid, StagBLArrayType system_type)
{
  PetscFunctionBegin;
  grid->system_type = system_type;
  PetscFunctionReturn(0);
}
