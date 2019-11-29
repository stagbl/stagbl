#include "stagbl/private/stagblgridimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLGridCreate(StagBLGrid *stagblgrid)
{
  PetscErrorCode ierr;

  ierr = PetscMalloc1(1,stagblgrid);CHKERRQ(ierr);
  ierr = PetscCalloc1(1,&(*stagblgrid)->ops);CHKERRQ(ierr);

  // Setting Type and calling creation routine hard-coded for now
  (*stagblgrid)->type = STAGBLGRIDPETSC;
  (*stagblgrid)->ops->create = StagBLGridCreate_PETSc; // Sets other ops
  ierr = ((*stagblgrid)->ops->create)(*stagblgrid);CHKERRQ(ierr);
  return 0;
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
  ierr = StagBLSystemCreate(grid,system);CHKERRQ(ierr);
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
