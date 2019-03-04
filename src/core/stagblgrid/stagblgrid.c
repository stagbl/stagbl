#include "stagbl/private/stagblgridimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLGridCreate(StagBLGrid *stagblgrid)
{
  StagBLErrorCode ierr;

  *stagblgrid = malloc(sizeof(struct _p_StagBLGrid));
  (*stagblgrid)->ops = calloc(1,sizeof(struct _p_StagBLGridOps));

  // Setting Type and calling creation routine hard-coded for now
  (*stagblgrid)->type = STAGBLGRIDPETSC;
  (*stagblgrid)->ops->create = StagBLGridCreate_PETSc; // Sets other ops
  ierr = ((*stagblgrid)->ops->create)(*stagblgrid);CHKERRQ(ierr);
  return 0;
}

StagBLErrorCode StagBLGridCreateCompatibleStagBLGrid(StagBLGrid grid,StagBLInt dof0,StagBLInt dof1,StagBLInt dof2,StagBLInt dof3,StagBLGrid *newgrid)
{
  StagBLErrorCode ierr;

  if (grid->ops->createcompatiblestagblgrid) {
    ierr = (grid->ops->createcompatiblestagblgrid)(grid,dof0,dof1,dof2,dof3,newgrid);CHKERRQ(ierr);
  } else {
    StagBLError(MPI_COMM_SELF,"StagBLCreateCompatibleStagBLGrid not implemented for this type");
  }
  return 0;
}

StagBLErrorCode StagBLGridCreateStagBLArray(StagBLGrid grid,StagBLArray *array)
{
  StagBLErrorCode ierr;

  if (grid->ops->createstagblarray) {
    ierr = (grid->ops->createstagblarray)(grid,array);CHKERRQ(ierr);
  } else {
    StagBLError(MPI_COMM_SELF,"StagBLGridCreateStagBLArray not implemented for this type");
  }
  return 0;
} 

StagBLErrorCode StagBLGridCreateStagBLSystem(StagBLGrid grid,StagBLSystem *system)
{
  StagBLErrorCode ierr;
  // TODO this is basically a placeholder, creating a generic system
  ierr = StagBLSystemCreate(grid,system);CHKERRQ(ierr);
  return 0;
} 

StagBLErrorCode StagBLGridDestroy(StagBLGrid *stagblgrid)
{
  StagBLErrorCode ierr;

  if (!*stagblgrid) return 0;
  if ((*stagblgrid)->ops->destroy) {
    ierr = ((*stagblgrid)->ops->destroy)(*stagblgrid);CHKERRQ(ierr);
  }
  free((*stagblgrid)->ops);
  free(*stagblgrid);
  *stagblgrid = NULL;
  return 0;
}
