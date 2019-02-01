#include "stagbl/private/stagblgridimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLGridCreate(StagBLGrid *stagblgrid)
{
  *stagblgrid = malloc(sizeof(struct _p_StagBLGrid));
  (*stagblgrid)->ops = calloc(1,sizeof(struct _p_StagBLGridOps));

  // Setting Type and calling creation routine hard-coded for now
  (*stagblgrid)->type = STAGBLGRIDPETSC;
  (*stagblgrid)->ops->create = StagBLGridCreate_PETSc; // Sets other ops
  ((*stagblgrid)->ops->create)(*stagblgrid);
  return 0;
}

StagBLErrorCode StagBLGridCreateCompatibleStagBLGrid(StagBLGrid grid,StagBLInt dof0,StagBLInt dof1,StagBLInt dof2,StagBLInt dof3,StagBLGrid *newgrid)
{
  if (grid->ops->createcompatiblestagblgrid) {
    (grid->ops->createcompatiblestagblgrid)(grid,dof0,dof1,dof2,dof3,newgrid);
  } else {
    StagBLError(MPI_COMM_SELF,"StagBLCreateCompatibleStagBLGrid not implemented for this type");
  }
  return 0;
}

StagBLErrorCode StagBLGridCreateStagBLArray(StagBLGrid grid,StagBLArray *array)
{
  if (grid->ops->createstagblarray) {
    (grid->ops->createstagblarray)(grid,array);
  } else {
    StagBLError(MPI_COMM_SELF,"StagBLGridCreateStagBLArray not implemented for this type");
  }
  return 0;
} 

StagBLErrorCode StagBLGridCreateStagBLSystem(StagBLGrid grid,StagBLSystem *system)
{
  StagBLErrorCode ierr;
  // TODO this is basically a placeholder, createing a generic system
  ierr = StagBLSystemCreate(system);CHKERRQ(ierr);
  return 0;
} 

StagBLErrorCode StagBLGridDestroy(StagBLGrid *stagblgrid)
{
  if ((*stagblgrid)->ops->destroy) {
    ((*stagblgrid)->ops->destroy)(*stagblgrid);
  }
  free((*stagblgrid)->ops);
  free(*stagblgrid);
  *stagblgrid = NULL;
  return 0;
}
