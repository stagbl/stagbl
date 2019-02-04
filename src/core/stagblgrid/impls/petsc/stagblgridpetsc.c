#include "stagbl/private/stagblgridimpl.h"
#include "stagblgridpetscimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLGridDestroy_PETSc(StagBLGrid stagblgrid)
{
  StagBLGrid_PETSc *data = (StagBLGrid_PETSc*) stagblgrid->data;
  if (data->dm) {
    DMDestroy(&data->dm);
  }
  free(stagblgrid->data);
  stagblgrid->data = NULL;
  return 0;
}

StagBLErrorCode StagBLGridCreateCompatibleStagBLGrid_PETSc(StagBLGrid grid,StagBLInt dof0,StagBLInt dof1,StagBLInt dof2,StagBLInt dof3,StagBLGrid *newgrid)
{
  StagBLGrid_PETSc *data,*dataNew;
  PetscErrorCode   ierr;
  DM               coordinateDM;

  data = (StagBLGrid_PETSc*) grid->data;
  StagBLGridCreate(newgrid);
  dataNew = (StagBLGrid_PETSc*) (*newgrid)->data;
  ierr = DMStagCreateCompatibleDMStag(data->dm,dof0,dof1,dof2,dof3,&(dataNew->dm));CHKERRQ(ierr);
  ierr = DMSetUp(dataNew->dm);CHKERRQ(ierr);

  /* Use the same coordinate DM (Different from PETSc's behavior) */
  ierr = DMGetCoordinateDM(data->dm,&coordinateDM);CHKERRQ(ierr);
  ierr = DMSetCoordinateDM(dataNew->dm,coordinateDM);CHKERRQ(ierr);
  return 0;
}

StagBLErrorCode StagBLGridCreateStagBLArray_PETSc(StagBLGrid grid,StagBLArray *array)
{
  StagBLErrorCode  ierr;

  ierr = StagBLArrayCreate(grid,array);CHKERRQ(ierr); // The default type is PETSc, so we do nothing special for now
  return 0;
}

StagBLErrorCode StagBLGridPETScGetDM(StagBLGrid stagblgrid,DM *dm)
{
  StagBLGrid_PETSc * const data = (StagBLGrid_PETSc*) stagblgrid->data;
  *dm = (data->dm);
  return 0;
}

StagBLErrorCode StagBLGridPETScGetDMPointer(StagBLGrid stagblgrid,DM **dm)
{
  StagBLGrid_PETSc * const data = (StagBLGrid_PETSc*) stagblgrid->data;
  *dm = &(data->dm);
  return 0;
}

StagBLErrorCode StagBLGridCreate_PETSc(StagBLGrid stagblgrid)
{
  StagBLGrid_PETSc *data;
  stagblgrid->data = (void*) malloc(sizeof(StagBLGrid_PETSc));
  data = (StagBLGrid_PETSc*) stagblgrid->data;
  data->dm = NULL;
  stagblgrid->ops->createcompatiblestagblgrid = StagBLGridCreateCompatibleStagBLGrid_PETSc;
  stagblgrid->ops->createstagblarray          = StagBLGridCreateStagBLArray_PETSc;
  stagblgrid->ops->destroy                    = StagBLGridDestroy_PETSc;
  return 0;
}
