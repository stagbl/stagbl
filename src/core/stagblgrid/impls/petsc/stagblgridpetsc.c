#include "stagbl/private/stagblgridimpl.h"
#include "stagblgridpetscimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLGridDestroy_PETSc(StagBLGrid stagblgrid)
{
  StagBLGrid_PETSc *data = (StagBLGrid_PETSc*) stagblgrid->data;
  if (data->dm) {
    DMDestroy(&data->dm);
  }
  free(stagblgrid->data);
  stagblgrid->data = NULL;
  return 0;
}

PetscErrorCode StagBLGridCreateCompatibleStagBLGrid_PETSc(StagBLGrid grid,PetscInt dof0,PetscInt dof1,PetscInt dof2,PetscInt dof3,StagBLGrid *newgrid)
{
  StagBLGrid_PETSc *data,*dataNew;
  PetscErrorCode   ierr;
  DM               coordinateDM;
  PetscInt         dof0_from, dof1_from, dof2_from, dof3_from;

  data = (StagBLGrid_PETSc*) grid->data;
  StagBLGridCreate(newgrid);
  dataNew = (StagBLGrid_PETSc*) (*newgrid)->data;
  ierr = DMStagCreateCompatibleDMStag(data->dm,dof0,dof1,dof2,dof3,&(dataNew->dm));CHKERRQ(ierr);
  ierr = DMSetUp(dataNew->dm);CHKERRQ(ierr);

  /* Use the same coordinate DM (Different from PETSc's behavior) */
  ierr = DMStagGetDOF(data->dm,&dof0_from,&dof1_from,&dof2_from,&dof3_from);CHKERRQ(ierr);
  if (  // Check that the same strata are active (non-zero dof)
      ((dof0 == 0) == (dof0_from ==0)) ||
      ((dof1 == 0) == (dof1_from ==0)) ||
      ((dof2 == 0) == (dof2_from ==0)) ||
      ((dof3 == 0) == (dof3_from ==0)) 
     ) {
    ierr = DMGetCoordinateDM(data->dm,&coordinateDM);CHKERRQ(ierr);
    ierr = DMSetCoordinateDM(dataNew->dm,coordinateDM);CHKERRQ(ierr);
  } else {
    SETERRQ(PetscObjectComm((PetscObject)data->dm),PETSC_ERR_SUP,"Coordinate DM creation not implemented when active strata change");
  }
  return 0;
}

PetscErrorCode StagBLGridCreateStagBLArray_PETSc(StagBLGrid grid,StagBLArray *array)
{
  PetscErrorCode  ierr;

  ierr = StagBLArrayCreate(grid,array);CHKERRQ(ierr); // The default type is PETSc, so we do nothing special for now
  return 0;
}

PetscErrorCode StagBLGridPETScGetDM(StagBLGrid stagblgrid,DM *dm)
{
  StagBLGrid_PETSc * const data = (StagBLGrid_PETSc*) stagblgrid->data;
  *dm = (data->dm);
  return 0;
}

PetscErrorCode StagBLGridPETScGetDMPointer(StagBLGrid stagblgrid,DM **dm)
{
  StagBLGrid_PETSc * const data = (StagBLGrid_PETSc*) stagblgrid->data;
  *dm = &(data->dm);
  return 0;
}

PetscErrorCode StagBLGridCreate_PETSc(StagBLGrid stagblgrid)
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
