#include "stagbl.h"

StagBLErrorCode StagBLGridCreateStokes2DBox(MPI_Comm comm, StagBLInt nx,StagBLInt ny,StagBLReal xmin, StagBLReal xmax, StagBLReal ymin, StagBLReal ymax,StagBLGrid *pgrid)
{
  // TODO this function assumes PETSc is included, and that the defaults types for things are PETSc
  PetscErrorCode ierr;
  DM *pdm;
  DM dmStokes;

  StagBLGridCreate(pgrid);
  StagBLGridPETScGetDMPointer(*pgrid,&pdm);
  ierr = DMStagCreate2d(
      comm,
      DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
      nx,ny,                                   /* Global element counts */
      PETSC_DECIDE,PETSC_DECIDE,               /* Determine parallel decomposition automatically */
      0,1,1,                                   /* dof: 0 per vertex, 1 per edge, 1 per face/element */
      DMSTAG_STENCIL_BOX,
      1,                                       /* elementwise stencil width */
      NULL,NULL,
      pdm);
  dmStokes = *pdm;
  ierr = DMSetFromOptions(dmStokes);CHKERRQ(ierr);
  ierr = DMSetUp(dmStokes);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesProduct(dmStokes,xmin,xmax,ymin,ymax,0.0,0.0);CHKERRQ(ierr);
  return 0;
}

