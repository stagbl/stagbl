#include "dump.h"

PetscErrorCode DumpParticles(Ctx ctx,PetscInt timestep)
{
  PetscErrorCode ierr;
  char           filename[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"particles_%.4D.xmf",timestep);CHKERRQ(ierr);
  ierr = DMSwarmViewXDMF(ctx->dm_particles,filename);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DumpTemperature(Ctx ctx, PetscInt timestep)
{
  PetscErrorCode ierr;
  PetscInt       dim;
  DM             dmTemp,daTemp;
  Vec            vecTemp,vecTempDa;
  PetscViewer    viewer;
  char           filename[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  ierr = StagBLArrayPETScGetGlobalVec(ctx->temperature_array,&vecTemp);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(ctx->temperature_grid,&dmTemp);CHKERRQ(ierr);
  ierr = DMGetDimension(dmTemp,&dim);CHKERRQ(ierr);
  if (dim != 2) StagBLError(PetscObjectComm((PetscObject)dmTemp),"Only implemented for 2D");

  /* Split to DMDA */
  ierr = DMStagVecSplitToDMDA(dmTemp,vecTemp,DMSTAG_DOWN_LEFT,0,&daTemp,&vecTempDa);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecTempDa,"Temperature");CHKERRQ(ierr);

  /* View */
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"out_temp_vertex_%.4D.vtr",timestep);CHKERRQ(ierr);
  ierr = PetscViewerVTKOpen(PetscObjectComm((PetscObject)daTemp),filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(vecTempDa,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* Destroy DMDAs and Vecs */
  ierr = VecDestroy(&vecTempDa);CHKERRQ(ierr);
  ierr = DMDestroy(&daTemp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
