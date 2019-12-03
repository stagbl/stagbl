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

PetscErrorCode DumpStokes(Ctx ctx,PetscInt timestep)
{
  PetscErrorCode ierr;
  DM             dmStokes,dmCoeff;
  DM             dmVelAvg;
  Vec            vecx;
  Vec            velAvg;
  DM             daVelAvg,daP,daEtaElement,daEtaCorner,daRho;
  Vec            vecVelAvg,vecP,vecEtaElement,vecEtaCorner,vecRho;
  Vec            vecCoeff,vecCoeffLocal;
  StagBLArray    x;

  PetscFunctionBeginUser;

  /* Use the "escape hatch" */
  x = ctx->stokes_array;
  ierr = StagBLArrayPETScGetGlobalVec(x,&vecx);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecx,"solution");CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(ctx->stokes_grid,&dmStokes);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(ctx->coefficient_grid,&dmCoeff);CHKERRQ(ierr);
  ierr = StagBLArrayPETScGetLocalVec(ctx->coefficient_array,&vecCoeffLocal);CHKERRQ(ierr);

  /* For convenience, create a new DM and Vec which will hold averaged velocities
     Note that this could also be accomplished with direct array access, using
     DMStagVecGetArrayDOF() and related functions */
  ierr = DMStagCreateCompatibleDMStag(dmStokes,0,0,2,0,&dmVelAvg); /* 2 dof per element */
  ierr = DMSetUp(dmVelAvg);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesProduct(dmVelAvg,0.0,ctx->xmax,0.0,ctx->ymax,0.0,0.0);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dmVelAvg,&velAvg);CHKERRQ(ierr);
  {
    PetscInt    ex,ey,startx,starty,nx,ny;
    Vec         stokesLocal,velAvgLocal;
    PetscInt    iVxLeft,iVxRight,iVyDown,iVyUp,iVxCenter,iVyCenter;
    PetscScalar ***arrStokes,***arrVelAvg;

    ierr = DMStagGetLocationSlot(dmStokes,DMSTAG_LEFT,   0,&iVxLeft  );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmStokes,DMSTAG_RIGHT,  0,&iVxRight );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmStokes,DMSTAG_DOWN,   0,&iVyDown  );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmStokes,DMSTAG_UP,     0,&iVyUp    );CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmVelAvg,DMSTAG_ELEMENT,0,&iVxCenter);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmVelAvg,DMSTAG_ELEMENT,1,&iVyCenter);CHKERRQ(ierr);
    ierr = DMStagGetCorners(dmVelAvg,&startx,&starty,NULL,&nx,&ny,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMGetLocalVector(dmStokes,&stokesLocal);CHKERRQ(ierr);
    ierr = DMGetLocalVector(dmVelAvg,&velAvgLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocal(dmStokes,vecx,INSERT_VALUES,stokesLocal);CHKERRQ(ierr);
    ierr = DMStagVecGetArrayRead(dmStokes,stokesLocal,&arrStokes);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(    dmVelAvg,velAvgLocal,&arrVelAvg);CHKERRQ(ierr);
    for (ey = starty; ey<starty+ny; ++ey) {
      for (ex = startx; ex<startx+nx; ++ex) {
        arrVelAvg[ey][ex][iVxCenter] = 0.5 * (arrStokes[ey][ex][iVxLeft] + arrStokes[ey][ex][iVxRight]);
        arrVelAvg[ey][ex][iVyCenter] = 0.5 * (arrStokes[ey][ex][iVyDown] + arrStokes[ey][ex][iVyUp]);
      }
    }
    ierr = DMStagVecRestoreArrayRead(dmStokes,stokesLocal,&arrStokes);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(   dmVelAvg, velAvgLocal,&arrVelAvg);CHKERRQ(ierr);
    ierr = DMLocalToGlobal(dmVelAvg,velAvgLocal,INSERT_VALUES,velAvg);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmStokes,&stokesLocal);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dmVelAvg,&velAvgLocal);CHKERRQ(ierr);
  }

  /* Create a global coefficient vector (otherwise not needed) */
  ierr = DMGetGlobalVector(dmCoeff,&vecCoeff);CHKERRQ(ierr);
  ierr = DMLocalToGlobal(dmCoeff,vecCoeffLocal,INSERT_VALUES,vecCoeff);CHKERRQ(ierr);

  /* Split to DMDAs */
  ierr = DMStagVecSplitToDMDA(dmStokes,vecx,DMSTAG_ELEMENT,0,&daP,&vecP);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecP,"p (scaled)");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(dmCoeff,vecCoeff,DMSTAG_DOWN_LEFT,0,&daEtaCorner,&vecEtaCorner);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecEtaCorner,"eta");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(dmCoeff,vecCoeff,DMSTAG_ELEMENT,0,&daEtaElement,&vecEtaElement);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecEtaElement,"eta");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(dmCoeff,vecCoeff,DMSTAG_DOWN_LEFT,1,&daRho,&vecRho);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecRho,"density");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(dmVelAvg,velAvg,DMSTAG_ELEMENT,-3,&daVelAvg,&vecVelAvg);CHKERRQ(ierr); /* note -3 : pad with zero */
  ierr = PetscObjectSetName((PetscObject)vecVelAvg,"Velocity (Averaged)");CHKERRQ(ierr);

  ierr = DMRestoreGlobalVector(dmCoeff,&vecCoeff);CHKERRQ(ierr);

  /* Dump element-based fields to a .vtr file */
  {
    PetscViewer viewer;
    char        filename[PETSC_MAX_PATH_LEN];

    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"out_element_%.4D.vtr",timestep);CHKERRQ(ierr);
    ierr = PetscViewerVTKOpen(PetscObjectComm((PetscObject)daVelAvg),filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(vecVelAvg,viewer);CHKERRQ(ierr);
    ierr = VecView(vecP,viewer);CHKERRQ(ierr);
    ierr = VecView(vecEtaElement,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Dump vertex-based fields to a second .vtr file */
  {
    PetscViewer viewer;
    char        filename[PETSC_MAX_PATH_LEN];

    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"out_vertex_%.4D.vtr",timestep);CHKERRQ(ierr);
    ierr = PetscViewerVTKOpen(PetscObjectComm((PetscObject)daEtaCorner),filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(vecEtaCorner,viewer);CHKERRQ(ierr);
    ierr = VecView(vecRho,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Edge-based fields could similarly be dumped */

  /* For testing, option to dump the solution to an ASCII file */
  if (timestep == 0) {
    PetscBool debug_ascii_dump = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,NULL,"-debug_ascii_dump",&debug_ascii_dump,NULL);CHKERRQ(ierr);
    if (debug_ascii_dump) {
      PetscViewer viewer;
      PetscViewerASCIIOpen(PetscObjectComm((PetscObject)daVelAvg),"x.matlabascii.txt",&viewer);CHKERRQ(ierr);
      PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
      ierr = VecView(vecx,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }
  }

  /* Destroy DMDAs and Vecs */
  ierr = VecDestroy(&vecVelAvg);CHKERRQ(ierr);
  ierr = VecDestroy(&vecP);CHKERRQ(ierr);
  ierr = VecDestroy(&vecEtaCorner);CHKERRQ(ierr);
  ierr = VecDestroy(&vecEtaElement);CHKERRQ(ierr);
  ierr = VecDestroy(&vecRho);CHKERRQ(ierr);
  ierr = VecDestroy(&velAvg);CHKERRQ(ierr);
  ierr = DMDestroy(&daVelAvg);CHKERRQ(ierr);
  ierr = DMDestroy(&daP);CHKERRQ(ierr);
  ierr = DMDestroy(&daEtaCorner);CHKERRQ(ierr);
  ierr = DMDestroy(&daEtaElement);CHKERRQ(ierr);
  ierr = DMDestroy(&daRho);CHKERRQ(ierr);
  ierr = DMDestroy(&dmVelAvg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DumpTemperature(Ctx ctx, PetscInt timestep)
{
  PetscErrorCode ierr;
  DM             dmTemp,daTemp;
  Vec            vecTemp,vecTempDa;
  PetscViewer    viewer;
  char           filename[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;

  /* Use the "escape hatch" */
  ierr = StagBLArrayPETScGetGlobalVec(ctx->temperature_array,&vecTemp);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(ctx->temperature_grid,&dmTemp);CHKERRQ(ierr);

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
