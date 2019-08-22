#include "particles.h"

PetscErrorCode InterpolateTemperatureToParticles(Ctx ctx)
{
  PetscErrorCode ierr;
  DM             dmTemp,dm_mpoint;
  Vec            temp,tempLocal;
  PetscInt       p,e,npoints;
  PetscInt       *mpfield_cell;
  PetscReal      *array_temperature;
  PetscScalar    ***arr;
  PetscInt       slotDownLeft,slotDownRight,slotUpLeft,slotUpRight;

  PetscFunctionBeginUser;
  dm_mpoint = ctx->dm_particles;
  ierr = StagBLGridPETScGetDM(ctx->temperature_grid,&dmTemp);CHKERRQ(ierr);
  ierr = StagBLArrayPETScGetGlobalVec(ctx->temperature_array,&temp);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dmTemp,&tempLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocal(dmTemp,temp,INSERT_VALUES,tempLocal);CHKERRQ(ierr);

  ierr = DMStagGetLocationSlot(dmTemp,DMSTAG_DOWN_LEFT,0,&slotDownLeft);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmTemp,DMSTAG_DOWN_RIGHT,0,&slotDownRight);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmTemp,DMSTAG_UP_RIGHT,0,&slotUpRight);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmTemp,DMSTAG_UP_LEFT,0,&slotUpLeft);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayRead(dmTemp,tempLocal,&arr);CHKERRQ(ierr);
  ierr = DMSwarmGetField(dm_mpoint,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);
  ierr = DMSwarmGetField(ctx->dm_particles,"Temperature",NULL,NULL,(void**)&array_temperature);CHKERRQ(ierr);
  ierr = DMSwarmGetLocalSize(dm_mpoint,&npoints);CHKERRQ(ierr);
  for (p=0; p<npoints; p++) {
    PetscInt    ind[2];

    e       = mpfield_cell[p];
    ierr = DMStagGetLocalElementGlobalIndices(dmTemp,e,ind);CHKERRQ(ierr);

    /* simply average corners (could do linear interp) */
    array_temperature[p] = 0.25*(
            arr[ind[1]][ind[0]][slotDownLeft]
          + arr[ind[1]][ind[0]][slotDownRight]
          + arr[ind[1]][ind[0]][slotUpLeft]
          + arr[ind[1]][ind[0]][slotUpRight]
        );
  }
  ierr = DMSwarmRestoreField(ctx->dm_particles,"Temperature",NULL,NULL,(void**)&array_temperature);CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(dm_mpoint,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(dmTemp,&tempLocal);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestTeleport(Ctx ctx)
{
  PetscErrorCode ierr;
  DM             dmCell;
  PetscInt       N[2],ex,ey,exloc,eyloc,n[2],start[2];
  PetscInt       p,npoints;
  PetscInt       *mpfield_cell;
  PetscReal      *mpfield_coor;
  PetscReal      **massBefore,**massAfter,**displacement;
  PetscReal      *array_temperature;

  PetscFunctionBeginUser;
  ierr = StagBLGridPETScGetDM(ctx->stokes_grid,&dmCell);CHKERRQ(ierr);
  ierr = DMSwarmGetLocalSize(ctx->dm_particles,&npoints);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(dmCell,&N[0],&N[1],NULL);CHKERRQ(ierr);
  ierr = DMStagGetCorners(dmCell,&start[0],&start[1],NULL,&n[0],&n[1],NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  /* These arrays hold mass before and after teleportation, and displacements
     at each boundary in y */
  ierr = PetscMalloc3(n[0],&massBefore,n[0],&massAfter,n[0],&displacement);CHKERRQ(ierr);
  for (eyloc=0; eyloc<n[0]; ++eyloc) {
    ierr = PetscCalloc3(n[1],&massBefore[eyloc],n[1],&massAfter[eyloc],n[1]+1,&displacement[eyloc]);CHKERRQ(ierr);
  }


  ierr = DMSwarmGetField(ctx->dm_particles,DMSwarmPICField_coor,NULL,NULL,(void**)&mpfield_coor);CHKERRQ(ierr);
  ierr = DMSwarmGetField(ctx->dm_particles,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);
  ierr = DMSwarmGetField(ctx->dm_particles,"Temperature",NULL,NULL,(void**)&array_temperature);CHKERRQ(ierr);

  /* Iterate over all particles.
     For each, if it's the third particle in its element, move it.
     Increment a count of particles before and after */
  for (p=0; p<npoints; p++) {
    PetscReal   *coor_p;
    PetscInt    ind[2],e;

    e       = mpfield_cell[p];
    coor_p  = &mpfield_coor[2*p];
    ierr = DMStagGetLocalElementGlobalIndices(dmCell,e,ind);CHKERRQ(ierr);
    eyloc = ind[1]-start[1];
    exloc = ind[0]-start[0];

    /* Accumulate */
    /* Note the order of the indices in the arrays (ex first), unlike DMStag */
    massBefore[exloc][eyloc] += 1.0;
    /* Teleport the third particle, in a given band */
    if (massBefore[exloc][eyloc] == 3.0 && ind[1] > 0.4*N[1] && ind[1] < N[1]-1) {
      coor_p[1] = ctx->ymax - 0.5 * (ctx->ymax-ctx->ymin)/N[1];
      array_temperature[p] = 0;
    } else {
      massAfter[exloc][eyloc] += 1.0;
    }
  }
  ierr = DMSwarmRestoreField(ctx->dm_particles,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(ctx->dm_particles,DMSwarmPICField_coor,NULL,NULL,(void**)&mpfield_coor);CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(ctx->dm_particles,"Temperature",NULL,NULL,(void**)&array_temperature);CHKERRQ(ierr);

  /* Iterate over the auxiliary arrays.
     Compute a local displacement for the top edge relative to the bottom edge,
     by seeing how much mass was removed. Use a columnwise MPI communicator to perform
     a prefix sum, to see how much the top and bottom of each element should
     be displaced */
  // TODO
  {
    MPI_Comm          comm_column;
    PetscMPIInt       rankx;
    const PetscScalar **arrCoordY;
    PetscInt          slotPrev,slotNext;
    const PetscMPIInt key = 0;


    ierr = DMStagGetProductCoordinateLocationSlot(dmCell,DMSTAG_LEFT,&slotPrev);CHKERRQ(ierr);
    ierr = DMStagGetProductCoordinateLocationSlot(dmCell,DMSTAG_RIGHT,&slotNext);CHKERRQ(ierr);
    ierr = DMStagGetProductCoordinateArraysRead(dmCell,NULL,&arrCoordY,NULL);CHKERRQ(ierr);
    ierr = DMStagGetRank(dmCell,&rankx,NULL,NULL);CHKERRQ(ierr);
    ierr = MPI_Comm_split(PetscObjectComm((PetscObject)dmCell),key,rankx,&comm_column);CHKERRQ(ierr);
    for (ex = start[0]; ex<start[0]+n[0]; ++ex) { /* Note loop bounds reverse from normal */
      displacement[exloc][0] = 0.0; /* locally, no displacement at bottom of column */
      for (ey = start[1]; ey<start[1]+n[1]; ++ey) {
        const PetscReal scale = 1.0 - (massAfter[exloc][eyloc]/massBefore[exloc][eyloc]);
        const PetscReal height = arrCoordY[ey][slotNext] - arrCoordY[ey][slotPrev];

        exloc = ex-start[0]; eyloc = ey-start[1];
        displacement[exloc][eyloc+1] = displacement[exloc][eyloc] - scale*height; /* note sign */
      }
    }
    // TODO MPI logic not in place - should break for more than one rank in Y direction
    ierr = MPI_Comm_free(&comm_column);CHKERRQ(ierr);
  }

  /* Iterate over particles again, displacing based on the boundaries of the
     (currently) containing cell */
  ierr = DMSwarmGetField(ctx->dm_particles,DMSwarmPICField_coor,NULL,NULL,(void**)&mpfield_coor);CHKERRQ(ierr);
  ierr = DMSwarmGetField(ctx->dm_particles,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);

  /* Iterate over all particles.
     For each, if it's the third particle in its element, move it.
     Increment a count of particles before and after */
  for (p=0; p<npoints; p++) {
    PetscReal   *coor_p;
    PetscInt    ind[2],e;

    e       = mpfield_cell[p];
    coor_p  = &mpfield_coor[2*p];
    ierr = DMStagGetLocalElementGlobalIndices(dmCell,e,ind);CHKERRQ(ierr);
    eyloc = ind[1]-start[1];
    exloc = ind[0]-start[0];

    /* Displace down */
    coor_p[1] += displacement[exloc][eyloc];
  }
  ierr = DMSwarmRestoreField(ctx->dm_particles,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(ctx->dm_particles,DMSwarmPICField_coor,NULL,NULL,(void**)&mpfield_coor);CHKERRQ(ierr);

  for (eyloc=0; eyloc<n[0]; ++eyloc) { /* eyloc abused for local element number */
    ierr = PetscFree3(massBefore[eyloc],massAfter[eyloc],displacement[eyloc]);CHKERRQ(ierr);
  }
  ierr = PetscFree3(massBefore,massAfter,displacement);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscReal RangeMod(PetscReal a,PetscReal xmin,PetscReal xmax) { PetscReal range = xmax-xmin; return xmin +PetscFmodReal(range+PetscFmodReal(a,range),range); }

PetscErrorCode MaterialPoint_AdvectRK1(Ctx ctx,Vec vp,PetscReal dt)
{
  PetscErrorCode    ierr;
  Vec               vp_l;
  const PetscScalar ***LA_vp;
  PetscInt          p,e,npoints;
  PetscInt          *mpfield_cell;
  PetscReal         *mpfield_coor;
  PetscScalar       **cArrX,**cArrY;
  PetscInt          iPrev,iNext,iCenter,iVxLeft,iVxRight,iVyDown,iVyUp,n[2];
  DMBoundaryType    boundaryTypes[2];
  PetscBool         periodicx;
  DM                dm_vp,dm_mpoint;

  PetscFunctionBeginUser;
  ierr =StagBLGridPETScGetDM(ctx->stokes_grid,&dm_vp);CHKERRQ(ierr);
  dm_mpoint = ctx->dm_particles;
  ierr = DMStagGetBoundaryTypes(dm_vp,&boundaryTypes[0],&boundaryTypes[1],NULL);CHKERRQ(ierr);
  if (boundaryTypes[1] == DM_BOUNDARY_PERIODIC) SETERRQ(PetscObjectComm((PetscObject)dm_vp),PETSC_ERR_SUP,"Periodic boundary conditions in y-direction not supported");
  periodicx = boundaryTypes[0] == DM_BOUNDARY_PERIODIC;

  ierr = DMStagGetProductCoordinateArraysRead(dm_vp,&cArrX,&cArrY,NULL);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dm_vp,DMSTAG_ELEMENT,&iCenter);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dm_vp,DMSTAG_LEFT,&iPrev);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dm_vp,DMSTAG_RIGHT,&iNext);CHKERRQ(ierr);

  ierr = DMStagGetLocationSlot(dm_vp,DMSTAG_LEFT,0,&iVxLeft);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm_vp,DMSTAG_RIGHT,0,&iVxRight);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm_vp,DMSTAG_DOWN,0,&iVyDown);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm_vp,DMSTAG_UP,0,&iVyUp);CHKERRQ(ierr);

  ierr = DMGetLocalVector(dm_vp,&vp_l);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm_vp,vp,INSERT_VALUES,vp_l);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm_vp,vp,INSERT_VALUES,vp_l);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayRead(dm_vp,vp_l,&LA_vp);CHKERRQ(ierr);

  ierr = DMStagGetLocalSizes(dm_vp,&n[0],&n[1],NULL);CHKERRQ(ierr);
  ierr = DMSwarmGetLocalSize(dm_mpoint,&npoints);CHKERRQ(ierr);
  ierr = DMSwarmGetField(dm_mpoint,DMSwarmPICField_coor,NULL,NULL,(void**)&mpfield_coor);CHKERRQ(ierr);
  ierr = DMSwarmGetField(dm_mpoint,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);
  for (p=0; p<npoints; p++) {
    PetscReal   *coor_p;
    PetscScalar vel_p[2],vLeft,vRight,vUp,vDown;
    PetscScalar x0[2],dx[2],xloc_p[2],xi_p[2];
    PetscInt    ind[2];

    e       = mpfield_cell[p];
    coor_p  = &mpfield_coor[2*p];
    ierr = DMStagGetLocalElementGlobalIndices(dm_vp,e,ind);CHKERRQ(ierr);

    /* compute local coordinates: (xp-x0)/dx = (xip+1)/2 */
    x0[0] = cArrX[ind[0]][iPrev];
    x0[1] = cArrY[ind[1]][iPrev];

    dx[0] = cArrX[ind[0]][iNext] - x0[0];
    dx[1] = cArrY[ind[1]][iNext] - x0[1];

    xloc_p[0] = (coor_p[0] - x0[0])/dx[0];
    xloc_p[1] = (coor_p[1] - x0[1])/dx[1];

    /* Checks (xi_p is only used for this, here) */
    xi_p[0] = 2.0 * xloc_p[0] -1.0;
    xi_p[1] = 2.0 * xloc_p[1] -1.0;
    if (PetscRealPart(xi_p[0]) < -1.0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"value (xi) too small %1.4e [e=%D]\n",(double)PetscRealPart(xi_p[0]),e);
    if (PetscRealPart(xi_p[0]) >  1.0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"value (xi) too large %1.4e [e=%D]\n",(double)PetscRealPart(xi_p[0]),e);
    if (PetscRealPart(xi_p[1]) < -1.0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"value (eta) too small %1.4e [e=%D]\n",(double)PetscRealPart(xi_p[1]),e);
    if (PetscRealPart(xi_p[1]) >  1.0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"value (eta) too large %1.4e [e=%D]\n",(double)PetscRealPart(xi_p[1]),e);

    /* interpolate velocity */
    vLeft  = LA_vp[ind[1]][ind[0]][iVxLeft];
    vRight = LA_vp[ind[1]][ind[0]][iVxRight];
    vUp    = LA_vp[ind[1]][ind[0]][iVyUp];
    vDown  = LA_vp[ind[1]][ind[0]][iVyDown];
    vel_p[0] = xloc_p[0]*vRight + (1.0-xloc_p[0])*vLeft;
    vel_p[1] = xloc_p[1]*vUp    + (1.0-xloc_p[1])*vDown;

    /* Update Coordinates */
    coor_p[0] += dt * PetscRealPart(vel_p[0]);
    coor_p[1] += dt * PetscRealPart(vel_p[1]);

    /* Wrap in periodic directions */
    if (periodicx) {
      coor_p[0] = RangeMod(coor_p[0],ctx->xmin,ctx->xmax);CHKERRQ(ierr);
    }
  }

  ierr = DMSwarmRestoreField(dm_mpoint,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(dm_mpoint,DMSwarmPICField_coor,NULL,NULL,(void**)&mpfield_coor);CHKERRQ(ierr);
  ierr = DMStagVecRestoreArrayRead(dm_vp,vp_l,&LA_vp);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm_vp,&vp_l);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
