#include "temperature.h"

PetscErrorCode InitializeTemperature(Ctx ctx)
{
  PetscErrorCode ierr;
  Vec            *pTemp;
  Vec            temp,tempLocal;
  DM             dmTemp;
  PetscInt       ex,ey,startx,starty,nx,ny,nExtrax,nExtray,slot,N[2];
  PetscScalar    ***arr;

  PetscFunctionBeginUser;
  if (!ctx->temperature_array) {
    ierr = StagBLGridCreateStagBLArray(ctx->temperature_grid,&ctx->temperature_array);CHKERRQ(ierr);
  }
  ierr = StagBLGridPETScGetDM(ctx->temperature_grid,&dmTemp);CHKERRQ(ierr);
  ierr = StagBLArrayPETScGetGlobalVecPointer(ctx->temperature_array,&pTemp);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dmTemp,pTemp);CHKERRQ(ierr);
  temp = *pTemp;

  ierr = DMGetLocalVector(dmTemp,&tempLocal);CHKERRQ(ierr);
  ierr = DMStagGetCorners(dmTemp,&startx,&starty,NULL,&nx,&ny,NULL,&nExtrax,&nExtray,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(dmTemp,&N[0],&N[1],NULL);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmTemp,DMSTAG_DOWN_LEFT,0,&slot);CHKERRQ(ierr);
  ierr = DMStagVecGetArray(dmTemp,tempLocal,&arr);CHKERRQ(ierr);
  /* A simple inital temperature field, for the Blankenbach benchmark */
  for (ey = starty; ey<starty+ny+nExtray; ++ey) {
    for (ex = startx; ex<startx+nx+nExtrax; ++ex) {
      if (ex == 0 || ex < (PetscInt) (0.1 * N[0])){
          arr[ey][ex][slot] = 1000;
      } else if (ex == N[0]-1 || ex > (PetscInt) (0.9 * N[0])){
          arr[ey][ex][slot] = 0;
      } else {
          arr[ey][ex][slot] = 500;
      }
    }
  }
  ierr = DMStagVecRestoreArray(dmTemp,tempLocal,&arr);CHKERRQ(ierr);
  ierr = DMLocalToGlobal(dmTemp,tempLocal,INSERT_VALUES,temp);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dmTemp,&tempLocal);CHKERRQ(ierr);

// TODO check mode here and set something coordinate-based (say give bouyant hot anomaly on the left!)
  PetscFunctionReturn(0);
}

/**
  * Create a system corresponding to an implicit (backwards) Euler step
  * for the energy equation.
  */
PetscErrorCode PopulateTemperatureSystem(Ctx ctx)
{
  PetscErrorCode  ierr;
  DM              dmTemp,dmStokes;
  PetscInt        N[2];
  PetscInt        ex,ey,startx,starty,nx,ny,nExtrax,nExtray;
  PetscBool       isLastx,isLasty,isFirstx,isFirsty;
  Mat             *pA;
  Vec             *pRhs;
  Mat             A;
  Vec             rhs;
  PetscReal       hx,hy;
  Vec             tempPrev,tempPrevLocal;
  PetscScalar     ***arr;
  PetscInt        slot;
  PetscScalar     ***arrStokes;
  PetscInt        slotVxLeft,slotVyDown;
  Vec             stokes,stokesLocal;

  PetscFunctionBeginUser;
  ierr = StagBLSystemPETScGetMatPointer(ctx->temperature_system,&pA);CHKERRQ(ierr);
  ierr = StagBLSystemPETScGetVecPointer(ctx->temperature_system,&pRhs);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(ctx->temperature_grid,&dmTemp);CHKERRQ(ierr);

  ierr = StagBLArrayPETScGetGlobalVec(ctx->temperature_array,&tempPrev);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dmTemp,&tempPrevLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocal(dmTemp,tempPrev,INSERT_VALUES,tempPrevLocal);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayRead(dmTemp,tempPrevLocal,&arr);CHKERRQ(ierr);
  slot = 0; // TODO lazy

  ierr = StagBLArrayPETScGetGlobalVec(ctx->stokes_array,&stokes);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(ctx->stokes_grid,&dmStokes);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmStokes,DMSTAG_LEFT,0,&slotVxLeft);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmStokes,DMSTAG_DOWN,0,&slotVyDown);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dmStokes,&stokesLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocal(dmStokes,stokes,INSERT_VALUES,stokesLocal);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayRead(dmStokes,stokesLocal,&arrStokes);CHKERRQ(ierr);

  if (!*pA) {
    ierr = DMCreateMatrix(dmTemp,pA);CHKERRQ(ierr);
  }
  A = *pA;
  if (!*pRhs) {
    ierr = DMCreateGlobalVector(dmTemp,pRhs);CHKERRQ(ierr);
  }
  rhs = *pRhs;
  ierr = DMStagGetCorners(dmTemp,&startx,&starty,NULL,&nx,&ny,NULL,&nExtrax,&nExtray,NULL);CHKERRQ(ierr);
  ierr = DMStagGetIsLastRank(dmTemp,&isLastx,&isLasty,NULL);CHKERRQ(ierr);
  ierr = DMStagGetIsFirstRank(dmTemp,&isFirstx,&isFirsty,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(dmTemp,&N[0],&N[1],NULL);CHKERRQ(ierr);
  hx = ctx->hxCharacteristic;
  hy = ctx->hyCharacteristic;
  if (ctx->uniform_grid) {
    hx = (ctx->xmax-ctx->xmin)/N[0];
    hy = (ctx->ymax-ctx->ymin)/N[1];
  } else StagBLError(PetscObjectComm((PetscObject)dmTemp),"Non-uniform grids not supported yet");

  /* Top boundary - fixed temperature */
  if (isLasty) {
    for (ex = startx; ex<startx+nx+nExtrax; ++ex) {
      DMStagStencil row;
      PetscScalar   valA,valRhs;

      const PetscInt ey = N[1]; /* The "extra" row */

      row.i = ex; row.j = ey; row.c = 0; row.loc = DMSTAG_DOWN_LEFT;
      valA = ctx->KTemp;
      ierr = DMStagMatSetValuesStencil(dmTemp,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
      valRhs = ctx->KTemp * ctx->temperature_top;
      ierr = DMStagVecSetValuesStencil(dmTemp,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  /* Bottom boundary - fixed temperature */
  if (isFirsty) {
    for (ex = startx; ex<startx+nx+nExtrax; ++ex) {
      DMStagStencil row;
      PetscScalar   valA,valRhs;

      const PetscInt ey = 0;

      row.i = ex; row.j = ey; row.c = 0; row.loc = DMSTAG_DOWN_LEFT;
      valA = ctx->KTemp;
      ierr = DMStagMatSetValuesStencil(dmTemp,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
      valRhs = ctx->KTemp * ctx->temperature_bottom;
      ierr = DMStagVecSetValuesStencil(dmTemp,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  /* Insulating right boundary */
  if (isLastx) {
    for (ey = PetscMax(1,starty); ey<starty+ny; ++ey) { /* Skip corners */
      DMStagStencil row,col[4];
      PetscScalar   valA[4],valRhs;
      PetscScalar   vx,vy,avx,avy,vxp,vyp,vym;

      const PetscInt ex = N[0]; /* the "extra" column of elements */

      row.i = ex; row.j = ey; row.c = 0; row.loc = DMSTAG_DOWN_LEFT;

      /* Diffusion terms */
      col[0].i = ex;   col[0].j = ey;   col[0].c = 0; col[0].loc=DMSTAG_DOWN_LEFT; valA[0] =  ctx->kappa*(1.0/(hx*hx) + 2.0/(hy*hy));
      col[1].i = ex-1; col[1].j = ey;   col[1].c = 0; col[1].loc=DMSTAG_DOWN_LEFT; valA[1] = -1.0*ctx->kappa/(hx*hx);
      /* Missing right */
      col[2].i = ex;   col[2].j = ey-1; col[2].c = 0; col[2].loc=DMSTAG_DOWN_LEFT; valA[2] = -1.0*ctx->kappa/(hy*hy);
      col[3].i = ex;   col[3].j = ey+1; col[3].c = 0; col[3].loc=DMSTAG_DOWN_LEFT; valA[3] = -1.0*ctx->kappa/(hy*hy);

      /* Additional terms from upwinding (note +=) */
      vx = 0.5*(arrStokes[ey][ex][slotVxLeft] + arrStokes[ey-1][ex][slotVxLeft]);
      avx = PetscAbsScalar(vx); vxp = 0.5*(vx+avx);
      vy = arrStokes[ey][ex-1][slotVyDown]; /* only one entry (assume free slip) */
      avy = PetscAbsScalar(vy); vyp = 0.5*(vy+avy); vym = 0.5*(vy-avy);
      valA[0] +=  avx/hx + avy/hy;
      valA[1] += -vxp/hx;
      /* Missing Right */ // TODO WRONG forgot something??
      valA[2] += -vyp/hy;
      valA[3] +=  vym/hy;

      /* Time Discretization Term */
      valA[0] += 1.0/ctx->dt;

      ierr = DMStagMatSetValuesStencil(dmTemp,A,1,&row,4,col,valA,INSERT_VALUES);CHKERRQ(ierr);
      valRhs = arr[ey][ex][slot]/ctx->dt; /* No heat source */

      ierr = DMStagVecSetValuesStencil(dmTemp,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  /* Insulating left boundary */
  if (isFirstx) {
    for (ey = PetscMax(1,starty); ey<starty+ny; ++ey) { /* Skip corners */
      /* Note: this is duplicated mostly from the interior points section */
      DMStagStencil row,col[4];
      PetscScalar   valA[4],valRhs;
      PetscScalar   vx,vy,avx,avy,vxm,vyp,vym;

      const PetscInt ex = 0;

      row.i = ex; row.j = ey; row.c = 0; row.loc = DMSTAG_DOWN_LEFT;

      /* Diffusion terms */
      col[0].i = ex;   col[0].j = ey;   col[0].c = 0; col[0].loc=DMSTAG_DOWN_LEFT; valA[0] =  ctx->kappa*(1.0/(hx*hx) + 2.0/(hy*hy));
      /* Missing left */
      col[1].i = ex+1; col[1].j = ey;   col[1].c = 0; col[1].loc=DMSTAG_DOWN_LEFT; valA[1] = -1.0*ctx->kappa/(hx*hx);
      col[2].i = ex;   col[2].j = ey-1; col[2].c = 0; col[2].loc=DMSTAG_DOWN_LEFT; valA[2] = -1.0*ctx->kappa/(hy*hy);
      col[3].i = ex;   col[3].j = ey+1; col[3].c = 0; col[3].loc=DMSTAG_DOWN_LEFT; valA[3] = -1.0*ctx->kappa/(hy*hy);

      /* Additional terms from upwinding (note +=) */
      vx = 0.5*(arrStokes[ey][ex][slotVxLeft] + arrStokes[ey-1][ex][slotVxLeft]);
      avx = PetscAbsScalar(vx); vxm = 0.5*(vx-avx);
      vy = arrStokes[ey][ex][slotVyDown]; /* only use one value  (assume free slip) */
      avy = PetscAbsScalar(vy); vyp = 0.5*(vy+avy); vym = 0.5*(vy-avy);
      valA[0] +=  avx/hx + avy/hy;
      /* Missing left */ // TODO WRONG forgot something??
      valA[1] +=  vxm/hx;
      valA[2] += -vyp/hy;
      valA[3] +=  vym/hy;

      /* Time Discretization Term */
      valA[0] += 1.0/ctx->dt;

      ierr = DMStagMatSetValuesStencil(dmTemp,A,1,&row,4,col,valA,INSERT_VALUES);CHKERRQ(ierr);
      valRhs = arr[ey][ex][slot]/ctx->dt; /* No heat source */

      ierr = DMStagVecSetValuesStencil(dmTemp,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  /* Loop over all interior points */
  for (ey = PetscMax(1,starty); ey<starty+ny; ++ey) {
    for (ex = PetscMax(1,startx); ex<startx+nx; ++ex) {
      /* Note: this is duplicated above, mostly, for the boundary cases */
      DMStagStencil row,col[5];
      PetscScalar   valA[5],valRhs;
      PetscScalar   vx,vy,avx,avy,vxp,vxm,vyp,vym;

      row.i = ex; row.j = ey; row.c = 0; row.loc = DMSTAG_DOWN_LEFT;

      /* Diffusion terms */
      col[0].i = ex;   col[0].j = ey;   col[0].c = 0; col[0].loc=DMSTAG_DOWN_LEFT; valA[0] =  ctx->kappa*(2.0/(hx*hx) + 2.0/(hy*hy));
      col[1].i = ex-1; col[1].j = ey;   col[1].c = 0; col[1].loc=DMSTAG_DOWN_LEFT; valA[1] = -1.0*ctx->kappa/(hx*hx);
      col[2].i = ex+1; col[2].j = ey;   col[2].c = 0; col[2].loc=DMSTAG_DOWN_LEFT; valA[2] = -1.0*ctx->kappa/(hx*hx);
      col[3].i = ex;   col[3].j = ey-1; col[3].c = 0; col[3].loc=DMSTAG_DOWN_LEFT; valA[3] = -1.0*ctx->kappa/(hy*hy);
      col[4].i = ex;   col[4].j = ey+1; col[4].c = 0; col[4].loc=DMSTAG_DOWN_LEFT; valA[4] = -1.0*ctx->kappa/(hy*hy);

      /* Additional terms from upwinding (note +=) */
      vx = 0.5*(arrStokes[ey][ex][slotVxLeft] + arrStokes[ey-1][ex][slotVxLeft]);
      avx = PetscAbsScalar(vx); vxp = 0.5*(vx+avx); vxm = 0.5*(vx-avx);
      vy = 0.5*(arrStokes[ey][ex][slotVyDown] + arrStokes[ey][ex-1][slotVyDown]);
      avy = PetscAbsScalar(vy); vyp = 0.5*(vy+avy); vym = 0.5*(vy-avy);
      valA[0] +=  avx/hx + avy/hy;
      valA[1] += -vxp/hx;
      valA[2] +=  vxm/hx;
      valA[3] += -vyp/hy;
      valA[4] +=  vym/hy;

      /* Time Discretization Term */
      valA[0] += 1.0/ctx->dt;

      ierr = DMStagMatSetValuesStencil(dmTemp,A,1,&row,5,col,valA,INSERT_VALUES);CHKERRQ(ierr);
      valRhs = arr[ey][ex][slot]/ctx->dt; /* No heat source */

      ierr = DMStagVecSetValuesStencil(dmTemp,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = DMStagVecRestoreArrayRead(dmStokes,stokesLocal,&arrStokes);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dmStokes,&stokesLocal);CHKERRQ(ierr);
  ierr = DMStagVecRestoreArrayRead(dmTemp,tempPrevLocal,&arr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dmTemp,&tempPrevLocal);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rhs);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode UpdateTemperature(Ctx ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = StagBLSolverSolve(ctx->temperature_solver,ctx->temperature_array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
