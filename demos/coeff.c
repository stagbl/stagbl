#include "coeff.h"

/* Coefficient/forcing Functions */

/* Constant */
static PetscScalar getRho_constant(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  STAGBL_UNUSED(x);
  STAGBL_UNUSED(y);
  STAGBL_UNUSED(z);
  return ctx->rho1;
}

static PetscScalar getEta_constant(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  STAGBL_UNUSED(x);
  STAGBL_UNUSED(y);
  STAGBL_UNUSED(z);
  return ctx->eta1;
}

/* Sinker */
static PetscScalar getRho_sinker(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  const PetscScalar zz = z/d - 0.5;
  return (xx*xx + yy*yy + zz*zz) > 0.3*0.3 ? ctx->rho1 : ctx->rho2;
}

static PetscScalar getEta_sinker(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  const PetscScalar zz = z/d - 0.5;
  return (xx*xx + yy*yy + zz*zz) > 0.3*0.3 ? ctx->eta1 : ctx->eta2;
}

/* Vertical layers */
static PetscScalar getRho_gerya72(void *ptr,PetscScalar x,PetscScalar y, PetscScalar z)
{
  Ctx ctx = (Ctx) ptr;
  STAGBL_UNUSED(z);
  if (x + 0.0*y < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->rho1;
  } else {
    return ctx->rho2;
  }
}

static PetscScalar getEta_gerya72(void *ptr,PetscScalar x,PetscScalar y, PetscScalar z)
{
  STAGBL_UNUSED(z);
  Ctx ctx = (Ctx) ptr;
  if (x  + 0.0*y < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->eta1;
  } else {
    return ctx->eta2;
  }
}

PetscErrorCode PopulateCoefficientData(Ctx ctx,const char* mode)
{
  PetscErrorCode ierr;
  PetscInt       dim;
  PetscInt       N[3];
  PetscInt       ex,ey,ez,startx,starty,startz,nx,ny,nz;
  PetscInt       iprev,icenter; // TODO change names
  PetscInt       slot_rho_downleft,slot_rho_backleft,slot_rho_backdown,slot_eta_element,slot_eta_downleft,slot_eta_backleft,slot_eta_backdown;
  DM             dmCoeff;
  Vec            *pcoeffLocal;
  Vec            coeffLocal;
  PetscReal      **cArrX,**cArrY,**cArrZ; // TODO names in this whole function are ugly
  PetscBool      flg;

  PetscFunctionBeginUser;
  flg = PETSC_FALSE;
  ierr = PetscStrcmp(mode,"gerya72",&flg);CHKERRQ(ierr);
  if (flg) {
    ctx->getEta = getEta_gerya72;
    ctx->getRho = getRho_gerya72;
  }
  if (!flg) {
    ierr = PetscStrcmp(mode,"sinker",&flg);CHKERRQ(ierr);
    if (flg) {
      ctx->getEta = getEta_sinker;
      ctx->getRho = getRho_sinker;
    }
  }
  if (!flg) {
    ierr = PetscStrcmp(mode,"blankenbach",&flg);CHKERRQ(ierr);
    if (flg) {
      ctx->getEta = getEta_constant;
      ctx->getRho = getRho_constant;
    }
  }
  if (!flg) {
    SETERRQ1(ctx->comm,PETSC_ERR_ARG_OUTOFRANGE,"Unrecognized mode %s",mode);
  }

  /* Pull out DM object */
  ierr = StagBLGridPETScGetDM(ctx->coefficient_grid,&dmCoeff);CHKERRQ(ierr);
  ierr = DMGetDimension(dmCoeff,&dim);CHKERRQ(ierr);

  /* If array doesnt exist, create it and pull out a local Vec. Otherwise, get the local Vec */
  if (!ctx->coefficient_array) {
    ierr = StagBLGridCreateStagBLArray(ctx->coefficient_grid,&ctx->coefficient_array);CHKERRQ(ierr);
    ierr = StagBLArrayPETScGetLocalVecPointer(ctx->coefficient_array,&pcoeffLocal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(dmCoeff,pcoeffLocal);CHKERRQ(ierr);
    coeffLocal = *pcoeffLocal;
  } else {
    ierr = StagBLArrayPETScGetLocalVec(ctx->coefficient_array,&coeffLocal);CHKERRQ(ierr);
  }

  ierr = DMStagGetGhostCorners(dmCoeff,&startx,&starty,&startz,&nx,&ny,&nz);CHKERRQ(ierr); /* Iterate over all local elements */
  ierr = DMStagGetGlobalSizes(dmCoeff,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateArraysRead(dmCoeff,&cArrX,&cArrY,&cArrZ);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dmCoeff,DMSTAG_ELEMENT,&icenter);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dmCoeff,DMSTAG_LEFT,&iprev);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_ELEMENT,  0,&slot_eta_element);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_DOWN_LEFT,0,&slot_eta_downleft);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_DOWN_LEFT,1,&slot_rho_downleft);CHKERRQ(ierr);
  if (dim == 2) {
    PetscScalar ***coeffArr;

    ierr = DMStagVecGetArray(dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);
    for (ey = starty; ey<starty+ny; ++ey) {
      for (ex = startx; ex<startx+nx; ++ex) {
        coeffArr[ey][ex][slot_eta_element]  = ctx->getEta(ctx,cArrX[ex][icenter],cArrY[ey][icenter],0.0);
        coeffArr[ey][ex][slot_eta_downleft] = ctx->getEta(ctx,cArrX[ex][iprev],  cArrY[ey][iprev],  0.0);
        coeffArr[ey][ex][slot_rho_downleft] = ctx->getRho(ctx,cArrX[ex][iprev],  cArrY[ey][iprev],  0.0);
      }
    }
    ierr = DMStagVecRestoreArray(dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);
  } else if (dim == 3) {
    PetscScalar ****coeffArr;

    ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_BACK_LEFT,0,&slot_eta_backleft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_BACK_LEFT,1,&slot_rho_backleft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_BACK_DOWN,0,&slot_eta_backdown);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_BACK_DOWN,1,&slot_rho_backdown);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);
    for (ez = startz; ez<startz+nz; ++ez) {
      for (ey = starty; ey<starty+ny; ++ey) {
        for (ex = startx; ex<startx+nx; ++ex) {
          coeffArr[ez][ey][ex][slot_eta_element]  = ctx->getEta(ctx,cArrX[ex][icenter],cArrY[ey][icenter],cArrZ[ez][icenter]);
          coeffArr[ez][ey][ex][slot_eta_downleft] = ctx->getEta(ctx,cArrX[ex][iprev],  cArrY[ey][iprev],  cArrZ[ez][icenter]);
          coeffArr[ez][ey][ex][slot_rho_downleft] = ctx->getRho(ctx,cArrX[ex][iprev],  cArrY[ey][iprev],  cArrZ[ez][icenter]);
          coeffArr[ez][ey][ex][slot_eta_backleft] = ctx->getEta(ctx,cArrX[ex][iprev],  cArrY[ey][icenter],cArrZ[ez][iprev]);
          coeffArr[ez][ey][ex][slot_rho_backleft] = ctx->getRho(ctx,cArrX[ex][iprev],  cArrY[ey][icenter],cArrZ[ez][iprev]);
          coeffArr[ez][ey][ex][slot_eta_backdown] = ctx->getEta(ctx,cArrX[ex][icenter],cArrY[ey][iprev],  cArrZ[ez][iprev]);
          coeffArr[ez][ey][ex][slot_rho_backdown] = ctx->getRho(ctx,cArrX[ex][icenter],cArrY[ey][iprev],  cArrZ[ez][iprev]);
        }
      }
    }
    ierr = DMStagVecRestoreArray(dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);
  } else SETERRQ1(PetscObjectComm((PetscObject)dmCoeff),PETSC_ERR_SUP,"Unsupported dimension %d",dim);
  ierr = DMStagRestoreProductCoordinateArraysRead(dmCoeff,&cArrX,&cArrY,&cArrZ);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
