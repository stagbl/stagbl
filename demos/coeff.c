#include "coeff.h"

/* Coefficient/forcing Functions */

/* Constant */
static PetscScalar getRho_constant(void *ptr,PetscScalar x, PetscScalar y) {
  STAGBL_UNUSED(x);
  STAGBL_UNUSED(y);
  Ctx ctx = (Ctx) ptr;
  return ctx->rho1;
}

static PetscScalar getEta_constant(void *ptr,PetscScalar x, PetscScalar y) {
  STAGBL_UNUSED(x);
  STAGBL_UNUSED(y);
  Ctx ctx = (Ctx) ptr;
  return ctx->eta1;
}

/* Sinker */
static PetscScalar getRho_sinker(void *ptr,PetscScalar x, PetscScalar y) {
  Ctx ctx = (Ctx) ptr;
  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  return (xx*xx + yy*yy) > 0.3*0.3 ? ctx->rho1 : ctx->rho2;
}

static PetscScalar getEta_sinker(void *ptr,PetscScalar x, PetscScalar y) {
  Ctx ctx = (Ctx) ptr;
  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  return (xx*xx + yy*yy) > 0.3*0.3 ? ctx->eta1 : ctx->eta2;
}

/* Vertical layers */
static PetscScalar getRho_gerya72(void *ptr,PetscScalar x,PetscScalar y)
{
  Ctx ctx = (Ctx) ptr;
  if (x + 0.0*y < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->rho1;
  } else {
    return ctx->rho2;
  }
}

static PetscScalar getEta_gerya72(void *ptr,PetscScalar x,PetscScalar y)
{
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
  PetscInt       N[2];
  PetscInt       ex,ey,startx,starty,nx,ny,ietaCorner,ietaElement,irho,iprev,icenter;
  DM             dmCoeff;
  Vec            *pcoeffLocal;
  Vec            coeffLocal;
  PetscReal      **cArrX,**cArrY;
  PetscReal      ***coeffArr;
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

  /* If array doesnt exist, create it an pull out a local Vec. Otherwise,
     zero the local Vec */
  if (!ctx->coefficient_array) {
    ierr = StagBLGridCreateStagBLArray(ctx->coefficient_grid,&ctx->coefficient_array);CHKERRQ(ierr);
    ierr = StagBLArrayPETScGetLocalVecPointer(ctx->coefficient_array,&pcoeffLocal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(dmCoeff,pcoeffLocal);CHKERRQ(ierr);
    coeffLocal = *pcoeffLocal;
  } else {
    ierr = StagBLArrayPETScGetLocalVec(ctx->coefficient_array,&coeffLocal);CHKERRQ(ierr);
  }

  ierr = DMStagGetGhostCorners(dmCoeff,&startx,&starty,NULL,&nx,&ny,NULL);CHKERRQ(ierr); /* Iterate over all local elements */
  ierr = DMStagGetGlobalSizes(dmCoeff,&N[0],&N[1],NULL);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_DOWN_LEFT,0,&ietaCorner);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_ELEMENT,0,&ietaElement);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_DOWN_LEFT,1,&irho);CHKERRQ(ierr);

  ierr = DMStagGetProductCoordinateArraysRead(dmCoeff,&cArrX,&cArrY,NULL);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dmCoeff,DMSTAG_ELEMENT,&icenter);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dmCoeff,DMSTAG_LEFT,&iprev);CHKERRQ(ierr);

  ierr = DMStagVecGetArray(dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);

  for (ey = starty; ey<starty+ny; ++ey) {
    for (ex = startx; ex<startx+nx; ++ex) {
      coeffArr[ey][ex][ietaElement] = ctx->getEta(ctx,cArrX[ex][icenter],cArrY[ey][icenter]);
      coeffArr[ey][ex][ietaCorner]  = ctx->getEta(ctx,cArrX[ex][iprev],cArrY[ey][iprev]);
      coeffArr[ey][ex][irho]        = ctx->getRho(ctx,cArrX[ex][iprev],cArrY[ey][iprev]);
    }
  }
  ierr = DMStagRestoreProductCoordinateArraysRead(dmCoeff,&cArrX,&cArrY,NULL);CHKERRQ(ierr);
  ierr = DMStagVecRestoreArray(dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

