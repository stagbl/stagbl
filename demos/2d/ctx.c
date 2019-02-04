#include "ctx.h"

StagBLErrorCode CtxCreate(MPI_Comm comm,Ctx *pctx)
{
  StagBLErrorCode ierr;
  Ctx             ctx;

  ierr = PetscMalloc1(1,pctx);CHKERRQ(ierr);
  ctx = *pctx;
  ctx->comm = comm;
  ctx->xmin = 0.0;
  ctx->xmax = 1e6;
  ctx->ymin = 0.0;
  ctx->ymax = 1.5e6;
  ctx->rho1 = 3200;
  ctx->rho2 = 3300;
  ctx->eta1 = 1e20;
  ctx->eta2 = 1e22;
  ctx->gy   = 10.0;
  return 0;
}

StagBLErrorCode CtxSetupFromGrid(Ctx ctx)
{
  StagBLErrorCode ierr;
  DM         dmStokes;
  StagBLInt  N[2];
  StagBLReal hxAvgInv;

  /* Escape hatch */
  ierr = StagBLGridPETScGetDM(ctx->stokesGrid,&dmStokes);CHKERRQ(ierr);

  ierr = DMStagGetGlobalSizes(dmStokes,&N[0],&N[1],NULL);CHKERRQ(ierr);
  ctx->hxCharacteristic = (ctx->xmax-ctx->xmin)/N[0];
  ctx->hyCharacteristic = (ctx->ymax-ctx->ymin)/N[1];
  ctx->etaCharacteristic = PetscMin(ctx->eta1,ctx->eta2);
  hxAvgInv = 2.0/(ctx->hxCharacteristic + ctx->hyCharacteristic);
  ctx->Kcont = ctx->etaCharacteristic*hxAvgInv;
  ctx->Kbound = ctx->etaCharacteristic*hxAvgInv*hxAvgInv;
  if (N[0] < 2) SETERRQ(ctx->comm,PETSC_ERR_SUP,"Not implemented for a single element in the x direction");
  ctx->pinx = 1; ctx->piny = 0;
  return 0;
}
