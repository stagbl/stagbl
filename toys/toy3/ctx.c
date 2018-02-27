#include "ctx.h" 

PetscErrorCode CreateCtx(Ctx *pctx)
{
  PetscErrorCode  ierr;
  Ctx             ctx;
  PetscReal       hxAvgInv;

  PetscFunctionBeginUser;
  *pctx = malloc(sizeof(CtxData));
  ctx = *pctx;
  ctx->stokesGrid = NULL;
  ctx->paramGrid = NULL;
  ctx->isoviscous = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-isoviscous",&ctx->isoviscous,NULL);CHKERRQ(ierr);
  ctx->M = 30; ctx->N= 20;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&ctx->M,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&ctx->N,NULL);CHKERRQ(ierr);
  ctx->xmin = 0.0;
  ctx->xmax = 1e6;
  ierr = PetscOptionsGetReal(NULL,NULL,"-xmax",&ctx->xmax,NULL);CHKERRQ(ierr);
  ctx->ymin = 0.0;
  ctx->ymax = 1.5e6;
  ierr = PetscOptionsGetReal(NULL,NULL,"-ymax",&ctx->ymax,NULL);CHKERRQ(ierr);
  ctx->rho1 = 3200;
  ctx->rho2 = 3300;
  ctx->eta1 = 1e20;
  ctx->eta2 = 1e22;
  if (ctx->isoviscous) {
    ctx->etaCharacteristic = 1e21;
  } else {
    ctx->etaCharacteristic = PetscMin(ctx->eta1,ctx->eta2);
  }
  ctx->hxCharacteristic = (ctx->xmax - ctx->xmin)/ctx->M; // assumes uniform
  ctx->hyCharacteristic = (ctx->ymax - ctx->ymin)/ctx->N; // assumes uniform (otherwise might be better as a min/max value)
  ctx->gx = 0.0;
  ctx->gy = 10.0;
  if (ctx->M > 2 && ctx->N > 3) {
    ctx->pFixx = 1; ctx->pFixy = 2; // To correspond with Chapter 7 solutions 
  } else {
    ctx->pFixx = 0; ctx->pFixy = 1;
  }
  hxAvgInv = 2.0/(ctx->hxCharacteristic + ctx->hyCharacteristic);
  ctx->Kcont         = ctx->etaCharacteristic*hxAvgInv;
  ctx->Kbound        = ctx->etaCharacteristic*hxAvgInv*hxAvgInv;
  PetscFunctionReturn(0);
}

PetscErrorCode DestroyCtx(Ctx *pctx)
{
  PetscErrorCode ierr;
  Ctx            ctx=*pctx;

  PetscFunctionBeginUser;
  if (ctx->stokesGrid) {ierr = DMDestroy (&ctx->stokesGrid);CHKERRQ(ierr);}
  if (ctx->paramGrid)  {ierr = DMDestroy (&ctx->paramGrid); CHKERRQ(ierr);}
  if (ctx->param)      {ierr = VecDestroy(&ctx->param);     CHKERRQ(ierr);}
  free(ctx);
  *pctx = NULL;
  PetscFunctionReturn(0);
}
