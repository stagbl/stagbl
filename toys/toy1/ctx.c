#include "ctx.h" 

/* =========================================================================== */
PetscErrorCode CreateCtx(Ctx *pctx)
{
  PetscErrorCode  ierr;
  Ctx             ctx;

  PetscFunctionBeginUser;
  *pctx = malloc(sizeof(CtxData));
  ctx = *pctx;
  ctx->isoviscous = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-isoviscous",&ctx->isoviscous,NULL);CHKERRQ(ierr);
  ctx->rho1 = 3200;
  ctx->rho2 = 3300;
  ctx->eta1 = 1e20;
  ctx->eta2 = 1e22;
  if (ctx->isoviscous) {
    ctx->etaCharacteristic = 1e21;
  } else {
    ctx->etaCharacteristic = PetscMin(ctx->eta1,ctx->eta2);
  }
  ctx->gx = 0.0;
  ctx->gy = 10.0;
  PetscFunctionReturn(0);
}

/* =========================================================================== */
PetscErrorCode DestroyCtx(Ctx *pctx)
{
  PetscErrorCode ierr;
  Ctx            ctx=*pctx;

  PetscFunctionBeginUser;
  if (ctx->rho)      {ierr = VecDestroy(&ctx->rho);CHKERRQ(ierr);}
  if (ctx->etaC)     {ierr = VecDestroy(&ctx->etaC);CHKERRQ(ierr);}
  if (ctx->etaCor)   {ierr = VecDestroy(&ctx->etaCor);CHKERRQ(ierr);}
  if (ctx->daC)      {ierr = DMDestroy(&ctx->daC);CHKERRQ(ierr);}
  if (ctx->daCor)    {ierr = DMDestroy(&ctx->daCor);CHKERRQ(ierr);}
  if (ctx->daX)      {ierr = DMDestroy(&ctx->daX);CHKERRQ(ierr);}
  if (ctx->daY)      {ierr = DMDestroy(&ctx->daY);CHKERRQ(ierr);}
  if (ctx->dmStokes) {ierr = DMDestroy(&ctx->dmStokes);CHKERRQ(ierr);}
  free(ctx);
  *pctx = NULL;
  PetscFunctionReturn(0);
}
