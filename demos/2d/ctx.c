#include "ctx.h"

PetscErrorCode CtxCreate(MPI_Comm comm,const char* mode,Ctx *pctx)
{
  PetscErrorCode ierr;
  Ctx             ctx;
  PetscBool       flg;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,pctx);CHKERRQ(ierr);
  ctx = *pctx;
  ierr = PetscStrcpy(ctx->mode,mode);CHKERRQ(ierr);
  ctx->comm = comm;
  flg = PETSC_FALSE;
  ierr = PetscStrcmp(mode,"gerya72",&flg);CHKERRQ(ierr);
  if (!flg) {
    ierr = PetscStrcmp(mode,"sinker",&flg);CHKERRQ(ierr);
  }
  if (flg) {
    ctx->xmin = 0.0;
    ctx->xmax = 1e6;
    ctx->ymin = 0.0;
    ctx->ymax = 1.5e6;
    ctx->rho1 = 3200;
    ctx->rho2 = 3300;
    ctx->eta1 = 1e20;
    ctx->eta2 = 1e22;
    ctx->gy   = 10.0;
  } else SETERRQ1(comm,PETSC_ERR_SUP,"Unrecognized mode %s",mode);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode CtxDestroy(Ctx *pctx)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscFree(*pctx);CHKERRQ(ierr);
  *pctx = NULL;
  PetscFunctionReturn(0);
}
