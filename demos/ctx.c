#include "ctx.h"

/**
  * Create application context and set parameters based on mode
  */
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

  /* Default Settings */
  ctx->totalTimesteps     = 3;
  ctx->dt                 = 1e13;
  ctx->temperature_top    = 0;
  ctx->temperature_bottom = 1000;
  ctx->kappa              = 1e-6;
  ctx->alpha              = 2.5e-5;
  ctx->cp                 = 1.25e3;
  ctx->KTemp              = 1.0; /* Note: not chosen in a princpled way, yet */

  ctx->compute_nusselt_number = PETSC_FALSE;

  /* Initialize Data */
  ctx->stokes_array       = NULL;
  ctx->coefficient_array  = NULL;
  ctx->temperature_array  = NULL;
  ctx->stokes_system      = NULL;
  ctx->stokes_solver      = NULL;
  ctx->temperature_system = NULL;
  ctx->temperature_solver = NULL;

  /* Additional, mode-dependent settings */
  flg = PETSC_FALSE;
  ierr = PetscStrcmp(mode,"gerya72",&flg);CHKERRQ(ierr);
  if (!flg) {
    ierr = PetscStrcmp(mode,"sinker",&flg);CHKERRQ(ierr);
  }
  if (!flg) {
    ierr = PetscStrcmp(mode,"sinker_box",&flg);CHKERRQ(ierr);
  }
  if (flg) {
    ctx->uniform_grid = PETSC_TRUE;
    ctx->boussinesq_forcing = PETSC_FALSE;
    ctx->xmin = 0.0;
    ctx->xmax = 1e6;
    ctx->ymin = 0.0;
    ctx->ymax = 1.5e6;
    ctx->zmin = 0.0;
    ctx->zmax = 1e6;
    ctx->rho1 = 3200;
    ctx->rho2 = 3300;
    ctx->eta1 = 1e20;
    ctx->eta2 = 1e22;
    ctx->eta_characteristic = ctx->eta1;
    ctx->gy   = 10.0; /* gravity points in positive y direction */
  } else {
    ierr = PetscStrcmp(mode,"blankenbach",&flg);CHKERRQ(ierr);
    if (flg) {
      ctx->uniform_grid = PETSC_TRUE;
      ctx->boussinesq_forcing = PETSC_TRUE;
      ctx->compute_nusselt_number = PETSC_TRUE;
      ctx->xmin = 0.0;
      ctx->xmax = 1e6;
      ctx->ymin = 0.0;
      ctx->ymax = 1e6;
      ctx->zmin = 0.0;
      ctx->zmax = 1e6;
      ctx->rho1 = 4000;
      ctx->rho2 = 4400; /* not used normally */
      ctx->eta1 = 2.5e19; /* case 1a */
      ctx->eta2 = 2.5e19; /* Not normally used */
      ctx->eta_characteristic = ctx->eta1;
      ctx->gy   = -10.0; /* gravity points in negative y direction */
    }
  }

  if (!flg) {
    SETERRQ1(comm,PETSC_ERR_SUP,"Unrecognized mode %s",mode);CHKERRQ(ierr);
  }

  /* Override settings from the command line */
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&ctx->dt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-nsteps",&ctx->totalTimesteps,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-kappa",&ctx->kappa,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-alpha",&ctx->alpha,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-eta1",&ctx->eta1,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-eta2",&ctx->eta2,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-eta_characteristic",&ctx->eta_characteristic,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-rho1",&ctx->rho1,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-rho2",&ctx->rho2,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-xmin",&ctx->xmin,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-xmax",&ctx->xmax,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-ymin",&ctx->ymin,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-ymax",&ctx->ymax,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-zmin",&ctx->zmin,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-zmax",&ctx->zmax,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-gy",&ctx->gy,NULL);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode CtxDestroy(Ctx *pctx)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = StagBLArrayDestroy(&(*pctx)->stokes_array);CHKERRQ(ierr);
  ierr = StagBLSystemDestroy(&(*pctx)->stokes_system);CHKERRQ(ierr);
  ierr = StagBLSolverDestroy(&(*pctx)->stokes_solver);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&(*pctx)->stokes_grid);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&(*pctx)->coefficient_grid);CHKERRQ(ierr);
  ierr = StagBLArrayDestroy(&(*pctx)->coefficient_array);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&(*pctx)->temperature_grid);CHKERRQ(ierr);
  ierr = StagBLArrayDestroy(&(*pctx)->temperature_array);CHKERRQ(ierr);
  ierr = StagBLSystemDestroy(&(*pctx)->temperature_system);CHKERRQ(ierr);
  ierr = StagBLSolverDestroy(&(*pctx)->temperature_solver);CHKERRQ(ierr);
  ierr = DMDestroy(&(*pctx)->dm_particles);CHKERRQ(ierr);
  ierr = PetscFree(*pctx);CHKERRQ(ierr);
  *pctx = NULL;
  PetscFunctionReturn(0);
}
