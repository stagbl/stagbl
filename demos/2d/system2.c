/*
Implement an alternate system of equations to solve.
This is intended to exactly match exercise 7.2 in Gerya's textbook,
"Introduction to Numerical Geodynamical Modelling" (1st ed., 2009),
and attendant MATLAB implementation. Note that the MATLAB implementation
in the solutions for the textbook appears to have an error in the interpolation
of velocities to elements, so only the "raw" velocities and pressures will
match that output.
*/
#include "system2.h"

PetscErrorCode CreateSystem2(const Ctx ctx,Mat *pA,Vec *pRhs)
{
  PetscErrorCode  ierr;
  DM              dmStokes,dmCoeff;
  PetscInt        ex,ey,startx,starty,nx,ny;
  Mat             A;
  Vec             rhs;
  PetscReal       hx,hy;
  PetscInt        N[2];
  const PetscBool pinPressure = PETSC_TRUE;
  Vec             coeffLocal;

  PetscFunctionBeginUser;

  /* Use the "escape hatch" */
  ierr = StagBLGridPETScGetDM(ctx->stokesGrid,&dmStokes);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(ctx->coeffGrid,&dmCoeff);CHKERRQ(ierr);

  ierr = DMCreateMatrix(dmStokes,pA);CHKERRQ(ierr);
  A = *pA;
  ierr = DMCreateGlobalVector(dmStokes,pRhs);CHKERRQ(ierr);
  rhs = *pRhs;
  ierr = DMStagGetCorners(dmStokes,&startx,&starty,NULL,&nx,&ny,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(dmStokes,&N[0],&N[1],NULL);CHKERRQ(ierr);
  hx = ctx->hxCharacteristic;
  hy = ctx->hyCharacteristic;
  ierr = DMGetLocalVector(dmCoeff,&coeffLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocal(dmCoeff,ctx->coeff,INSERT_VALUES,coeffLocal);CHKERRQ(ierr);

  /* Loop over all elements at once */
  for (ey = starty; ey<starty+ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
    for (ex = startx; ex<startx+nx; ++ex) {

      /* Vy terms */
      if (ey == N[1]-1) {
        /* Top boundary velocity Dirichlet */
        DMStagStencil     row;
        PetscScalar       valRhs;
        const PetscScalar valA = ctx->Kbound;
        row.i = ex; row.j = ey; row.loc = DMSTAG_UP; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      if (ey == 0) {
        /* Bottom boundary velocity Dirichlet */
        DMStagStencil     row;
        PetscScalar       valRhs;
        const PetscScalar valA = ctx->Kbound;
        row.i = ex; row.j = ey; row.loc = DMSTAG_DOWN; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        PetscInt      nEntries;
        DMStagStencil row,col[11];
        PetscScalar   valRhs,valA[11];
        if (ex == 0) {
          /* Left boundary y velocity - zero gradient */
          nEntries = 2;
          row.i = ex; row.j = ey; row.loc = DMSTAG_DOWN; row.c = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DMSTAG_DOWN; col[0].c  = 0; valA[0] =  ctx->Kbound;
          col[1].i = ex+1; col[1].j  = ey  ; col[1].loc  = DMSTAG_DOWN; col[1].c  = 0; valA[1] = -ctx->Kbound;
          valRhs  = 0.0;
        } else if (ex == N[0]-1) {
          /* Right boundary y velocity - zero gradient */
          nEntries = 2;
          row.i = ex; row.j = ey; row.loc = DMSTAG_DOWN; row.c = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DMSTAG_DOWN; col[0].c  = 0; valA[0] =  ctx->Kbound;
          col[1].i = ex-1; col[1].j  = ey  ; col[1].loc  = DMSTAG_DOWN; col[1].c  = 0; valA[1] = -ctx->Kbound;
          valRhs  = 0.0;
        } else  {
          /* y velocity - interior */
          DMStagStencil rhoPoint[2];
          PetscScalar   rho[2];
          DMStagStencil etaPoint[4];
          PetscScalar   eta[4],etaLeft,etaRight,etaUp,etaDown;

          /* get rho values  and compute rhs value*/
          rhoPoint[0].i = ex; rhoPoint[0].j = ey; rhoPoint[0].loc = DMSTAG_DOWN_LEFT;  rhoPoint[0].c = 1;
          rhoPoint[1].i = ex; rhoPoint[1].j = ey; rhoPoint[1].loc = DMSTAG_DOWN_RIGHT; rhoPoint[1].c = 1;
          ierr = DMStagVecGetValuesStencil(dmCoeff,coeffLocal,2,rhoPoint,rho);CHKERRQ(ierr);
          valRhs = -ctx->gy * 0.5 * (rho[0] + rho[1]);

          /* Get eta values */
          etaPoint[0].i = ex; etaPoint[0].j = ey;   etaPoint[0].loc = DMSTAG_DOWN_LEFT;  etaPoint[0].c = 0; /* Left  */
          etaPoint[1].i = ex; etaPoint[1].j = ey;   etaPoint[1].loc = DMSTAG_DOWN_RIGHT; etaPoint[1].c = 0; /* Right */
          etaPoint[2].i = ex; etaPoint[2].j = ey+1; etaPoint[2].loc = DMSTAG_ELEMENT;    etaPoint[2].c = 0; /* Up    */
          etaPoint[3].i = ex; etaPoint[3].j = ey-1; etaPoint[3].loc = DMSTAG_ELEMENT;    etaPoint[3].c = 0; /* Down  */
          ierr = DMStagVecGetValuesStencil(dmCoeff,coeffLocal,4,etaPoint,eta);CHKERRQ(ierr);
          etaLeft = eta[0]; etaRight = eta[1]; etaUp = eta[2]; etaDown = eta[3];

          nEntries = 11;

          row.i    = ex  ; row.j     = ey  ; row.loc     = DMSTAG_DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DMSTAG_DOWN;     col[0].c  = 0; valA[0]  = -2.0 * (etaDown + etaUp) / (hy*hy) - (etaLeft + etaRight) /(hx*hx);
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DMSTAG_DOWN;     col[1].c  = 0; valA[1]  =  2.0 * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DMSTAG_DOWN;     col[2].c  = 0; valA[2]  =  2.0 * etaUp    / (hy*hy);
          col[3].i = ex-1; col[3].j  = ey  ; col[3].loc  = DMSTAG_DOWN;     col[3].c  = 0; valA[3]  =        etaLeft  / (hx*hx);
          col[4].i = ex+1; col[4].j  = ey  ; col[4].loc  = DMSTAG_DOWN;     col[4].c  = 0; valA[4]  =        etaRight / (hx*hx);
          col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = DMSTAG_LEFT;     col[5].c  = 0; valA[5]  =        etaLeft  / (hx*hy); /* down left x edge */
          col[6].i = ex  ; col[6].j  = ey-1; col[6].loc  = DMSTAG_RIGHT;    col[6].c  = 0; valA[6]  = -      etaRight / (hx*hy); /* down right x edge */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = DMSTAG_LEFT;     col[7].c  = 0; valA[7]  = -      etaLeft  / (hx*hy); /* up left x edge */
          col[8].i = ex  ; col[8].j  = ey  ; col[8].loc  = DMSTAG_RIGHT;    col[8].c  = 0; valA[8]  =        etaRight / (hx*hy); /* up right x edge */
          col[9].i = ex  ; col[9].j  = ey-1; col[9].loc  = DMSTAG_ELEMENT;  col[9].c  = 0; valA[9]  =  ctx->Kcont / hy;
          col[10].i = ex ; col[10].j = ey  ; col[10].loc = DMSTAG_ELEMENT; col[10].c  = 0; valA[10] = -ctx->Kcont / hy;

        }
        ierr = DMStagMatSetValuesStencil(dmStokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      if (ex == N[0]-1) {
        /* Right Boundary velocity Dirichlet */
        /* Redundant in the corner */
        DMStagStencil row;
        PetscScalar   valRhs;

        const PetscScalar valA = ctx->Kbound;
        row.i = ex; row.j = ey; row.loc = DMSTAG_RIGHT; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }
      if (ex == 0) {
        /* Left velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = ctx->Kbound;
        row.i = ex; row.j = ey; row.loc = DMSTAG_LEFT; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        PetscInt nEntries;
        DMStagStencil row,col[11];
        PetscScalar   valRhs,valA[11];

        if (ey == 0) {
          /* Bottom boundary x velocity stencil (with zero vel deriv) */
          nEntries = 2;
          row.i = ex; row.j = ey; row.loc = DMSTAG_LEFT; row.c = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DMSTAG_LEFT; col[0].c  = 0; valA[0] =  ctx->Kbound;
          col[1].i = ex  ; col[1].j  = ey+1; col[1].loc  = DMSTAG_LEFT; col[1].c  = 0; valA[1] = -ctx->Kbound;
          valRhs = 0.0;
        } else if (ey == N[1]-1) {
          /* Top boundary x velocity stencil (with zero vel deriv) */
          nEntries = 2;
          row.i = ex; row.j = ey; row.loc = DMSTAG_LEFT; row.c = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DMSTAG_LEFT; col[0].c  = 0; valA[0] =  ctx->Kbound;
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DMSTAG_LEFT; col[1].c  = 0; valA[1] = -ctx->Kbound;
          valRhs = 0.0;
        } else {
          /* U_x interior equation */
          DMStagStencil etaPoint[4];
          PetscScalar eta[4],etaLeft,etaRight,etaUp,etaDown;


          /* Get eta values */
          etaPoint[0].i = ex-1; etaPoint[0].j = ey; etaPoint[0].loc = DMSTAG_ELEMENT;   etaPoint[0].c = 0; /* Left  */
          etaPoint[1].i = ex;   etaPoint[1].j = ey; etaPoint[1].loc = DMSTAG_ELEMENT;   etaPoint[1].c = 0; /* Right */
          etaPoint[2].i = ex;   etaPoint[2].j = ey; etaPoint[2].loc = DMSTAG_UP_LEFT;   etaPoint[2].c = 0; /* Up    */
          etaPoint[3].i = ex;   etaPoint[3].j = ey; etaPoint[3].loc = DMSTAG_DOWN_LEFT; etaPoint[3].c = 0; /* Down  */
          ierr = DMStagVecGetValuesStencil(dmCoeff,coeffLocal,4,etaPoint,eta);CHKERRQ(ierr);
          etaLeft = eta[0]; etaRight = eta[1]; etaUp = eta[2]; etaDown = eta[3];

          nEntries = 11;
          row.i     = ex  ; row.j     = ey  ; row.loc     = DMSTAG_LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = DMSTAG_LEFT;    col[0].c   = 0; valA[0]  = -2.0 * (etaLeft + etaRight) / (hx*hx) -(etaUp + etaDown) / (hy*hy);
          col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc  = DMSTAG_LEFT;    col[1].c   = 0; valA[1]  =        etaDown  / (hy*hy);
          col[2].i  = ex  ; col[2].j  = ey+1; col[2].loc  = DMSTAG_LEFT;    col[2].c   = 0; valA[2]  =        etaUp    / (hy*hy);
          col[3].i  = ex-1; col[3].j  = ey  ; col[3].loc  = DMSTAG_LEFT;    col[3].c   = 0; valA[3]  =  2.0 * etaLeft  / (hx*hx);
          col[4].i  = ex+1; col[4].j  = ey  ; col[4].loc  = DMSTAG_LEFT;    col[4].c   = 0; valA[4]  =  2.0 * etaRight / (hx*hx);
          col[5].i  = ex-1; col[5].j  = ey  ; col[5].loc  = DMSTAG_DOWN;    col[5].c   = 0; valA[5]  =        etaDown  / (hx*hy); /* down left */
          col[6].i  = ex  ; col[6].j  = ey  ; col[6].loc  = DMSTAG_DOWN;    col[6].c   = 0; valA[6]  = -      etaDown  / (hx*hy); /* down right */
          col[7].i  = ex-1; col[7].j  = ey  ; col[7].loc  = DMSTAG_UP;      col[7].c   = 0; valA[7]  = -      etaUp    / (hx*hy); /* up left */
          col[8].i  = ex  ; col[8].j  = ey  ; col[8].loc  = DMSTAG_UP;      col[8].c   = 0; valA[8]  =        etaUp    / (hx*hy); /* up right */
          col[9].i  = ex-1; col[9].j  = ey  ; col[9].loc  = DMSTAG_ELEMENT; col[9].c   = 0; valA[9]  =  ctx->Kcont / hx;
          col[10].i = ex  ; col[10].j = ey  ; col[10].loc = DMSTAG_ELEMENT; col[10].c  = 0; valA[10] = -ctx->Kcont / hx;
          valRhs = 0.0;
        }
        ierr = DMStagMatSetValuesStencil(dmStokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      /* P terms */
      if (pinPressure && ex == ctx->pinx && ey == ctx->piny) { /* Pin a pressure node to zero, if requested */
        DMStagStencil row;
        PetscScalar valA,valRhs;
        row.i = ex; row.j = ey; row.loc = DMSTAG_ELEMENT; row.c = 0;
        valA = ctx->Kbound;
        ierr = DMStagMatSetValuesStencil(dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        PetscInt nEntries;
        DMStagStencil row,col[5];
        PetscScalar   valA[5],valRhs;

        row.i    = ex; row.j    = ey; row.loc    = DMSTAG_ELEMENT; row.c    = 0;

        if (ex == 0 && (ey == 0 || ey == N[1]-1)) {
          nEntries = 2;
          col[0].i = ex  ; col[0].j = ey; col[0].loc = DMSTAG_ELEMENT; col[0].c = 0; valA[0] =  ctx->Kbound;
          col[1].i = ex+1; col[1].j = ey; col[1].loc = DMSTAG_ELEMENT; col[1].c = 0; valA[1] = -ctx->Kbound;
        } else if (ex == N[0]-1 && (ey == 0 || ey == N[1]-1)) {
          nEntries = 2;
          col[0].i = ex  ; col[0].j = ey; col[0].loc = DMSTAG_ELEMENT; col[0].c = 0; valA[0] =  ctx->Kbound;
          col[1].i = ex-1; col[1].j = ey; col[1].loc = DMSTAG_ELEMENT; col[1].c = 0; valA[1] = -ctx->Kbound;
        } else {
          nEntries = 5;
          col[0].i = ex; col[0].j = ey; col[0].loc = DMSTAG_LEFT;    col[0].c = 0; valA[0] = -ctx->Kcont / hx;
          col[1].i = ex; col[1].j = ey; col[1].loc = DMSTAG_RIGHT;   col[1].c = 0; valA[1] =  ctx->Kcont / hx;
          col[2].i = ex; col[2].j = ey; col[2].loc = DMSTAG_DOWN;    col[2].c = 0; valA[2] = -ctx->Kcont / hy;
          col[3].i = ex; col[3].j = ey; col[3].loc = DMSTAG_UP;      col[3].c = 0; valA[3] =  ctx->Kcont / hy;
          col[4] = row;                                                            valA[4] = 0.0; // Explicit zero for direct solve only
        }
        valRhs = 0.0;
        ierr = DMStagMatSetValuesStencil(dmStokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = DMRestoreLocalVector(dmCoeff,&coeffLocal);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rhs);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
