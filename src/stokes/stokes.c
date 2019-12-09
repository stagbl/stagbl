#include "stagbl.h"

PetscErrorCode StagBLStokesParametersCreate(StagBLStokesParameters *parameters)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,parameters);CHKERRQ(ierr); /* Zero all fields */
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLStokesParametersDestroy(StagBLStokesParameters *parameters)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(*parameters);CHKERRQ(ierr);
  *parameters = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLGridCreateStokes2DBox(MPI_Comm comm, PetscInt nx,PetscInt ny,PetscScalar xmin, PetscScalar xmax, PetscScalar ymin, PetscScalar ymax,StagBLGrid *pgrid)
{
  PetscErrorCode ierr;
  DM *pdm;
  DM dm_stokes;

  PetscFunctionBegin;
  StagBLGridCreate(pgrid);
  StagBLGridPETScGetDMPointer(*pgrid,&pdm);
  ierr = DMStagCreate2d(
      comm,
      DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
      nx,ny,                                   /* Global element counts */
      PETSC_DECIDE,PETSC_DECIDE,               /* Determine parallel decomposition automatically */
      0,1,1,                                   /* dof: 0 per vertex, 1 per edge, 1 per face/element */
      DMSTAG_STENCIL_BOX,
      1,                                       /* elementwise stencil width */
      NULL,NULL,
      pdm);
  dm_stokes = *pdm;
  ierr = DMSetFromOptions(dm_stokes);CHKERRQ(ierr);
  ierr = DMSetUp(dm_stokes);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesProduct(dm_stokes,xmin,xmax,ymin,ymax,0.0,0.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLGridCreateStokes3DBox(MPI_Comm comm, PetscInt nx,PetscInt ny,PetscInt nz,PetscScalar xmin, PetscScalar xmax, PetscScalar ymin, PetscScalar ymax, PetscScalar zmin, PetscScalar zmax, StagBLGrid *pgrid)
{
  PetscErrorCode ierr;
  DM *pdm;
  DM dm_stokes;

  PetscFunctionBegin;
  StagBLGridCreate(pgrid);
  StagBLGridPETScGetDMPointer(*pgrid,&pdm);
  ierr = DMStagCreate3d(
      comm,
      DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
      nx,ny,nz,                                /* Global element counts */
      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,  /* Determine parallel decomposition automatically */
      0,0,1,1,                                 /* dof: 1 per face, 1 per element */
      DMSTAG_STENCIL_BOX,
      1,                                       /* elementwise stencil width */
      NULL,NULL,NULL,
      pdm);
  dm_stokes = *pdm;
  ierr = DMSetFromOptions(dm_stokes);CHKERRQ(ierr);
  ierr = DMSetUp(dm_stokes);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesProduct(dm_stokes,xmin,xmax,ymin,ymax,zmin,zmax);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateSystem_2D_FreeSlip(StagBLStokesParameters,StagBLSystem);
static PetscErrorCode CreateSystem_3D_FreeSlip(StagBLStokesParameters,StagBLSystem);

/**
  * A general function which creates Stokes StagBLSystem objects. It accepts
  * a struct containing all relevant parameters.
  */
PetscErrorCode StagBLCreateStokesSystem(StagBLStokesParameters parameters, StagBLSystem *system)
{
  PetscErrorCode ierr;
  DM             dm_stokes;
  PetscInt       dim;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(parameters->stokes_grid,&dm_stokes);CHKERRQ(ierr);
  ierr = DMGetDimension(dm_stokes,&dim);CHKERRQ(ierr);
  ierr = StagBLGridCreateStagBLSystem(parameters->stokes_grid,system);CHKERRQ(ierr);
  switch (dim) {
    case 2:
    ierr = CreateSystem_2D_FreeSlip(parameters,*system);CHKERRQ(ierr);
  break;
    case 3:
    ierr = CreateSystem_3D_FreeSlip(parameters,*system);CHKERRQ(ierr);
    break;
    default: StagBLError1(PetscObjectComm((PetscObject)dm_stokes),"Unsupported dimension %D",dim);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateSystem_2D_FreeSlip(StagBLStokesParameters parameters,StagBLSystem system)
{
  PetscErrorCode  ierr;
  DM              dm_stokes,dm_coefficient;
  PetscInt        N[2];
  PetscInt        ex,ey,startx,starty,nx,ny;
  Mat             *pA;
  Vec             *pRhs;
  Mat             A;
  Vec             rhs;
  PetscReal       hx,hy,dv,hxAvgInv,Kcont,Kbound;
  PetscInt        pinx,piny;
  const PetscBool pin_pressure = PETSC_TRUE;
  Vec             coeff_local;
  StagBLGrid      coefficient_grid;
  DM              dm_temperature;
  Vec             temperature,temperature_local;
  PetscScalar     ***arr_temperature;
  PetscInt        slot_temperature_downleft,slot_temperature_downright;

  PetscFunctionBeginUser;
  ierr = StagBLSystemPETScGetMatPointer(system,&pA);CHKERRQ(ierr);
  ierr = StagBLSystemPETScGetVecPointer(system,&pRhs);CHKERRQ(ierr);
  if (!parameters->coefficient_array) StagBLError(PETSC_COMM_SELF,"coefficient_array field not set in StagBLStokesParameters argument");
  ierr = StagBLArrayGetStagBLGrid(parameters->coefficient_array,&coefficient_grid);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(parameters->stokes_grid,&dm_stokes);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(coefficient_grid,&dm_coefficient);CHKERRQ(ierr);
  ierr = StagBLArrayPETScGetLocalVec(parameters->coefficient_array,&coeff_local);CHKERRQ(ierr);

  /* Compute some parameters */
  ierr = DMStagGetGlobalSizes(dm_stokes,&N[0],&N[1],NULL);CHKERRQ(ierr);
  if (N[0] < 2 || N[1] < 2) StagBLError(PetscObjectComm((PetscObject)dm_stokes),"Stokes system construction not implemented for a single element in any direction");
  if (parameters->uniform_grid) {
    hx = (parameters->xmax-parameters->xmin)/N[0];
    hy = (parameters->ymax-parameters->ymin)/N[1];
    dv = hx*hy;
  } else StagBLError(PetscObjectComm((PetscObject)dm_stokes),"Non-uniform grids not supported yet");
  hxAvgInv = 2.0/(hx + hy);
  Kcont = parameters->eta_characteristic*hxAvgInv;
  Kbound = parameters->eta_characteristic*hxAvgInv*hxAvgInv;
  pinx = 1; piny = 0; /* Depends on the assertion above that we have 2 or more elemnets in the x direction */

  ierr = DMCreateMatrix(dm_stokes,pA);CHKERRQ(ierr);
  A = *pA;
  ierr = DMCreateGlobalVector(dm_stokes,pRhs);CHKERRQ(ierr);
  rhs = *pRhs;
  ierr = DMStagGetCorners(dm_stokes,&startx,&starty,NULL,&nx,&ny,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  /* If using Boussinesq forcing from a temperature field, get access */
  if (parameters->boussinesq_forcing) {
    ierr = StagBLGridPETScGetDM(parameters->temperature_grid,&dm_temperature);CHKERRQ(ierr);
    ierr = StagBLArrayPETScGetGlobalVec(parameters->temperature_array,&temperature);CHKERRQ(ierr);
    ierr = DMGetLocalVector(dm_temperature,&temperature_local);CHKERRQ(ierr);
    ierr = DMGlobalToLocal(dm_temperature,temperature,INSERT_VALUES,temperature_local);CHKERRQ(ierr);
    ierr = DMStagVecGetArrayRead(dm_temperature,temperature_local,&arr_temperature);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_DOWN_LEFT,  0,&slot_temperature_downleft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_DOWN_RIGHT, 0,&slot_temperature_downright);CHKERRQ(ierr);
  }

  /* Loop over all local elements. */
  for (ey = starty; ey<starty+ny; ++ey) {
    for (ex = startx; ex<startx+nx; ++ex) {

      if (ey == N[1]-1) {
        /* Top boundary velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   val_rhs;
        const PetscScalar val_A = Kbound;
        row.i = ex; row.j = ey; row.loc = DMSTAG_UP; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
        val_rhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      if (ey == 0) {
        /* Bottom boundary velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   val_rhs;
        const PetscScalar val_A = Kbound;
        row.i = ex; row.j = ey; row.loc = DMSTAG_DOWN; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
        val_rhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        /* Y-momentum equation : (u_xx + u_yy) - p_y = f^y : includes non-zero forcing */
        PetscInt      nEntries;
        DMStagStencil row,col[11];
        PetscScalar   val_A[11];
        DMStagStencil rhoPoint[2];
        PetscScalar   rho[2],val_rhs;
        DMStagStencil eta_point[4];
        PetscScalar   eta[4],eta_left,eta_right,eta_up,eta_down;

        /* get rho values  and compute rhs value*/
        rhoPoint[0].i = ex; rhoPoint[0].j = ey; rhoPoint[0].loc = DMSTAG_DOWN_LEFT;  rhoPoint[0].c = 1;
        rhoPoint[1].i = ex; rhoPoint[1].j = ey; rhoPoint[1].loc = DMSTAG_DOWN_RIGHT; rhoPoint[1].c = 1;
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,2,rhoPoint,rho);CHKERRQ(ierr);
        val_rhs = -parameters->gy * dv * 0.5 * (rho[0] + rho[1]);

        /* get rho values  */
        rhoPoint[0].i = ex; rhoPoint[0].j = ey; rhoPoint[0].loc = DMSTAG_DOWN_LEFT;  rhoPoint[0].c = 1;
        rhoPoint[1].i = ex; rhoPoint[1].j = ey; rhoPoint[1].loc = DMSTAG_DOWN_RIGHT; rhoPoint[1].c = 1;
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,2,rhoPoint,rho);CHKERRQ(ierr);

        /* Compute forcing */
        {
          const PetscReal rho_avg = 0.5 * (rho[0] + rho[1]);

          /* Note Boussinesq forcing has the opposite sign */
          if (parameters->boussinesq_forcing) {
            const PetscScalar temperature_average = 0.5 * (arr_temperature[ey][ex][slot_temperature_downleft] + arr_temperature[ey][ex][slot_temperature_downright]);
            val_rhs = parameters->alpha * rho_avg * parameters->gy * dv * temperature_average;
          } else {
            val_rhs = -parameters->gy * dv * rho_avg;
          }
        }

        /* Get eta values */
        eta_point[0].i = ex; eta_point[0].j = ey;   eta_point[0].loc = DMSTAG_DOWN_LEFT;  eta_point[0].c = 0; /* Left  */
        eta_point[1].i = ex; eta_point[1].j = ey;   eta_point[1].loc = DMSTAG_DOWN_RIGHT; eta_point[1].c = 0; /* Right */
        eta_point[2].i = ex; eta_point[2].j = ey;   eta_point[2].loc = DMSTAG_ELEMENT;    eta_point[2].c = 0; /* Up    */
        eta_point[3].i = ex; eta_point[3].j = ey-1; eta_point[3].loc = DMSTAG_ELEMENT;    eta_point[3].c = 0; /* Down  */
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,4,eta_point,eta);CHKERRQ(ierr);
        eta_left = eta[0]; eta_right = eta[1]; eta_up = eta[2]; eta_down = eta[3];

        if (ex == 0) {
          /* Left boundary y velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j    = ey  ; row.loc    = DMSTAG_DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j = ey  ; col[0].loc = DMSTAG_DOWN;     col[0].c  = 0; val_A[0]  = -2.0 * dv * (eta_down + eta_up) / (hy*hy) - dv * (eta_right) /(hx*hx);
          col[1].i = ex  ; col[1].j = ey-1; col[1].loc = DMSTAG_DOWN;     col[1].c  = 0; val_A[1]  =  2.0 * dv * eta_down  / (hy*hy);
          col[2].i = ex  ; col[2].j = ey+1; col[2].loc = DMSTAG_DOWN;     col[2].c  = 0; val_A[2]  =  2.0 * dv * eta_up    / (hy*hy);
          /* No left entry */
          col[3].i = ex+1; col[3].j = ey  ; col[3].loc = DMSTAG_DOWN;     col[3].c  = 0; val_A[3]  =        dv * eta_right / (hx*hx);
          col[4].i = ex  ; col[4].j = ey-1; col[4].loc = DMSTAG_LEFT;     col[4].c  = 0; val_A[4]  =        dv * eta_left  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j = ey-1; col[5].loc = DMSTAG_RIGHT;    col[5].c  = 0; val_A[5]  = -1.0 * dv * eta_right / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j = ey  ; col[6].loc = DMSTAG_LEFT;     col[6].c  = 0; val_A[6]  = -1.0 * dv * eta_left  / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j = ey  ; col[7].loc = DMSTAG_RIGHT;    col[7].c  = 0; val_A[7]  =        dv * eta_right / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j = ey-1; col[8].loc = DMSTAG_ELEMENT;  col[8].c  = 0; val_A[8]  =        Kcont * dv / hy;
          col[9].i = ex  ; col[9].j = ey  ; col[9].loc = DMSTAG_ELEMENT;  col[9].c  = 0; val_A[9]  = -1.0 * Kcont * dv / hy;
        } else if (ex == N[0]-1) {
          /* Right boundary y velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j    = ey  ; row.loc    = DMSTAG_DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j = ey  ; col[0].loc = DMSTAG_DOWN;     col[0].c  = 0; val_A[0]  = -2.0 * dv * (eta_down + eta_up) / (hy*hy) - dv * (eta_left) /(hx*hx);
          col[1].i = ex  ; col[1].j = ey-1; col[1].loc = DMSTAG_DOWN;     col[1].c  = 0; val_A[1]  =  2.0 * dv * eta_down  / (hy*hy);
          col[2].i = ex  ; col[2].j = ey+1; col[2].loc = DMSTAG_DOWN;     col[2].c  = 0; val_A[2]  =  2.0 * dv * eta_up    / (hy*hy);
          col[3].i = ex-1; col[3].j = ey  ; col[3].loc = DMSTAG_DOWN;     col[3].c  = 0; val_A[3]  =        dv * eta_left  / (hx*hx);
          /* No right element */
          col[4].i = ex  ; col[4].j = ey-1; col[4].loc = DMSTAG_LEFT;     col[4].c  = 0; val_A[4]  =        dv * eta_left  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j = ey-1; col[5].loc = DMSTAG_RIGHT;    col[5].c  = 0; val_A[5]  = -1.0 * dv * eta_right / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j = ey  ; col[6].loc = DMSTAG_LEFT;     col[6].c  = 0; val_A[7]  = -1.0 * dv * eta_left  / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j = ey  ; col[7].loc = DMSTAG_RIGHT;    col[7].c  = 0; val_A[7]  =        dv * eta_right / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j = ey-1; col[8].loc = DMSTAG_ELEMENT;  col[8].c  = 0; val_A[8]  =        Kcont * dv / hy;
          col[9].i = ex  ; col[9].j = ey  ; col[9].loc = DMSTAG_ELEMENT;  col[9].c  = 0; val_A[9]  = -1.0 * Kcont * dv / hy;
        } else {
          /* U_y interior equation */
          nEntries = 11;
          row.i    = ex  ; row.j     = ey  ; row.loc     = DMSTAG_DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DMSTAG_DOWN;     col[0].c  = 0; val_A[0]  = -2.0 * dv * (eta_down + eta_up) / (hy*hy) - dv * (eta_left + eta_right) /(hx*hx);
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DMSTAG_DOWN;     col[1].c  = 0; val_A[1]  =  2.0 * dv * eta_down  / (hy*hy);
          col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DMSTAG_DOWN;     col[2].c  = 0; val_A[2]  =  2.0 * dv * eta_up    / (hy*hy);
          col[3].i = ex-1; col[3].j  = ey  ; col[3].loc  = DMSTAG_DOWN;     col[3].c  = 0; val_A[3]  =        dv * eta_left  / (hx*hx);
          col[4].i = ex+1; col[4].j  = ey  ; col[4].loc  = DMSTAG_DOWN;     col[4].c  = 0; val_A[4]  =        dv * eta_right / (hx*hx);
          col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = DMSTAG_LEFT;     col[5].c  = 0; val_A[5]  =        dv * eta_left  / (hx*hy); /* down left x edge */
          col[6].i = ex  ; col[6].j  = ey-1; col[6].loc  = DMSTAG_RIGHT;    col[6].c  = 0; val_A[6]  = -1.0 * dv * eta_right / (hx*hy); /* down right x edge */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = DMSTAG_LEFT;     col[7].c  = 0; val_A[7]  = -1.0 * dv * eta_left  / (hx*hy); /* up left x edge */
          col[8].i = ex  ; col[8].j  = ey  ; col[8].loc  = DMSTAG_RIGHT;    col[8].c  = 0; val_A[8]  =        dv * eta_right / (hx*hy); /* up right x edge */
          col[9].i = ex  ; col[9].j  = ey-1; col[9].loc  = DMSTAG_ELEMENT;  col[9].c  = 0; val_A[9]  =        Kcont * dv / hy;
          col[10].i = ex ; col[10].j = ey  ; col[10].loc = DMSTAG_ELEMENT; col[10].c  = 0; val_A[10] = -1.0 * Kcont * dv / hy;
        }

        /* Insert Y-momentum entries */
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,nEntries,col,val_A,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      if (ex == N[0]-1) {
        /* Right Boundary velocity Dirichlet */
        /* Redundant in the corner */
        DMStagStencil row;
        PetscScalar   val_rhs;

        const PetscScalar val_A = Kbound;
        row.i = ex; row.j = ey; row.loc = DMSTAG_RIGHT; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
        val_rhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
      }
      if (ex == 0) {
        /* Left velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   val_rhs;
        const PetscScalar val_A = Kbound;
        row.i = ex; row.j = ey; row.loc = DMSTAG_LEFT; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
        val_rhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
        PetscInt nEntries;
        DMStagStencil row,col[11];
        PetscScalar   val_rhs,val_A[11];
        DMStagStencil eta_point[4];
        PetscScalar eta[4],eta_left,eta_right,eta_up,eta_down;

        /* Get eta values */
        eta_point[0].i = ex-1; eta_point[0].j = ey; eta_point[0].loc = DMSTAG_ELEMENT;   eta_point[0].c = 0; /* Left  */
        eta_point[1].i = ex;   eta_point[1].j = ey; eta_point[1].loc = DMSTAG_ELEMENT;   eta_point[1].c = 0; /* Right */
        eta_point[2].i = ex;   eta_point[2].j = ey; eta_point[2].loc = DMSTAG_UP_LEFT;   eta_point[2].c = 0; /* Up    */
        eta_point[3].i = ex;   eta_point[3].j = ey; eta_point[3].loc = DMSTAG_DOWN_LEFT; eta_point[3].c = 0; /* Down  */
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,4,eta_point,eta);CHKERRQ(ierr);
        eta_left = eta[0]; eta_right = eta[1]; eta_up = eta[2]; eta_down = eta[3];

        if (ey == 0) {
          /* Bottom boundary x velocity stencil (with zero vel deriv) */
          nEntries = 10;
          row.i     = ex  ; row.j     = ey  ; row.loc    = DMSTAG_LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc = DMSTAG_LEFT;    col[0].c   = 0; val_A[0]  = -2.0 * dv * (eta_left + eta_right) / (hx*hx) - dv * (eta_up) / (hy*hy);
          /* Missing element below */
          col[1].i = ex  ; col[1].j  = ey+1; col[1].loc  = DMSTAG_LEFT;    col[1].c   = 0; val_A[1]  =        dv * eta_up    / (hy*hy);
          col[2].i = ex-1; col[2].j  = ey  ; col[2].loc  = DMSTAG_LEFT;    col[2].c   = 0; val_A[2]  =  2.0 * dv * eta_left  / (hx*hx);
          col[3].i = ex+1; col[3].j  = ey  ; col[3].loc  = DMSTAG_LEFT;    col[3].c   = 0; val_A[3]  =  2.0 * dv * eta_right / (hx*hx);
          col[4].i = ex-1; col[4].j  = ey  ; col[4].loc  = DMSTAG_DOWN;    col[4].c   = 0; val_A[4]  =        dv * eta_down  / (hx*hy); /* down left */
          col[5].i = ex  ; col[5].j  = ey  ; col[5].loc  = DMSTAG_DOWN;    col[5].c   = 0; val_A[5]  = -1.0 * dv * eta_down  / (hx*hy); /* down right */
          col[6].i = ex-1; col[6].j  = ey  ; col[6].loc  = DMSTAG_UP;      col[6].c   = 0; val_A[6]  = -1.0 * dv * eta_up    / (hx*hy); /* up left */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = DMSTAG_UP;      col[7].c   = 0; val_A[7]  =        dv * eta_up    / (hx*hy); /* up right */
          col[8].i = ex-1; col[8].j  = ey  ; col[8].loc  = DMSTAG_ELEMENT; col[8].c   = 0; val_A[8]  =        Kcont * dv / hx;
          col[9].i = ex  ; col[9].j  = ey  ; col[9].loc  = DMSTAG_ELEMENT; col[9].c   = 0; val_A[9]  = -1.0 * Kcont * dv / hx;
          val_rhs = 0.0;
        } else if (ey == N[1]-1) {
          /* Top boundary x velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j     = ey  ; row.loc     = DMSTAG_LEFT;    row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DMSTAG_LEFT;    col[0].c  = 0; val_A[0]  = -2.0 * dv * (eta_left + eta_right) / (hx*hx) - dv * (eta_down) / (hy*hy);
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DMSTAG_LEFT;    col[1].c  = 0; val_A[1]  =        dv * eta_down  / (hy*hy);
          /* Missing element above */
          col[2].i = ex-1; col[2].j  = ey  ; col[2].loc  = DMSTAG_LEFT;    col[2].c  = 0; val_A[2]  =  2.0 * dv * eta_left  / (hx*hx);
          col[3].i = ex+1; col[3].j  = ey  ; col[3].loc  = DMSTAG_LEFT;    col[3].c  = 0; val_A[3]  =  2.0 * dv * eta_right / (hx*hx);
          col[4].i = ex-1; col[4].j  = ey  ; col[4].loc  = DMSTAG_DOWN;    col[4].c  = 0; val_A[4]  =        dv * eta_down  / (hx*hy); /* down left */
          col[5].i = ex  ; col[5].j  = ey  ; col[5].loc  = DMSTAG_DOWN;    col[5].c  = 0; val_A[5]  = -1.0 * dv * eta_down  / (hx*hy); /* down right */
          col[6].i = ex-1; col[6].j  = ey  ; col[6].loc  = DMSTAG_UP;      col[6].c  = 0; val_A[6]  = -1.0 * dv * eta_up    / (hx*hy); /* up left */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = DMSTAG_UP;      col[7].c  = 0; val_A[7]  =        dv * eta_up    / (hx*hy); /* up right */
          col[8].i = ex-1; col[8].j  = ey  ; col[8].loc  = DMSTAG_ELEMENT; col[8].c  = 0; val_A[8]  =        Kcont * dv / hx;
          col[9].i = ex  ; col[9].j  = ey  ; col[9].loc  = DMSTAG_ELEMENT; col[9].c  = 0; val_A[9]  = -1.0 * Kcont * dv / hx;
          val_rhs = 0.0;
        } else {
          /* U_x interior equation */
          nEntries = 11;
          row.i     = ex  ; row.j     = ey  ; row.loc     = DMSTAG_LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = DMSTAG_LEFT;    col[0].c   = 0; val_A[0]  = -2.0 * dv * (eta_left + eta_right) / (hx*hx) - dv * (eta_up + eta_down) / (hy*hy);
          col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc  = DMSTAG_LEFT;    col[1].c   = 0; val_A[1]  =        dv * eta_down  / (hy*hy);
          col[2].i  = ex  ; col[2].j  = ey+1; col[2].loc  = DMSTAG_LEFT;    col[2].c   = 0; val_A[2]  =        dv * eta_up    / (hy*hy);
          col[3].i  = ex-1; col[3].j  = ey  ; col[3].loc  = DMSTAG_LEFT;    col[3].c   = 0; val_A[3]  =  2.0 * dv * eta_left  / (hx*hx);
          col[4].i  = ex+1; col[4].j  = ey  ; col[4].loc  = DMSTAG_LEFT;    col[4].c   = 0; val_A[4]  =  2.0 * dv * eta_right / (hx*hx);
          col[5].i  = ex-1; col[5].j  = ey  ; col[5].loc  = DMSTAG_DOWN;    col[5].c   = 0; val_A[5]  =        dv * eta_down  / (hx*hy); /* down left */
          col[6].i  = ex  ; col[6].j  = ey  ; col[6].loc  = DMSTAG_DOWN;    col[6].c   = 0; val_A[6]  = -1.0 * dv * eta_down  / (hx*hy); /* down right */
          col[7].i  = ex-1; col[7].j  = ey  ; col[7].loc  = DMSTAG_UP;      col[7].c   = 0; val_A[7]  = -1.0 * dv * eta_up    / (hx*hy); /* up left */
          col[8].i  = ex  ; col[8].j  = ey  ; col[8].loc  = DMSTAG_UP;      col[8].c   = 0; val_A[8]  =        dv * eta_up    / (hx*hy); /* up right */
          col[9].i  = ex-1; col[9].j  = ey  ; col[9].loc  = DMSTAG_ELEMENT; col[9].c   = 0; val_A[9]  =        Kcont * dv / hx;
          col[10].i = ex  ; col[10].j = ey  ; col[10].loc = DMSTAG_ELEMENT; col[10].c  = 0; val_A[10] = -1.0 * Kcont * dv / hx;
          val_rhs = 0.0;
        }
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,nEntries,col,val_A,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      /* P equation : u_x + v_y = 0
         Note that this includes an explicit zero on the diagonal. This is only needed for
         direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace)

         Note: the scaling by dv is not chosen in a principled way and is likely sub-optimal */
      if (pin_pressure && ex == pinx && ey == piny) { /* Pin a pressure node to zero, if requested */
        DMStagStencil row;
        PetscScalar val_A,val_rhs;
        row.i = ex; row.j = ey; row.loc = DMSTAG_ELEMENT; row.c = 0;
        val_A = Kbound;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
        val_rhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        DMStagStencil row,col[5];
        PetscScalar   val_A[5],val_rhs;

        row.i    = ex; row.j    = ey; row.loc    = DMSTAG_ELEMENT; row.c    = 0;
        col[0].i = ex; col[0].j = ey; col[0].loc = DMSTAG_LEFT;    col[0].c = 0; val_A[0] = - 1.0 * Kcont * dv / hx;
        col[1].i = ex; col[1].j = ey; col[1].loc = DMSTAG_RIGHT;   col[1].c = 0; val_A[1] =         Kcont * dv / hx;
        col[2].i = ex; col[2].j = ey; col[2].loc = DMSTAG_DOWN;    col[2].c = 0; val_A[2] = - 1.0 * Kcont * dv / hy;
        col[3].i = ex; col[3].j = ey; col[3].loc = DMSTAG_UP;      col[3].c = 0; val_A[3] =         Kcont * dv / hy;
        col[4] = row;                                                            val_A[4] = 0.0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,5,col,val_A,INSERT_VALUES);CHKERRQ(ierr);
        val_rhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  if (parameters->boussinesq_forcing) {
    ierr = DMStagVecRestoreArrayRead(dm_temperature,temperature_local,&arr_temperature);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm_temperature,&temperature_local);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rhs);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateSystem_3D_FreeSlip(StagBLStokesParameters parameters,StagBLSystem system)
{
  PetscErrorCode  ierr;
  DM              dm_stokes,dm_coefficient;
  PetscInt        N[3];
  PetscInt        ex,ey,ez,startx,starty,startz,nx,ny,nz;
  Mat             *pA;
  Vec             *pRhs;
  Mat             A;
  PetscReal       hx,hy,hz,dv,hxAvgInv,Kcont,Kbound;
  PetscInt        pinx,piny,pinz;
  const PetscBool pin_pressure = PETSC_TRUE;
  StagBLGrid      coefficient_grid;
  DM              dm_temperature;
  Vec             coeff_local,rhs,temperature,temperature_local;
  PetscScalar     ****arr_temperature;
  PetscInt        slot_temperature_backdownleft,slot_temperature_frontdownleft,slot_temperature_backdownright,slot_temperature_frontdownright;

  PetscFunctionBeginUser;
  ierr = StagBLSystemPETScGetMatPointer(system,&pA);CHKERRQ(ierr);
  ierr = StagBLSystemPETScGetVecPointer(system,&pRhs);CHKERRQ(ierr);
  if (!parameters->coefficient_array) StagBLError(PETSC_COMM_SELF,"coefficient_array field not set in StagBLStokesParameters argument");
  ierr = StagBLArrayGetStagBLGrid(parameters->coefficient_array,&coefficient_grid);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(parameters->stokes_grid,&dm_stokes);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(coefficient_grid,&dm_coefficient);CHKERRQ(ierr);
  ierr = StagBLArrayPETScGetLocalVec(parameters->coefficient_array,&coeff_local);CHKERRQ(ierr);

  /* Compute some parameters */
  ierr = DMStagGetGlobalSizes(dm_stokes,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
  if (parameters->uniform_grid) {
    hx = (parameters->xmax - parameters->xmin)/N[0];
    hy = (parameters->ymax - parameters->ymin)/N[1];
    hz = (parameters->zmax - parameters->zmin)/N[2];
    dv = hx*hy*hz;;
  } else StagBLError(PetscObjectComm((PetscObject)dm_stokes),"Non-uniform grids not supported yet");
  hxAvgInv = 3.0/(hx + hy + hz);
  Kcont  = parameters->eta_characteristic*hxAvgInv;
  Kbound = parameters->eta_characteristic*hxAvgInv*hxAvgInv;
  if (N[0] < 2) SETERRQ(PetscObjectComm((PetscObject)dm_stokes),PETSC_ERR_SUP,"Not implemented for a single element in the x direction");
  pinx = 1; piny = 0; pinz = 0;

  ierr = DMCreateMatrix(dm_stokes,pA);CHKERRQ(ierr);
  A = *pA;
  ierr = DMCreateGlobalVector(dm_stokes,pRhs);CHKERRQ(ierr);
  rhs = *pRhs;
  ierr = DMStagGetCorners(dm_stokes,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);

  /* If using Boussinesq forcing from a temperature field, get access */
  if (parameters->boussinesq_forcing) {
    ierr = StagBLGridPETScGetDM(parameters->temperature_grid,&dm_temperature);CHKERRQ(ierr);
    ierr = StagBLArrayPETScGetGlobalVec(parameters->temperature_array,&temperature);CHKERRQ(ierr);
    ierr = DMGetLocalVector(dm_temperature,&temperature_local);CHKERRQ(ierr);
    ierr = DMGlobalToLocal(dm_temperature,temperature,INSERT_VALUES,temperature_local);CHKERRQ(ierr);
    ierr = DMStagVecGetArrayRead(dm_temperature,temperature_local,&arr_temperature);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_BACK_DOWN_LEFT,  0,&slot_temperature_backdownleft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_BACK_DOWN_RIGHT, 0,&slot_temperature_backdownright);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_FRONT,           0,&slot_temperature_frontdownleft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_FRONT,           0,&slot_temperature_frontdownright);CHKERRQ(ierr);
  }

  /* Loop over all local elements. */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ey = starty; ey < starty + ny; ++ey) {
      for (ex = startx; ex < startx + nx; ++ex) {
        STAGBL_UNUSED(pin_pressure);
        StagBLError(PETSC_COMM_WORLD,"Not implemented");
      }
    }
  }

  if (parameters->boussinesq_forcing) {
    ierr = DMStagVecRestoreArrayRead(dm_temperature,temperature_local,&arr_temperature);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm_temperature,&temperature_local);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rhs);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
