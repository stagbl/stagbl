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
        PetscScalar val_rhs;
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
        DMStagStencil rho_point[2];
        PetscScalar   rho[2],val_rhs;
        DMStagStencil eta_point[4];
        PetscScalar   eta[4],eta_left,eta_right,eta_up,eta_down;

        /* get rho values  (note .c = 1) */
        rho_point[0].i = ex; rho_point[0].j = ey; rho_point[0].loc = DMSTAG_DOWN_LEFT;  rho_point[0].c = 1;
        rho_point[1].i = ex; rho_point[1].j = ey; rho_point[1].loc = DMSTAG_DOWN_RIGHT; rho_point[1].c = 1;
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,2,rho_point,rho);CHKERRQ(ierr);

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
        DMStagStencil     row;
        const PetscScalar val_rhs = 0.0;
        const PetscScalar val_A = Kbound;

        row.i = ex; row.j = ey; row.loc = DMSTAG_RIGHT; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
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
  if (N[0] < 2 || N[1] < 2 || N[2] < 2) StagBLError(PetscObjectComm((PetscObject)dm_stokes),"Stokes system construction  implemented for a single element in any direction");
  if (parameters->uniform_grid) {
    hx = (parameters->xmax - parameters->xmin)/N[0];
    hy = (parameters->ymax - parameters->ymin)/N[1];
    hz = (parameters->zmax - parameters->zmin)/N[2];
    dv = hx*hy*hz;
  } else StagBLError(PetscObjectComm((PetscObject)dm_stokes),"Non-uniform grids not supported yet");
  hxAvgInv = 3.0/(hx + hy + hz);
  Kcont  = parameters->eta_characteristic*hxAvgInv;
  Kbound = parameters->eta_characteristic*hxAvgInv*hxAvgInv;
  pinx = 1; piny = 0; pinz = 0; /* Depends on assertion above that there are at least two element in the x direction */

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
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_FRONT_DOWN_LEFT, 0,&slot_temperature_frontdownleft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_FRONT_DOWN_RIGHT,0,&slot_temperature_frontdownright);CHKERRQ(ierr);
  }

  /* Loop over all local elements.

     For each element, fill 4-7 rows of the matrix, corresponding to
     - the pressure degree of freedom (dof), centered on the element
     - the 3 velocity dofs on left, bottom, and back faces of the element
     - velocity dof on the right, top, and front faces of the element (only on domain boundaries)

   */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ey = starty; ey < starty + ny; ++ey) {
      for (ex = startx; ex < startx + nx; ++ex) {
        const PetscBool left_boundary   = ex == 0;
        const PetscBool right_boundary  = ex == N[0]-1;
        const PetscBool bottom_boundary = ey == 0;
        const PetscBool top_boundary    = ey == N[1]-1;
        const PetscBool back_boundary   = ez == 0;
        const PetscBool front_boundary  = ez == N[2]-1;

        /* Note that below, we depend on the check above that there is never one
           element (globally) in a given direction.  Thus, for example, an
           element is never both on the left and right boundary */

        /* X-faces - right boundary */
        if (right_boundary) {
          /* Right x-velocity Dirichlet */
          DMStagStencil     row;
          const PetscScalar val_rhs = 0.0;
          const PetscScalar val_A = Kbound;

          row.i = ex; row.j = ey; row.k = ez; row.loc = DMSTAG_RIGHT; row.c = 0;
          ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
          ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* X faces - left*/
        {
          DMStagStencil row;

          row.i = ex; row.j = ey; row.k = ez; row.loc = DMSTAG_LEFT; row.c = 0;

          if (left_boundary) {
            /* Left x-velocity Dirichlet */
            const PetscScalar val_rhs = 0.0;
            const PetscScalar val_A = Kbound;

            ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
            ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
          } else {
            /* X-momentum equation */
            PetscInt      count;
            DMStagStencil col[17];
            PetscScalar   val_rhs,val_A[17];
            DMStagStencil eta_point[6];
            PetscScalar   eta[6],eta_left,eta_right,eta_up,eta_down,eta_back,eta_front; /* relative to the left face */

            /* Get eta values */
            eta_point[0].i = ex-1; eta_point[0].j = ey; eta_point[0].k = ez; eta_point[0].loc = DMSTAG_ELEMENT;    eta_point[0].c = 0; /* Left  */
            eta_point[1].i = ex;   eta_point[1].j = ey; eta_point[1].k = ez; eta_point[1].loc = DMSTAG_ELEMENT;    eta_point[1].c = 0; /* Right */
            eta_point[2].i = ex;   eta_point[2].j = ey; eta_point[2].k = ez; eta_point[2].loc = DMSTAG_UP_LEFT;    eta_point[2].c = 0; /* Up    */
            eta_point[3].i = ex;   eta_point[3].j = ey; eta_point[3].k = ez; eta_point[3].loc = DMSTAG_DOWN_LEFT;  eta_point[3].c = 0; /* Down  */
            eta_point[4].i = ex;   eta_point[4].j = ey; eta_point[4].k = ez; eta_point[4].loc = DMSTAG_BACK_LEFT;  eta_point[4].c = 0; /* Back  */
            eta_point[5].i = ex;   eta_point[5].j = ey; eta_point[5].k = ez; eta_point[5].loc = DMSTAG_FRONT_LEFT; eta_point[5].c = 0; /* Front  */
            ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,6,eta_point,eta);CHKERRQ(ierr);
            eta_left = eta[0]; eta_right = eta[1]; eta_up = eta[2]; eta_down = eta[3]; eta_back = eta[4]; eta_front = eta[5];

            count = 0;

            col[count] = row;
            val_A[count] = -2.0 * dv * (eta_left + eta_right) / (hx*hx);
            if (!top_boundary)    val_A[count] += -1.0 * dv * eta_up    / (hy*hy);
            if (!bottom_boundary) val_A[count] += -1.0 * dv * eta_down  / (hy*hy);
            if (!back_boundary)   val_A[count] += -1.0 * dv * eta_back  / (hz*hz);
            if (!front_boundary)  val_A[count] += -1.0 * dv * eta_front / (hz*hz);
            ++count;

            col[count].i = ex-1; col[count].j = ey; col[count].k = ez; col[count].loc = DMSTAG_LEFT; col[count].c = 0;
            val_A[count] = 2.0 * dv * eta_left  / (hx*hx); ++count;
            col[count].i = ex+1; col[count].j = ey; col[count].k = ez; col[count].loc = DMSTAG_LEFT; col[count].c = 0;
            val_A[count] = 2.0 * dv * eta_right  / (hx*hx); ++count;
            if (!bottom_boundary) {
              col[count].i = ex; col[count].j = ey-1; col[count].k = ez; col[count].loc = DMSTAG_LEFT; col[count].c = 0;
              val_A[count] = dv * eta_down / (hy*hy); ++count;
            }
            if (!top_boundary) {
              col[count].i = ex; col[count].j = ey+1; col[count].k = ez; col[count].loc = DMSTAG_LEFT; col[count].c = 0;
              val_A[count] = dv * eta_up / (hy*hy); ++count;
            }
            if (!back_boundary) {
              col[count].i = ex; col[count].j = ey; col[count].k = ez-1; col[count].loc = DMSTAG_LEFT; col[count].c = 0;
              val_A[count] = dv * eta_back / (hz*hz); ++count;
            }
            if (!front_boundary) {
              col[count].i = ex; col[count].j = ey; col[count].k = ez+1; col[count].loc = DMSTAG_LEFT; col[count].c = 0;
              val_A[count] = dv * eta_front / (hz*hz); ++count;
            }

            col[count].i  = ex-1; col[count].j  = ey; col[count].k = ez; col[count].loc  = DMSTAG_DOWN;  col[count].c = 0;
            val_A[count]  =        dv * eta_down  / (hx*hy); ++count; /* down left */
            col[count].i  = ex  ; col[count].j  = ey; col[count].k = ez; col[count].loc  = DMSTAG_DOWN;  col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_down  / (hx*hy); ++count; /* down right */

            col[count].i  = ex-1; col[count].j  = ey; col[count].k = ez; col[count].loc  = DMSTAG_UP;    col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_up    / (hx*hy); ++count; /* up left */
            col[count].i  = ex  ; col[count].j  = ey; col[count].k = ez; col[count].loc  = DMSTAG_UP;    col[count].c = 0;
            val_A[count]  =        dv * eta_up    / (hx*hy); ++count; /* up right */

            col[count].i  = ex-1; col[count].j  = ey; col[count].k = ez; col[count].loc  = DMSTAG_BACK;  col[count].c = 0;
            val_A[count]  =        dv * eta_back  / (hx*hz); ++count; /* back left */
            col[count].i  = ex  ; col[count].j  = ey; col[count].k = ez; col[count].loc  = DMSTAG_BACK;  col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_back  / (hx*hz); ++count; /* back right */

            col[count].i  = ex-1; col[count].j  = ey; col[count].k = ez; col[count].loc  = DMSTAG_FRONT; col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_front / (hx*hz); ++count; /* front left */
            col[count].i  = ex  ; col[count].j  = ey; col[count].k = ez; col[count].loc  = DMSTAG_FRONT; col[count].c = 0;
            val_A[count]  =        dv * eta_front / (hx*hz); ++count; /* front right */

            col[count].i = ex-1; col[count].j = ey; col[count].k = ez; col[count].loc = DMSTAG_ELEMENT; col[count].c  = 0;
            val_A[count] = Kcont * dv / hx; ++count;
            col[count].i = ex;   col[count].j = ey; col[count].k = ez; col[count].loc = DMSTAG_ELEMENT; col[count].c  = 0;
            val_A[count] = -1.0 * Kcont * dv / hx; ++count;

            ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,count,col,val_A,INSERT_VALUES);CHKERRQ(ierr);
            val_rhs = 0.0; /* No forcing */
            ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
          }
        }

        /* Y faces - top boundary */
        if (top_boundary) {
          /* Top y-velocity Dirichlet */
          DMStagStencil     row;
          const PetscScalar val_rhs = 0.0;
          const PetscScalar val_A = Kbound;

          row.i = ex; row.j = ey; row.k = ez; row.loc = DMSTAG_UP; row.c = 0;
          ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
          ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Y faces - down */
        {
          DMStagStencil row;

          row.i = ex; row.j = ey; row.k = ez; row.loc = DMSTAG_DOWN; row.c = 0;

          if (bottom_boundary) {
            /* Bottom y-velocity Dirichlet */
            const PetscScalar val_rhs = 0.0;
            const PetscScalar val_A = Kbound;

            ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
            ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
          } else {
            /* Y-momentum equation (including non-zero forcing) */
            PetscInt      count;
            DMStagStencil col[17];
            PetscScalar   val_rhs,val_A[17];
            DMStagStencil eta_point[6],rho_point[4];
            PetscScalar   eta[6],rho[4],eta_left,eta_right,eta_up,eta_down,eta_back,eta_front; /* relative to the bottom face */

            /* get rho values  (note .c = 1) */
            // TODO we have rho at perhaps strange points (edges not corners).. don't see how it really matters when getting off a grid and averaging, though. This is kinda wacky.. why not just have rho at the faces??)
            rho_point[0].i = ex; rho_point[0].j = ey; rho_point[0].k = ez; rho_point[0].loc = DMSTAG_DOWN_LEFT;  rho_point[0].c = 1;
            rho_point[1].i = ex; rho_point[1].j = ey; rho_point[1].k = ez; rho_point[1].loc = DMSTAG_DOWN_RIGHT; rho_point[1].c = 1;
            rho_point[2].i = ex; rho_point[2].j = ey; rho_point[2].k = ez; rho_point[2].loc = DMSTAG_BACK_DOWN;  rho_point[2].c = 1;
            rho_point[3].i = ex; rho_point[3].j = ey; rho_point[3].k = ez; rho_point[3].loc = DMSTAG_FRONT_DOWN; rho_point[3].c = 1;
            ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,4,rho_point,rho);CHKERRQ(ierr);

            /* Compute forcing */
            {
              const PetscReal rho_avg = 0.25 * (rho[0] + rho[1] + rho[2] + rho[3]);

              /* Note Boussinesq forcing has the opposite sign */
              if (parameters->boussinesq_forcing) {
                const PetscScalar temperature_average = 0.25 * (
                      arr_temperature[ez][ey][ex][slot_temperature_backdownleft]
                    + arr_temperature[ez][ey][ex][slot_temperature_backdownright]
                    + arr_temperature[ez][ey][ex][slot_temperature_frontdownleft]
                    + arr_temperature[ez][ey][ex][slot_temperature_frontdownright]
                    );
                val_rhs = parameters->alpha * rho_avg * parameters->gy * dv * temperature_average;
              } else {
                val_rhs = -parameters->gy * dv * rho_avg; // TODO this sign convention is crazy, just pick something consistent (fix in 2D stokes function as well)
              }
            }

            /* Get eta values */
            eta_point[0].i = ex; eta_point[0].j = ey;   eta_point[0].k = ez; eta_point[0].loc = DMSTAG_DOWN_LEFT;  eta_point[0].c = 0; /* Left  */
            eta_point[1].i = ex; eta_point[1].j = ey;   eta_point[1].k = ez; eta_point[1].loc = DMSTAG_DOWN_RIGHT; eta_point[1].c = 0; /* Right */
            eta_point[2].i = ex; eta_point[2].j = ey;   eta_point[2].k = ez; eta_point[2].loc = DMSTAG_ELEMENT;    eta_point[2].c = 0; /* Up    */
            eta_point[3].i = ex; eta_point[3].j = ey-1; eta_point[3].k = ez; eta_point[3].loc = DMSTAG_ELEMENT;    eta_point[3].c = 0; /* Down  */ // TODO sucks that up and down are flipped (in all places like this in this file)
            eta_point[4].i = ex; eta_point[4].j = ey;   eta_point[4].k = ez; eta_point[4].loc = DMSTAG_BACK_DOWN;  eta_point[4].c = 0; /* Back  */
            eta_point[5].i = ex; eta_point[5].j = ey;   eta_point[5].k = ez; eta_point[5].loc = DMSTAG_FRONT_DOWN; eta_point[5].c = 0; /* Front  */
            ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,6,eta_point,eta);CHKERRQ(ierr);
            eta_left = eta[0]; eta_right = eta[1]; eta_up = eta[2]; eta_down = eta[3]; eta_back = eta[4]; eta_front = eta[5];

            count = 0;

            col[count] = row;
            val_A[count] = -2.0 * dv * (eta_up + eta_down) / (hy*hy);
            if (!left_boundary)  val_A[count] += -1.0 * dv * eta_left  / (hx*hx);
            if (!right_boundary) val_A[count] += -1.0 * dv * eta_right / (hx*hx);
            if (!back_boundary)  val_A[count] += -1.0 * dv * eta_back  / (hz*hz);
            if (!front_boundary) val_A[count] += -1.0 * dv * eta_front / (hz*hz);
            ++count;

            col[count].i = ex; col[count].j = ey-1; col[count].k = ez; col[count].loc = DMSTAG_DOWN; col[count].c = 0;
            val_A[count] = 2.0 * dv * eta_down / (hy*hy); ++count;
            col[count].i = ex; col[count].j = ey+1; col[count].k = ez; col[count].loc = DMSTAG_DOWN; col[count].c = 0;
            val_A[count] = 2.0 * dv * eta_up   / (hy*hy); ++count;

            if (!left_boundary) {
              col[count].i = ex-1; col[count].j = ey; col[count].k = ez; col[count].loc = DMSTAG_DOWN; col[count].c = 0;
              val_A[count] = dv * eta_left / (hx*hx); ++count;
            }
            if (!right_boundary) {
              col[count].i = ex+1; col[count].j = ey; col[count].k = ez; col[count].loc = DMSTAG_DOWN; col[count].c = 0;
              val_A[count] = dv * eta_right / (hx*hx); ++count;
            }
            if (!back_boundary) {
              col[count].i = ex; col[count].j = ey; col[count].k = ez-1; col[count].loc = DMSTAG_DOWN; col[count].c = 0;
              val_A[count] = dv * eta_back / (hz*hz); ++count;
            }
            if (!front_boundary) {
              col[count].i = ex; col[count].j = ey; col[count].k = ez+1; col[count].loc = DMSTAG_DOWN; col[count].c = 0;
              val_A[count] = dv * eta_front / (hz*hz); ++count;
            }

            col[count].i  = ex; col[count].j  = ey-1; col[count].k = ez; col[count].loc = DMSTAG_LEFT;  col[count].c = 0;
            val_A[count]  =        dv * eta_left  / (hx*hy); ++count; /* down left*/
            col[count].i  = ex; col[count].j  = ey;   col[count].k = ez; col[count].loc = DMSTAG_LEFT;  col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_left  / (hx*hy); ++count; /* up left*/

            col[count].i  = ex; col[count].j  = ey-1; col[count].k = ez; col[count].loc = DMSTAG_RIGHT; col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_right / (hx*hy); ++count; /* down right*/
            col[count].i  = ex; col[count].j  = ey;   col[count].k = ez; col[count].loc = DMSTAG_RIGHT; col[count].c = 0;
            val_A[count]  =        dv * eta_right / (hx*hy); ++count; /* up right*/

            col[count].i  = ex; col[count].j  = ey-1; col[count].k = ez; col[count].loc  = DMSTAG_BACK;  col[count].c = 0;
            val_A[count]  =        dv * eta_back  / (hy*hz); ++count; /* back down */
            col[count].i  = ex; col[count].j  = ey;   col[count].k = ez; col[count].loc  = DMSTAG_BACK;  col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_back  / (hy*hz); ++count;/* back up */

            col[count].i  = ex; col[count].j  = ey-1; col[count].k = ez; col[count].loc  = DMSTAG_FRONT; col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_front / (hy*hz); ++count; /* front down */
            col[count].i  = ex; col[count].j  = ey;   col[count].k = ez; col[count].loc  = DMSTAG_FRONT; col[count].c = 0;
            val_A[count]  =        dv * eta_front / (hy*hz); ++count;/* front up */

            col[count].i = ex; col[count].j = ey-1; col[count].k = ez; col[count].loc = DMSTAG_ELEMENT; col[count].c  = 0;
            val_A[count] = Kcont * dv / hy; ++count;
            col[count].i = ex; col[count].j = ey;   col[count].k = ez; col[count].loc = DMSTAG_ELEMENT; col[count].c  = 0;
            val_A[count] = -1.0 * Kcont * dv / hy; ++count;

            ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,count,col,val_A,INSERT_VALUES);CHKERRQ(ierr);
            ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
          }
        }

        if (front_boundary) {
          /* Front z-velocity Dirichlet */
          DMStagStencil     row;
          const PetscScalar val_rhs = 0.0;
          const PetscScalar val_A = Kbound;

          row.i = ex; row.j = ey; row.k = ez; row.loc = DMSTAG_FRONT; row.c = 0;
          ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
          ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Z faces - back */
        {
          DMStagStencil row;

          row.i = ex; row.j = ey; row.k = ez; row.loc = DMSTAG_BACK; row.c = 0;

          if (back_boundary) {
            /* Back z-velocity Dirichlet */
            const PetscScalar val_rhs = 0.0;
            const PetscScalar val_A = Kbound;

            ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
            ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
          } else {
            /* Z-momentum equation */
            PetscInt      count;
            DMStagStencil col[17];
            PetscScalar   val_rhs,val_A[17];
            DMStagStencil eta_point[6];
            PetscScalar   eta[6],eta_left,eta_right,eta_up,eta_down,eta_back,eta_front; /* relative to the back face */

            /* Get eta values */
            eta_point[0].i = ex; eta_point[0].j = ey; eta_point[0].k = ez;   eta_point[0].loc = DMSTAG_BACK_LEFT;  eta_point[0].c = 0; /* Left  */
            eta_point[1].i = ex; eta_point[1].j = ey; eta_point[1].k = ez;   eta_point[1].loc = DMSTAG_BACK_RIGHT; eta_point[1].c = 0; /* Right */
            eta_point[2].i = ex; eta_point[2].j = ey; eta_point[2].k = ez;   eta_point[2].loc = DMSTAG_BACK_UP;    eta_point[2].c = 0; /* Up    */
            eta_point[3].i = ex; eta_point[3].j = ey; eta_point[3].k = ez;   eta_point[3].loc = DMSTAG_BACK_DOWN;  eta_point[3].c = 0; /* Down  */
            eta_point[4].i = ex; eta_point[4].j = ey; eta_point[4].k = ez-1; eta_point[4].loc = DMSTAG_ELEMENT;    eta_point[4].c = 0; /* Back  */
            eta_point[5].i = ex; eta_point[5].j = ey; eta_point[5].k = ez;   eta_point[5].loc = DMSTAG_ELEMENT;    eta_point[5].c = 0; /* Front  */
            ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,6,eta_point,eta);CHKERRQ(ierr);
            eta_left = eta[0]; eta_right = eta[1]; eta_up = eta[2]; eta_down = eta[3]; eta_back = eta[4]; eta_front = eta[5];

            count = 0;

            col[count] = row;
            val_A[count] = -2.0 * dv * (eta_back + eta_front) / (hz*hz);
            if (!left_boundary)   val_A[count] += -1.0 * dv * eta_left  / (hx*hx);
            if (!right_boundary)  val_A[count] += -1.0 * dv * eta_right / (hx*hx);
            if (!top_boundary)    val_A[count] += -1.0 * dv * eta_up    / (hy*hy);
            if (!bottom_boundary) val_A[count] += -1.0 * dv * eta_down  / (hy*hy);
            ++count;

            col[count].i = ex; col[count].j = ey; col[count].k = ez-1; col[count].loc = DMSTAG_BACK; col[count].c = 0;
            val_A[count] = 2.0 * dv * eta_back  / (hz*hz); ++count;
            col[count].i = ex; col[count].j = ey; col[count].k = ez+1; col[count].loc = DMSTAG_BACK; col[count].c = 0;
            val_A[count] = 2.0 * dv * eta_front / (hz*hz); ++count;

            if (!left_boundary) {
              col[count].i = ex-1; col[count].j = ey; col[count].k = ez; col[count].loc = DMSTAG_BACK; col[count].c = 0;
              val_A[count] = dv * eta_left / (hx*hx); ++count;
            }
            if (!right_boundary) {
              col[count].i = ex+1; col[count].j = ey; col[count].k = ez; col[count].loc = DMSTAG_BACK; col[count].c = 0;
              val_A[count] = dv * eta_right / (hx*hx); ++count;
            }
            if (!bottom_boundary) {
              col[count].i = ex; col[count].j = ey-1; col[count].k = ez; col[count].loc = DMSTAG_BACK; col[count].c = 0;
              val_A[count] = dv * eta_down / (hy*hy); ++count;
            }
            if (!top_boundary) {
              col[count].i = ex; col[count].j = ey+1; col[count].k = ez; col[count].loc = DMSTAG_BACK; col[count].c = 0;
              val_A[count] = dv * eta_up  / (hy*hy); ++count;
            }

            col[count].i  = ex; col[count].j  = ey; col[count].k = ez-1; col[count].loc = DMSTAG_LEFT; col[count].c = 0;
            val_A[count]  =        dv * eta_left  / (hx*hz); ++count; /* back left*/
            col[count].i  = ex; col[count].j  = ey; col[count].k = ez;   col[count].loc = DMSTAG_LEFT; col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_left  / (hx*hz); ++count; /* front left*/

            col[count].i  = ex; col[count].j  = ey; col[count].k = ez-1; col[count].loc = DMSTAG_RIGHT; col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_right / (hx*hz); ++count; /* back right */
            col[count].i  = ex; col[count].j  = ey; col[count].k = ez;   col[count].loc = DMSTAG_RIGHT; col[count].c = 0;
            val_A[count]  =        dv * eta_right / (hx*hz); ++count; /* front right*/

            col[count].i  = ex; col[count].j  = ey; col[count].k = ez-1; col[count].loc = DMSTAG_DOWN; col[count].c = 0;
            val_A[count]  =        dv * eta_down  / (hy*hz); ++count; /* back down */
            col[count].i  = ex; col[count].j  = ey; col[count].k = ez;   col[count].loc = DMSTAG_DOWN; col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_down  / (hy*hz); ++count; /* back down */

            col[count].i  = ex; col[count].j  = ey; col[count].k = ez-1; col[count].loc = DMSTAG_UP; col[count].c = 0;
            val_A[count]  = -1.0 * dv * eta_up    / (hy*hz); ++count; /* back up */
            col[count].i  = ex; col[count].j  = ey; col[count].k = ez;   col[count].loc = DMSTAG_UP; col[count].c = 0;
            val_A[count]  =        dv * eta_up    / (hy*hz); ++count; /* back up */

            col[count].i = ex; col[count].j = ey; col[count].k = ez-1; col[count].loc = DMSTAG_ELEMENT; col[count].c  = 0;
            val_A[count] = Kcont * dv / hz; ++count;
            col[count].i = ex; col[count].j = ey;   col[count].k = ez; col[count].loc = DMSTAG_ELEMENT; col[count].c  = 0;
            val_A[count] = -1.0 * Kcont * dv / hz; ++count;

            ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,count,col,val_A,INSERT_VALUES);CHKERRQ(ierr);
            val_rhs = 0.0; /* No forcing */
            ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
          }
        }

        /* Elements */
        {
          DMStagStencil row;

          row.i = ex; row.j = ey; row.k = ez; row.loc = DMSTAG_ELEMENT; row.c = 0;

          if (pin_pressure && ex == pinx && ey == piny && ez == pinz) {
            /* Pin a pressure node to zero, if requested */
            const PetscScalar val_A = Kbound;
            const PetscScalar val_rhs = 0.0;

            ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&val_A,INSERT_VALUES);CHKERRQ(ierr);
            ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
          } else {
            /* Continuity equation */
            /* Note that this includes an explicit zero on the diagonal. This is only needed for
               some direct solvers (not required if using an iterative solver and setting a constant-pressure nullspace) */
            /* Note: the scaling by dv is not chosen in a principled way and is likely sub-optimal */
            DMStagStencil     col[7];
            PetscScalar       val_A[7];
            const PetscScalar val_rhs = 0.0;

            col[0].i = ex; col[0].j = ey; col[0].k = ez; col[0].loc = DMSTAG_LEFT;    col[0].c = 0; val_A[0] = - 1.0 * Kcont * dv / hx;
            col[1].i = ex; col[1].j = ey; col[1].k = ez; col[1].loc = DMSTAG_RIGHT;   col[1].c = 0; val_A[1] =         Kcont * dv / hx;
            col[2].i = ex; col[2].j = ey; col[2].k = ez; col[2].loc = DMSTAG_DOWN;    col[2].c = 0; val_A[2] = - 1.0 * Kcont * dv / hy;
            col[3].i = ex; col[3].j = ey; col[3].k = ez; col[3].loc = DMSTAG_UP;      col[3].c = 0; val_A[3] =         Kcont * dv / hy;
            col[4].i = ex; col[4].j = ey; col[4].k = ez; col[4].loc = DMSTAG_BACK;    col[4].c = 0; val_A[4] = - 1.0 * Kcont * dv / hz;
            col[5].i = ex; col[5].j = ey; col[5].k = ez; col[5].loc = DMSTAG_FRONT;   col[5].c = 0; val_A[5] =         Kcont * dv / hz;
            col[6] = row;                                                                           val_A[6] = 0.0;
            ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,7,col,val_A,INSERT_VALUES);CHKERRQ(ierr);
            ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&val_rhs,INSERT_VALUES);CHKERRQ(ierr);
          }
        }
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

  // TODO DEBUG
  {
    PetscBool dump = PETSC_FALSE;

    ierr = PetscOptionsGetBool(NULL,NULL,"-stokes_mat_binary_dump",&dump,NULL);CHKERRQ(ierr);
    if (dump) {
      ierr = MatView(A,PETSC_VIEWER_BINARY_WORLD);CHKERRQ(ierr);
    }
  }
  {
    PetscBool dump = PETSC_FALSE;

    ierr = PetscOptionsGetBool(NULL,NULL,"-stokes_rhs_binary_dump",&dump,NULL);CHKERRQ(ierr);
    if (dump) {
      ierr = VecView(rhs,PETSC_VIEWER_BINARY_WORLD);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}
