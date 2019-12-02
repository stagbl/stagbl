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

// TODO temp as we refactor
static PetscErrorCode CreateSystem_Temp(StagBLStokesParameters parameters,StagBLSystem system);

/**
  * A general function which creates Stokes StagBLSystem objects. It accepts
  * a struct containing all relevant parameters.
  */
PetscErrorCode StagBLCreateStokesSystem(StagBLStokesParameters parameters, StagBLSystem *system)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  // TODO this is temporary, as we refactor
  ierr = StagBLGridCreateStagBLSystem(parameters->stokes_grid,system);CHKERRQ(ierr);
  ierr = CreateSystem_Temp(parameters,*system);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// TODO this is temporary as we refactor

/* Shorter, more convenient names for DMStagLocation entries */
// TODO get rid of this
#define DOWN_LEFT  DMSTAG_DOWN_LEFT
#define DOWN       DMSTAG_DOWN
#define DOWN_RIGHT DMSTAG_DOWN_RIGHT
#define LEFT       DMSTAG_LEFT
#define ELEMENT    DMSTAG_ELEMENT
#define RIGHT      DMSTAG_RIGHT
#define UP_LEFT    DMSTAG_UP_LEFT
#define UP         DMSTAG_UP
#define UP_RIGHT   DMSTAG_UP_RIGHT

static PetscErrorCode CreateSystem_Temp(StagBLStokesParameters parameters,StagBLSystem system)
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
  PetscInt        slotTempDownLeft,slotTempDownRight,slotTempUpLeft,slotTempUpRight;

  PetscFunctionBeginUser;
  ierr = StagBLSystemPETScGetMatPointer(system,&pA);CHKERRQ(ierr);
  ierr = StagBLSystemPETScGetVecPointer(system,&pRhs);CHKERRQ(ierr);
  ierr = StagBLArrayGetStagBLGrid(parameters->coefficient_array,&coefficient_grid);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(parameters->stokes_grid,&dm_stokes);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(coefficient_grid,&dm_coefficient);CHKERRQ(ierr);
  ierr = StagBLArrayPETScGetLocalVec(parameters->coefficient_array,&coeff_local);CHKERRQ(ierr);

  /* Compute some parameters */
  ierr = DMStagGetGlobalSizes(dm_stokes,&N[0],&N[1],NULL);CHKERRQ(ierr);
  if (parameters->uniform_grid) {
    hx = (parameters->xmax-parameters->xmin)/N[0];
    hy = (parameters->ymax-parameters->ymin)/N[1];
    dv = hx*hy;
  } else StagBLError(PetscObjectComm((PetscObject)dm_stokes),"Non-uniform grids not supported yet");
  hxAvgInv = 2.0/(hx + hy);
  Kcont = parameters->eta_characteristic*hxAvgInv;
  Kbound = parameters->eta_characteristic*hxAvgInv*hxAvgInv;
  if (N[0] < 2) SETERRQ(PetscObjectComm((PetscObject)dm_stokes),PETSC_ERR_SUP,"Not implemented for a single element in the x direction");
  pinx = 1; piny = 0;

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
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_DOWN_LEFT,0,&slotTempDownLeft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_DOWN_RIGHT,0,&slotTempDownRight);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_UP_LEFT,0,&slotTempUpLeft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_UP_RIGHT,0,&slotTempUpRight);CHKERRQ(ierr);
  }

  /* Loop over all local elements. Note that it may be more efficient in real
     applications to loop over each boundary separately */
  for (ey = starty; ey<starty+ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
    for (ex = startx; ex<startx+nx; ++ex) {

      if (ey == N[1]-1) {
        /* Top boundary velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = Kbound;
        row.i = ex; row.j = ey; row.loc = UP; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      if (ey == 0) {
        /* Bottom boundary velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = Kbound;
        row.i = ex; row.j = ey; row.loc = DOWN; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        /* Y-momentum equation : (u_xx + u_yy) - p_y = f^y : includes non-zero forcing */
        PetscInt      nEntries;
        DMStagStencil row,col[11];
        PetscScalar   valA[11];
        DMStagStencil rhoPoint[2];
        PetscScalar   rho[2],valRhs;
        DMStagStencil etaPoint[4];
        PetscScalar   eta[4],etaLeft,etaRight,etaUp,etaDown;

        /* get rho values  and compute rhs value*/
        rhoPoint[0].i = ex; rhoPoint[0].j = ey; rhoPoint[0].loc = DOWN_LEFT;  rhoPoint[0].c = 1;
        rhoPoint[1].i = ex; rhoPoint[1].j = ey; rhoPoint[1].loc = DOWN_RIGHT; rhoPoint[1].c = 1;
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,2,rhoPoint,rho);CHKERRQ(ierr);
        valRhs = -parameters->gy * dv * 0.5 * (rho[0] + rho[1]);

        /* get rho values  */
        rhoPoint[0].i = ex; rhoPoint[0].j = ey; rhoPoint[0].loc = DOWN_LEFT;  rhoPoint[0].c = 1;
        rhoPoint[1].i = ex; rhoPoint[1].j = ey; rhoPoint[1].loc = DOWN_RIGHT; rhoPoint[1].c = 1;
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,2,rhoPoint,rho);CHKERRQ(ierr);

        /* Compute forcing */
        {
          const PetscReal rho_avg = 0.5 * (rho[0] + rho[1]);

          /* Note Boussinesq forcing has the opposite sign */
          if (parameters->boussinesq_forcing) {
            valRhs = parameters->alpha * rho_avg * parameters->gy * dv * arr_temperature[ey][ex][slotTempDownLeft];
          } else {
            valRhs = -parameters->gy * dv * rho_avg;
          }
        }

        /* Get eta values */
        etaPoint[0].i = ex; etaPoint[0].j = ey;   etaPoint[0].loc = DOWN_LEFT;  etaPoint[0].c = 0; /* Left  */
        etaPoint[1].i = ex; etaPoint[1].j = ey;   etaPoint[1].loc = DOWN_RIGHT; etaPoint[1].c = 0; /* Right */
        etaPoint[2].i = ex; etaPoint[2].j = ey;   etaPoint[2].loc = ELEMENT;    etaPoint[2].c = 0; /* Up    */
        etaPoint[3].i = ex; etaPoint[3].j = ey-1; etaPoint[3].loc = ELEMENT;    etaPoint[3].c = 0; /* Down  */
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,4,etaPoint,eta);CHKERRQ(ierr);
        etaLeft = eta[0]; etaRight = eta[1]; etaUp = eta[2]; etaDown = eta[3];

        if (ex == 0) {
          /* Left boundary y velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j    = ey  ; row.loc    = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j = ey  ; col[0].loc = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * dv * (etaDown + etaUp) / (hy*hy) - dv * (etaRight) /(hx*hx);
          col[1].i = ex  ; col[1].j = ey-1; col[1].loc = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * dv * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j = ey+1; col[2].loc = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * dv * etaUp    / (hy*hy);
          /* No left entry */
          col[3].i = ex+1; col[3].j = ey  ; col[3].loc = DOWN;     col[3].c  = 0; valA[3]  =        dv * etaRight / (hx*hx);
          col[4].i = ex  ; col[4].j = ey-1; col[4].loc = LEFT;     col[4].c  = 0; valA[4]  =        dv * etaLeft  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j = ey-1; col[5].loc = RIGHT;    col[5].c  = 0; valA[5]  = -1.0 * dv * etaRight / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j = ey  ; col[6].loc = LEFT;     col[6].c  = 0; valA[6]  = -1.0 * dv * etaLeft  / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j = ey  ; col[7].loc = RIGHT;    col[7].c  = 0; valA[7]  =        dv * etaRight / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j = ey-1; col[8].loc = ELEMENT;  col[8].c  = 0; valA[8]  =        Kcont * dv / hy;
          col[9].i = ex  ; col[9].j = ey  ; col[9].loc = ELEMENT;  col[9].c  = 0; valA[9]  = -1.0 * Kcont * dv / hy;
        } else if (ex == N[0]-1) {
          /* Right boundary y velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j    = ey  ; row.loc    = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j = ey  ; col[0].loc = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * dv * (etaDown + etaUp) / (hy*hy) - dv * (etaLeft) /(hx*hx);
          col[1].i = ex  ; col[1].j = ey-1; col[1].loc = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * dv * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j = ey+1; col[2].loc = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * dv * etaUp    / (hy*hy);
          col[3].i = ex-1; col[3].j = ey  ; col[3].loc = DOWN;     col[3].c  = 0; valA[3]  =        dv * etaLeft  / (hx*hx);
          /* No right element */
          col[4].i = ex  ; col[4].j = ey-1; col[4].loc = LEFT;     col[4].c  = 0; valA[4]  =        dv * etaLeft  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j = ey-1; col[5].loc = RIGHT;    col[5].c  = 0; valA[5]  = -1.0 * dv * etaRight / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j = ey  ; col[6].loc = LEFT;     col[6].c  = 0; valA[7]  = -1.0 * dv * etaLeft  / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j = ey  ; col[7].loc = RIGHT;    col[7].c  = 0; valA[7]  =        dv * etaRight / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j = ey-1; col[8].loc = ELEMENT;  col[8].c  = 0; valA[8]  =        Kcont * dv / hy;
          col[9].i = ex  ; col[9].j = ey  ; col[9].loc = ELEMENT;  col[9].c  = 0; valA[9]  = -1.0 * Kcont * dv / hy;
        } else {
          /* U_y interior equation */
          nEntries = 11;
          row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * dv * (etaDown + etaUp) / (hy*hy) - dv * (etaLeft + etaRight) /(hx*hx);
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * dv * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * dv * etaUp    / (hy*hy);
          col[3].i = ex-1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0; valA[3]  =        dv * etaLeft  / (hx*hx);
          col[4].i = ex+1; col[4].j  = ey  ; col[4].loc  = DOWN;     col[4].c  = 0; valA[4]  =        dv * etaRight / (hx*hx);
          col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = LEFT;     col[5].c  = 0; valA[5]  =        dv * etaLeft  / (hx*hy); /* down left x edge */
          col[6].i = ex  ; col[6].j  = ey-1; col[6].loc  = RIGHT;    col[6].c  = 0; valA[6]  = -1.0 * dv * etaRight / (hx*hy); /* down right x edge */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = LEFT;     col[7].c  = 0; valA[7]  = -1.0 * dv * etaLeft  / (hx*hy); /* up left x edge */
          col[8].i = ex  ; col[8].j  = ey  ; col[8].loc  = RIGHT;    col[8].c  = 0; valA[8]  =        dv * etaRight / (hx*hy); /* up right x edge */
          col[9].i = ex  ; col[9].j  = ey-1; col[9].loc  = ELEMENT;  col[9].c  = 0; valA[9]  =        Kcont * dv / hy;
          col[10].i = ex ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0; valA[10] = -1.0 * Kcont * dv / hy;
        }

        /* Insert Y-momentum entries */
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      if (ex == N[0]-1) {
        /* Right Boundary velocity Dirichlet */
        /* Redundant in the corner */
        DMStagStencil row;
        PetscScalar   valRhs;

        const PetscScalar valA = Kbound;
        row.i = ex; row.j = ey; row.loc = RIGHT; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }
      if (ex == 0) {
        /* Left velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = Kbound;
        row.i = ex; row.j = ey; row.loc = LEFT; row.c = 0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
        PetscInt nEntries;
        DMStagStencil row,col[11];
        PetscScalar   valRhs,valA[11];
        DMStagStencil etaPoint[4];
        PetscScalar eta[4],etaLeft,etaRight,etaUp,etaDown;

        /* Get eta values */
        etaPoint[0].i = ex-1; etaPoint[0].j = ey; etaPoint[0].loc = ELEMENT;   etaPoint[0].c = 0; /* Left  */
        etaPoint[1].i = ex;   etaPoint[1].j = ey; etaPoint[1].loc = ELEMENT;   etaPoint[1].c = 0; /* Right */
        etaPoint[2].i = ex;   etaPoint[2].j = ey; etaPoint[2].loc = UP_LEFT;   etaPoint[2].c = 0; /* Up    */
        etaPoint[3].i = ex;   etaPoint[3].j = ey; etaPoint[3].loc = DOWN_LEFT; etaPoint[3].c = 0; /* Down  */
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeff_local,4,etaPoint,eta);CHKERRQ(ierr);
        etaLeft = eta[0]; etaRight = eta[1]; etaUp = eta[2]; etaDown = eta[3];

        if (ey == 0) {
          /* Bottom boundary x velocity stencil (with zero vel deriv) */
          nEntries = 10;
          row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0; valA[0]  = -2.0 * dv * (etaLeft + etaRight) / (hx*hx) - dv * (etaUp) / (hy*hy);
          /* Missing element below */
          col[1].i = ex  ; col[1].j  = ey+1; col[1].loc  = LEFT;    col[1].c   = 0; valA[1]  =        dv * etaUp    / (hy*hy);
          col[2].i = ex-1; col[2].j  = ey  ; col[2].loc  = LEFT;    col[2].c   = 0; valA[2]  =  2.0 * dv * etaLeft  / (hx*hx);
          col[3].i = ex+1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0; valA[3]  =  2.0 * dv * etaRight / (hx*hx);
          col[4].i = ex-1; col[4].j  = ey  ; col[4].loc  = DOWN;    col[4].c   = 0; valA[4]  =        dv * etaDown  / (hx*hy); /* down left */
          col[5].i = ex  ; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0; valA[5]  = -1.0 * dv * etaDown  / (hx*hy); /* down right */
          col[6].i = ex-1; col[6].j  = ey  ; col[6].loc  = UP;      col[6].c   = 0; valA[6]  = -1.0 * dv * etaUp    / (hx*hy); /* up left */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  =        dv * etaUp    / (hx*hy); /* up right */
          col[8].i = ex-1; col[8].j  = ey  ; col[8].loc  = ELEMENT; col[8].c   = 0; valA[8]  =        Kcont * dv / hx;
          col[9].i = ex  ; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0; valA[9]  = -1.0 * Kcont * dv / hx;
          valRhs = 0.0;
        } else if (ey == N[1]-1) {
          /* Top boundary x velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c  = 0; valA[0]  = -2.0 * dv * (etaLeft + etaRight) / (hx*hx) - dv * (etaDown) / (hy*hy);
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = LEFT;    col[1].c  = 0; valA[1]  =        dv * etaDown  / (hy*hy);
          /* Missing element above */
          col[2].i = ex-1; col[2].j  = ey  ; col[2].loc  = LEFT;    col[2].c  = 0; valA[2]  =  2.0 * dv * etaLeft  / (hx*hx);
          col[3].i = ex+1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c  = 0; valA[3]  =  2.0 * dv * etaRight / (hx*hx);
          col[4].i = ex-1; col[4].j  = ey  ; col[4].loc  = DOWN;    col[4].c  = 0; valA[4]  =        dv * etaDown  / (hx*hy); /* down left */
          col[5].i = ex  ; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c  = 0; valA[5]  = -1.0 * dv * etaDown  / (hx*hy); /* down right */
          col[6].i = ex-1; col[6].j  = ey  ; col[6].loc  = UP;      col[6].c  = 0; valA[6]  = -1.0 * dv * etaUp    / (hx*hy); /* up left */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c  = 0; valA[7]  =        dv * etaUp    / (hx*hy); /* up right */
          col[8].i = ex-1; col[8].j  = ey  ; col[8].loc  = ELEMENT; col[8].c  = 0; valA[8]  =        Kcont * dv / hx;
          col[9].i = ex  ; col[9].j  = ey  ; col[9].loc = ELEMENT;  col[9].c  = 0; valA[9]  = -1.0 * Kcont * dv / hx;
          valRhs = 0.0;
        } else {
          /* U_x interior equation */
          nEntries = 11;
          row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0; valA[0]  = -2.0 * dv * (etaLeft + etaRight) / (hx*hx) - dv * (etaUp + etaDown) / (hy*hy);
          col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc  = LEFT;    col[1].c   = 0; valA[1]  =        dv * etaDown  / (hy*hy);
          col[2].i  = ex  ; col[2].j  = ey+1; col[2].loc  = LEFT;    col[2].c   = 0; valA[2]  =        dv * etaUp    / (hy*hy);
          col[3].i  = ex-1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0; valA[3]  =  2.0 * dv * etaLeft  / (hx*hx);
          col[4].i  = ex+1; col[4].j  = ey  ; col[4].loc  = LEFT;    col[4].c   = 0; valA[4]  =  2.0 * dv * etaRight / (hx*hx);
          col[5].i  = ex-1; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0; valA[5]  =        dv * etaDown  / (hx*hy); /* down left */
          col[6].i  = ex  ; col[6].j  = ey  ; col[6].loc  = DOWN;    col[6].c   = 0; valA[6]  = -1.0 * dv * etaDown  / (hx*hy); /* down right */
          col[7].i  = ex-1; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  = -1.0 * dv * etaUp    / (hx*hy); /* up left */
          col[8].i  = ex  ; col[8].j  = ey  ; col[8].loc  = UP;      col[8].c   = 0; valA[8]  =        dv * etaUp    / (hx*hy); /* up right */
          col[9].i  = ex-1; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0; valA[9]  =        Kcont * dv / hx;
          col[10].i = ex  ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0; valA[10] = -1.0 * Kcont * dv / hx;
          valRhs = 0.0;
        }
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      /* P equation : u_x + v_y = 0
         Note that this includes an explicit zero on the diagonal. This is only needed for
         direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace)

         Note: the scaling by dv is not chosen in a principled way and is likely sub-optimal */
      if (pin_pressure && ex == pinx && ey == piny) { /* Pin a pressure node to zero, if requested */
        DMStagStencil row;
        PetscScalar valA,valRhs;
        row.i = ex; row.j = ey; row.loc = ELEMENT; row.c = 0;
        valA = Kbound;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        DMStagStencil row,col[5];
        PetscScalar   valA[5],valRhs;

        row.i    = ex; row.j    = ey; row.loc    = ELEMENT; row.c    = 0;
        col[0].i = ex; col[0].j = ey; col[0].loc = LEFT;    col[0].c = 0; valA[0] = - 1.0 * Kcont * dv / hx;
        col[1].i = ex; col[1].j = ey; col[1].loc = RIGHT;   col[1].c = 0; valA[1] =         Kcont * dv / hx;
        col[2].i = ex; col[2].j = ey; col[2].loc = DOWN;    col[2].c = 0; valA[2] = - 1.0 * Kcont * dv / hy;
        col[3].i = ex; col[3].j = ey; col[3].loc = UP;      col[3].c = 0; valA[3] =         Kcont * dv / hy;
        col[4] = row;                                                     valA[4] = 0.0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,5,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
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
