#include "stagbl.h"

PetscErrorCode StagBLGridCreateStokes2DBox(MPI_Comm comm, PetscInt nx,PetscInt ny,PetscScalar xmin, PetscScalar xmax, PetscScalar ymin, PetscScalar ymax,StagBLGrid *pgrid)
{
  // TODO this function assumes PETSc is included, and that the defaults types for things are PETSc
  PetscErrorCode ierr;
  DM *pdm;
  DM dm_stokes;

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
  return 0;
}

PetscErrorCode StagBLCreateSimpleStokesSystem(StagBLGrid stokes_grid, StagBLArray coefficient_array, StagBLSystem *system)
{
  PetscErrorCode         ierr;
  StagBLStokesParameters parameters;

  PetscFunctionBegin;
  /* Check that the grids are compatible and have the correct numbers of dof */
  // TODO

  /* Pass to the general function */
  parameters.coefficient_array = coefficient_array;
  parameters.stokes_grid = stokes_grid;
  // TODO this is all temp - this needs to be computed from the grid and coordinate arrays
  parameters.xmin = 0; // TODO temp, rather just have coordinate arrays for now (later could optimize)
  parameters.xmax = 1e6;
  parameters.ymin = 0.0;
  parameters.ymax = 1.5e6;
  parameters.gy   = 10.0;
    parameters.eta1 = 1e20; // TODO should just have eta_characteristic field which you can supply, otherwise supply from coefficient_array
    parameters.eta2 = 1e22;

  ierr = StagBLCreateStokesSystem(stokes_grid,parameters,system);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// TODO temp as we refactor
static PetscErrorCode CreateSystem_Temp(StagBLStokesParameters *parameters,StagBLSystem system);

/**
  * A general function which creates Stokes StagBLSystem objects. It accepts
  * a struct containing all relevant parameters.
  */
PetscErrorCode StagBLCreateStokesSystem(StagBLGrid stokes_grid, StagBLStokesParameters parameters, StagBLSystem *system)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  // TODO this is temporary, as we refactor
  ierr = StagBLGridCreateStagBLSystem(stokes_grid,system);CHKERRQ(ierr);
  ierr = CreateSystem_Temp(&parameters,*system);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// TODO this is temporary as we refactor

/* Shorter, more convenient names for DMStagLocation entries */
#define DOWN_LEFT  DMSTAG_DOWN_LEFT
#define DOWN       DMSTAG_DOWN
#define DOWN_RIGHT DMSTAG_DOWN_RIGHT
#define LEFT       DMSTAG_LEFT
#define ELEMENT    DMSTAG_ELEMENT
#define RIGHT      DMSTAG_RIGHT
#define UP_LEFT    DMSTAG_UP_LEFT
#define UP         DMSTAG_UP
#define UP_RIGHT   DMSTAG_UP_RIGHT

static PetscErrorCode CreateSystem_Temp(StagBLStokesParameters *parameters,StagBLSystem system)
{
  PetscErrorCode  ierr;
  DM              dm_stokes,dm_coefficient;
  PetscInt        N[2];
  PetscInt        ex,ey,startx,starty,nx,ny;
  Mat             *pA;
  Vec             *pRhs;
  Mat             A;
  Vec             rhs;
  PetscReal       hx,hy,hxAvgInv,Kcont,Kbound,etaCharacteristic;
  PetscInt        pinx,piny;
  const PetscBool pinPressure = PETSC_TRUE;
  Vec             coeffLocal;
  StagBLGrid      coefficient_grid;

  PetscFunctionBeginUser;

  /* Use the "escape hatch" */
  ierr = StagBLSystemPETScGetMatPointer(system,&pA);CHKERRQ(ierr);
  ierr = StagBLSystemPETScGetVecPointer(system,&pRhs);CHKERRQ(ierr);
  ierr = StagBLArrayGetStagBLGrid(parameters->coefficient_array,&coefficient_grid);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(parameters->stokes_grid,&dm_stokes);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(coefficient_grid,&dm_coefficient);CHKERRQ(ierr);
  ierr = StagBLArrayPETScGetLocalVec(parameters->coefficient_array,&coeffLocal);CHKERRQ(ierr);

  /* Compute some parameters */
  ierr = DMStagGetGlobalSizes(dm_stokes,&N[0],&N[1],NULL);CHKERRQ(ierr);
  hx = (parameters->xmax-parameters->xmin)/N[0];
  hy = (parameters->ymax-parameters->ymin)/N[1];
  etaCharacteristic = PetscMin(parameters->eta1,parameters->eta2);
  hxAvgInv = 2.0/(hx + hy);
  Kcont = etaCharacteristic*hxAvgInv;
  Kbound = etaCharacteristic*hxAvgInv*hxAvgInv;
  if (N[0] < 2) SETERRQ(PetscObjectComm((PetscObject)dm_stokes),PETSC_ERR_SUP,"Not implemented for a single element in the x direction");
  pinx = 1; piny = 0;

  ierr = DMCreateMatrix(dm_stokes,pA);CHKERRQ(ierr);
  A = *pA;
  ierr = DMCreateGlobalVector(dm_stokes,pRhs);CHKERRQ(ierr);
  rhs = *pRhs;
  ierr = DMStagGetCorners(dm_stokes,&startx,&starty,NULL,&nx,&ny,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(dm_stokes,&N[0],&N[1],NULL);CHKERRQ(ierr);

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
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeffLocal,2,rhoPoint,rho);CHKERRQ(ierr);
        valRhs = -parameters->gy * 0.5 * (rho[0] + rho[1]);

        /* Get eta values */
        etaPoint[0].i = ex; etaPoint[0].j = ey;   etaPoint[0].loc = DOWN_LEFT;  etaPoint[0].c = 0; /* Left  */
        etaPoint[1].i = ex; etaPoint[1].j = ey;   etaPoint[1].loc = DOWN_RIGHT; etaPoint[1].c = 0; /* Right */
        etaPoint[2].i = ex; etaPoint[2].j = ey;   etaPoint[2].loc = ELEMENT;    etaPoint[2].c = 0; /* Up    */
        etaPoint[3].i = ex; etaPoint[3].j = ey-1; etaPoint[3].loc = ELEMENT;    etaPoint[3].c = 0; /* Down  */
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeffLocal,4,etaPoint,eta);CHKERRQ(ierr);
        etaLeft = eta[0]; etaRight = eta[1]; etaUp = eta[2]; etaDown = eta[3];

        if (ex == 0) {
          /* Left boundary y velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j    = ey  ; row.loc    = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j = ey  ; col[0].loc = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * (etaDown + etaUp) / (hy*hy) - (etaRight) /(hx*hx);
          col[1].i = ex  ; col[1].j = ey-1; col[1].loc = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j = ey+1; col[2].loc = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * etaUp    / (hy*hy);
          /* No left entry */
          col[3].i = ex+1; col[3].j = ey  ; col[3].loc = DOWN;     col[3].c  = 0; valA[3]  =        etaRight / (hx*hx);
          col[4].i = ex  ; col[4].j = ey-1; col[4].loc = LEFT;     col[4].c  = 0; valA[4]  =        etaLeft  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j = ey-1; col[5].loc = RIGHT;    col[5].c  = 0; valA[5]  = -      etaRight / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j = ey  ; col[6].loc = LEFT;     col[6].c  = 0; valA[6]  = -      etaLeft  / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j = ey  ; col[7].loc = RIGHT;    col[7].c  = 0; valA[7]  =        etaRight / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j = ey-1; col[8].loc = ELEMENT;  col[8].c  = 0; valA[8]  =  Kcont / hy;
          col[9].i = ex  ; col[9].j = ey  ; col[9].loc = ELEMENT;  col[9].c  = 0; valA[9]  = -Kcont / hy;
        } else if (ex == N[0]-1) {
          /* Right boundary y velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j    = ey  ; row.loc    = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j = ey  ; col[0].loc = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * (etaDown + etaUp) / (hy*hy) - (etaLeft) /(hx*hx );
          col[1].i = ex  ; col[1].j = ey-1; col[1].loc = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j = ey+1; col[2].loc = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * etaUp    / (hy*hy);
          col[3].i = ex-1; col[3].j = ey  ; col[3].loc = DOWN;     col[3].c  = 0; valA[3]  =        etaLeft  / (hx*hx);
          /* No right element */
          col[4].i = ex  ; col[4].j = ey-1; col[4].loc = LEFT;     col[4].c  = 0; valA[4]  =        etaLeft  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j = ey-1; col[5].loc = RIGHT;    col[5].c  = 0; valA[5]  = -      etaRight / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j = ey  ; col[6].loc = LEFT;     col[6].c  = 0; valA[7]  = -      etaLeft  / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j = ey  ; col[7].loc = RIGHT;    col[7].c  = 0; valA[7]  =        etaRight / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j = ey-1; col[8].loc = ELEMENT;  col[8].c  = 0; valA[8]  =  Kcont / hy;
          col[9].i = ex  ; col[9].j = ey  ; col[9].loc = ELEMENT;  col[9].c  = 0; valA[9]  = -Kcont / hy;
        } else {
          /* U_y interior equation */
          nEntries = 11;
          row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * (etaDown + etaUp) / (hy*hy) - (etaLeft + etaRight) /(hx*hx);
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * etaUp    / (hy*hy);
          col[3].i = ex-1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0; valA[3]  =        etaLeft  / (hx*hx);
          col[4].i = ex+1; col[4].j  = ey  ; col[4].loc  = DOWN;     col[4].c  = 0; valA[4]  =        etaRight / (hx*hx);
          col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = LEFT;     col[5].c  = 0; valA[5]  =        etaLeft  / (hx*hy); /* down left x edge */
          col[6].i = ex  ; col[6].j  = ey-1; col[6].loc  = RIGHT;    col[6].c  = 0; valA[6]  = -      etaRight / (hx*hy); /* down right x edge */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = LEFT;     col[7].c  = 0; valA[7]  = -      etaLeft  / (hx*hy); /* up left x edge */
          col[8].i = ex  ; col[8].j  = ey  ; col[8].loc  = RIGHT;    col[8].c  = 0; valA[8]  =        etaRight / (hx*hy); /* up right x edge */
          col[9].i = ex  ; col[9].j  = ey-1; col[9].loc  = ELEMENT;  col[9].c  = 0; valA[9]  =  Kcont / hy;
          col[10].i = ex ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0; valA[10] = -Kcont / hy;
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
        ierr = DMStagVecGetValuesStencil(dm_coefficient,coeffLocal,4,etaPoint,eta);CHKERRQ(ierr);
        etaLeft = eta[0]; etaRight = eta[1]; etaUp = eta[2]; etaDown = eta[3];

        if (ey == 0) {
          /* Bottom boundary x velocity stencil (with zero vel deriv) */
          nEntries = 10;
          row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0; valA[0]  = -2.0 * (etaLeft + etaRight) / (hx*hx) -(etaUp) / (hy*hy);
          /* Missing element below */
          col[1].i  = ex  ; col[1].j  = ey+1; col[1].loc  = LEFT;    col[1].c   = 0; valA[1]  =        etaUp    / (hy*hy);
          col[2].i  = ex-1; col[2].j  = ey  ; col[2].loc  = LEFT;    col[2].c   = 0; valA[2]  =  2.0 * etaLeft  / (hx*hx);
          col[3].i  = ex+1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0; valA[3]  =  2.0 * etaRight / (hx*hx);
          col[4].i  = ex-1; col[4].j  = ey  ; col[4].loc  = DOWN;    col[4].c   = 0; valA[4]  =        etaDown  / (hx*hy); /* down left */
          col[5].i  = ex  ; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0; valA[5]  = -      etaDown  / (hx*hy); /* down right */
          col[6].i  = ex-1; col[6].j  = ey  ; col[6].loc  = UP;      col[6].c   = 0; valA[6]  = -      etaUp    / (hx*hy); /* up left */
          col[7].i  = ex  ; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  =        etaUp    / (hx*hy); /* up right */
          col[8].i  = ex-1; col[8].j  = ey  ; col[8].loc  = ELEMENT; col[8].c   = 0; valA[8]  =  Kcont / hx;
          col[9].i = ex   ; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0; valA[9]  = -Kcont / hx;
          valRhs = 0.0;
        } else if (ey == N[1]-1) {
          /* Top boundary x velocity stencil */
          nEntries = 10;
          row.i     = ex  ; row.j     = ey  ; row.loc    = LEFT;    row.c     = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0]  = -2.0 * (etaLeft + etaRight) / (hx*hx) -(etaDown) / (hy*hy);
          col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc = LEFT;    col[1].c  = 0; valA[1]  =        etaDown  / (hy*hy);
          /* Missing element above */
          col[2].i  = ex-1; col[2].j  = ey  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2]  =  2.0 * etaLeft  / (hx*hx);
          col[3].i  = ex+1; col[3].j  = ey  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3]  =  2.0 * etaRight / (hx*hx);
          col[4].i  = ex-1; col[4].j  = ey  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4]  =        etaDown  / (hx*hy); /* down left */
          col[5].i  = ex  ; col[5].j  = ey  ; col[5].loc = DOWN;    col[5].c  = 0; valA[5]  = -      etaDown  / (hx*hy); /* down right */
          col[6].i  = ex-1; col[6].j  = ey  ; col[6].loc = UP;      col[6].c  = 0; valA[6]  = -      etaUp    / (hx*hy); /* up left */
          col[7].i  = ex  ; col[7].j  = ey  ; col[7].loc = UP;      col[7].c  = 0; valA[7]  =        etaUp    / (hx*hy); /* up right */
          col[8].i  = ex-1; col[8].j  = ey  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8]  =  Kcont / hx;
          col[9].i  = ex  ; col[9].j  = ey  ; col[9].loc = ELEMENT; col[9].c  = 0; valA[9]  = -Kcont / hx;
          valRhs = 0.0;
        } else {
          /* U_x interior equation */
          nEntries = 11;
          row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0; valA[0]  = -2.0 * (etaLeft + etaRight) / (hx*hx) -(etaUp + etaDown) / (hy*hy);
          col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc  = LEFT;    col[1].c   = 0; valA[1]  =        etaDown  / (hy*hy);
          col[2].i  = ex  ; col[2].j  = ey+1; col[2].loc  = LEFT;    col[2].c   = 0; valA[2]  =        etaUp    / (hy*hy);
          col[3].i  = ex-1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0; valA[3]  =  2.0 * etaLeft  / (hx*hx);
          col[4].i  = ex+1; col[4].j  = ey  ; col[4].loc  = LEFT;    col[4].c   = 0; valA[4]  =  2.0 * etaRight / (hx*hx);
          col[5].i  = ex-1; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0; valA[5]  =        etaDown  / (hx*hy); /* down left */
          col[6].i  = ex  ; col[6].j  = ey  ; col[6].loc  = DOWN;    col[6].c   = 0; valA[6]  = -      etaDown  / (hx*hy); /* down right */
          col[7].i  = ex-1; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  = -      etaUp    / (hx*hy); /* up left */
          col[8].i  = ex  ; col[8].j  = ey  ; col[8].loc  = UP;      col[8].c   = 0; valA[8]  =        etaUp    / (hx*hy); /* up right */
          col[9].i  = ex-1; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0; valA[9]  =  Kcont / hx;
          col[10].i = ex  ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0; valA[10] = -Kcont / hx;
          valRhs = 0.0;
        }
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      /* P equation : u_x + v_y = 0
         Note that this includes an explicit zero on the diagonal. This is only needed for
         direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace) */
      if (pinPressure && ex == pinx && ey == piny) { /* Pin a pressure node to zero, if requested */
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
        col[0].i = ex; col[0].j = ey; col[0].loc = LEFT;    col[0].c = 0; valA[0] = -Kcont / hx;
        col[1].i = ex; col[1].j = ey; col[1].loc = RIGHT;   col[1].c = 0; valA[1] =  Kcont / hx;
        col[2].i = ex; col[2].j = ey; col[2].loc = DOWN;    col[2].c = 0; valA[2] = -Kcont / hy;
        col[3].i = ex; col[3].j = ey; col[3].loc = UP;      col[3].c = 0; valA[3] =  Kcont / hy;
        col[4] = row;                                                     valA[4] = 0.0;
        ierr = DMStagMatSetValuesStencil(dm_stokes,A,1,&row,5,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(dm_stokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rhs);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
