#include "stagbl.h"
#include <stdio.h>
#include <petsc.h> // Note that we still have work to do guarding for the non-PETSc case (probably define STAGBL_HAVE_PETSC in a configured include file eventually)

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

/* An application context */
typedef struct {
  MPI_Comm    comm;
  DM          dmStokes,dmCoeff;
  Vec         coeff;
  PetscReal   xmax,ymax,xmin,ymin,hxCharacteristic,hyCharacteristic;
  PetscScalar eta1,eta2,rho1,rho2,gy,Kbound,Kcont,etaCharacteristic;
  PetscBool   hdfOutput; /* as opposed to binary */
} CtxData;
typedef CtxData* Ctx;

/* Helper functions */
static PetscErrorCode PopulateCoefficientData(Ctx);
static PetscErrorCode CreateSystem(const Ctx,Mat*,Vec*);
static PetscErrorCode DumpSolution(Ctx,Vec);

/* Coefficient/forcing Functions */
static PetscReal getRho(Ctx ctx,PetscReal x, PetscReal y) {
  const PetscReal d = ctx->xmax-ctx->xmin;
  const PetscReal xx = x/d - 0.5;
  const PetscReal yy = y/d - 0.5;
  return (xx*xx + yy*yy) > 0.3*0.3 ? ctx->rho1 : ctx->rho2;
}
static PetscReal getEta(Ctx ctx,PetscReal x, PetscReal y) {
  const PetscReal d = ctx->xmax-ctx->xmin;
  const PetscReal xx = x/d - 0.5;
  const PetscReal yy = y/d - 0.5;
  return (xx*xx + yy*yy) > 0.3*0.3 ? ctx->eta1 : ctx->eta2;
}

int main(int argc, char** argv)
{
  int                rank,size;
  StagBLGrid         grid;
  StagBLArray        x,b;
  StagBLOperator     A;
  StagBLLinearSolver solver;
  MPI_Comm           comm;

  Ctx            ctx;

  PetscErrorCode ierr;
  DM             *pdm;
  Vec            *pvecx,*pvecb;
  Vec            vecx,vecb;
  Mat            *pmatA;
  Mat            matA;
  KSP            *pksp;
  KSP            ksp;

  /* Initialize MPI and print a message */
  MPI_Init(&argc,&argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  if (rank == 0) {
    printf("=== StagBLDemo2d ===\n");
    printf("%d ranks\n",size);
  }
  MPI_Barrier(comm);

  /* Initialize StagBL (which will initialize PETSc if needbe) */
  StagBLInitialize(argc,argv,comm);

  /* Populate application context */
  ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);
  ctx->comm = PETSC_COMM_WORLD;
  ctx->xmin = 0.0;
  ctx->xmax = 1e6;
  ctx->ymin = 0.0;
  ctx->ymax = 1.5e6;
  ctx->rho1 = 3200;
  ctx->rho2 = 3300;
  ctx->eta1 = 1e20;
  ctx->eta2 = 1e22;
  ctx->gy    = 10.0;

  /* Create a Grid */
  StagBLGridCreate(&grid);
  StagBLGridPETScGetDMPointer(grid,&pdm);
  ierr = DMStagCreate2d(
      comm,
      DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
      30,20,                                   /* Global element counts */
      PETSC_DECIDE,PETSC_DECIDE,               /* Determine parallel decomposition automatically */
      0,1,1,                                   /* dof: 0 per vertex, 1 per edge, 1 per face/element */
      DMSTAG_GHOST_STENCIL_BOX,
      1,                                       /* elementwise stencil width */
      NULL,NULL,
      pdm);
  ctx->dmStokes = *pdm;
  ierr = DMSetFromOptions(ctx->dmStokes);CHKERRQ(ierr);
  ierr = DMSetUp(ctx->dmStokes);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesProduct(ctx->dmStokes,0.0,ctx->xmax,0.0,ctx->ymax,0.0,0.0);CHKERRQ(ierr);
  ierr = DMStagCreateCompatibleDMStag(ctx->dmStokes,2,0,1,0,&ctx->dmCoeff);CHKERRQ(ierr);
  ierr = DMSetUp(ctx->dmCoeff);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesProduct(ctx->dmCoeff,0.0,ctx->xmax,0.0,ctx->ymax,0.0,0.0);CHKERRQ(ierr);

  ierr = DMView(ctx->dmStokes,0);

  /* Get scaling constants, knowing grid spacing */
  {
    PetscInt N[2];
    PetscReal hxAvgInv;
    ierr = DMStagGetGlobalSizes(ctx->dmStokes,&N[0],&N[1],NULL);CHKERRQ(ierr);
    ctx->hxCharacteristic = (ctx->xmax-ctx->xmin)/N[0];
    ctx->hyCharacteristic = (ctx->ymax-ctx->ymin)/N[1];
    ctx->etaCharacteristic = PetscMin(ctx->eta1,ctx->eta2);
    hxAvgInv = 2.0/(ctx->hxCharacteristic + ctx->hyCharacteristic);
    ctx->Kcont = ctx->etaCharacteristic*hxAvgInv;
    ctx->Kbound = ctx->etaCharacteristic*hxAvgInv*hxAvgInv;
  }

  /* Populate coefficient data */
  ierr = PopulateCoefficientData(ctx);CHKERRQ(ierr);

  /* Create a system */
  StagBLOperatorCreate(&A);
  StagBLArrayCreate(&x);
  StagBLArrayCreate(&b);

  StagBLArrayPETScGetVecPointer(x,&pvecx);
  ierr = DMCreateGlobalVector(ctx->dmStokes,pvecx);
  vecx = *pvecx;

  StagBLArrayPETScGetVecPointer(b,&pvecb);
  StagBLOperatorPETScGetMatPointer(A,&pmatA);
  ierr = CreateSystem(ctx,pmatA,pvecb);CHKERRQ(ierr);
  matA = *pmatA;
  vecb = *pvecb;

  /* Solve the system (you will likely want to specify a solver from the command line) */
  StagBLLinearSolverCreate(&solver);
  StagBLLinearSolverPETScGetKSPPointer(solver,&pksp);

  ierr = VecDuplicate(vecb,pvecx);CHKERRQ(ierr);
  ierr = KSPCreate(ctx->comm,pksp);CHKERRQ(ierr);
  ksp = *pksp;
  ierr = KSPSetOperators(ksp,matA,matA);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,vecb,vecx);CHKERRQ(ierr);
  {
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
    if (reason < 0) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,"Linear solve failed");CHKERRQ(ierr);
  }

  /* Dump solution by converting to DMDAs and dumping */
  ierr = DumpSolution(ctx,vecx);CHKERRQ(ierr);

  /* Free data */
  StagBLArrayDestroy(&x);
  StagBLArrayDestroy(&b);
  StagBLOperatorDestroy(&A);
  StagBLLinearSolverDestroy(&solver);
  StagBLGridDestroy(&grid);

  StagBLFinalize();

  return 0;
}

static PetscErrorCode CreateSystem(const Ctx ctx,Mat *pA,Vec *pRhs)
{
  PetscErrorCode ierr;
  PetscInt       N[2];
  PetscInt       ex,ey,startx,starty,nx,ny;
  Mat            A;
  Vec            rhs;
  PetscReal      hx,hy;
  const PetscBool pinPressure = PETSC_TRUE;
  Vec            coeffLocal;

  PetscFunctionBeginUser;
  ierr = DMCreateMatrix(ctx->dmStokes,pA);CHKERRQ(ierr);
  A = *pA;
  ierr = DMCreateGlobalVector(ctx->dmStokes,pRhs);CHKERRQ(ierr);
  rhs = *pRhs;
  ierr = DMStagGetCorners(ctx->dmStokes,&startx,&starty,NULL,&nx,&ny,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(ctx->dmStokes,&N[0],&N[1],NULL);CHKERRQ(ierr);
  hx = ctx->hxCharacteristic;
  hy = ctx->hyCharacteristic;
  ierr = DMGetLocalVector(ctx->dmCoeff,&coeffLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->dmCoeff,ctx->coeff,INSERT_VALUES,coeffLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->dmCoeff,ctx->coeff,INSERT_VALUES,coeffLocal);CHKERRQ(ierr);

  /* Loop over all local elements. Note that it may be more efficient in real
     applications to loop over each boundary separately */
  for (ey = starty; ey<starty+ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
    for (ex = startx; ex<startx+nx; ++ex) {

      if (ey == N[1]-1) {
        /* Top boundary velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = ctx->Kbound;
        row.i = ex; row.j = ey; row.loc = UP; row.c = 0;
        ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      if (ey == 0) {
        /* Bottom boundary velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = ctx->Kbound;
        row.i = ex; row.j = ey; row.loc = DOWN; row.c = 0;
        ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
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
        ierr = DMStagVecGetValuesStencil(ctx->dmCoeff,coeffLocal,2,rhoPoint,rho);CHKERRQ(ierr);
        valRhs = -ctx->gy * 0.5 * (rho[0] + rho[1]);

        /* Get eta values */
        etaPoint[0].i = ex; etaPoint[0].j = ey;   etaPoint[0].loc = DOWN_LEFT;  etaPoint[0].c = 0; /* Left  */
        etaPoint[1].i = ex; etaPoint[1].j = ey;   etaPoint[1].loc = DOWN_RIGHT; etaPoint[1].c = 0; /* Right */
        etaPoint[2].i = ex; etaPoint[2].j = ey+1; etaPoint[2].loc = ELEMENT;    etaPoint[2].c = 0; /* Up    */
        etaPoint[3].i = ex; etaPoint[3].j = ey-1; etaPoint[3].loc = ELEMENT;    etaPoint[3].c = 0; /* Down  */
        ierr = DMStagVecGetValuesStencil(ctx->dmCoeff,coeffLocal,4,etaPoint,eta);CHKERRQ(ierr);
        etaLeft = eta[0]; etaRight = eta[1]; etaUp = eta[2]; etaDown = eta[3];

        if (ex == 0) {
          /* Left boundary y velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * (etaDown + etaUp) / (hy*hy) - (etaRight) /(hx*hx);
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * etaUp    / (hy*hy);
          /* No left entry */
          col[3].i = ex+1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0; valA[3]  =        etaRight / (hx*hx);
          col[4].i = ex  ; col[4].j  = ey-1; col[4].loc  = LEFT;     col[4].c  = 0; valA[4]  =        etaLeft  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = RIGHT;    col[5].c  = 0; valA[5]  = -      etaRight / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j  = ey  ; col[6].loc  = LEFT;     col[6].c  = 0; valA[6]  =        etaRight / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = RIGHT;    col[7].c  = 0; valA[7]  = -      etaRight / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j  = ey-1; col[8].loc  = ELEMENT;  col[8].c  = 0; valA[8]  =  ctx->Kcont / hy;
          col[9].i = ex  ; col[9].j = ey   ; col[9].loc = ELEMENT;   col[9].c  = 0; valA[9]  = -ctx->Kcont / hy;
        } else if (ex == N[0]-1) {
          /* Right boundary y velocity stencil */
          nEntries = 10;
          row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * (etaDown + etaUp) / (hy*hy) - (etaLeft) /(hx*hx );
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * etaUp    / (hy*hy);
          col[3].i = ex-1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0; valA[3]  =        etaLeft  / (hx*hx);
          /* No right element */
          col[4].i = ex  ; col[4].j  = ey-1; col[4].loc  = LEFT;     col[4].c  = 0; valA[4]  =        etaLeft  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = RIGHT;    col[5].c  = 0; valA[5]  = -      etaRight / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j  = ey  ; col[6].loc  = LEFT;     col[6].c  = 0; valA[7]  =        etaRight / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = RIGHT;    col[7].c  = 0; valA[7]  = -      etaRight / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j  = ey-1; col[8].loc  = ELEMENT;  col[8].c  = 0; valA[8]  =  ctx->Kcont / hy;
          col[9].i = ex  ; col[9].j = ey   ; col[9].loc = ELEMENT;   col[9].c  = 0; valA[9]  = -ctx->Kcont / hy;
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
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = LEFT;     col[7].c  = 0; valA[7]  =        etaRight / (hx*hy); /* up left x edge */
          col[8].i = ex  ; col[8].j  = ey  ; col[8].loc  = RIGHT;    col[8].c  = 0; valA[8]  = -      etaRight / (hx*hy); /* up right x edge */
          col[9].i = ex  ; col[9].j  = ey-1; col[9].loc  = ELEMENT;  col[9].c  = 0; valA[9]  =  ctx->Kcont / hy;
          col[10].i = ex ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0; valA[10] = -ctx->Kcont / hy;
        }

        /* Insert Y-momentum entries */
        ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      if (ex == N[0]-1) {
        /* Right Boundary velocity Dirichlet */
        /* Redundant in the corner */
        DMStagStencil row;
        PetscScalar   valRhs;

        const PetscScalar valA = ctx->Kbound;
        row.i = ex; row.j = ey; row.loc = RIGHT; row.c = 0;
        ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }
      if (ex == 0) {
        /* Left velocity Dirichlet */
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = ctx->Kbound;
        row.i = ex; row.j = ey; row.loc = LEFT; row.c = 0;
        ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
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
        ierr = DMStagVecGetValuesStencil(ctx->dmCoeff,coeffLocal,4,etaPoint,eta);CHKERRQ(ierr);
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
          col[6].i  = ex-1; col[6].j  = ey  ; col[6].loc  = UP;      col[6].c   = 0; valA[6]  =        etaUp    / (hx*hy); /* up left */
          col[7].i  = ex  ; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  = -      etaUp    / (hx*hy); /* up right */
          col[8].i  = ex-1; col[8].j  = ey  ; col[8].loc  = ELEMENT; col[8].c   = 0; valA[8]  =  ctx->Kcont / hx;
          col[9].i = ex   ; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0; valA[9]  = -ctx->Kcont / hx;
          valRhs = 0.0;
        } else if (ey == N[1]-1) {
          /* Top boundary x velocity stencil */
          nEntries = 10;
          row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0; valA[0]  = -2.0 * (etaLeft + etaRight) / (hx*hx) -(etaDown) / (hy*hy);
          col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc  = LEFT;    col[1].c   = 0; valA[1]  =        etaDown  / (hy*hy);
          /* Missing element above */
          col[2].i  = ex-1; col[2].j  = ey  ; col[2].loc  = LEFT;    col[2].c   = 0; valA[2]  =  2.0 * etaLeft  / (hx*hx);
          col[3].i  = ex+1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0; valA[3]  =  2.0 * etaRight / (hx*hx);
          col[4].i  = ex-1; col[4].j  = ey  ; col[4].loc  = DOWN;    col[4].c   = 0; valA[4]  =        etaDown  / (hx*hy); /* down left */
          col[5].i  = ex  ; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0; valA[5]  = -      etaDown  / (hx*hy); /* down right */
          col[6].i  = ex-1; col[6].j  = ey  ; col[6].loc  = UP;      col[6].c   = 0; valA[6]  =        etaUp    / (hx*hy); /* up left */
          col[7].i  = ex  ; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  = -      etaUp    / (hx*hy); /* up right */
          col[8].i  = ex-1; col[8].j  = ey  ; col[8].loc  = ELEMENT; col[8].c   = 0; valA[8]  =  ctx->Kcont / hx;
          col[9].i = ex   ; col[9].j  = ey   ; col[9].loc = ELEMENT;  col[9].c  = 0; valA[9]  = -ctx->Kcont / hx;
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
          col[7].i  = ex-1; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  =        etaUp    / (hx*hy); /* up left */
          col[8].i  = ex  ; col[8].j  = ey  ; col[8].loc  = UP;      col[8].c   = 0; valA[8]  = -      etaUp    / (hx*hy); /* up right */
          col[9].i  = ex-1; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0; valA[9]  =  ctx->Kcont / hx;
          col[10].i = ex  ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0; valA[10] = -ctx->Kcont / hx;
          valRhs = 0.0;
        }
        ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      /* P equation : u_x + v_y = 0
         Note that this includes an explicit zero on the diagonal. This is only needed for
         direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace) */
      if (pinPressure && ex == 0 && ey == 0) { /* Pin the first pressure node to zero, if requested */
        DMStagStencil row;
        PetscScalar valA,valRhs;
        row.i = ex; row.j = ey; row.loc = ELEMENT; row.c = 0;
        valA = ctx->Kbound;
        ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        DMStagStencil row,col[5];
        PetscScalar   valA[5],valRhs;

        row.i    = ex; row.j    = ey; row.loc    = ELEMENT; row.c    = 0;
        col[0].i = ex; col[0].j = ey; col[0].loc = LEFT;    col[0].c = 0; valA[0] = -ctx->Kcont / hx;
        col[1].i = ex; col[1].j = ey; col[1].loc = RIGHT;   col[1].c = 0; valA[1] =  ctx->Kcont / hx;
        col[2].i = ex; col[2].j = ey; col[2].loc = DOWN;    col[2].c = 0; valA[2] = -ctx->Kcont / hy;
        col[3].i = ex; col[3].j = ey; col[3].loc = UP;      col[3].c = 0; valA[3] =  ctx->Kcont / hy;
        col[4] = row;                                                     valA[4] = 0.0;
        ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,5,col,valA,INSERT_VALUES);CHKERRQ(ierr);
        valRhs = 0.0;
        ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = DMRestoreLocalVector(ctx->dmCoeff,&coeffLocal);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rhs);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Here, we demonstrate getting coordinates from a vector by using DMStagStencil.
This would usually be done with direct array access, though. */
static PetscErrorCode PopulateCoefficientData(Ctx ctx)
{
  PetscErrorCode ierr;
  PetscInt       N[2],nDummy[2];
  PetscInt       ex,ey,startx,starty,nx,ny,ietaCorner,ietaElement,irho,iprev,icenter;
  Vec            coeffLocal;
  PetscReal      **cArrX,**cArrY;
  PetscScalar    ***coeffArr;

  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(ctx->dmCoeff,&coeffLocal);CHKERRQ(ierr);
  ierr = DMStagGetCorners(ctx->dmCoeff,&startx,&starty,NULL,&nx,&ny,NULL,&nDummy[0],&nDummy[1],NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(ctx->dmCoeff,&N[0],&N[1],NULL);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(ctx->dmCoeff,DMSTAG_DOWN_LEFT,0,&ietaCorner);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(ctx->dmCoeff,DMSTAG_ELEMENT,0,&ietaElement);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(ctx->dmCoeff,DMSTAG_DOWN_LEFT,1,&irho);CHKERRQ(ierr);

  ierr = DMStagGet1DCoordinateArraysDOFRead(ctx->dmCoeff,&cArrX,&cArrY,NULL);CHKERRQ(ierr);
  ierr = DMStagGet1DCoordinateLocationSlot(ctx->dmCoeff,DMSTAG_ELEMENT,&icenter);CHKERRQ(ierr);
  ierr = DMStagGet1DCoordinateLocationSlot(ctx->dmCoeff,DMSTAG_LEFT,&iprev);CHKERRQ(ierr);

  ierr = DMStagVecGetArrayDOF(ctx->dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);

  for (ey = starty; ey<starty+ny+nDummy[1]; ++ey) {
    for (ex = startx; ex<startx+nx+nDummy[0]; ++ex) {
      coeffArr[ey][ex][ietaElement] = getEta(ctx,cArrX[ex][icenter],cArrY[ey][icenter]); // Note dummy value filled here, needlessly
      coeffArr[ey][ex][ietaCorner]  = getEta(ctx,cArrX[ex][iprev],cArrY[ey][iprev]);
      coeffArr[ey][ex][irho]        = getRho(ctx,cArrX[ex][iprev],cArrY[ey][iprev]);
    }
  }
  ierr = DMStagRestore1DCoordinateArraysDOFRead(ctx->dmCoeff,&cArrX,&cArrY,NULL);CHKERRQ(ierr);
  ierr = DMStagVecRestoreArrayDOF(ctx->dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->dmCoeff,&ctx->coeff);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ctx->dmCoeff,coeffLocal,INSERT_VALUES,ctx->coeff);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ctx->dmCoeff,coeffLocal,INSERT_VALUES,ctx->coeff);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ctx->dmCoeff,&coeffLocal);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DumpSolution(Ctx ctx,Vec x)
{
  PetscErrorCode ierr;
  DM             dmVelAvg;
  Vec            velAvg;
  DM             daVelAvg,daP,daEtaElement,daEtaCorner,daRho;
  Vec            vecVelAvg,vecP,vecEtaElement,vecEtaCorner,vecRho;

  PetscFunctionBeginUser;

  /* For convenience, create a new DM and Vec which will hold averaged velocities
     Note that this could also be accomplished with direct array access, using
     DMStagVecGetArrayDOF() and related functions */
  ierr = DMStagCreateCompatibleDMStag(ctx->dmStokes,0,0,2,0,&dmVelAvg); /* 2 dof per element */
  ierr = DMSetUp(dmVelAvg);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesExplicit(dmVelAvg,0.0,ctx->xmax,0.0,ctx->ymax,0.0,0.0);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dmVelAvg,&velAvg);CHKERRQ(ierr);
  {
    PetscInt ex,ey,startx,starty,nx,ny;
    Vec      stokesLocal;
    ierr = DMGetLocalVector(ctx->dmStokes,&stokesLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(ctx->dmStokes,x,INSERT_VALUES,stokesLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(ctx->dmStokes,x,INSERT_VALUES,stokesLocal);CHKERRQ(ierr);
    ierr = DMStagGetCorners(dmVelAvg,&startx,&starty,NULL,&nx,&ny,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    for (ey = starty; ey<starty+ny; ++ey) {
      for (ex = startx; ex<startx+nx; ++ex) {
        DMStagStencil from[4],to[2];
        PetscScalar   valFrom[4],valTo[2];
        from[0].i = ex; from[0].j = ey; from[0].loc = UP;    from[0].c = 0;
        from[1].i = ex; from[1].j = ey; from[1].loc = DOWN;  from[1].c = 0;
        from[2].i = ex; from[2].j = ey; from[2].loc = LEFT;  from[2].c = 0;
        from[3].i = ex; from[3].j = ey; from[3].loc = RIGHT; from[3].c = 0;
        ierr = DMStagVecGetValuesStencil(ctx->dmStokes,stokesLocal,4,from,valFrom);CHKERRQ(ierr);
        to[0].i = ex; to[0].j = ey; to[0].loc = ELEMENT;    to[0].c = 0; valTo[0] = 0.5 * (valFrom[2] + valFrom[3]);
        to[1].i = ex; to[1].j = ey; to[1].loc = ELEMENT;    to[1].c = 1; valTo[1] = 0.5 * (valFrom[0] + valFrom[1]);
        ierr = DMStagVecSetValuesStencil(dmVelAvg,velAvg,2,to,valTo,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    ierr = VecAssemblyBegin(velAvg);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(velAvg);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(ctx->dmStokes,&stokesLocal);CHKERRQ(ierr);
  }

  ierr = DMStagVecSplitToDMDA(ctx->dmStokes,x,DMSTAG_ELEMENT,0,&daP,&vecP);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecP,"p (scaled)");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(ctx->dmCoeff,ctx->coeff,DMSTAG_DOWN_LEFT,0,&daEtaCorner,&vecEtaCorner);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecEtaCorner,"eta");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(ctx->dmCoeff,ctx->coeff,DMSTAG_ELEMENT,0,&daEtaElement,&vecEtaElement);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecEtaElement,"eta");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(ctx->dmCoeff,ctx->coeff,DMSTAG_DOWN_LEFT,1,&daRho,&vecRho);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecRho,"density");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(dmVelAvg,velAvg,DMSTAG_ELEMENT,-3,&daVelAvg,&vecVelAvg);CHKERRQ(ierr); /* note -3 : pad with zero */
  ierr = PetscObjectSetName((PetscObject)vecVelAvg,"Velocity (Averaged)");CHKERRQ(ierr);

  /* Dump element-based fields to a .vtr file */
  {
    PetscViewer viewer;
    ierr = PetscViewerVTKOpen(PetscObjectComm((PetscObject)daVelAvg),"out_element.vtr",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(vecVelAvg,viewer);CHKERRQ(ierr);
    ierr = VecView(vecP,viewer);CHKERRQ(ierr);
    ierr = VecView(vecEtaElement,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Dump vertex-based fields to a second .vtr file */
  {
    PetscViewer viewer;
    ierr = PetscViewerVTKOpen(PetscObjectComm((PetscObject)daEtaCorner),"out_vertex.vtr",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(vecEtaCorner,viewer);CHKERRQ(ierr);
    ierr = VecView(vecRho,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Edge-based fields could similarly be dumped */

  /* Destroy DMDAs and Vecs */
  ierr = VecDestroy(&vecVelAvg);CHKERRQ(ierr);
  ierr = VecDestroy(&vecP);CHKERRQ(ierr);
  ierr = VecDestroy(&vecEtaCorner);CHKERRQ(ierr);
  ierr = VecDestroy(&vecEtaElement);CHKERRQ(ierr);
  ierr = VecDestroy(&vecRho);CHKERRQ(ierr);
  ierr = DMDestroy(&daVelAvg);CHKERRQ(ierr);
  ierr = DMDestroy(&daP);CHKERRQ(ierr);
  ierr = DMDestroy(&daEtaCorner);CHKERRQ(ierr);
  ierr = DMDestroy(&daEtaElement);CHKERRQ(ierr);
  ierr = DMDestroy(&daRho);CHKERRQ(ierr);
  ierr = DMDestroy(&dmVelAvg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
