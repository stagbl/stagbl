static const char *help = "StagBLDemo3D: Demonstrate features and usage of StagBL, in 3 dimensions, with simple geodynamic box model setups\n\n";

#include <stagbl.h>
#include <stdio.h>
#include <petsc.h>

// Note: This is not a complete demo yet.
// It is mainly here to test 3d grid functionality and will be completely changed.

/* Shorter, more convenient names for DMStagLocation entries */
#define BACK_DOWN_LEFT   DMSTAG_BACK_DOWN_LEFT
#define BACK_DOWN        DMSTAG_BACK_DOWN
#define BACK_DOWN_RIGHT  DMSTAG_BACK_DOWN_RIGHT
#define BACK_LEFT        DMSTAG_BACK_LEFT
#define BACK             DMSTAG_BACK
#define BACK_RIGHT       DMSTAG_BACK_RIGHT
#define BACK_UP_LEFT     DMSTAG_BACK_UP_LEFT
#define BACK_UP          DMSTAG_BACK_UP
#define BACK_UP_RIGHT    DMSTAG_BACK_UP_RIGHT
#define DOWN_LEFT        DMSTAG_DOWN_LEFT
#define DOWN             DMSTAG_DOWN
#define DOWN_RIGHT       DMSTAG_DOWN_RIGHT
#define LEFT             DMSTAG_LEFT
#define ELEMENT          DMSTAG_ELEMENT
#define RIGHT            DMSTAG_RIGHT
#define UP_LEFT          DMSTAG_UP_LEFT
#define UP               DMSTAG_UP
#define UP_RIGHT         DMSTAG_UP_RIGHT
#define FRONT_DOWN_LEFT  DMSTAG_FRONT_DOWN_LEFT
#define FRONT_DOWN       DMSTAG_FRONT_DOWN
#define FRONT_DOWN_RIGHT DMSTAG_FRONT_DOWN_RIGHT
#define FRONT_LEFT       DMSTAG_FRONT_LEFT
#define FRONT            DMSTAG_FRONT
#define FRONT_RIGHT      DMSTAG_FRONT_RIGHT
#define FRONT_UP_LEFT    DMSTAG_FRONT_UP_LEFT
#define FRONT_UP         DMSTAG_FRONT_UP
#define FRONT_UP_RIGHT   DMSTAG_FRONT_UP_RIGHT

/* An application context */
typedef struct {
  MPI_Comm    comm;
  DM          dmStokes,dmCoeff;
  Vec         coeff;
  PetscReal   xmax,ymax,zmax,xmin,ymin,zmin,hxCharacteristic,hyCharacteristic,hzCharacteristic;
  PetscScalar eta1,eta2,rho1,rho2,gy,Kbound,Kcont,etaCharacteristic;
} CtxData;
typedef CtxData* Ctx;

/* Helper functions */
static PetscErrorCode PopulateCoefficientData(Ctx);
static PetscErrorCode CreateSystem(const Ctx,Mat*,Vec*);
static PetscErrorCode DumpSolution(Ctx,Vec);

/* Coefficient/forcing Functions */

static PetscReal getRho(Ctx ctx,PetscReal x, PetscReal y, PetscReal z) {
  const PetscReal d = ctx->xmax-ctx->xmin;
  const PetscReal xx = x/d - 0.5;
  const PetscReal yy = y/d - 0.5;
  const PetscReal zz = z/d - 0.5;
  return (xx*xx + yy*yy + zz*zz) > 0.3*0.3*0.3 ? ctx->rho1 : ctx->rho2;
}
static PetscReal getEta(Ctx ctx,PetscReal x, PetscReal y, PetscReal z) {
  const PetscReal d = ctx->xmax-ctx->xmin;
  const PetscReal xx = x/d - 0.5;
  const PetscReal yy = y/d - 0.5;
  const PetscReal zz = z/d - 0.5;
  return (xx*xx + yy*yy + zz*zz) > 0.3*0.3*0.3 ? ctx->eta1 : ctx->eta2;
}

int main(int argc, char** argv)
{
  int          rank,size;
  StagBLGrid   grid;
  StagBLArray  x;
  StagBLSystem system;
  StagBLSolver solver;
  MPI_Comm     comm;
  PetscInt     Nx,Ny,Nz;

  Ctx          ctx;

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
    printf("=== StagBLDemo3d ===\n");
    printf("%d ranks\n",size);
  }
  MPI_Barrier(comm);

  /* Initialize StagBL (which will initialize PETSc if needbe) */
  ierr = StagBLInitialize(argc,argv,help,comm);CHKERRQ(ierr);

  /* Populate application context */
  ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);
  ctx->comm = PETSC_COMM_WORLD;
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
  ctx->gy   = 10.0;

  // Create a Grid
  Nx = 20; Ny = 40; Nz = 20;
  {
    PetscBool flg;
    ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&Nx,&flg);CHKERRQ(ierr);
    if (flg) {
      Ny = Nz = Nx;
    }
    ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&Ny,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-nz",&Nz,NULL);CHKERRQ(ierr);
  }
  ierr = StagBLGridCreate(&grid);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDMPointer(grid,&pdm);CHKERRQ(ierr);
  ierr = DMStagCreate3d(
      comm,
      DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
      Nx,Ny,Nz,                                /* Global element counts */
      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,  /* Determine parallel decomposition automatically */
      0,0,1,1,                                 /* dof: 1 dof on each face and 3-cell */
      DMSTAG_STENCIL_BOX,
      1,                                       /* elementwise stencil width */
      NULL,NULL,NULL,
      pdm);
  ctx->dmStokes = *pdm;
  ierr = DMSetFromOptions(ctx->dmStokes);CHKERRQ(ierr);
  ierr = DMSetUp(ctx->dmStokes);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesProduct(ctx->dmStokes,0.0,ctx->xmax,0.0,ctx->ymax,0.0,ctx->zmax);CHKERRQ(ierr);
  ierr = DMStagCreateCompatibleDMStag(ctx->dmStokes,1,0,0,2,&ctx->dmCoeff);CHKERRQ(ierr); /* 1 dof per vertex, 2 per element */
  ierr = DMSetUp(ctx->dmCoeff);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesProduct(ctx->dmCoeff,0.0,ctx->xmax,0.0,ctx->ymax,0.0,ctx->zmax);CHKERRQ(ierr);

  /* Get scaling constants, knowing grid spacing */
  {
    PetscInt N[3];
    PetscReal hxAvgInv;
    ierr = DMStagGetGlobalSizes(ctx->dmStokes,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
    ctx->hxCharacteristic = (ctx->xmax-ctx->xmin)/N[0];
    ctx->hyCharacteristic = (ctx->ymax-ctx->ymin)/N[1];
    ctx->hzCharacteristic = (ctx->zmax-ctx->zmin)/N[2];
    ctx->etaCharacteristic = PetscMin(ctx->eta1,ctx->eta2);
    hxAvgInv = 2.0/(ctx->hxCharacteristic + ctx->hyCharacteristic);
    ctx->Kcont = ctx->etaCharacteristic*hxAvgInv;
    ctx->Kbound = ctx->etaCharacteristic*hxAvgInv*hxAvgInv;
  }

  /* Populate coefficient data */
  ierr = PopulateCoefficientData(ctx);CHKERRQ(ierr);

  // Create a system
  ierr = StagBLGridCreateStagBLSystem(grid,&system);CHKERRQ(ierr);
  ierr = StagBLGridCreateStagBLArray(grid,&x);CHKERRQ(ierr);

  ierr = StagBLArrayPETScGetGlobalVecPointer(x,&pvecx);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->dmStokes,pvecx);
  vecx = *pvecx;

  ierr = StagBLSystemPETScGetVecPointer(system,&pvecb);CHKERRQ(ierr);
  ierr = StagBLSystemPETScGetMatPointer(system,&pmatA);CHKERRQ(ierr);
  ierr = CreateSystem(ctx,pmatA,pvecb);CHKERRQ(ierr);
  matA = *pmatA;
  vecb = *pvecb;

  // Solve the system (you will likely want to choose a solver from the command line)
  ierr = StagBLSolverCreate(system,&solver);CHKERRQ(ierr);
  ierr = StagBLSolverPETScGetKSPPointer(solver,&pksp);CHKERRQ(ierr);

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
  ierr = VecDestroy(pvecx);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->coeff);CHKERRQ(ierr);
  ierr = StagBLArrayDestroy(&x);CHKERRQ(ierr);
  ierr = StagBLSystemDestroy(&system);CHKERRQ(ierr);
  ierr = StagBLSolverDestroy(&solver);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&grid);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->dmCoeff);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);

  ierr = StagBLFinalize();

  return ierr;
}

// Note: this is not  properly scaled
static PetscErrorCode CreateSystem(const Ctx ctx,Mat *pA,Vec *pRhs)
{
  PetscErrorCode ierr;
  PetscInt       N[3];
  PetscInt       ex,ey,ez,startx,starty,startz,nx,ny,nz;
  Mat            A;
  Vec            rhs;
  PetscReal      hx,hy,hz;
  const PetscBool pinPressure = PETSC_TRUE;
  Vec            coeffLocal;

  PetscFunctionBeginUser;
  ierr = DMCreateMatrix(ctx->dmStokes,pA);CHKERRQ(ierr);
  A = *pA;
  ierr = DMCreateGlobalVector(ctx->dmStokes,pRhs);CHKERRQ(ierr);
  rhs = *pRhs;
  ierr = DMStagGetCorners(ctx->dmStokes,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(ctx->dmStokes,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
  hx = ctx->hxCharacteristic;
  hy = ctx->hyCharacteristic;
  hz = ctx->hzCharacteristic;
  ierr = DMGetLocalVector(ctx->dmCoeff,&coeffLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ctx->dmCoeff,ctx->coeff,INSERT_VALUES,coeffLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ctx->dmCoeff,ctx->coeff,INSERT_VALUES,coeffLocal);CHKERRQ(ierr);

  // TODO the rest of this function is a placeholder, from DMStag tutorial 3. To be replaced by an actual stokes sytem.
  for (ez = startz; ez<startz+nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
    for (ey = starty; ey<starty+ny; ++ey) {
      for (ex = startx; ex<startx+nx; ++ex) {

        if (ex == N[0]-1) {
          /* Right Boundary velocity Dirichlet */
          DMStagStencil row;
          PetscScalar   valRhs;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = RIGHT; row.c = 0;
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        }
        if (ey == N[1]-1) {
          /* Top boundary velocity Dirichlet */
          DMStagStencil row;
          PetscScalar   valRhs;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = UP; row.c = 0;
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        }
        if (ez == N[2]-1) {
          /* Top boundary velocity Dirichlet */
          DMStagStencil row;
          PetscScalar   valRhs;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = FRONT; row.c = 0;
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Equation on left face of this element */
        if (ex == 0) {
          /* Left velocity Dirichlet */
          DMStagStencil row;
          PetscScalar   valRhs;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = LEFT; row.c = 0;
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          /* X-momentum interior equation : (u_xx + u_yy + u_zz) - p_x = f^x */
          DMStagStencil row,col[9];
          PetscScalar   valA[9],valRhs;
          PetscInt      nEntries;

          row.i = ex; row.j = ey; row.k = ez; row.loc = LEFT; row.c = 0;
          if (ey == 0) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -1.0 / (hz*hz);
              /* Missing down term */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Missing back term */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -1.0 / (hz*hz);
              /* Missing down term */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Missing front term */
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              /* Missing down term */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
            }
          } else if (ey == N[1]-1) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Missing up term */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Missing back entry */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Missing up term */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Missing front term */
              col[5].i = ex-1; col[5].j = ey  ;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hx;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hx;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Missing up term */
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
            }
          } else if (ez == 0) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            /* Missing back term */
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
          } else if (ez == N[2]-1) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            /* Missing front term */
            col[6].i = ex-1; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hx;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hx;
          } else {
            nEntries = 9;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = LEFT;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = LEFT;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = LEFT;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = LEFT;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = LEFT;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = LEFT;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez+1; col[6].loc = LEFT;    col[6].c  = 0; valA[6] =  1.0 / (hz*hz);
            col[7].i = ex-1; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] =  1.0 / hx;
            col[8].i = ex  ; col[8].j = ey  ;  col[8].k = ez  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8] = -1.0 / hx;
          }
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Equation on bottom face of this element */
        if (ey == 0) {
          /* Bottom boundary velocity Dirichlet */
          DMStagStencil row;
          PetscScalar   valRhs;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = DOWN; row.c = 0;
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          /* Y-momentum equation, (v_xx + v_yy + v_zz) - p_y = f^y */
          DMStagStencil row,col[9],rhoPoint[2];
          PetscScalar   valA[9],rho[2],valRhs;
          PetscInt      nEntries;

          /* get rho values  and compute rhs value*/
          rhoPoint[0].i = ex; rhoPoint[0].j = ey;   rhoPoint[0].k = ez; rhoPoint[0].loc = ELEMENT; rhoPoint[0].c = 1;
          rhoPoint[1].i = ex; rhoPoint[1].j = ey-1; rhoPoint[1].k = ez; rhoPoint[1].loc = ELEMENT; rhoPoint[1].c = 1;
          ierr = DMStagVecGetValuesStencil(ctx->dmCoeff,coeffLocal,2,rhoPoint,rho);CHKERRQ(ierr);
          valRhs = -ctx->gy * 0.5 * (rho[0] + rho[1]);

          row.i    = ex  ; row.j    = ey  ;  row.k    = ez  ; row.loc    = DOWN;    row.c     = 0;
          if (ex ==0) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Back term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Front term missing */
              col[5].i = ex  ; col[4].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[5].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
            }
          } else if (ex == N[0]-1) {
            if (ez == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              /* Back term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else if (ez == N[2]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              /* Front term missing */
              col[5].i = ex  ; col[5].j = ey-1;  col[5].k = ez  ; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hy;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hy;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
            }
          } else if (ez == 0) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            /* Back term missing */
            col[5].i = ez  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
          } else if (ez == N[2]-1) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -1.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            /* Front term missing */
            col[6].i = ex  ; col[6].j = ey-1;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hy;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hy;
          } else {
            nEntries = 9;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = DOWN;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = DOWN;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = DOWN;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = DOWN;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = DOWN;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = DOWN;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez+1; col[6].loc = DOWN;    col[6].c  = 0; valA[6] =  1.0 / (hz*hz);
            col[7].i = ex  ; col[7].j = ey-1;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] =  1.0 / hy;
            col[8].i = ex  ; col[8].j = ey  ;  col[8].k = ez  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8] = -1.0 / hy;
          }
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Equation on back face of this element */
        if (ez == 0) {
          /* Back boundary velocity Dirichlet */
          DMStagStencil row;
          PetscScalar   valRhs;
          const PetscScalar valA = 1.0;
          row.i = ex; row.j = ey; row.k = ez; row.loc = BACK; row.c = 0;
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          /* Z-momentum equation, (w_xx + w_yy + w_zz) - p_z = f^z */
          DMStagStencil row,col[9];
          PetscScalar   valA[9],valRhs;
          PetscInt      nEntries;

          row.i = ex; row.j = ey; row.k = ez; row.loc = BACK; row.c = 0;
          if (ex == 0) {
            if (ey == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              /* Down term missing */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[2].i = ex+1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else if (ey == N[1]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Up term missing */
              /* Left term missing */
              col[2].i = ex+1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              /* Left term missing */
              col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
            }
          } else if (ex == N[0]-1) {
            if (ey == 0) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              /* Down term missing */
              col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
              /* Right term missing */
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else if (ey == N[1]-1) {
              nEntries = 7;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              /* Up term missing */
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[3].i = ex  ; col[3].j = ey  ;  col[3].k = ez-1; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hz*hz);
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez+1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = ELEMENT; col[5].c  = 0; valA[5] =  1.0 / hz;
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez  ; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] = -1.0 / hz;
            } else {
              nEntries = 8;
              col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -1.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
              col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
              col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
              col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
              /* Right term missing */
              col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
              col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
              col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
              col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
            }
          } else if (ey == 0) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
            /* Down term missing */
            col[1].i = ex  ; col[1].j = ey+1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
            col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
          } else if (ey == N[1]-1) {
            nEntries = 8;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -1.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            /* Up term missing */
            col[2].i = ex-1; col[2].j = ey  ;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hx*hx);
            col[3].i = ex+1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex  ; col[4].j = ey  ;  col[4].k = ez-1; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hz*hz);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez+1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez-1; col[6].loc = ELEMENT; col[6].c  = 0; valA[6] =  1.0 / hz;
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez  ; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] = -1.0 / hz;
          } else {
            nEntries = 9;
            col[0].i = ex  ; col[0].j = ey  ;  col[0].k = ez  ; col[0].loc = BACK;    col[0].c  = 0; valA[0] = -2.0 / (hx*hx) + -2.0 / (hy*hy) -2.0 / (hz*hz);
            col[1].i = ex  ; col[1].j = ey-1;  col[1].k = ez  ; col[1].loc = BACK;    col[1].c  = 0; valA[1] =  1.0 / (hy*hy);
            col[2].i = ex  ; col[2].j = ey+1;  col[2].k = ez  ; col[2].loc = BACK;    col[2].c  = 0; valA[2] =  1.0 / (hy*hy);
            col[3].i = ex-1; col[3].j = ey  ;  col[3].k = ez  ; col[3].loc = BACK;    col[3].c  = 0; valA[3] =  1.0 / (hx*hx);
            col[4].i = ex+1; col[4].j = ey  ;  col[4].k = ez  ; col[4].loc = BACK;    col[4].c  = 0; valA[4] =  1.0 / (hx*hx);
            col[5].i = ex  ; col[5].j = ey  ;  col[5].k = ez-1; col[5].loc = BACK;    col[5].c  = 0; valA[5] =  1.0 / (hz*hz);
            col[6].i = ex  ; col[6].j = ey  ;  col[6].k = ez+1; col[6].loc = BACK;    col[6].c  = 0; valA[6] =  1.0 / (hz*hz);
            col[7].i = ex  ; col[7].j = ey  ;  col[7].k = ez-1; col[7].loc = ELEMENT; col[7].c  = 0; valA[7] =  1.0 / hz;
            col[8].i = ex  ; col[8].j = ey  ;  col[8].k = ez  ; col[8].loc = ELEMENT; col[8].c  = 0; valA[8] = -1.0 / hz;
          }
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,nEntries,col,valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* P equation : u_x + v_y + w_z = g
           Note that this includes an explicit zero on the diagonal. This is only needed for
           direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace) */
        if (pinPressure && ex == 0 && ey == 0 && ez == 0) { /* Pin the first pressure node, if requested */
          DMStagStencil row;
          PetscScalar valA,valRhs;
          row.i = ex; row.j = ey; row.k = ez; row.loc  = ELEMENT; row.c = 0;
          valA = 1.0;
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,1,&row,&valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          DMStagStencil row,col[7];
          PetscScalar   valA[7],valRhs;

          row.i    = ex; row.j    = ey; row.k    = ez; row.loc    = ELEMENT; row.c    = 0;
          col[0].i = ex; col[0].j = ey; col[0].k = ez; col[0].loc = LEFT;    col[0].c = 0; valA[0] = -1.0 / hx;
          col[1].i = ex; col[1].j = ey; col[1].k = ez; col[1].loc = RIGHT;   col[1].c = 0; valA[1] =  1.0 / hx;
          col[2].i = ex; col[2].j = ey; col[2].k = ez; col[2].loc = DOWN;    col[2].c = 0; valA[2] = -1.0 / hy;
          col[3].i = ex; col[3].j = ey; col[3].k = ez; col[3].loc = UP;      col[3].c = 0; valA[3] =  1.0 / hy;
          col[4].i = ex; col[4].j = ey; col[4].k = ez; col[4].loc = BACK;    col[4].c = 0; valA[4] = -1.0 / hz;
          col[5].i = ex; col[5].j = ey; col[5].k = ez; col[5].loc = FRONT;   col[5].c = 0; valA[5] =  1.0 / hz;
          col[6]   = row;                                                                  valA[6] =  0.0;
          ierr = DMStagMatSetValuesStencil(ctx->dmStokes,A,1,&row,7,col,valA,INSERT_VALUES);CHKERRQ(ierr);
          valRhs = 0.0;
          ierr = DMStagVecSetValuesStencil(ctx->dmStokes,rhs,1,&row,&valRhs,INSERT_VALUES);CHKERRQ(ierr);
        }
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
  PetscInt       N[3],nExtra[3];
  PetscInt       ex,ey,ez,startx,starty,startz,nx,ny,nz,ietaCorner,ietaElement,irho,iprev,icenter;
  Vec            coeffLocal;
  PetscReal      **cArrX,**cArrY,**cArrZ;
  PetscScalar    ****coeffArr;

  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(ctx->dmCoeff,&coeffLocal);CHKERRQ(ierr);
  ierr = DMStagGetCorners(ctx->dmCoeff,&startx,&starty,&startz,&nx,&ny,&nz,&nExtra[0],&nExtra[1],&nExtra[2]);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(ctx->dmCoeff,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(ctx->dmCoeff,DMSTAG_BACK_DOWN_LEFT,0,&ietaCorner);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(ctx->dmCoeff,DMSTAG_ELEMENT,0,&ietaElement);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(ctx->dmCoeff,DMSTAG_ELEMENT,1,&irho);CHKERRQ(ierr);

  ierr = DMStagGetProductCoordinateArraysRead(ctx->dmCoeff,&cArrX,&cArrY,&cArrZ);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(ctx->dmCoeff,DMSTAG_ELEMENT,&icenter);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(ctx->dmCoeff,DMSTAG_LEFT,&iprev);CHKERRQ(ierr);

  ierr = DMStagVecGetArray(ctx->dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);

  for (ez = startz; ez<startz+nz+nExtra[2]; ++ez) {
    for (ey = starty; ey<starty+ny+nExtra[1]; ++ey) {
      for (ex = startx; ex<startx+nx+nExtra[0]; ++ex) {
        coeffArr[ez][ey][ex][ietaElement] = getEta(ctx,cArrX[ex][icenter],cArrY[ey][icenter],cArrZ[ez][icenter]); // Note dummy value filled here, needlessly
        coeffArr[ez][ey][ex][ietaCorner]  = getEta(ctx,cArrX[ex][iprev],  cArrY[ey][iprev],cArrZ[ez][iprev]);
        coeffArr[ez][ey][ex][irho]        = getRho(ctx,cArrX[ex][icenter],cArrY[ey][icenter],cArrZ[ez][icenter]);
      }
    }
  }
  ierr = DMStagRestoreProductCoordinateArraysRead(ctx->dmCoeff,&cArrX,&cArrY,&cArrZ);CHKERRQ(ierr);
  ierr = DMStagVecRestoreArray(ctx->dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);
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
     DMStagVecGetArray() and related functions */
  ierr = DMStagCreateCompatibleDMStag(ctx->dmStokes,0,0,0,3,&dmVelAvg); /* 2 dof per element */
  ierr = DMSetUp(dmVelAvg);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesProduct(dmVelAvg,0.0,ctx->xmax,0.0,ctx->ymax,0.0,ctx->zmax);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dmVelAvg,&velAvg);CHKERRQ(ierr);
  {
    PetscInt ex,ey,ez,startx,starty,startz,nx,ny,nz;
    Vec      stokesLocal;
    ierr = DMGetLocalVector(ctx->dmStokes,&stokesLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(ctx->dmStokes,x,INSERT_VALUES,stokesLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(ctx->dmStokes,x,INSERT_VALUES,stokesLocal);CHKERRQ(ierr);
    ierr = DMStagGetCorners(dmVelAvg,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);
    for (ez = startz; ez<startz+nz; ++ez) {
      for (ey = starty; ey<starty+ny; ++ey) {
        for (ex = startx; ex<startx+nx; ++ex) {
          DMStagStencil from[6],to[3];
          PetscScalar   valFrom[6],valTo[3];
          from[0].i = ex; from[0].j = ey; from[0].k = ez; from[0].loc = UP;    from[0].c = 0;
          from[1].i = ex; from[1].j = ey; from[1].k = ez; from[1].loc = DOWN;  from[1].c = 0;
          from[2].i = ex; from[2].j = ey; from[2].k = ez; from[2].loc = LEFT;  from[2].c = 0;
          from[3].i = ex; from[3].j = ey; from[3].k = ez; from[3].loc = RIGHT; from[3].c = 0;
          from[4].i = ex; from[4].j = ey; from[4].k = ez; from[4].loc = BACK;  from[4].c = 0;
          from[5].i = ex; from[5].j = ey; from[5].k = ez; from[5].loc = FRONT; from[5].c = 0;
          ierr = DMStagVecGetValuesStencil(ctx->dmStokes,stokesLocal,6,from,valFrom);CHKERRQ(ierr);
          to[0].i = ex; to[0].j = ey; to[0].k = ez; to[0].loc = ELEMENT;    to[0].c = 0; valTo[0] = 0.5 * (valFrom[2] + valFrom[3]);
          to[1].i = ex; to[1].j = ey; to[1].k = ez; to[1].loc = ELEMENT;    to[1].c = 1; valTo[1] = 0.5 * (valFrom[0] + valFrom[1]);
          to[2].i = ex; to[2].j = ey; to[2].k = ez; to[2].loc = ELEMENT;    to[2].c = 2; valTo[2] = 0.5 * (valFrom[4] + valFrom[5]);
          ierr = DMStagVecSetValuesStencil(dmVelAvg,velAvg,3,to,valTo,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
    ierr = VecAssemblyBegin(velAvg);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(velAvg);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(ctx->dmStokes,&stokesLocal);CHKERRQ(ierr);
  }

  ierr = DMStagVecSplitToDMDA(ctx->dmStokes,x,DMSTAG_ELEMENT,0,&daP,&vecP);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecP,"p (scaled)");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(ctx->dmCoeff,ctx->coeff,DMSTAG_BACK_DOWN_LEFT,0, &daEtaCorner, &vecEtaCorner);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecEtaCorner,"eta");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(ctx->dmCoeff,ctx->coeff,DMSTAG_ELEMENT,  0, &daEtaElement,&vecEtaElement);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecEtaElement,"eta");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(ctx->dmCoeff,ctx->coeff,DMSTAG_ELEMENT,  1, &daRho,       &vecRho);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecRho,"density");CHKERRQ(ierr);
  ierr = DMStagVecSplitToDMDA(dmVelAvg,    velAvg,    DMSTAG_ELEMENT,  -3,&daVelAvg,    &vecVelAvg);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)vecVelAvg,"Velocity (Averaged)");CHKERRQ(ierr);

  /* Dump element-based fields to a .vtr file */
  {
    PetscViewer viewer;
    ierr = PetscViewerVTKOpen(PetscObjectComm((PetscObject)daVelAvg),"out_element.vtr",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(vecVelAvg,viewer);CHKERRQ(ierr);
    ierr = VecView(vecP,viewer);CHKERRQ(ierr);
    ierr = VecView(vecEtaElement,viewer);CHKERRQ(ierr);
    ierr = VecView(vecRho,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Dump vertex-based fields to a second .vtr file */
  {
    PetscViewer viewer;
    ierr = PetscViewerVTKOpen(PetscObjectComm((PetscObject)daEtaCorner),"out_vertex.vtr",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(vecEtaCorner,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Edge-based fields could similarly be dumped */

  /* Destroy DMDAs and Vecs */
  ierr = VecDestroy(&velAvg);CHKERRQ(ierr);
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
