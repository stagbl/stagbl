#include "args.h"
#include "ctx.h"
#include "dump.h"
#include "system.h"
#include "system2.h"
#include <stagbl.h>
#include <stdio.h>
#include <petsc.h> // TODO REMOVE THIS INCLUDE
#include <mpi.h>

/* Helper functions */
static PetscErrorCode PopulateCoefficientData(Ctx);

/* Coefficient/forcing Functions */

/* Sinker */
static StagBLReal getRho_sinker(void *ptr,StagBLReal x, StagBLReal y) {
  Ctx ctx = (Ctx) ptr;
  const StagBLReal d = ctx->xmax-ctx->xmin;
  const StagBLReal xx = x/d - 0.5;
  const StagBLReal yy = y/d - 0.5;
  return (xx*xx + yy*yy) > 0.3*0.3 ? ctx->rho1 : ctx->rho2;
}
static StagBLReal getEta_sinker(void *ptr,StagBLReal x, StagBLReal y) {
  Ctx ctx = (Ctx) ptr;
  const StagBLReal d = ctx->xmax-ctx->xmin;
  const StagBLReal xx = x/d - 0.5;
  const StagBLReal yy = y/d - 0.5;
  return (xx*xx + yy*yy) > 0.3*0.3 ? ctx->eta1 : ctx->eta2;
}
/* Vertical layers */
StagBLReal getRho_gerya72(void *ptr,StagBLReal x,StagBLReal y)
{
  Ctx ctx = (Ctx) ptr;
  if (x + 0.0*y < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->rho1;
  } else {
    return ctx->rho2;
  }
}
StagBLReal getEta_gerya72(void *ptr,StagBLReal x,StagBLReal y)
{
  Ctx ctx = (Ctx) ptr;
  if (x  + 0.0*y < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->eta1;
  } else {
    return ctx->eta2;
  }
}

int main(int argc, char** argv)
{
  StagBLErrorCode ierr;
  int             rank,size;
  StagBLArray     x;
  StagBLSystem    system;
  StagBLSolver    solver;
  MPI_Comm        comm;
  StagBLInt       systemtype,structure;
  Ctx             ctx;

  // TODO remove all this petsc-dependence
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

  /* Accept argument for system type and coefficient structure */
  ierr = GetIntArg("-system",1,&systemtype);CHKERRQ(ierr);
  ierr = GetIntArg("-structure",1,&structure);CHKERRQ(ierr);

  /* Populate application context */
  ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);
  ctx->comm = comm;
  ctx->xmin = 0.0;
  ctx->xmax = 1e6;
  ctx->ymin = 0.0;
  ctx->ymax = 1.5e6;
  ctx->rho1 = 3200;
  ctx->rho2 = 3300;
  ctx->eta1 = 1e20;
  ctx->eta2 = 1e22;
  ctx->gy   = 10.0;

  switch (structure) {
    case 1:
      ctx->getEta = getEta_gerya72;
      ctx->getRho = getRho_gerya72;
      break;
    case 2:
      ctx->getEta = getEta_sinker;
      ctx->getRho = getRho_sinker;
      break;
    default: SETERRQ1(comm,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported viscosity structure %d",structure);
  }

  /* Create a Grid
     We call a helper function from "stokes" to create an interlaced p-v
     grid of the correct size. Note that the choice of whether to create a single
     grid with both pressure and velocity, or two grids, must be made early, 
     as this affects data layout. This single-grid choice is appropriate for
     monolithic mulitigrid or direct solution, whereas two grids would be appropriated
     for segregated or Approximate block factorization based solvers. */
  ierr = StagBLGridCreateStokes2DBox(comm,30,20,0.0,ctx->xmax,0.0,ctx->ymax,&ctx->stokesGrid);CHKERRQ(ierr);

  /* Get scaling constants and node to pin, knowing grid dimensions */
  {
    DM dmStokes;
    StagBLInt N[2];
    StagBLReal hxAvgInv;
    // TODO remove this escape hatch logic from main..
    ierr = StagBLGridPETScGetDM(ctx->stokesGrid,&dmStokes);CHKERRQ(ierr);
    ierr = DMStagGetGlobalSizes(dmStokes,&N[0],&N[1],NULL);CHKERRQ(ierr);
    ctx->hxCharacteristic = (ctx->xmax-ctx->xmin)/N[0];
    ctx->hyCharacteristic = (ctx->ymax-ctx->ymin)/N[1];
    ctx->etaCharacteristic = PetscMin(ctx->eta1,ctx->eta2);
    hxAvgInv = 2.0/(ctx->hxCharacteristic + ctx->hyCharacteristic);
    ctx->Kcont = ctx->etaCharacteristic*hxAvgInv;
    ctx->Kbound = ctx->etaCharacteristic*hxAvgInv*hxAvgInv;
    if (N[0] < 2) SETERRQ(comm,PETSC_ERR_SUP,"Not implemented for a single element in the x direction");
    ctx->pinx = 1; ctx->piny = 0;
  }

  /* Create another, compatible grid to represent coefficients */
 {
   const StagBLInt dofPerVertex  = 2;
   const StagBLInt dofPerElement = 1;
   ierr = StagBLGridCreateCompatibleStagBLGrid(ctx->stokesGrid,dofPerVertex,0,dofPerElement,0,&ctx->coeffGrid);CHKERRQ(ierr);
 }

 {
   DM dmCoeff;
  ierr = StagBLGridPETScGetDM(ctx->coeffGrid,&dmCoeff);CHKERRQ(ierr);
  // TODO by default set this same coordinate DM
  ierr = DMStagSetUniformCoordinatesProduct(dmCoeff,0.0,ctx->xmax,0.0,ctx->ymax,0.0,0.0);CHKERRQ(ierr);
 }

  /* Coefficient data */
  ierr = PopulateCoefficientData(ctx);CHKERRQ(ierr);

  /* Create a system */
  ierr = StagBLGridCreateStagBLSystem(ctx->stokesGrid,&system);CHKERRQ(ierr);

  ierr = StagBLGridCreateStagBLArray(ctx->stokesGrid,&x);CHKERRQ(ierr);

  // TODO move this escape hatching in to system building calls
  ierr = StagBLSystemPETScGetMatPointer(system,&pmatA);CHKERRQ(ierr);
  ierr = StagBLSystemPETScGetVecPointer(system,&pvecb);CHKERRQ(ierr);

  {
    DM dmStokes;
    ierr = StagBLGridPETScGetDM(ctx->stokesGrid,&dmStokes);CHKERRQ(ierr);
    // TODO make sure this escape hatching is gone
    ierr = StagBLArrayPETScGetGlobalVecPointer(x,&pvecx);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dmStokes,pvecx);CHKERRQ(ierr);
  }
  vecx = *pvecx;
  ierr = PetscObjectSetName((PetscObject)vecx,"solution");CHKERRQ(ierr);

  if (systemtype == 1) {
    ierr = CreateSystem(ctx,pmatA,pvecb);CHKERRQ(ierr);
  } else if (systemtype == 2) {
    ierr = CreateSystem2(ctx,pmatA,pvecb);CHKERRQ(ierr);
  } else SETERRQ1(ctx->comm,PETSC_ERR_SUP,"Unsupported system type %D",system);
  matA = *pmatA;
  vecb = *pvecb;

  // TODO we call a function to create the default solver for the created StagBLSystem function,
  // set any additional options we'd like, and set up

  // TODO StagBLSystemCreateStagBLSolver(stokesSystem,&solver);
  // TODO StagBLSolverSolve(solver,...);

  /* Solve the system (you will likely want to specify a solver from the command line) */
  StagBLSolverCreate(&solver);
  StagBLSolverPETScGetKSPPointer(solver,&pksp);

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
  ierr = StagBLArrayDestroy(&x);CHKERRQ(ierr);
  ierr = StagBLSystemDestroy(&system);CHKERRQ(ierr);
  ierr = StagBLSolverDestroy(&solver);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&ctx->stokesGrid);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&ctx->coeffGrid);CHKERRQ(ierr);

  StagBLFinalize();

  return 0;
}

static PetscErrorCode PopulateCoefficientData(Ctx ctx)
{
  PetscErrorCode ierr;
  PetscInt       N[2];
  PetscInt       ex,ey,startx,starty,nx,ny,ietaCorner,ietaElement,irho,iprev,icenter;
  DM             dmCoeff;
  Vec            *pcoeffLocal;
  Vec            coeffLocal;
  PetscReal      **cArrX,**cArrY;
  PetscReal      ***coeffArr;

  PetscFunctionBeginUser;

  ierr = StagBLGridCreateStagBLArray(ctx->coeffGrid,&ctx->coeffArray);CHKERRQ(ierr);

  /* Escape Hatch */
  ierr = StagBLGridPETScGetDM(ctx->coeffGrid,&dmCoeff);CHKERRQ(ierr);
  ierr = StagBLArrayPETScGetLocalVecPointer(ctx->coeffArray,&pcoeffLocal);CHKERRQ(ierr);

  ierr = DMCreateLocalVector(dmCoeff,pcoeffLocal);CHKERRQ(ierr);
  coeffLocal = *pcoeffLocal;
  ierr = DMStagGetGhostCorners(dmCoeff,&startx,&starty,NULL,&nx,&ny,NULL);CHKERRQ(ierr); /* Iterate over all local elements */
  ierr = DMStagGetGlobalSizes(dmCoeff,&N[0],&N[1],NULL);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_DOWN_LEFT,0,&ietaCorner);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_ELEMENT,0,&ietaElement);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dmCoeff,DMSTAG_DOWN_LEFT,1,&irho);CHKERRQ(ierr);

  ierr = DMStagGet1dCoordinateArraysDOFRead(dmCoeff,&cArrX,&cArrY,NULL);CHKERRQ(ierr);
  ierr = DMStagGet1dCoordinateLocationSlot(dmCoeff,DMSTAG_ELEMENT,&icenter);CHKERRQ(ierr);
  ierr = DMStagGet1dCoordinateLocationSlot(dmCoeff,DMSTAG_LEFT,&iprev);CHKERRQ(ierr);

  ierr = DMStagVecGetArrayDOF(dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);

  for (ey = starty; ey<starty+ny; ++ey) {
    for (ex = startx; ex<startx+nx; ++ex) {
      coeffArr[ey][ex][ietaElement] = ctx->getEta(ctx,cArrX[ex][icenter],cArrY[ey][icenter]);
      coeffArr[ey][ex][ietaCorner]  = ctx->getEta(ctx,cArrX[ex][iprev],cArrY[ey][iprev]);
      coeffArr[ey][ex][irho]        = ctx->getRho(ctx,cArrX[ex][iprev],cArrY[ey][iprev]);
    }
  }
  ierr = DMStagRestore1dCoordinateArraysDOFRead(dmCoeff,&cArrX,&cArrY,NULL);CHKERRQ(ierr);
  ierr = DMStagVecRestoreArrayDOF(dmCoeff,coeffLocal,&coeffArr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
