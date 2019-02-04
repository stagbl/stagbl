#include "args.h"
#include "coeff.h"
#include "ctx.h"
#include "dump.h"
#include "system.h"
#include "system2.h"
#include <stagbl.h>
#include <stdio.h>
#include <petsc.h> // TODO REMOVE THIS INCLUDE (and make sure it's not included indirectly! Should FAIL if you use PETSc directly here)
#include <mpi.h>

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

  /* Populate application context (sets parameters) */
  ierr = CtxCreate(comm,&ctx);CHKERRQ(ierr);

  /* Create a Grid
     We call a helper function from "stokes" to create an interlaced p-v
     grid of the correct size. Note that the choice of whether to create a single
     grid with both pressure and velocity, or two grids, must be made early, 
     as this affects data layout. This single-grid choice is appropriate for
     monolithic mulitigrid or direct solution, whereas two grids would be appropriated
     for segregated or Approximate block factorization based solvers. */
  ierr = StagBLGridCreateStokes2DBox(comm,30,20,0.0,ctx->xmax,0.0,ctx->ymax,&ctx->stokesGrid);CHKERRQ(ierr);

  /* Get scaling constants and node to pin, knowing grid dimensions */
  ierr = CtxSetupFromGrid(ctx);CHKERRQ(ierr);

  /* Create another, compatible grid to represent coefficients */
 {
   const StagBLInt dofPerVertex  = 2;
   const StagBLInt dofPerElement = 1;
   ierr = StagBLGridCreateCompatibleStagBLGrid(ctx->stokesGrid,dofPerVertex,0,dofPerElement,0,&ctx->coeffGrid);CHKERRQ(ierr);
 }

 {
   DM dmCoeff;
  ierr = StagBLGridPETScGetDM(ctx->coeffGrid,&dmCoeff);CHKERRQ(ierr);
  // TODO by default set this same coordinate DM (confirm that ref count is incremented..)
  ierr = DMStagSetUniformCoordinatesProduct(dmCoeff,0.0,ctx->xmax,0.0,ctx->ymax,0.0,0.0);CHKERRQ(ierr);
 }

  /* Coefficient data */
  ierr = PopulateCoefficientData(ctx,structure);CHKERRQ(ierr);

  /* Create a system */
  ierr = StagBLGridCreateStagBLSystem(ctx->stokesGrid,&system);CHKERRQ(ierr);

  if (systemtype == 1) {
    ierr = CreateSystem(ctx,system);CHKERRQ(ierr);
  } else if (systemtype == 2) {
    ierr = CreateSystem2(ctx,system);CHKERRQ(ierr);
  } else StagBLError(ctx->comm,"Unsupported system type");

  /* Solve the system (you will likely want to specify a solver from the command line,
     e.g. -pc_type lu -pc_factor_mat_solver_type umfpack) */
  ierr = StagBLSystemCreateStagBLSolver(system,&solver);CHKERRQ(ierr);
  ierr = StagBLGridCreateStagBLArray(ctx->stokesGrid,&x);CHKERRQ(ierr);
  ierr = StagBLSolverSolve(solver,x);CHKERRQ(ierr);

  /* Dump solution by converting to DMDAs and dumping */
  ierr = DumpSolution(ctx,x);CHKERRQ(ierr);

  /* Free data */
  ierr = StagBLArrayDestroy(&x);CHKERRQ(ierr);
  ierr = StagBLSystemDestroy(&system);CHKERRQ(ierr);
  ierr = StagBLSolverDestroy(&solver);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&ctx->stokesGrid);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&ctx->coeffGrid);CHKERRQ(ierr);

  StagBLFinalize();
  return 0;
}

