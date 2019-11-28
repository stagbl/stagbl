static const char *help = "StagBLDemo2D: Demonstrate features and usage of StagBL, in 2 dimensions, with simple geodynamic box model setups\n\n";

#include "args.h"
#include "coeff.h"
#include "ctx.h"
#include "dump.h"
#include "system.h"
#include "system2.h"
#include <stagbl.h>
#include <stdio.h>
#include <mpi.h>

/* Note: This demonstration code is written using the PETSc library
    to provide various control, options processing, printing, and
    other operations. This is convenient, since StagBL itself relies
    on PETSc and we can thus reuse the configuration. For those
    developing new codes or integrating StagBL into an existing code,
    it is important to realize that an application may use other
    tools for these operations, while still using StagBL for the
    discretization and solver "base layer". */

int main(int argc, char** argv)

{
  PetscErrorCode ierr;
  int             rank,size;
  StagBLArray     x;
  StagBLSystem    system;
  StagBLSolver    solver;
  MPI_Comm        comm;
  char            mode[1024];
  PetscInt       systemtype;
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

  /* Initialize StagBL (which will initialize PETSc, and MPI if not initialized) */
  ierr = StagBLInitialize(argc,argv,help,comm);CHKERRQ(ierr);

  /* Accept an argument for the "mode". StagBLDemo2D is not intended
     to be a full application, but rather a demonstration, test case,
     and mini-app. As such, it is assumed that a small number of
     simple setups (often known benchmarks) will be used without
     major modification, and most parameters are set by choosing a
     mode. Other options can of course modify some of these */
  ierr = GetStringArg("-mode","gerya72",sizeof(mode),mode);CHKERRQ(ierr);

  /* Accept argument for system type and coefficient structure */
  // TODO : make a special gerya72_repro mode which uses the other system...
  systemtype = 1;

  /* Populate application context (Create with a given mode, then set up) */
  ierr = CtxCreate(comm,mode,&ctx);CHKERRQ(ierr);

  /* Create a Grid
     We call a helper function to create an interlaced p-v grid of the correct size.
     Note that the choice of whether to create a single grid with both pressure and
     velocity, or two grids, must be made early, as this affects data layout. This
     single-grid choice is appropriate for monolithic mulitigrid or direct solution,
     whereas two grids would be appropriated for segregated or Approximate block
     factorization-based solvers. */
  ierr = StagBLGridCreateStokes2DBox(comm,30,20,0.0,ctx->xmax,0.0,ctx->ymax,&ctx->stokesGrid);CHKERRQ(ierr);

  /* Get scaling constants and node to pin, knowing grid dimensions */
  ierr = CtxSetupFromGrid(ctx);CHKERRQ(ierr);

  /* Create another, compatible grid to represent coefficients */
 {
   const PetscInt dofPerVertex  = 2;
   const PetscInt dofPerElement = 1;
   ierr = StagBLGridCreateCompatibleStagBLGrid(ctx->stokesGrid,dofPerVertex,0,dofPerElement,0,&ctx->coeffGrid);CHKERRQ(ierr);
 }

  /* Coefficient data */
  ierr = PopulateCoefficientData(ctx,mode);CHKERRQ(ierr);

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
