static const char *help = "StagBLDemo2D: Demonstrate features and usage of StagBL, in 2 dimensions, with simple geodynamic box model setups\n\n";

#include "args.h"
#include "coeff.h"
#include "ctx.h"
#include "dump.h"
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
  PetscErrorCode         ierr;
  int                    rank,size;
  StagBLArray            x;
  StagBLSystem           system;
  StagBLSolver           solver;
  MPI_Comm               comm;
  char                   mode[1024];
  Ctx                    ctx;
  StagBLStokesParameters parameters;

  /* Initialize MPI and print a message

     This is not required: one can simply call StagBLInitialize,
     which will also initialize MPI and allow one to work
     with PETSC_COMM_WORLD. We include this logic here to
     demonstrate how one can work with StagBL as a library within
     a larger application which already uses MPI.  */
  MPI_Init(&argc,&argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  if (rank == 0) {
    printf("=== StagBLDemo2d ===\n");
    printf("%d ranks\n",size);
  }
  fflush(stdout);
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

  /* Populate application context (Create with a given mode, then set up) */
  ierr = CtxCreate(comm,mode,&ctx);CHKERRQ(ierr);

  /* Create a Grid
     We call a helper function to create an interlaced p-v grid of the correct size.
     Note that the choice of whether to create a single grid with both pressure and
     velocity, or two grids, must be made early, as this affects data layout. This
     single-grid choice is appropriate for monolithic mulitigrid or direct solution,
     whereas two grids would be appropriated for segregated or Approximate block
     factorization-based solvers. */
  ierr = StagBLGridCreateStokes2DBox(comm,30,20,0.0,ctx->xmax,0.0,ctx->ymax,&ctx->stokes_grid);CHKERRQ(ierr);

  /* Get scaling constants and node to pin, knowing grid dimensions */
  ierr = CtxSetupFromGrid(ctx);CHKERRQ(ierr);

  /* Create another, compatible grid to represent coefficients */
  {
    const PetscInt dofPerVertex  = 2;
    const PetscInt dofPerElement = 1;
    ierr = StagBLGridCreateCompatibleStagBLGrid(ctx->stokes_grid,dofPerVertex,0,dofPerElement,0,&ctx->coefficient_grid);CHKERRQ(ierr);
  }

  /* Coefficient data in an application-determined way */
  ierr = PopulateCoefficientData(ctx,mode);CHKERRQ(ierr);

  /* Create a simple Stokes system by directly populating some fields
     of a struct and passing to a StagBL function */
  ierr = StagBLStokesParametersCreate(&parameters);CHKERRQ(ierr);
  parameters->coefficient_array  = ctx->coefficient_array;
  parameters->stokes_grid        = ctx->stokes_grid;
  parameters->uniform_grid       = PETSC_TRUE;
  parameters->xmin               = 0;
  parameters->xmax               = 1e6;
  parameters->ymin               = 0.0;
  parameters->ymax               = 1.5e6;
  parameters->gy                 = 10.0;
  parameters->eta_characteristic = 1e20; /* A minimum viscosity */
  ierr = StagBLCreateStokesSystem(parameters,&system);CHKERRQ(ierr);
  ierr = StagBLStokesParametersDestroy(&parameters);CHKERRQ(ierr);

  /* Solve the system  */
  ierr = StagBLSystemCreateStagBLSolver(system,&solver);CHKERRQ(ierr);
  ierr = StagBLGridCreateStagBLArray(ctx->stokes_grid,&x);CHKERRQ(ierr);
  ierr = StagBLSolverSolve(solver,x);CHKERRQ(ierr);

  /* Dump solution by converting to DMDAs and dumping */
  ierr = DumpSolution(ctx,x);CHKERRQ(ierr);

  /* Free data */
  ierr = StagBLArrayDestroy(&x);CHKERRQ(ierr);
  ierr = StagBLArrayDestroy(&ctx->coefficient_array);CHKERRQ(ierr);
  ierr = StagBLSystemDestroy(&system);CHKERRQ(ierr);
  ierr = StagBLSolverDestroy(&solver);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&ctx->stokes_grid);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&ctx->coefficient_grid);CHKERRQ(ierr);
  ierr = CtxDestroy(&ctx);CHKERRQ(ierr);

  /* Finalize StagBL (which includes finalizing PETSc) */
  ierr = StagBLFinalize();
  return ierr;
}
