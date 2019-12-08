static const char *help = "StagBLDemo3D: Demonstrate features and usage of StagBL, in 3 dimensions, with simple geodynamic box model setups\n\n";

#include "args.h"
#include "coeff.h"
#include "ctx.h"
#include "dump.h"
#include <stagbl.h>
#include <mpi.h>
#include <stdio.h>

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
  int                    rank;
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
  MPI_Comm_rank(comm,&rank);
  if (rank == 0) {
    printf("=== StagBLDemo3d ===\n");
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
  ierr = GetStringArg("-mode","sinker",sizeof(ctx->mode),mode);CHKERRQ(ierr);
  // TODO sinker currently doesn't use z coord, so it's a cylinder!

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Mode: %s\n",mode);CHKERRQ(ierr);

  /* Populate application context (Create with a given mode, then set up) */
  ierr = CtxCreate(comm,mode,&ctx);CHKERRQ(ierr);

  /* Create a Grid
     We call a helper function to create an interlaced p-v grid of the correct size.
     Note that the choice of whether to create a single grid with both pressure and
     velocity, or two grids, must be made early, as this affects data layout. This
     single-grid choice is appropriate for monolithic mulitigrid or direct solution,
     whereas two grids would be appropriated for segregated or Approximate block
     factorization-based solvers. */
  ierr = StagBLGridCreateStokes3DBox(comm,30,20,20,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,ctx->zmin,ctx->zmax,&ctx->stokes_grid);CHKERRQ(ierr);

  /* Create another, compatible grid to represent coefficients */
  {
    const PetscInt dofPerEdge  = 2;
    const PetscInt dofPerElement = 1;
    ierr = StagBLGridCreateCompatibleStagBLGrid(ctx->stokes_grid,0,dofPerEdge,0,dofPerElement,&ctx->coefficient_grid);CHKERRQ(ierr);
  }

  /* Populate Coefficient data */
  ierr = PopulateCoefficientData(ctx,mode);CHKERRQ(ierr);

  /* Create parameters for a Stokes system by directly populating some fields
     of a struct and passing to a StagBL function  */
  ierr = StagBLStokesParametersCreate(&parameters);CHKERRQ(ierr);
  parameters->coefficient_array  = ctx->coefficient_array;
  parameters->stokes_grid        = ctx->stokes_grid;
  parameters->temperature_grid   = ctx->temperature_grid;
  parameters->temperature_array  = ctx->temperature_array;
  parameters->uniform_grid       = ctx->uniform_grid;
  parameters->xmin               = ctx->xmin;
  parameters->xmax               = ctx->xmax;
  parameters->ymin               = ctx->ymin;
  parameters->ymax               = ctx->ymax;
  parameters->zmin               = ctx->zmin;
  parameters->zmax               = ctx->zmax;
  parameters->gy                 = ctx->gy;
  parameters->alpha              = ctx->alpha;
  parameters->eta_characteristic = ctx->eta_characteristic;
  parameters->boussinesq_forcing = ctx->boussinesq_forcing;

  /* Create the Stokes system */
  ierr = StagBLGridCreateStagBLArray(ctx->stokes_grid,&ctx->stokes_array);CHKERRQ(ierr);
  ierr = StagBLCreateStokesSystem(parameters,&ctx->stokes_system);CHKERRQ(ierr);
#if 0
  ierr = StagBLSystemCreateStagBLSolver(ctx->stokes_system,&ctx->stokes_solver);CHKERRQ(ierr);

  /* Solve the system  */
  ierr = StagBLSolverSolve(ctx->stokes_solver,ctx->stokes_array);CHKERRQ(ierr);
  ierr = StagBLSystemDestroy(&ctx->stokes_system);CHKERRQ(ierr);
  ierr = StagBLSolverDestroy(&ctx->stokes_solver);CHKERRQ(ierr);
#endif

  /* Output Stokes data to file */
  ierr = DumpStokes(ctx,/* timestep */ 0);CHKERRQ(ierr);

  /* Free data */
  ierr = StagBLStokesParametersDestroy(&parameters);CHKERRQ(ierr);
  ierr = CtxDestroy(&ctx);CHKERRQ(ierr);

  /* Finalize StagBL (which includes finalizing PETSc) */
  ierr = StagBLFinalize();
  return ierr;
}
