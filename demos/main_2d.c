static const char *help = "StagBLDemo2D: Demonstrate features and usage of StagBL, in 2 dimensions, with simple geodynamic box model setups\n\n";

#include "args.h"
#include "coeff.h"
#include "ctx.h"
#include "dump.h"
#include "temperature.h"
#include "particles.h"
#include "nusselt.h"
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
  PetscInt               timestep;
  PetscReal              t;
  PetscBool              test_magmatism;

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
    printf("=== StagBLDemo2d ===\n");
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
  ierr = GetStringArg("-mode","gerya72",sizeof(ctx->mode),mode);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Mode: %s\n",mode);CHKERRQ(ierr);

  // TODO clean up later
  test_magmatism = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-test_magmatism",&test_magmatism,NULL);CHKERRQ(ierr);

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

  /* Create another, compatible grid to represent coefficients */
  {
    const PetscInt dofPerVertex  = 2;
    const PetscInt dofPerElement = 1;
    ierr = StagBLGridCreateCompatibleStagBLGrid(ctx->stokes_grid,dofPerVertex,0,dofPerElement,0,&ctx->coefficient_grid);CHKERRQ(ierr);
  }

  /* Create another, compatible grid for the temperature field */
  ierr = StagBLGridCreateCompatibleStagBLGrid(ctx->stokes_grid,1,0,0,0,&ctx->temperature_grid);CHKERRQ(ierr);

  /* Coefficient data in an application-determined way */
  ierr = PopulateCoefficientData(ctx,mode);CHKERRQ(ierr);

  /* Create a Temperature system */
  ierr = StagBLGridCreateStagBLSystem(ctx->temperature_grid,&ctx->temperature_system);CHKERRQ(ierr);
  ierr = StagBLSystemCreateStagBLSolver(ctx->temperature_system,&ctx->temperature_solver);CHKERRQ(ierr);
  ierr = StagBLGridCreateStagBLArray(ctx->temperature_grid,&ctx->temperature_array);CHKERRQ(ierr);

  /* Create a System of Particles using PETSc DMSwarm object

     Note that StagBL does not force you to use any particular particle system,
     but this is a convenient one, as we have added code allow interaction
     between DMSwarm and DMStag, which was written as part of the StagBL project.
  */
  ierr = CreateParticleSystem(ctx);CHKERRQ(ierr);


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
  parameters->gy                 = ctx->gy;
  parameters->alpha              = ctx->alpha;
  parameters->eta_characteristic = 1e20; /* A minimum viscosity */ // TODO wrong
  parameters->boussinesq_forcing = ctx->boussinesq_forcing;

  // TODO print Rayleigh number?

  /* Main solver loop */
  timestep = 0;
  t = 0.0;
  do { /* Always take the 0th timestep */
    if (rank == 0) {
      printf("Timestep %d, starts at t=%g (~%g Myr)\n",timestep,(double)t,(double)t/3.154e13);
    }
    ierr = MPI_Barrier(comm);CHKERRQ(ierr);

    /* initialize or update temperature field */
    if (timestep == 0) {
      ierr = InitializeTemperature(ctx);CHKERRQ(ierr);
    } else {
      ierr = PopulateTemperatureSystem(ctx);CHKERRQ(ierr);
      ierr = UpdateTemperature(ctx);CHKERRQ(ierr);
    }

    /* Interpolate temperature from grid to particles */
    ierr = InterpolateTemperatureToParticles(ctx);CHKERRQ(ierr);

    /* (Re-)Populate Coefficient data */
    ierr = PopulateCoefficientData(ctx,mode);CHKERRQ(ierr);

    /* Create the Stokes system */
    /* Solve the system  */
    ierr = StagBLCreateStokesSystem(parameters,&ctx->stokes_system);CHKERRQ(ierr);
    ierr = StagBLSystemCreateStagBLSolver(ctx->stokes_system,&ctx->stokes_solver);CHKERRQ(ierr);
    ierr = StagBLGridCreateStagBLArray(ctx->stokes_grid,&ctx->stokes_array);CHKERRQ(ierr);
    ierr = StagBLSolverSolve(ctx->stokes_solver,ctx->stokes_array);CHKERRQ(ierr);
    ierr = StagBLSystemDestroy(&ctx->stokes_system);CHKERRQ(ierr);
    ierr = StagBLSolverDestroy(&ctx->stokes_solver);CHKERRQ(ierr); // TODO this sucks - we want to keep the solver and upate the system..

    /* Analyze temperature field */
    if (ctx->compute_nusselt_number) {
      PetscScalar nu;
      ierr = ComputeNusseltNumber(ctx,&nu);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Nusselt Number: %g\n",(double)nu);CHKERRQ(ierr);
    }

    /* Output stokes data, temperature data, and particle data to files */
    ierr = DumpStokes(ctx,timestep);CHKERRQ(ierr);
    ierr = DumpTemperature(ctx,timestep);CHKERRQ(ierr);
    ierr = DumpParticles(ctx,timestep);CHKERRQ(ierr);

    /* Advect and Migrate */
    {
      Vec x;

      ierr = StagBLArrayPETScGetGlobalVec(ctx->stokes_array,&x);CHKERRQ(ierr);
      ierr = MaterialPoint_AdvectRK1(ctx,x,ctx->dt);CHKERRQ(ierr);
      ierr = DMSwarmMigrate(ctx->dm_particles,PETSC_TRUE);CHKERRQ(ierr);
    }

    /* Test particle "teleportation" for use with magmatic processes */
    if (test_magmatism) {
      ierr = TestMagmatism(ctx);CHKERRQ(ierr);
      ierr = DMSwarmMigrate(ctx->dm_particles,PETSC_TRUE);CHKERRQ(ierr);
    }

    /* Update timestepping state */
    ++timestep;
    t += ctx->dt;
  } while (timestep <= ctx->totalTimesteps);

  /* Free data */
  ierr = StagBLStokesParametersDestroy(&parameters);CHKERRQ(ierr);
  ierr = CtxDestroy(&ctx);CHKERRQ(ierr);

  /* Finalize StagBL (which includes finalizing PETSc) */
  ierr = StagBLFinalize();
  return ierr;
}
