static const char *help = "StagBLDemo2D: Demonstrate features and usage of StagBL, in 2 dimensions, with simple geodynamic box model setups\n\n";

#include "args.h"
#include "coeff.h"
#include "ctx.h"
#include "dump.h"
#include "temperature.h"
#include "particles.h"
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
  PetscBool              test_teleport;

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
  test_teleport = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-test_teleport",&test_teleport,NULL);CHKERRQ(ierr);

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
  // TODO this is too long for main and needs to move to a file
  {
    DM dmStokes;
    PetscInt particlesPerElementPerDim = 3,n[2],nel;

    ierr = PetscOptionsGetInt(NULL,NULL,"-p",&particlesPerElementPerDim,NULL);CHKERRQ(ierr);
    ierr = DMCreate(PETSC_COMM_WORLD,&ctx->dm_particles);CHKERRQ(ierr);
    ierr = DMSetType(ctx->dm_particles,DMSWARM);CHKERRQ(ierr);
    ierr = DMSetDimension(ctx->dm_particles,2);CHKERRQ(ierr);
    ierr = DMSwarmSetType(ctx->dm_particles,DMSWARM_PIC);CHKERRQ(ierr);
    ierr = StagBLGridPETScGetDM(ctx->stokes_grid,&dmStokes);CHKERRQ(ierr);
    ierr = DMStagGetLocalSizes(dmStokes,&n[0],&n[1],NULL);CHKERRQ(ierr);
    nel = n[0]*n[1];
    ierr = DMSwarmSetCellDM(ctx->dm_particles,dmStokes);CHKERRQ(ierr);
    // TODO actually carry things on the particles
    //ierr = DMSwarmRegisterPetscDatatypeField(ctx->dm_particles,"rho",1,PETSC_REAL);CHKERRQ(ierr);
    //ierr = DMSwarmRegisterPetscDatatypeField(ctx->dm_particles,"eta",1,PETSC_REAL);CHKERRQ(ierr);
    //ierr = DMSwarmRegisterPetscDatatypeField(ctx->dm_particles,"Temperature",1,PETSC_REAL);CHKERRQ(ierr);
    ierr = DMSwarmFinalizeFieldRegister(ctx->dm_particles);CHKERRQ(ierr);
    ierr = DMSwarmSetLocalSizes(ctx->dm_particles,nel*particlesPerElementPerDim,100);CHKERRQ(ierr);

    /* Define initial point positions */
    ierr = DMSwarmInsertPointsUsingCellDM(ctx->dm_particles,DMSWARMPIC_LAYOUT_REGULAR,particlesPerElementPerDim);CHKERRQ(ierr);

    /* Set properties for particles.  Each particle has a viscosity and density
       (in some applications one would use material ids ) */
    {
      PetscReal   *array_x,jiggle,dh;
      //PetscReal   *array_e,*array_r;
      PetscInt    npoints,p,N[2];
      PetscRandom r;
      PetscMPIInt rank;

      /* Displace particles randomly */
      jiggle = 0.1;
      ierr = PetscOptionsGetReal(NULL,NULL,"-jiggle",&jiggle,NULL);CHKERRQ(ierr);
      ierr = DMStagGetGlobalSizes(dmStokes,&N[0],&N[1],NULL);CHKERRQ(ierr);
      dh = PetscMin((ctx->xmax-ctx->xmin)/N[0],(ctx->ymax-ctx->ymin)/N[1]);
      if (jiggle > 0.0) {
        ierr = PetscRandomCreate(PETSC_COMM_SELF,&r);CHKERRQ(ierr);
        ierr = PetscRandomSetInterval(r,-jiggle*dh,jiggle*dh);CHKERRQ(ierr);
        ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dmStokes),&rank);CHKERRQ(ierr);
        ierr = PetscRandomSetSeed(r,(unsigned long)rank);CHKERRQ(ierr);
        ierr = PetscRandomSeed(r);CHKERRQ(ierr);
      }

      /* Compute coefficient values on particles */
      ierr = DMSwarmGetField(ctx->dm_particles,DMSwarmPICField_coor,NULL,NULL,(void**)&array_x);CHKERRQ(ierr);
      // TODO actually store data on particles
      //ierr = DMSwarmGetField(ctx->dm_particles,"eta",               NULL,NULL,(void**)&array_e);CHKERRQ(ierr);
      //ierr = DMSwarmGetField(ctx->dm_particles,"rho",               NULL,NULL,(void**)&array_r);CHKERRQ(ierr);
      ierr = DMSwarmGetLocalSize(ctx->dm_particles,&npoints);CHKERRQ(ierr);
      for (p = 0; p < npoints; p++) {
        PetscReal x_p[2];

        if (jiggle > 0.0) {
          PetscReal rr[2];
          ierr = PetscRandomGetValueReal(r,&rr[0]);CHKERRQ(ierr);
          ierr = PetscRandomGetValueReal(r,&rr[1]);CHKERRQ(ierr);
          array_x[2*p + 0] += rr[0];
          array_x[2*p + 1] += rr[1];
        }

        /* Get the coordinates of point p */
        x_p[0] = array_x[2*p + 0];
        x_p[1] = array_x[2*p + 1];

        /* Call functions to compute eta and rho at that location */
        // TODO actually carry things on the particles
        //array_e[p] = getEta(ctx,x_p[0]); // no y dependence
        //array_r[p] = getRho(ctx,x_p[0]); // no y dependence

      }
      // TODO actually store data on particles
      //ierr = DMSwarmRestoreField(ctx->dm_particles,"rho",NULL,NULL,(void**)&array_r);CHKERRQ(ierr);
      //ierr = DMSwarmRestoreField(ctx->dm_particles,"eta",NULL,NULL,(void**)&array_e);CHKERRQ(ierr);
      ierr = DMSwarmRestoreField(ctx->dm_particles,DMSwarmPICField_coor,NULL,NULL,(void**)&array_x);CHKERRQ(ierr);

      if (jiggle > 0.0) {
        ierr = PetscRandomDestroy(&r);CHKERRQ(ierr);
      }
    }

    /* Migrate - this is important since we may have "jiggled" particles
       into different cells, or off the local subdomain */
    ierr = DMSwarmMigrate(ctx->dm_particles,PETSC_TRUE);CHKERRQ(ierr);

    /* View DMSwarm object */
    ierr = DMView(ctx->dm_particles,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }


  /* Create parameters for a Stokes system by directly populating some fields
     of a struct and passing to a StagBL function  */
  ierr = StagBLStokesParametersCreate(&parameters);CHKERRQ(ierr);
  parameters->coefficient_array  = ctx->coefficient_array;
  parameters->stokes_grid        = ctx->stokes_grid;
  parameters->temperature_grid   = ctx->temperature_grid;
  parameters->temperature_array  = ctx->temperature_array;
  parameters->uniform_grid       = ctx->uniform_grid;
  parameters->xmin               = 0;
  parameters->xmax               = 1e6;
  parameters->ymin               = 0.0;
  parameters->ymax               = 1.5e6;
  parameters->gy                 = 10.0;
  parameters->alpha              = ctx->alpha;
  parameters->eta_characteristic = 1e20; /* A minimum viscosity */
  parameters->boussinesq_forcing = ctx->boussinesq_forcing;

  /* Main solver loop */
  timestep = 0;
  t = 0.0;
  do { /* Always take the 0th timestep */
    if (rank == 0) {
      printf("Timestep %d, starts at t=%g (~%g Myr)\n",timestep,(double)t,(double)t/3.154e13);
    }
    ierr = MPI_Barrier(comm);CHKERRQ(ierr);

    /* Initialize or Update Temperature Field */
    if (timestep == 0) {
      ierr = InitializeTemperature(ctx);CHKERRQ(ierr);
    } else {
      ierr = PopulateTemperatureSystem(ctx);CHKERRQ(ierr);
      ierr = UpdateTemperature(ctx);CHKERRQ(ierr);
    }

    /* Dump Temperature field */
    ierr = DumpTemperature(ctx,timestep);CHKERRQ(ierr);

    /* Dump Particles */
    // TODO put in dump.c
    {
      char filename[PETSC_MAX_PATH_LEN];

      ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"particles_%.4D.xmf",timestep);CHKERRQ(ierr);
      ierr = DMSwarmViewXDMF(ctx->dm_particles,filename);CHKERRQ(ierr);
    }

    /* Interpolate temperature from grid to particles */
    // TODO did this work befor?
    //ierr = InterpolateTemperatureToParticles(ctx);CHKERRQ(ierr);

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

    /* Dump Stokes solution */
    ierr = DumpStokes(ctx,timestep);CHKERRQ(ierr);

    /* Advect and Migrate */
    {
      Vec x;

      ierr = StagBLArrayPETScGetGlobalVec(ctx->stokes_array,&x);CHKERRQ(ierr);
      ierr = MaterialPoint_AdvectRK1(ctx,x,ctx->dt);CHKERRQ(ierr);
      ierr = DMSwarmMigrate(ctx->dm_particles,PETSC_TRUE);CHKERRQ(ierr);
    }

    /* Test particle teleportation */
    if (test_teleport) {
      ierr = TestTeleport(ctx);CHKERRQ(ierr);
      ierr = DMSwarmMigrate(ctx->dm_particles,PETSC_TRUE);CHKERRQ(ierr);
    }

    ++timestep;
    t += ctx->dt;
  } while (timestep <= ctx->totalTimesteps);

  /* Free data */
  ierr = StagBLStokesParametersDestroy(&parameters);CHKERRQ(ierr); // TODO put in ctx as well?
  ierr = CtxDestroy(&ctx);CHKERRQ(ierr);

  /* Finalize StagBL (which includes finalizing PETSc) */
  ierr = StagBLFinalize();
  return ierr;
}
