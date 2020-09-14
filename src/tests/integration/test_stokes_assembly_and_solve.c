#include <stagbl.h>

int main(int argc, char** argv)
{
  PetscErrorCode         ierr;
  StagBLGrid             stokes_grid,coefficient_grid;
  StagBLArray            stokes_array,coefficient_array;
  StagBLStokesParameters parameters;
  StagBLSystem           stokes_system;
  StagBLSolver           stokes_solver;
  PetscReal              xmin,xmax,ymin,ymax,zmin,zmax;

  ierr = StagBLInitialize(argc,argv,NULL,MPI_COMM_NULL);CHKERRQ(ierr);

  xmin = 0.0; xmax = 3.0;
  ymin = 0.0; ymax = 3.0;
  zmin = 0.0; zmax = 3.0;
  ierr = StagBLGridCreateStokes3DBox(PETSC_COMM_WORLD,3,3,3,xmin,xmax,ymin,ymax,zmin,zmax,&stokes_grid);CHKERRQ(ierr);
  ierr = StagBLGridCreateCompatibleStagBLGrid(stokes_grid,0,2,0,1,&coefficient_grid);CHKERRQ(ierr);

  // FIXME this depends on non-DM PETSc! (just wrap in a helper)
  {
    DM  dm_coefficients;
    Vec *p_coefficients_local;

    ierr = StagBLGridCreateStagBLArray(coefficient_grid,&coefficient_array);CHKERRQ(ierr);
    ierr = StagBLArrayPETScGetLocalVecPointer(coefficient_array,&p_coefficients_local);CHKERRQ(ierr);
    ierr = StagBLGridPETScGetDM(coefficient_grid,&dm_coefficients);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(dm_coefficients,p_coefficients_local);CHKERRQ(ierr);
    ierr = VecSet(*p_coefficients_local,1.0);CHKERRQ(ierr);
  }

  ierr = StagBLStokesParametersCreate(&parameters);CHKERRQ(ierr);
  parameters->coefficient_array  = coefficient_array;
  parameters->stokes_grid        = stokes_grid;
  parameters->temperature_grid   = NULL;
  parameters->temperature_array  = NULL;
  parameters->uniform_grid       = PETSC_TRUE;
  parameters->xmin               = xmin;
  parameters->xmax               = xmax;
  parameters->ymin               = ymin;
  parameters->ymax               = ymax;
  parameters->zmin               = zmin;
  parameters->zmax               = zmax;
  parameters->gy                 = 1.0;
  parameters->eta_characteristic = 1.0;

  ierr = StagBLGridCreateStagBLArray(stokes_grid,&stokes_array);CHKERRQ(ierr);
  ierr = StagBLCreateStokesSystem(parameters,&stokes_system);CHKERRQ(ierr);
  ierr = StagBLSystemCreateStagBLSolver(stokes_system,&stokes_solver);CHKERRQ(ierr);

  ierr = StagBLSolverSolve(stokes_solver,stokes_array);CHKERRQ(ierr);

  // FIXME visualize or otherwise interpret result here

  ierr = StagBLDumpStokes(parameters,stokes_array,0);CHKERRQ(ierr);

  ierr = StagBLArrayDestroy(&stokes_array);CHKERRQ(ierr);
  ierr = StagBLSystemDestroy(&stokes_system);CHKERRQ(ierr);
  ierr = StagBLSolverDestroy(&stokes_solver);CHKERRQ(ierr);
  ierr = StagBLStokesParametersDestroy(&parameters);CHKERRQ(ierr);
  ierr = StagBLFinalize();
  return ierr;
}
