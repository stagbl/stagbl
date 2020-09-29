#include <stagbl.h>

int main(int argc, char** argv)
{
  PetscErrorCode         ierr;
  StagBLGrid             stokes_grid,coefficient_grid;
  StagBLArray            stokes_array,coefficient_array;
  StagBLStokesParameters parameters;
  StagBLSystem           stokes_system;
  PetscReal              xmin,xmax,ymin,ymax,zmin,zmax;
  PetscInt               dim,Nx,Ny,Nz;
  PetscBool              use_simple;

  ierr = StagBLInitialize(argc,argv,NULL,MPI_COMM_NULL);CHKERRQ(ierr);

  xmin = 0.0; xmax = 3.0;
  ymin = 0.0; ymax = 3.0;
  zmin = 0.0; zmax = 3.0;
  dim = 2;
  ierr = PetscOptionsGetInt(NULL,NULL,"-dim",&dim,NULL);CHKERRQ(ierr);
  Nx = 5;
  Ny = 5;
  Nz = 5;
  if (dim == 2) {
    ierr = StagBLGridCreateStokes2DBox(PETSC_COMM_WORLD,Nx,Ny,xmin,xmax,ymin,ymax,&stokes_grid);CHKERRQ(ierr);
  } else if (dim == 3) {
    ierr = StagBLGridCreateStokes3DBox(PETSC_COMM_WORLD,Nx,Ny,Nz,xmin,xmax,ymin,ymax,zmin,zmax,&stokes_grid);CHKERRQ(ierr);
  } else StagBLError1(PETSC_COMM_WORLD,"unsupported dimension %D",dim);
  use_simple = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-simple",&use_simple,NULL);CHKERRQ(ierr);
  if (use_simple) {
    ierr = StagBLGridSetArrayType(stokes_grid,STAGBLARRAYSIMPLE);CHKERRQ(ierr);
    ierr = StagBLGridSetSystemType(stokes_grid,STAGBLSYSTEMSIMPLE);CHKERRQ(ierr);
  } else {
    ierr = StagBLGridSetArrayType(stokes_grid,STAGBLARRAYPETSC);CHKERRQ(ierr);
    ierr = StagBLGridSetSystemType(stokes_grid,STAGBLSYSTEMPETSC);CHKERRQ(ierr);
  }
  if (dim == 2) {
    ierr = StagBLGridCreateCompatibleStagBLGrid(stokes_grid,2,0,1,0,&coefficient_grid);CHKERRQ(ierr);
  } else if (dim == 3) {
    ierr = StagBLGridCreateCompatibleStagBLGrid(stokes_grid,0,2,0,1,&coefficient_grid);CHKERRQ(ierr);
  } else StagBLError1(PETSC_COMM_WORLD,"unsupported dimension %D",dim);
  ierr = StagBLGridCreateStagBLArray(coefficient_grid,&coefficient_array);CHKERRQ(ierr);

  ierr = StagBLArraySetLocalConstant(coefficient_array,1.0);CHKERRQ(ierr);

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

  ierr = StagBLSystemSolve(stokes_system,stokes_array);CHKERRQ(ierr);

  ierr = StagBLDumpStokes(parameters,stokes_array,0);CHKERRQ(ierr);

  ierr = StagBLArrayDestroy(&stokes_array);CHKERRQ(ierr);
  ierr = StagBLSystemDestroy(&stokes_system);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&stokes_grid);CHKERRQ(ierr);
  ierr = StagBLArrayDestroy(&coefficient_array);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&coefficient_grid);CHKERRQ(ierr);
  ierr = StagBLStokesParametersDestroy(&parameters);CHKERRQ(ierr);
  ierr = StagBLFinalize();
  return ierr;
}
