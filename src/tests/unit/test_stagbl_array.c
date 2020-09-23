#include <stagbl.h>

int main(int argc, char** argv)
{
  PetscErrorCode  ierr;
  StagBLGrid      grid;
  StagBLArray     array;
  StagBLArrayType array_type;
  PetscInt        N[3];
  PetscReal       xmin,xmax,ymin,ymax,zmin,zmax;

  ierr = StagBLInitialize(argc,argv,NULL,MPI_COMM_NULL);CHKERRQ(ierr);

  xmin = 0.0; xmax = 3.0;
  ymin = 0.0; ymax = 3.0;
  zmin = 0.0; zmax = 3.0;
  N[0] = 2;
  N[1] = 3;
  N[2] = 2;
  ierr = StagBLGridCreateStokes3DBox(PETSC_COMM_WORLD,N[0],N[1],N[2],xmin,xmax,ymin,ymax,zmin,zmax,&grid);CHKERRQ(ierr);

  array_type = "simple";
  //array_type = "petsc"; // TODO control with flag
  ierr = StagBLGridSetArrayType(grid,array_type);CHKERRQ(ierr);
  ierr = StagBLGridCreateStagBLArray(grid,&array);CHKERRQ(ierr);

  ierr = StagBLArrayPrint(array);CHKERRQ(ierr);  // Should print nothing

  ierr = StagBLArraySetLocalConstant(array, 3.0);CHKERRQ(ierr);

  ierr = StagBLArrayLocalToGlobal(array);CHKERRQ(ierr);

  ierr = StagBLArrayGlobalToLocal(array);CHKERRQ(ierr);

  ierr = StagBLArrayPrint(array);CHKERRQ(ierr);


  ierr = StagBLArrayDestroy(&array);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&grid);CHKERRQ(ierr);
  ierr = StagBLFinalize();
  return ierr;
}
