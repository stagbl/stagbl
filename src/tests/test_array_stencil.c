#include <stagbl.h>

int main(int argc, char** argv)
{
  PetscErrorCode  ierr;
  StagBLGrid      grid;
  StagBLArray     array;
  StagBLArrayType array_type;
  PetscInt        N[3];
  PetscReal       xmin,xmax,ymin,ymax,zmin,zmax;
  PetscBool       use_simple;

  ierr = StagBLInitialize(argc,argv,NULL,MPI_COMM_NULL);CHKERRQ(ierr);

  xmin = 0.0; xmax = 3.0;
  ymin = 0.0; ymax = 3.0;
  zmin = 0.0; zmax = 3.0;
  N[0] = 2;
  N[1] = 3;
  N[2] = 2;
  ierr = StagBLGridCreateStokes3DBox(PETSC_COMM_WORLD,N[0],N[1],N[2],xmin,xmax,ymin,ymax,zmin,zmax,&grid);CHKERRQ(ierr);

  use_simple = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-simple",&use_simple,NULL);CHKERRQ(ierr);
  if (use_simple) {
    array_type = STAGBLARRAYSIMPLE;
  } else {
    array_type = STAGBLARRAYPETSC;
  }
  ierr = StagBLGridSetArrayType(grid,array_type);CHKERRQ(ierr);
  ierr = StagBLGridCreateStagBLArray(grid,&array);CHKERRQ(ierr);

  ierr = StagBLArrayPrint(array);CHKERRQ(ierr);  // Should print nothing

  ierr = StagBLArraySetLocalConstant(array, 4.567);CHKERRQ(ierr);

  ierr = StagBLArrayLocalToGlobal(array);CHKERRQ(ierr);

  {
    DMStagStencil  rows[2];
    const PetscInt nrows = 2;
    PetscScalar    values[2];

    rows[0].i = 1;
    rows[0].j = 2;
    rows[0].k = 0;
    rows[0].c = 0;
    rows[0].loc = DMSTAG_FRONT;

    rows[1].i   = 1;
    rows[1].j   = 2;
    rows[1].k   = 1;
    rows[1].c   = 0;
    rows[1].loc = DMSTAG_BACK;


    values[0] = 1.2345;
    values[0] = 2.2345;

    ierr = StagBLArraySetLocalValuesStencil(array,nrows,rows,values);CHKERRQ(ierr);
  }

  ierr = StagBLArrayGlobalToLocal(array);CHKERRQ(ierr);

  {
    DMStagStencil  rows[2];
    const PetscInt nrows = 2;
    PetscScalar    values[2];

    rows[0].i = 1;
    rows[0].j = 2;
    rows[0].k = 0;
    rows[0].c = 0;
    rows[0].loc = DMSTAG_FRONT;

    rows[1].i   = 1;
    rows[1].j   = 2;
    rows[1].k   = 1;
    rows[1].c   = 0;
    rows[1].loc = DMSTAG_BACK;

    ierr = StagBLArrayGetLocalValuesStencil(array,nrows,rows,values);CHKERRQ(ierr);
    {
      PetscMPIInt rank;

      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
      for (PetscInt i=0; i<nrows; ++i) {
        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] value[%d]: %g\n",rank,i,values[i]);
      }
      ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);CHKERRQ(ierr);
    }
  }

  ierr = StagBLArrayPrint(array);CHKERRQ(ierr);

  ierr = StagBLArrayDestroy(&array);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&grid);CHKERRQ(ierr);
  ierr = StagBLFinalize();
  return ierr;
}
