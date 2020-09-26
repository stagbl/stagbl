#include <stagbl.h>

int main(int argc, char** argv)
{
  PetscErrorCode         ierr;
  PetscBool              use_simple;
  StagBLGrid             grid;
  StagBLSystem           system;
  StagBLSystemType       system_type;

  ierr = StagBLInitialize(argc,argv,NULL,MPI_COMM_NULL);CHKERRQ(ierr);

  ierr = StagBLGridCreateStokes3DBox(PETSC_COMM_WORLD,3,3,3,0.0,1.0,0.0,1.0,0.0,1.0,&grid);CHKERRQ(ierr);

  use_simple = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-simple",&use_simple,NULL);CHKERRQ(ierr);
  if (use_simple) {
    system_type = STAGBLSYSTEMSIMPLE;
  } else {
    system_type = STAGBLSYSTEMPETSC;
  }
  ierr = StagBLGridSetSystemType(grid,system_type);CHKERRQ(ierr);
  ierr = StagBLGridCreateStagBLSystem(grid,&system);CHKERRQ(ierr);

  {
    const PetscInt nrows = 2;
    const PetscInt ncols = 3;
    DMStagStencil  rows[2];
    DMStagStencil  cols[3];
    PetscScalar    values[6];
    PetscScalar    rhs_values[3];

    rows[0].i   = 1;
    rows[0].j   = 2;
    rows[0].k   = 0;
    rows[0].c   = 0;
    rows[0].loc = DMSTAG_FRONT;

    rows[1].i   = 1;
    rows[1].j   = 2;
    rows[1].k   = 1;
    rows[1].c   = 0;
    rows[1].loc = DMSTAG_BACK;

    cols[0].i   = 1;
    cols[0].j   = 2;
    cols[0].k   = 0;
    cols[0].c   = 0;
    cols[0].loc = DMSTAG_ELEMENT;

    cols[1].i   = 1;
    cols[1].j   = 2;
    cols[1].k   = 0;
    cols[1].c   = 0;
    cols[1].loc = DMSTAG_LEFT;

    cols[2].i   = 1;
    cols[2].j   = 0;
    cols[2].k   = 1;
    cols[2].c   = 0;
    cols[2].loc = DMSTAG_ELEMENT;

    values[0]   = 1.2345;
    values[1]   = 1.2345;
    values[2]   = 1.2345;
    values[3]   = 1.2345;
    values[4]   = 1.2345;
    values[5]   = 1.2345;

    rhs_values[0] = -5.4321;
    rhs_values[1] = -6.4321;
    rhs_values[2] = -7.4321;

    ierr = StagBLSystemRHSSetValuesStencil(system,nrows,rows,rhs_values);CHKERRQ(ierr);
    ierr = StagBLSystemOperatorSetValuesStencil(system,nrows,rows,ncols,cols,values);CHKERRQ(ierr);
  }

  ierr = StagBLSystemDestroy(&system);CHKERRQ(ierr);
  ierr = StagBLGridDestroy(&grid);CHKERRQ(ierr);
  ierr = StagBLFinalize();
  return ierr;
}
