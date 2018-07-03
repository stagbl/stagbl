static char help[] = "Perform some standard operations to compare DMStag and a collection of DMDAs\n";

/* Proceeds by allowing the user to select one of a set of operations and whether to
   do it with DMStag or with (a collection of) DMDA(s) */

#include <petscdm.h>
#include <petscdmstag.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>

/* Supply one of these with -test  */
typedef enum {
  STAGTEST=0,
  DATEST1=1,
  DATEST2=2
} Test;

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DM             dm,dms[8];
  PetscInt       dof0,dof1,dof2,dof3,stencilWidth,Nx,Ny,Nz,i;
  PetscLogStage  creationStage,mainStage,destructionStage;
  PetscInt       test;

  /* Initialize and obtain parameters */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscLogStageRegister("Creation",&creationStage);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Main",&mainStage);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Destruction",&destructionStage);CHKERRQ(ierr);
  dof0 = 1;
  dof1 = 1;
  dof2 = 1;
  dof3 = 1;
  stencilWidth = 1;
  Nx = 30;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nx",&Nx,NULL);CHKERRQ(ierr);
  Ny = Nx;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Ny",&Ny,NULL);CHKERRQ(ierr);
  Nz = Nx;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nz",&Nz,NULL);CHKERRQ(ierr);
  test = STAGTEST;
  ierr = PetscOptionsGetInt(NULL,NULL,"-test",&test,NULL);CHKERRQ(ierr);

  /* STAGE : creation */
  // Note: we are only concerned with benchmarking well-balanced 3d problems for now
  ierr = PetscLogStagePush(creationStage);CHKERRQ(ierr);
  if (test == STAGTEST) {
    ierr = DMStagCreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,dof1,dof2,dof3,DMSTAG_GHOST_STENCIL_BOX,stencilWidth,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
  ierr = DMSetFromOptions(dm);CHKERRQ(ierr); // Use sparingly (for debugging)! Tests should compare the same functionality, remember.
  ierr = DMSetUp(dm);CHKERRQ(ierr);
  } else if (test == DATEST1) {
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0+3*dof1+3*dof2+dof3,stencilWidth,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
  ierr = DMSetFromOptions(dm);CHKERRQ(ierr); // Use sparingly (for debugging)! Tests should compare the same functionality, remember.
  ierr = DMSetUp(dm);CHKERRQ(ierr);
  } else if (test == DATEST2) {
    i = 0;
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,stencilWidth,NULL,NULL,NULL,&dms[i]);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dms[i]);CHKERRQ(ierr); // Use sparingly (for debugging)! Tests should compare the same functionality, remember.
    ierr = DMSetUp(dms[i]);CHKERRQ(ierr);
    for (i=1; i<4; ++i) {
      ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof1,stencilWidth,NULL,NULL,NULL,&dms[i]);CHKERRQ(ierr);
      ierr = DMSetFromOptions(dms[i]);CHKERRQ(ierr); // Use sparingly (for debugging)! Tests should compare the same functionality, remember.
      ierr = DMSetUp(dms[i]);CHKERRQ(ierr);
    }
    for (; i<7; ++i) {
      ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof2,stencilWidth,NULL,NULL,NULL,&dms[i]);CHKERRQ(ierr);
      ierr = DMSetFromOptions(dms[i]);CHKERRQ(ierr); // Use sparingly (for debugging)! Tests should compare the same functionality, remember.
      ierr = DMSetUp(dms[i]);CHKERRQ(ierr);
    }
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,stencilWidth,NULL,NULL,NULL,&dms[i]);CHKERRQ(ierr);
    ierr = DMSetFromOptions(dms[i]);CHKERRQ(ierr); // Use sparingly (for debugging)! Tests should compare the same functionality, remember.
    ierr = DMSetUp(dms[i]);CHKERRQ(ierr);
    ierr = DMCompositeCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
    for (i=0; i<8; ++i) {
      ierr = DMCompositeAddDM(dm,dms[i]);CHKERRQ(ierr);
    }
  } else SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported test %D",test);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* STAGE: main operation */
  ierr = PetscLogStagePush(mainStage);CHKERRQ(ierr);
  {
    Vec x;
    ierr = DMCreateGlobalVector(dm,&x);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);
  }
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* STAGE: destruction */
  ierr = PetscLogStagePush(destructionStage);CHKERRQ(ierr);
  if (test == DATEST2) {
    for (i=0; i<8; ++i) {
      ierr = DMDestroy(&dms[i]);CHKERRQ(ierr);
    }
  } else {
    ierr = DMDestroy(&dm);CHKERRQ(ierr);
  }
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* Finalize */
  ierr = PetscFinalize();
  return ierr;
}
