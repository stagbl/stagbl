static char help[] = "Perform some standard operations to compare DMStag and a collection of DMDAs\n";

/* Proceeds by allowing the user to select one of a set of operations and whether to
   do it with DMStag or with (a collection of) DMDA(s)

   Note that, since this is a performance test, there is no output to stdout, by default
   (Use -log_summary).

   Note that we don't call DMSetFromOptions() here, so DMStag/DMDA-specific command
   line options won't work.  */

#include <petscdm.h>
#include <petscdmstag.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>

/* Supply one of these with -test  */
typedef enum {
  STAGTEST = 0, /* Use DMStag */
  DATEST1  = 1, /* Use DMDA with uniform blocksize */
  DATEST2  = 2  /* Use 8 DMDAs */
} Test;

/* Supply one of these with -op */
typedef enum {
  OP_NONE              = 0,
  OP_CREATE_GLOBAL_VEC = 1,
  OP_LOCAL_TO_GLOBAL   = 2,
  OP_GLOBAL_TO_LOCAL   = 3,
} Op;

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DM             dm,dms[8];
  PetscInt       dof0,dof1,dof2,dof3,stencilWidth,Nx,Ny,Nz,i;
  PetscLogStage  creationStage,mainStage,destructionStage;
  PetscInt       test,op;

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
  Nx = 100;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nx",&Nx,NULL);CHKERRQ(ierr);
  Ny = Nx;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Ny",&Ny,NULL);CHKERRQ(ierr);
  Nz = Nx;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nz",&Nz,NULL);CHKERRQ(ierr);
  test = STAGTEST;
  ierr = PetscOptionsGetInt(NULL,NULL,"-test",&test,NULL);CHKERRQ(ierr);
  op = OP_NONE;
  ierr = PetscOptionsGetInt(NULL,NULL,"-op",&op,NULL);CHKERRQ(ierr);

  /* STAGE : creation */
  ierr = PetscLogStagePush(creationStage);CHKERRQ(ierr);
  if (test == STAGTEST) {
    ierr = DMStagCreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,dof1,dof2,dof3,DMSTAG_STENCIL_BOX,stencilWidth,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
    ierr = DMSetUp(dm);CHKERRQ(ierr);
  } else if (test == DATEST1) {
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0+3*dof1+3*dof2+dof3,stencilWidth,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
    ierr = DMSetUp(dm);CHKERRQ(ierr);
  } else if (test == DATEST2) {
    i = 0;
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,stencilWidth,NULL,NULL,NULL,&dms[i]);CHKERRQ(ierr);
    ierr = DMSetUp(dms[i]);CHKERRQ(ierr);
    for (i=1; i<4; ++i) {
      ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof1,stencilWidth,NULL,NULL,NULL,&dms[i]);CHKERRQ(ierr);
      ierr = DMSetUp(dms[i]);CHKERRQ(ierr);
    }
    for (; i<7; ++i) {
      ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof2,stencilWidth,NULL,NULL,NULL,&dms[i]);CHKERRQ(ierr);
      ierr = DMSetUp(dms[i]);CHKERRQ(ierr);
    }
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,stencilWidth,NULL,NULL,NULL,&dms[i]);CHKERRQ(ierr);
    ierr = DMSetUp(dms[i]);CHKERRQ(ierr);
    ierr = DMCompositeCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
    for (i=0; i<8; ++i) {
      ierr = DMCompositeAddDM(dm,dms[i]);CHKERRQ(ierr);
    }
  } else SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported test %D",test);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* STAGE: main operation */
  /* Note: the log stage does NOT necessarily cover everything here,
           e.g. if timing scatters, we don't time creating the vectors */
  if (op == OP_NONE) {
    ierr = PetscLogStagePush(mainStage);CHKERRQ(ierr);
    /* Empty */
    ierr = PetscLogStagePop();CHKERRQ(ierr);
  } else if (op == OP_CREATE_GLOBAL_VEC) {
    Vec x;
    ierr = PetscLogStagePush(mainStage);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm,&x);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = PetscLogStagePop();CHKERRQ(ierr);
  } else if (op == OP_LOCAL_TO_GLOBAL) {
      Vec local,global;
      ierr = DMCreateGlobalVector(dm,&global);CHKERRQ(ierr);
      ierr = DMCreateLocalVector(dm,&local);CHKERRQ(ierr);
      ierr = PetscLogStagePush(mainStage);CHKERRQ(ierr);
      ierr = DMLocalToGlobalBegin(dm,local,INSERT_VALUES,global);CHKERRQ(ierr);
      ierr = DMLocalToGlobalEnd(dm,local,INSERT_VALUES,global);CHKERRQ(ierr);
      ierr = PetscLogStagePop();CHKERRQ(ierr);
      ierr = VecDestroy(&global);CHKERRQ(ierr);
      ierr = VecDestroy(&local);CHKERRQ(ierr);
  } else if (op == OP_GLOBAL_TO_LOCAL) {
      Vec local,global;
      ierr = DMCreateGlobalVector(dm,&global);CHKERRQ(ierr);
      ierr = DMCreateLocalVector(dm,&local);CHKERRQ(ierr);
      ierr = PetscLogStagePush(mainStage);CHKERRQ(ierr);
      ierr = DMGlobalToLocalBegin(dm,global,INSERT_VALUES,local);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dm,global,INSERT_VALUES,local);CHKERRQ(ierr);
      ierr = PetscLogStagePop();CHKERRQ(ierr);
      ierr = VecDestroy(&global);CHKERRQ(ierr);
      ierr = VecDestroy(&local);CHKERRQ(ierr);
  } else SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported op %D",op);CHKERRQ(ierr);

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
