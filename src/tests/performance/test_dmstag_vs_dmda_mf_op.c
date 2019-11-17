static char help[] = "For a toy matrix-free operator, compare DMStag and a collection of DMDAs\n";

/* Note that, since this is a performance test, there is no output to stdout, by default
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

#define NUM_DAS 8

PetscErrorCode ApplyStag(DM,Vec,Vec);
PetscErrorCode ApplyDA1(DM,Vec,Vec);
//PetscErrorCode ApplyDA2(DM,Vec,Vec);

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DM             dm,dms[4];
  PetscInt       dof0,dof1,dof2,dof3,stencilWidth,Nx,Ny,Nz,i;
  PetscLogStage  creationStage,mainStage,destructionStage;
  PetscInt       test;
  Vec            in,out;

  /* Initialize and obtain parameters */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscLogStageRegister("Creation",&creationStage);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Main",&mainStage);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Destruction",&destructionStage);CHKERRQ(ierr);
  dof0 = 0;
  dof1 = 0;
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

  /* STAGE : creation */
  ierr = PetscLogStagePush(creationStage);CHKERRQ(ierr);
  if (test == STAGTEST) {
    ierr = DMStagCreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,dof1,dof2,dof3,DMSTAG_STENCIL_BOX,stencilWidth,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
    ierr = DMSetUp(dm);CHKERRQ(ierr);
  } else if (test == DATEST1) {
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0+3*dof1+3*dof2+dof3,stencilWidth,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
    ierr = DMSetUp(dm);CHKERRQ(ierr);
  } else if (test == DATEST2) {
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,stencilWidth,NULL,NULL,NULL,&dms[0]);CHKERRQ(ierr);
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx+1,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,stencilWidth,NULL,NULL,NULL,&dms[1]);CHKERRQ(ierr);
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny+1,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,stencilWidth,NULL,NULL,NULL,&dms[2]);CHKERRQ(ierr);
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz+1,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,stencilWidth,NULL,NULL,NULL,&dms[3]);CHKERRQ(ierr);
    for (i=0; i<NUM_DAS; ++i) {
      ierr = DMSetUp(dms[i]);CHKERRQ(ierr);
    }
    ierr = DMCompositeCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
    for (i=0; i<NUM_DAS; ++i) {
      ierr = DMCompositeAddDM(dm,dms[i]);CHKERRQ(ierr);
    }
  } else SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported test %D",test);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* STAGE: main operation */
  /* Note that the log stage does NOT include creating and destroying the vectors
     and that the usual global<-->local transfers are also not included  (both of
     these should be individually testable in test_dmstag_vs_dmda.c) */
  ierr = DMCreateLocalVector(dm,&in);CHKERRQ(ierr);
  ierr = VecDuplicate(in,&out);CHKERRQ(ierr);
  ierr = PetscLogStagePush(mainStage);CHKERRQ(ierr);
  switch (test) {
    case STAGTEST :
      ierr = ApplyStag(dm,in,out);CHKERRQ(ierr);
      break;
    case DATEST1 :
      ierr = ApplyDA1(dm,in,out);CHKERRQ(ierr);
      break;
#if 0
    case DATEST2 :
      ierr = ApplyDA2(dm,in,out);CHKERRQ(ierr);
      break;
#endif
    default: SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported test %D",test);CHKERRQ(ierr);
  }
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  ierr = VecDestroy(&in);CHKERRQ(ierr);
  ierr = VecDestroy(&out);CHKERRQ(ierr);

  /* STAGE: destruction */
  ierr = PetscLogStagePush(destructionStage);CHKERRQ(ierr);
  if (test == DATEST2) {
    for (i=0; i<NUM_DAS; ++i) {
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

PetscErrorCode ApplyStag(DM dm,Vec in,Vec out)
{
  PetscErrorCode ierr;
  PetscInt       start[3],end[3],n[3],N[3],i,j,k,p,vx,vy,vz,d;
  PetscScalar    ****arrOut,****arrIn;

  PetscFunctionBeginUser;
  ierr = DMStagGetCorners(dm,&start[0],&start[1],&start[2],&n[0],&n[1],&n[2],NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(dm,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
  /* for convenience, only process internal elements */
  for (d=0; d<3; ++d) end[d] = start[d]+n[d] == N[d] ? N[d]-1 : start[d]+n[d];
  for (d=0; d<3; ++d) start[d] = start[d] == 0 ? 1 : start[d];
  ierr = DMStagGetLocationSlot(dm,DMSTAG_ELEMENT,0,&p);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm,DMSTAG_LEFT,0,&vx);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm,DMSTAG_DOWN,0,&vy);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm,DMSTAG_BACK,0,&vz);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayDOFRead(dm,in,&arrIn);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayDOF(dm,out,&arrOut);CHKERRQ(ierr);
  for (k=start[2]; k<end[2]; ++k) {
    for (j=start[1]; j<end[1]; ++j) {
      for (i=start[0]; i<end[0]; ++i) {
        arrOut[k][j][i][vy] = arrIn[k][j][i][p] + arrIn[k][j+1][i][p];
        arrOut[k][j][i][vx] = arrIn[k][j][i][p] + arrIn[k][j][i+1][p];
        arrOut[k][j][i][p] =
          arrIn[k][j+1][i][vy] + arrIn[k][j+1][i][vy] +
          arrIn[k][j][i][vx]   + arrIn[k][j][i+1][vx] +
          arrIn[k+1][j][i][vz] + arrIn[k+1][j][i][vz];
        arrOut[k][j][i][vz] = arrIn[k][j][i][p] + arrIn[k+1][j][i][p];
      }
    }
  }
  ierr = DMStagVecRestoreArrayDOFRead(dm,in,&arrIn);CHKERRQ(ierr);
  ierr = DMStagVecRestoreArrayDOF(dm,in,&arrOut);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ApplyDA1(DM dm,Vec in,Vec out)
{
  PetscErrorCode ierr;
  PetscInt       start[3],end[3],n[3],N[3],i,j,k,d,p,vx,vy,vz;
  PetscScalar    ****arrOut,****arrIn;

  PetscFunctionBeginUser;
  ierr = DMDAGetCorners(dm,&start[0],&start[1],&start[2],&n[0],&n[1],&n[2]);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dm,NULL,&N[0],&N[1],&N[2],NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  /* for convenience, only process internal elements */
  for (d=0; d<3; ++d) end[d] = start[d]+n[d] == N[d] ? N[d]-1 : start[d]+n[d];
  for (d=0; d<3; ++d) start[d] = start[d] == 0 ? 1 : start[d];
  vy = 0;
  vx = 1;
  p  = 2;
  vz = 3;
  ierr = DMDAVecGetArrayDOFRead(dm,in,&arrIn);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(dm,in,&arrOut);CHKERRQ(ierr);
  for (i=start[0]; i<end[0]; ++i) {
    for (j=start[1]; j<end[1]; ++j) {
      for (k=start[2]; k<end[2]; ++k) {
        arrOut[k][j][i][vy] = arrIn[k][j][i][p] + arrIn[k][j+1][i][p];
        arrOut[k][j][i][vx] = arrIn[k][j][i][p] + arrIn[k][j][i+1][p];
        arrOut[k][j][i][p] =
          arrIn[k][j+1][i][vy] + arrIn[k][j+1][i][vy] +
          arrIn[k][j][i][vx]   + arrIn[k][j][i+1][vx] +
          arrIn[k+1][j][i][vz] + arrIn[k+1][j][i][vz];
        arrOut[k][j][i][vz] = arrIn[k][j][i][p] + arrIn[k+1][j][i][p];
      }
    }
  }
  ierr = DMDAVecRestoreArrayDOFRead(dm,in,&arrIn);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(dm,in,&arrOut);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// TODO ApplyDA2 using DMComposite interface
