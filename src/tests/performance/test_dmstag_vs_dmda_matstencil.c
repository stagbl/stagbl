static char help[] = "For stencil-based matrix assembly, compare DMStag and a collection of DMDAs\n";
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

#define NUM_DAS 4

PetscErrorCode AssembleStag(DM,Mat);
PetscErrorCode AssembleDA1(DM,Mat);
//PetscErrorCode AssembleDA2(DM,Vec,Vec);

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DM             dm,dms[4];
  PetscInt       dof0,dof1,dof2,dof3,stencilWidth,Nx,Ny,Nz,i;
  PetscLogStage  creationStage,mainStage,destructionStage;
  PetscInt       test;
  Mat            A;

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
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx+1,Ny+1,Nz+1,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0+3*dof1+3*dof2+dof3,stencilWidth,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
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

  /* Overestimate */
  ierr = DMCreateMatrix(dm,&A);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(A,13,NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,13,NULL,13,NULL);CHKERRQ(ierr);

  ierr = PetscLogStagePush(mainStage);CHKERRQ(ierr);
  switch (test) {
    case STAGTEST :
      ierr = AssembleStag(dm,A);CHKERRQ(ierr);
      break;
    case DATEST1 :
      ierr = AssembleDA1(dm,A);CHKERRQ(ierr);
      break;
#if 0
    case DATEST2 :
      ierr = AssembleDA2(dm,in,out);CHKERRQ(ierr);
      break;
#endif
    default: SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported test %D",test);CHKERRQ(ierr);
  }
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

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

PetscErrorCode AssembleStag(DM dm,Mat A)
{
  PetscErrorCode ierr;
  PetscInt       start[3],end[3],n[3],N[3],i,j,k,d;
  DMStagStencil  row,col[7];
  PetscScalar    v[7];

  PetscFunctionBeginUser;
  ierr = DMStagGetCorners(dm,&start[0],&start[1],&start[2],&n[0],&n[1],&n[2],NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(dm,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
  /* for convenience, only process internal elements */
  for (d=0; d<3; ++d) end[d]   = start[d]+n[d] == N[d] ? N[d]-2 : start[d]+n[d];
  for (d=0; d<3; ++d) start[d] = start[d] == 0 ? 1 : start[d];
  row.c = 0;
  col[0].c = 0;
  col[1].c = 0;
  col[2].c = 0;
  col[3].c = 0;
  col[4].c = 0;
  col[5].c = 0;
  col[6].c = 0;
  for (k=start[2]; k<end[2]; ++k) {
    for (j=start[1]; j<end[1]; ++j) {
      for (i=start[0]; i<end[0]; ++i) {

        /* vy */
        row.i = i; row.j = j; row.k = k; row.loc = DMSTAG_DOWN;
        col[0].i = i; col[0].j = j  ; col[0].k = k; col[0].loc = DMSTAG_ELEMENT; v[0] = 1.0;
        col[1].i = i; col[1].j = j+1; col[1].k = k; col[1].loc = DMSTAG_ELEMENT; v[1] = 1.0;
        ierr = DMStagMatSetValuesStencil(dm,A,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);

        /* vx */
        row.loc = DMSTAG_LEFT;
        col[0].i = i  ; col[0].j = j; col[0].k = k; col[0].loc = DMSTAG_ELEMENT; v[0] = 1.0;
        col[1].i = i+1; col[1].j = j; col[1].k = k; col[1].loc = DMSTAG_ELEMENT; v[1] = 1.0;
        ierr = DMStagMatSetValuesStencil(dm,A,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);

        /* p */
        row.loc = DMSTAG_ELEMENT;
        col[0].i = i  ; col[0].j = j  ; col[0].k = k  ; col[0].loc = DMSTAG_LEFT;    v[0] = 1.0;
        col[1].i = i+1; col[1].j = j  ; col[1].k = k  ; col[1].loc = DMSTAG_LEFT;    v[1] = 1.0;
        col[2].i = i  ; col[2].j = j  ; col[2].k = k  ; col[2].loc = DMSTAG_DOWN;    v[2] = 1.0;
        col[3].i = i  ; col[3].j = j+1; col[3].k = k  ; col[3].loc = DMSTAG_DOWN;    v[3] = 1.0;
        col[4].i = i  ; col[4].j = j  ; col[4].k = k  ; col[4].loc = DMSTAG_BACK;    v[4] = 1.0;
        col[5].i = i  ; col[5].j = j  ; col[5].k = k+1; col[5].loc = DMSTAG_BACK;    v[5] = 1.0;
        col[6].i = i  ; col[6].j = j  ; col[6].k = k  ; col[6].loc = DMSTAG_ELEMENT; v[6] = 1.0;
        ierr = DMStagMatSetValuesStencil(dm,A,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);

        /* vz */
        row.loc = DMSTAG_BACK;
        col[0].i = i; col[0].j = j; col[0].k = k  ; col[0].loc = DMSTAG_ELEMENT; v[0] = 1.0;
        col[1].i = i; col[1].j = j; col[1].k = k+1; col[1].loc = DMSTAG_ELEMENT; v[1] = 1.0;
        ierr = DMStagMatSetValuesStencil(dm,A,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode AssembleDA1(DM dm,Mat A)
{
  PetscErrorCode ierr;
  PetscInt       start[3],end[3],n[3],N[3],i,j,k,d,p,vx,vy,vz;
  MatStencil     row,col[7];
  PetscScalar    v[7];

  PetscFunctionBeginUser;
  ierr = DMDAGetCorners(dm,&start[0],&start[1],&start[2],&n[0],&n[1],&n[2]);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dm,NULL,&N[0],&N[1],&N[2],NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  /* for convenience, only process internal elements */
  for (d=0; d<3; ++d) end[d]   = start[d]+n[d] == N[d] ? N[d]-2 : start[d]+n[d]-1;
  for (d=0; d<3; ++d) start[d] = start[d] == 0 ? 1 : start[d];
  vy = 0;
  vx = 1;
  p  = 2;
  vz = 3;
  for (k=start[2]; k<end[2]; ++k) {
    for (j=start[1]; j<end[1]; ++j) {
      for (i=start[0]; i<end[0]; ++i) {
        /* vy */
        row.i    = i; row.j    = j  ; row.k    = k; row.c    = vy;
        col[0].i = i; col[0].j = j  ; col[0].k = k; col[0].c = p; v[0] = 1.0;
        col[1].i = i; col[1].j = j+1; col[1].k = k; col[1].c = p; v[1] = 1.0;
        ierr = MatSetValuesStencil(A,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);

        /* vx */
        row.c = vx;
        col[0].i = i  ; col[0].j = j; col[0].k = k; col[0].c = p; v[0] = 1.0;
        col[1].i = i+1; col[1].j = j; col[1].k = k; col[1].c = p; v[1] = 1.0;
        ierr = MatSetValuesStencil(A,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);

        /* p */
        row.c = p;
        col[0].i = i  ; col[0].j = j  ; col[0].k = k  ; col[0].c = vx; v[0] = 1.0;
        col[1].i = i+1; col[1].j = j  ; col[1].k = k  ; col[1].c = vx; v[1] = 1.0;
        col[2].i = i  ; col[2].j = j  ; col[2].k = k  ; col[2].c = vy; v[2] = 1.0;
        col[3].i = i  ; col[3].j = j+1; col[3].k = k  ; col[3].c = vy; v[3] = 1.0;
        col[4].i = i  ; col[4].j = j  ; col[4].k = k  ; col[4].c = vz; v[4] = 1.0;
        col[5].i = i  ; col[5].j = j  ; col[5].k = k+1; col[5].c = vz; v[5] = 1.0;
        col[6].i = i  ; col[6].j = j  ; col[6].k = k  ; col[6].c = p ; v[6] = 1.0;
        ierr = MatSetValuesStencil(A,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);

        /* vz */
        row.c = vz;
        col[0].i = i; col[0].j = j; col[0].k = k  ; col[0].c = p; v[0] = 1.0;
        col[1].i = i; col[1].j = j; col[1].k = k+1; col[1].c = p; v[1] = 1.0;
        ierr = MatSetValuesStencil(A,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// TODO AssembleDA2 using DMComposite interface
