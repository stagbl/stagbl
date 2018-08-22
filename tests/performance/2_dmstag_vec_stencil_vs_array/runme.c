static char help[] = "Compare performance for stencil- and array-based vector entry access, using DMStag";

/* Note: for the control case, with DMDA, the command line options are different!
         Thus, use -Nx instead of -stag_grid_x etc. */

#include <petscdm.h>
#include <petscdmstag.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>

/* Supply one of these with -test  */
typedef enum {
  TEST_ARRAY_DMDA = 0, /* Control */
  TEST_ARRAY      = 1, /* Typical array use */
  TEST_STENCIL    = 2,
  TEST_STENCIL2   = 3, /* Probably the most typical stencil use */
  TEST_STENCIL3   = 4,
} Test;

PetscErrorCode TestArrayDMDA(Vec);
PetscErrorCode TestArray(Vec);
PetscErrorCode TestStencil(Vec);
PetscErrorCode TestStencil2(Vec);
PetscErrorCode TestStencil3(Vec);

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DM             dm;
  Vec            vec;
  PetscInt       dof0,dof1,dof2,dof3,stencilWidth,Nx,Ny,Nz;
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
  Nx = 100;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nx",&Nx,NULL);CHKERRQ(ierr);
  Ny = Nx;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Ny",&Ny,NULL);CHKERRQ(ierr);
  Nz = Nx;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nz",&Nz,NULL);CHKERRQ(ierr);
  test = TEST_STENCIL;
  ierr = PetscOptionsGetInt(NULL,NULL,"-test",&test,NULL);CHKERRQ(ierr);


  /***** STAGE : creation ******************************************************/
  ierr = PetscLogStagePush(creationStage);CHKERRQ(ierr);
  if (test == TEST_ARRAY_DMDA) {
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0+3*dof1+3*dof2+dof3,stencilWidth,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
  } else {
    ierr = DMStagCreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,Nx,Ny,Nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof0,dof1,dof2,dof3,DMSTAG_STENCIL_BOX,stencilWidth,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
  }
  ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
  ierr = DMSetUp(dm);CHKERRQ(ierr);
  ierr = DMGetGlobalVector(dm,&vec);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  /*****************************************************************************/

  if (test != TEST_ARRAY_DMDA) {
    ierr = DMStagGetDOF(dm,&dof0,&dof1,&dof2,&dof3);CHKERRQ(ierr);
    if (dof0 != 1 || dof1 != 1 || dof2 != 1 || dof3 != 1) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Dof on all strata must be 1");
  }

  /* *** STAGE: main operation *************************************************/
  ierr = PetscLogStagePush(mainStage);CHKERRQ(ierr);
  switch (test) {
    case TEST_ARRAY_DMDA:
      ierr = TestArrayDMDA(vec); CHKERRQ(ierr);
      break;
    case TEST_ARRAY:
      ierr = TestArray(vec); CHKERRQ(ierr);
      break;
    case TEST_STENCIL:
      ierr = TestStencil(vec); CHKERRQ(ierr);
      break;
    case TEST_STENCIL2:
      ierr = TestStencil2(vec); CHKERRQ(ierr);
      break;
    case TEST_STENCIL3:
      ierr = TestStencil3(vec); CHKERRQ(ierr);
      break;
    default: SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported test %D",test);CHKERRQ(ierr);
  }
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  /*****************************************************************************/


  /***** STAGE: destruction ****************************************************/
  ierr = PetscLogStagePush(destructionStage);CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(dm,&vec);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  /*****************************************************************************/

  /* Finalize */
  ierr = PetscFinalize();
  return ierr;
}

PetscErrorCode TestArrayDMDA(Vec vec)
{
  PetscErrorCode ierr;
  DM             dm;
  Vec            vecLocal;
  PetscScalar    ****arr;
  PetscInt       startx,starty,startz,nx,ny,nz,i,j,k,s;

  PetscFunctionBeginUser;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&vecLocal);CHKERRQ(ierr);
  ierr = DMDAGetCorners(dm,&startx,&starty,&startz,&nx,&ny,&nz);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(dm,vecLocal,&arr);CHKERRQ(ierr);
  for (k=startz; k<startz + nz; ++k) {
    for (j=starty; j<starty + ny; ++j) {
      for (i=startx; i<startx + nx; ++i) {
        for (s=0; s<8; ++s) {
          arr[k][j][i][s] = i + j + k + s; /* No lookup for s! */
        }
      }
    }
  }
  ierr = DMDAVecRestoreArrayDOF(dm,vecLocal,&arr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,vecLocal,INSERT_VALUES,vec);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,vecLocal,INSERT_VALUES,vec);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&vecLocal);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestArray(Vec vec)
{
  PetscErrorCode        ierr;
  DM                    dm;
  Vec                   vecLocal;
  PetscScalar           ****arr;
  PetscInt              startx,starty,startz,nx,ny,nz,i,j,k,s;
  DMStagStencilLocation loc[8];
  PetscInt              slot[8];

  PetscFunctionBeginUser;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&vecLocal);CHKERRQ(ierr);
  ierr = DMStagGetCorners(dm,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);
  loc[0] = DMSTAG_BACK_DOWN_LEFT;
  loc[1] = DMSTAG_BACK_DOWN;
  loc[2] = DMSTAG_BACK_LEFT;
  loc[3] = DMSTAG_BACK;
  loc[4] = DMSTAG_DOWN_LEFT;
  loc[5] = DMSTAG_DOWN;
  loc[6] = DMSTAG_LEFT;
  loc[7] = DMSTAG_ELEMENT;
  for (s=0; s<8; ++s) {
    ierr = DMStagGetLocationSlot(dm,loc[s],0,&slot[s]);CHKERRQ(ierr);
  }
  ierr = DMStagVecGetArrayDOF(dm,vecLocal,&arr);CHKERRQ(ierr);
  for (k=startz; k<startz + nz; ++k) {
    for (j=starty; j<starty + ny; ++j) {
      for (i=startx; i<startx + nx; ++i) {
        for (s=0; s<8; ++s) {
            arr[k][j][i][slot[s]] = i + j + k + s;
        }
      }
    }
  }
  ierr = DMStagVecRestoreArrayDOF(dm,vecLocal,&arr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,vecLocal,INSERT_VALUES,vec);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,vecLocal,INSERT_VALUES,vec);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&vecLocal);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestStencil(Vec vec)
{
  PetscErrorCode        ierr;
  DM                    dm;
  PetscInt              startx,starty,startz,nx,ny,nz,i,j,k,s,d;
  DMStagStencilLocation loc[8];
  PetscInt              dof[8];
  DMStagStencil         pos;

  PetscFunctionBeginUser;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMStagGetCorners(dm,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);
  loc[0] = DMSTAG_BACK_DOWN_LEFT;
  loc[1] = DMSTAG_BACK_DOWN;
  loc[2] = DMSTAG_BACK_LEFT;
  loc[3] = DMSTAG_BACK;
  loc[4] = DMSTAG_DOWN_LEFT;
  loc[5] = DMSTAG_DOWN;
  loc[6] = DMSTAG_LEFT;
  loc[7] = DMSTAG_ELEMENT;
  for (s=0; s<8; ++s) {
    ierr = DMStagGetLocationDOF(dm,loc[s],&dof[s]);CHKERRQ(ierr);
  }
  for (k=startz; k<startz + nz; ++k) {
    pos.k = k;
    for (j=starty; j<starty + ny; ++j) {
      pos.j = j;
      for (i=startx; i<startx + nx; ++i) {
        pos.i = i;
        for (s=0; s<8; ++s) {
          pos.loc = loc[s];
          for (d=0; d<dof[s]; ++d) {
            pos.c = d;
            PetscScalar val = i+j+k+s+d;
            ierr = DMStagVecSetValuesStencil(dm,vec,1,&pos,&val,INSERT_VALUES);CHKERRQ(ierr);
          }
        }
      }
    }
  }
  ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if 1
PetscErrorCode TestStencil2(Vec vec)
{
  PetscErrorCode ierr;
  DM             dm;
  PetscInt       startx,starty,startz,nx,ny,nz,i,j,k;
  DMStagStencil  pos[8]; /* Not appropriate if you change dof! */
  PetscScalar    val[8];

  PetscFunctionBeginUser;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMStagGetCorners(dm,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);
  for (k=startz; k<startz + nz; ++k) {
    for (j=starty; j<starty + ny; ++j) {
      for (i=startx; i<startx + nx; ++i) {
        pos[0].i = i; pos[0].j = j; pos[0].k = k; pos[0].loc = DMSTAG_BACK_DOWN_LEFT; pos[0].c = 0; val[0] = i+j+k;
        pos[1].i = i; pos[1].j = j; pos[1].k = k; pos[1].loc = DMSTAG_BACK_DOWN;      pos[1].c = 0; val[1] = i+j+k+1;
        pos[2].i = i; pos[2].j = j; pos[2].k = k; pos[2].loc = DMSTAG_BACK_LEFT;      pos[2].c = 0; val[2] = i+j+k+2;
        pos[3].i = i; pos[3].j = j; pos[3].k = k; pos[3].loc = DMSTAG_BACK;           pos[3].c = 0; val[3] = i+j+k+3;
        pos[4].i = i; pos[4].j = j; pos[4].k = k; pos[4].loc = DMSTAG_DOWN_LEFT;      pos[4].c = 0; val[4] = i+j+k+4;
        pos[5].i = i; pos[5].j = j; pos[5].k = k; pos[5].loc = DMSTAG_DOWN;           pos[5].c = 0; val[5] = i+j+k+5;
        pos[6].i = i; pos[6].j = j; pos[6].k = k; pos[6].loc = DMSTAG_LEFT;           pos[6].c = 0; val[6] = i+j+k+6;
        pos[7].i = i; pos[7].j = j; pos[7].k = k; pos[7].loc = DMSTAG_ELEMENT;        pos[7].c = 0; val[7] = i+j+k+7;
        ierr = DMStagVecSetValuesStencil(dm,vec,8,pos,val,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#else
/* A marginally-faster variant */
PetscErrorCode TestStencil2(Vec vec)
{
  PetscErrorCode ierr;
  DM             dm;
  PetscInt       startx,starty,startz,nx,ny,nz,i,j,k;
  DMStagStencil  pos[8]; /* Not appropriate if you change dof! */
  PetscScalar    val[8];

  PetscFunctionBeginUser;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMStagGetCorners(dm,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);
  pos[0].loc = DMSTAG_BACK_DOWN_LEFT; pos[0].c = 0;
  pos[1].loc = DMSTAG_BACK_DOWN;      pos[1].c = 0;
  pos[2].loc = DMSTAG_BACK_LEFT;      pos[2].c = 0;
  pos[3].loc = DMSTAG_BACK;           pos[3].c = 0;
  pos[4].loc = DMSTAG_DOWN_LEFT;      pos[4].c = 0;
  pos[5].loc = DMSTAG_DOWN;           pos[5].c = 0;
  pos[6].loc = DMSTAG_LEFT;           pos[6].c = 0;
  pos[7].loc = DMSTAG_ELEMENT;        pos[7].c = 0;
  for (k=startz; k<startz + nz; ++k) {
    pos[0].k = k;
    pos[1].k = k;
    pos[2].k = k;
    pos[3].k = k;
    pos[4].k = k;
    pos[5].k = k;
    pos[6].k = k;
    pos[7].k = k;
    for (j=starty; j<starty + ny; ++j) {
      pos[0].j = j;
      pos[1].j = j;
      pos[2].j = j;
      pos[3].j = j;
      pos[4].j = j;
      pos[5].j = j;
      pos[6].j = j;
      pos[7].j = j;
      for (i=startx; i<startx + nx; ++i) {
        pos[0].i = i;  val[0] = i+j+k;
        pos[1].i = i;  val[1] = i+j+k+1;
        pos[2].i = i;  val[2] = i+j+k+2;
        pos[3].i = i;  val[3] = i+j+k+3;
        pos[4].i = i;  val[4] = i+j+k+4;
        pos[5].i = i;  val[5] = i+j+k+5;
        pos[6].i = i;  val[6] = i+j+k+6;
        pos[7].i = i;  val[7] = i+j+k+7;
        ierr = DMStagVecSetValuesStencil(dm,vec,8,pos,val,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif

PetscErrorCode TestStencil3(Vec vec)
{
  PetscErrorCode ierr;
  DM             dm;
  PetscInt       startx,starty,startz,nx,ny,nz,i,j,k,nEntries;
  DMStagStencil  *pos;
  PetscScalar    *val;

  PetscFunctionBeginUser;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMStagGetCorners(dm,&startx,&starty,&startz,&nx,&ny,&nz,NULL,NULL,NULL);CHKERRQ(ierr);
  nEntries = 8*nx*ny*nz;
  ierr = PetscMalloc1(nEntries,&pos);CHKERRQ(ierr);
  ierr = PetscMalloc1(nEntries,&val);CHKERRQ(ierr);
  for (k=startz; k<startz + nz; ++k) {
    for (j=starty; j<starty + ny; ++j) {
      for (i=startx; i<startx + nx; ++i) {
        const PetscInt offs = (nx*ny*k + nx*j + i)*8;
        pos[offs + 0].i = i; pos[offs + 0].j = j; pos[offs + 0].k = k; pos[offs + 0].loc = DMSTAG_BACK_DOWN_LEFT; pos[offs + 0].c = 0; val[offs + 0] = i+j+k;
        pos[offs + 1].i = i; pos[offs + 1].j = j; pos[offs + 1].k = k; pos[offs + 1].loc = DMSTAG_BACK_DOWN;      pos[offs + 1].c = 0; val[offs + 1] = i+j+k+1;
        pos[offs + 2].i = i; pos[offs + 2].j = j; pos[offs + 2].k = k; pos[offs + 2].loc = DMSTAG_BACK_LEFT;      pos[offs + 2].c = 0; val[offs + 2] = i+j+k+2;
        pos[offs + 3].i = i; pos[offs + 3].j = j; pos[offs + 3].k = k; pos[offs + 3].loc = DMSTAG_BACK;           pos[offs + 3].c = 0; val[offs + 3] = i+j+k+3;
        pos[offs + 4].i = i; pos[offs + 4].j = j; pos[offs + 4].k = k; pos[offs + 4].loc = DMSTAG_DOWN_LEFT;      pos[offs + 4].c = 0; val[offs + 4] = i+j+k+4;
        pos[offs + 5].i = i; pos[offs + 5].j = j; pos[offs + 5].k = k; pos[offs + 5].loc = DMSTAG_DOWN;           pos[offs + 5].c = 0; val[offs + 5] = i+j+k+5;
        pos[offs + 6].i = i; pos[offs + 6].j = j; pos[offs + 6].k = k; pos[offs + 6].loc = DMSTAG_LEFT;           pos[offs + 6].c = 0; val[offs + 6] = i+j+k+6;
        pos[offs + 7].i = i; pos[offs + 7].j = j; pos[offs + 7].k = k; pos[offs + 7].loc = DMSTAG_ELEMENT;        pos[offs + 7].c = 0; val[offs + 7] = i+j+k+7;
      }
    }
  }
  ierr = DMStagVecSetValuesStencil(dm,vec,nEntries,pos,val,INSERT_VALUES);CHKERRQ(ierr);
  ierr = PetscFree(pos);CHKERRQ(ierr);
  ierr = PetscFree(val);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
