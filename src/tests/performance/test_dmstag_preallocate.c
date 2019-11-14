static char help[] = "Solve simple assembled Stokes linear solve with DMStag, including MatPreallocator usage\n";

/* Note: this could be extended to compare with a DMDA-based appraoch, as
   a performance comparison. See PETSc SNES tutorial ex30  */

/* Note: the additional code wrt MatPreallocate is likely unnecessary as of
   PETSc 3.12 */

#include <petscdm.h>
#include <petscdmstag.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>
#include <petscksp.h>

#define DOWN_LEFT  DMSTAG_DOWN_LEFT
#define DOWN       DMSTAG_DOWN
#define DOWN_RIGHT DMSTAG_DOWN_RIGHT
#define LEFT       DMSTAG_LEFT
#define ELEMENT    DMSTAG_ELEMENT
#define RIGHT      DMSTAG_RIGHT
#define UP_LEFT    DMSTAG_UP_LEFT
#define UP         DMSTAG_UP
#define UP_RIGHT   DMSTAG_UP_RIGHT

typedef enum {
  STAGTEST = 0,
  DATEST   = 1,
} Test;

static PetscErrorCode CreateSystem_Stag(DM,Mat*,Vec*,PetscBool);
static PetscErrorCode MatGetPreallocator(Mat A,Mat *preallocator);
static PetscErrorCode MatPreallocatePhaseBegin(Mat A,Mat *preallocator);
static PetscErrorCode MatPreallocatePhaseEnd(Mat A);

#define RHO(x,y) PetscSinScalar(PETSC_PI * y) * PetscCosScalar(PETSC_PI * x)
#define ETA(x,y) x > 0.5 ? 1.0e6 : 1.0
#define GY 1.0

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DM             dm,dm_vel,dm_p;
  PetscInt       Nx,Ny;
  PetscInt       test;
  KSP            ksp;
  Mat            A;
  Vec            b,x;
  PetscBool      preallocate;
  PetscLogStage  creation_stage,assembly_stage,solve_stage,destruction_stage;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscLogStageRegister("Creation",&creation_stage);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Assembly",&assembly_stage);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Solve",&solve_stage);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Destruction",&destruction_stage);CHKERRQ(ierr);

  Nx = 17;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Nx",&Nx,NULL);CHKERRQ(ierr);
  Ny = Nx;
  ierr = PetscOptionsGetInt(NULL,NULL,"-Ny",&Ny,NULL);CHKERRQ(ierr);
  test = STAGTEST;
  ierr = PetscOptionsGetInt(NULL,NULL,"-test",&test,NULL);CHKERRQ(ierr);
  preallocate = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-preallocate",&preallocate,NULL);CHKERRQ(ierr);

  ierr = PetscLogStagePush(creation_stage);CHKERRQ(ierr);
  if (test == STAGTEST) {
    ierr = DMStagCreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,0,1,1,DMSTAG_STENCIL_BOX,1,NULL,NULL,&dm);CHKERRQ(ierr);
    ierr = DMSetUp(dm);CHKERRQ(ierr);
    ierr = DMStagSetUniformCoordinatesProduct(dm, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0); CHKERRQ(ierr);
  } else if (test == DATEST) {
    /* Note: this takes the "dummy variable" approach, adding unused pressure dof around the boundary */
    ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,Nx+1,Ny+1,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&dm_p);CHKERRQ(ierr);
    ierr = DMSetUp(dm_p);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates(dm_p, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0); CHKERRQ(ierr);
    ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,Nx+1,Ny+1,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&dm_vel);CHKERRQ(ierr);
    ierr = DMSetUp(dm_vel);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates(dm_vel, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0); CHKERRQ(ierr);
    ierr = DMCompositeCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
    ierr = DMCompositeAddDM(dm,dm_vel);CHKERRQ(ierr);
    ierr = DMCompositeAddDM(dm,dm_p);CHKERRQ(ierr);
  } else SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported test %D",test);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  ierr = PetscLogStagePush(assembly_stage);CHKERRQ(ierr);
  switch (test) {
    case STAGTEST:
      ierr = CreateSystem_Stag(dm, &A, &b, preallocate); CHKERRQ(ierr);
      break;
    default: SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unsupported test %D",test);CHKERRQ(ierr);
  }
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  ierr = PetscLogStagePush(solve_stage);CHKERRQ(ierr);
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  ierr = PetscLogStagePush(destruction_stage);CHKERRQ(ierr);
  if (test == DATEST) {
      ierr = DMDestroy(&dm_vel);CHKERRQ(ierr);
      ierr = DMDestroy(&dm_p);CHKERRQ(ierr);
    }
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* Finalize */
  ierr = PetscFinalize();
  return ierr;
}

static PetscErrorCode CreateSystem_Stag(DM dm, Mat *pA, Vec *pRhs, PetscBool preallocate)
{
  PetscErrorCode ierr;
  PetscInt       Nx, Ny;                         // global variables
  PetscInt       ex, ey, startx, starty, nx, ny; // local variables
  Mat            A, preallocator = NULL;
  Vec            rhs;
  PetscScalar    hx, hy, Kbound, Kcont;
  PetscBool      pinPressure = PETSC_TRUE;

  PetscFunctionBeginUser;

  //ierr = DMSetMatrixPreallocateOnly(dm, PETSC_TRUE); CHKERRQ(ierr);

  /* Create stiffness matrix A and rhs vector */
  ierr = DMCreateMatrix      (dm, pA  ); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm, pRhs); CHKERRQ(ierr);

  /* Assign pointers and other variables */
  A   = *pA;
  rhs = *pRhs;

  ierr = DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);CHKERRQ(ierr);
  hx = 1.0 / (Nx+1);
  hy = 1.0 / (Ny+1);
  Kbound = (Nx+1)*(Ny+1);
  Kcont = (Nx+1);

  /* Get local domain */
  ierr = DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

  if (preallocate) {
    ierr = MatPreallocatePhaseBegin(*pA,&preallocator);CHKERRQ(ierr);


    /* Get non-zero pattern - Loop over all local elements
       This should be put into a routine that is called twice: first for non-zero structure, second for setting values*/
    PetscInt      nEntries;
    PetscScalar   vv[11] = {0.0};
    DMStagStencil row,col[11];

    for (ey = starty; ey<starty+ny; ++ey) {
      for (ex = startx; ex<startx+nx; ++ex) {

        /* Top boundary velocity Dirichlet */
        if (ey == Ny-1) {
          nEntries = 1;
          row.i    = ex; row.j = ey; row.loc = UP; row.c = 0;
          col[0]   = row;

          ierr = DMStagMatSetValuesStencil(dm,preallocator,1,&row,nEntries,col,vv,INSERT_VALUES); CHKERRQ(ierr);
        }

        /* Bottom boundary velocity Dirichlet */
        if (ey == 0) {
          nEntries = 1;
          row.i    = ex; row.j = ey; row.loc = DOWN; row.c = 0;
          col[0]   = row;

          ierr = DMStagMatSetValuesStencil(dm,preallocator,1,&row,nEntries,col,vv,INSERT_VALUES); CHKERRQ(ierr);
        }

        /* Y-momentum equation : (u_xx + u_yy) - p_y = f^y : includes non-zero forcing */
        else {
          /* Left boundary y velocity stencil */
          if (ex == 0) {
            nEntries = 10;
            row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
            col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0;
            col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0;
            col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0;
            /* No left entry */
            col[3].i = ex+1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0;
            col[4].i = ex  ; col[4].j  = ey-1; col[4].loc  = LEFT;     col[4].c  = 0;
            col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = RIGHT;    col[5].c  = 0;
            col[6].i = ex  ; col[6].j  = ey  ; col[6].loc  = LEFT;     col[6].c  = 0;
            col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = RIGHT;    col[7].c  = 0;
            col[8].i = ex  ; col[8].j  = ey-1; col[8].loc  = ELEMENT;  col[8].c  = 0;
            col[9].i = ex  ; col[9].j  = ey  ; col[9].loc = ELEMENT;   col[9].c  = 0;
          }

          /* Right boundary y velocity stencil */
          else if (ex == Nx-1) {
            nEntries = 10;
            row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
            col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0;
            col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0;
            col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0;
            col[3].i = ex-1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0;
            /* No right element */
            col[4].i = ex  ; col[4].j  = ey-1; col[4].loc  = LEFT;     col[4].c  = 0;
            col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = RIGHT;    col[5].c  = 0;
            col[6].i = ex  ; col[6].j  = ey  ; col[6].loc  = LEFT;     col[6].c  = 0;
            col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = RIGHT;    col[7].c  = 0;
            col[8].i = ex  ; col[8].j  = ey-1; col[8].loc  = ELEMENT;  col[8].c  = 0;
            col[9].i = ex  ; col[9].j = ey   ; col[9].loc = ELEMENT;   col[9].c  = 0;
          }

          /* U_y interior equation */
          else {
            nEntries = 11;
            row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
            col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0;
            col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0;
            col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0;
            col[3].i = ex-1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0;
            col[4].i = ex+1; col[4].j  = ey  ; col[4].loc  = DOWN;     col[4].c  = 0;
            col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = LEFT;     col[5].c  = 0;
            col[6].i = ex  ; col[6].j  = ey-1; col[6].loc  = RIGHT;    col[6].c  = 0;
            col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = LEFT;     col[7].c  = 0;
            col[8].i = ex  ; col[8].j  = ey  ; col[8].loc  = RIGHT;    col[8].c  = 0;
            col[9].i = ex  ; col[9].j  = ey-1; col[9].loc  = ELEMENT;  col[9].c  = 0;
            col[10].i = ex ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0;
          }

          /* Insert Y-momentum entries */
          ierr = DMStagMatSetValuesStencil(dm,preallocator,1,&row,nEntries,col,vv,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Right Boundary velocity Dirichlet */
        if (ex == Nx-1) {
          nEntries = 1;
          row.i    = ex; row.j = ey; row.loc = RIGHT; row.c = 0;
          col[0]   = row;
          ierr = DMStagMatSetValuesStencil(dm,preallocator,1,&row,nEntries,col,vv,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* Left velocity Dirichlet */
        if (ex == 0) {
          nEntries = 1;
          row.i = ex; row.j = ey; row.loc = LEFT; row.c = 0;
          col[0]   = row;
          ierr = DMStagMatSetValuesStencil(dm,preallocator,1,&row,nEntries,col,vv,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
        else {

          /* Bottom boundary x velocity stencil (with zero vel deriv) */
          if (ey == 0) {
            nEntries = 10;
            row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
            col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0;
            /* Missing element below */
            col[1].i  = ex  ; col[1].j  = ey+1; col[1].loc  = LEFT;    col[1].c   = 0;
            col[2].i  = ex-1; col[2].j  = ey  ; col[2].loc  = LEFT;    col[2].c   = 0;
            col[3].i  = ex+1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0;
            col[4].i  = ex-1; col[4].j  = ey  ; col[4].loc  = DOWN;    col[4].c   = 0;
            col[5].i  = ex  ; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0;
            col[6].i  = ex-1; col[6].j  = ey  ; col[6].loc  = UP;      col[6].c   = 0;
            col[7].i  = ex  ; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0;
            col[8].i  = ex-1; col[8].j  = ey  ; col[8].loc  = ELEMENT; col[8].c   = 0;
            col[9].i = ex   ; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0;
          }

          /* Top boundary x velocity stencil */
          else if (ey == Ny-1) {
            nEntries = 10;
            row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
            col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0;
            col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc  = LEFT;    col[1].c   = 0;
            /* Missing element above */
            col[2].i  = ex-1; col[2].j  = ey  ; col[2].loc  = LEFT;    col[2].c   = 0;
            col[3].i  = ex+1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0;
            col[4].i  = ex-1; col[4].j  = ey  ; col[4].loc  = DOWN;    col[4].c   = 0;
            col[5].i  = ex  ; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0;
            col[6].i  = ex-1; col[6].j  = ey  ; col[6].loc  = UP;      col[6].c   = 0;
            col[7].i  = ex  ; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0;
            col[8].i  = ex-1; col[8].j  = ey  ; col[8].loc  = ELEMENT; col[8].c   = 0;
            col[9].i = ex   ; col[9].j  = ey   ; col[9].loc = ELEMENT;  col[9].c  = 0;
          }

          /* U_x interior equation */
          else {
            nEntries = 11;
            row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
            col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0;
            col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc  = LEFT;    col[1].c   = 0;
            col[2].i  = ex  ; col[2].j  = ey+1; col[2].loc  = LEFT;    col[2].c   = 0;
            col[3].i  = ex-1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0;
            col[4].i  = ex+1; col[4].j  = ey  ; col[4].loc  = LEFT;    col[4].c   = 0;
            col[5].i  = ex-1; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0;
            col[6].i  = ex  ; col[6].j  = ey  ; col[6].loc  = DOWN;    col[6].c   = 0;
            col[7].i  = ex-1; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0;
            col[8].i  = ex  ; col[8].j  = ey  ; col[8].loc  = UP;      col[8].c   = 0;
            col[9].i  = ex-1; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0;
            col[10].i = ex  ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0;
          }

          /* Insert X-momentum entries */
          ierr = DMStagMatSetValuesStencil(dm,preallocator,1,&row,nEntries,col,vv,INSERT_VALUES);CHKERRQ(ierr);
        }

        /* P equation : u_x + v_y = 0 */
        if (pinPressure && ex == 0 && ey == 0) {
          nEntries = 1;
          row.i = ex; row.j = ey; row.loc = ELEMENT; row.c = 0;
          col[0]   = row;
          ierr = DMStagMatSetValuesStencil(dm,preallocator,1,&row,nEntries,col,vv,INSERT_VALUES);CHKERRQ(ierr);
        }
        else {
          nEntries = 5;
          row.i = ex; row.j = ey; row.loc = ELEMENT; row.c = 0;
          col[0].i = ex; col[0].j = ey; col[0].loc = LEFT;    col[0].c = 0;
          col[1].i = ex; col[1].j = ey; col[1].loc = RIGHT;   col[1].c = 0;
          col[2].i = ex; col[2].j = ey; col[2].loc = DOWN;    col[2].c = 0;
          col[3].i = ex; col[3].j = ey; col[3].loc = UP;      col[3].c = 0;
          col[4] = row;
          ierr = DMStagMatSetValuesStencil(dm,preallocator,1,&row,nEntries,col,vv,INSERT_VALUES);CHKERRQ(ierr);
        }


      }
    }

    /* Push the non-zero pattern defined within preallocator into A */
    ierr = MatPreallocatePhaseEnd(*pA);CHKERRQ(ierr);
  }

  /* View preallocated struct of A */
  //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  /* Insert non-zero values in A and assemble rhs - Loop over all local elements */
  for (ey = starty; ey<starty+ny; ++ey) {
    for (ex = startx; ex<startx+nx; ++ex) {
      PetscScalar x_left,x_center,x_right,y_down,y_center,y_up;

      /* Use regular grid of known size to directly compute coordinates */
      x_left = ex * hx;
      x_right = (ex+1) * hx;
      x_center = (x_left + x_right) * 0.5;
      y_down = ey * hy;
      y_up = (ey+1) * hy;
      y_center = (y_down + y_up) * 0.5;

      /* Top boundary velocity Dirichlet */
      if (ey == Ny-1) {
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar   valA = Kbound;

        row.i = ex; row.j = ey; row.loc = UP; row.c = 0;
        valRhs = 0.0;

        ierr   = DMStagMatSetValuesStencil(dm,A  ,1,&row,1,&row,&valA  ,INSERT_VALUES);CHKERRQ(ierr);
        ierr   = DMStagVecSetValuesStencil(dm,rhs,1,&row,       &valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      /* Bottom boundary velocity Dirichlet */
      if (ey == 0) {
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = Kbound;

        row.i = ex; row.j = ey; row.loc = DOWN; row.c = 0;
        valRhs = 0.0;

        ierr = DMStagMatSetValuesStencil(dm,A  ,1,&row,1,&row,&valA  ,INSERT_VALUES);CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm,rhs,1,&row,       &valRhs,INSERT_VALUES);CHKERRQ(ierr);
      }

      /* Y-momentum equation : (u_xx + u_yy) - p_y = f^y : includes non-zero forcing */
      else {
        PetscInt      nEntries;
        DMStagStencil row,col[11];
        PetscScalar   valA[11];
        PetscScalar   rho[2],valRhs;
        PetscScalar   etaLeft,etaRight,etaUp,etaDown;

        /* Get rho values  and compute rhs value*/
        rho[0] = RHO(x_left,y_down);
        rho[1] = RHO(x_right,y_down);
        valRhs = - GY * 0.5 * (rho[0] + rho[1]);

        /* Get eta values */
        etaLeft  = ETA(x_left,y_down);
        etaRight = ETA(x_right,y_down);
        etaUp    = ETA(x_center,y_center);
        etaDown  = ETA(x_center,y_center-hy);

        /* Left boundary y velocity stencil */
        if (ex == 0) {
          nEntries = 10;
          row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * (etaDown + etaUp) / (hy*hy) - (etaRight) /(hx*hx);
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * etaUp    / (hy*hy);
          /* No left entry */
          col[3].i = ex+1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0; valA[3]  =        etaRight / (hx*hx);
          col[4].i = ex  ; col[4].j  = ey-1; col[4].loc  = LEFT;     col[4].c  = 0; valA[4]  =        etaLeft  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = RIGHT;    col[5].c  = 0; valA[5]  = -      etaRight / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j  = ey  ; col[6].loc  = LEFT;     col[6].c  = 0; valA[6]  = -      etaLeft  / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = RIGHT;    col[7].c  = 0; valA[7]  =        etaRight / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j  = ey-1; col[8].loc  = ELEMENT;  col[8].c  = 0; valA[8]  =  Kcont / hy;
          col[9].i = ex  ; col[9].j = ey   ; col[9].loc = ELEMENT;   col[9].c  = 0; valA[9]  = -Kcont / hy;
        }

        /* Right boundary y velocity stencil */
        else if (ex == Nx-1) {
          nEntries = 10;
          row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * (etaDown + etaUp) / (hy*hy) - (etaLeft) /(hx*hx );
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * etaUp    / (hy*hy);
          col[3].i = ex-1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0; valA[3]  =        etaLeft  / (hx*hx);
          /* No right element */
          col[4].i = ex  ; col[4].j  = ey-1; col[4].loc  = LEFT;     col[4].c  = 0; valA[4]  =        etaLeft  / (hx*hy); /* down left x edge */
          col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = RIGHT;    col[5].c  = 0; valA[5]  = -      etaRight / (hx*hy); /* down right x edge */
          col[6].i = ex  ; col[6].j  = ey  ; col[6].loc  = LEFT;     col[6].c  = 0; valA[7]  = -      etaLeft  / (hx*hy); /* up left x edge */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = RIGHT;    col[7].c  = 0; valA[7]  =        etaRight / (hx*hy); /* up right x edge */
          col[8].i = ex  ; col[8].j  = ey-1; col[8].loc  = ELEMENT;  col[8].c  = 0; valA[8]  =  Kcont / hy;
          col[9].i = ex  ; col[9].j = ey   ; col[9].loc = ELEMENT;   col[9].c  = 0; valA[9]  = -Kcont / hy;
        }

        /* U_y interior equation */
        else {
          nEntries = 11;
          row.i    = ex  ; row.j     = ey  ; row.loc     = DOWN;     row.c     = 0;
          col[0].i = ex  ; col[0].j  = ey  ; col[0].loc  = DOWN;     col[0].c  = 0; valA[0]  = -2.0 * (etaDown + etaUp) / (hy*hy) - (etaLeft + etaRight) /(hx*hx);
          col[1].i = ex  ; col[1].j  = ey-1; col[1].loc  = DOWN;     col[1].c  = 0; valA[1]  =  2.0 * etaDown  / (hy*hy);
          col[2].i = ex  ; col[2].j  = ey+1; col[2].loc  = DOWN;     col[2].c  = 0; valA[2]  =  2.0 * etaUp    / (hy*hy);
          col[3].i = ex-1; col[3].j  = ey  ; col[3].loc  = DOWN;     col[3].c  = 0; valA[3]  =        etaLeft  / (hx*hx);
          col[4].i = ex+1; col[4].j  = ey  ; col[4].loc  = DOWN;     col[4].c  = 0; valA[4]  =        etaRight / (hx*hx);
          col[5].i = ex  ; col[5].j  = ey-1; col[5].loc  = LEFT;     col[5].c  = 0; valA[5]  =        etaLeft  / (hx*hy); /* down left x edge */
          col[6].i = ex  ; col[6].j  = ey-1; col[6].loc  = RIGHT;    col[6].c  = 0; valA[6]  = -      etaRight / (hx*hy); /* down right x edge */
          col[7].i = ex  ; col[7].j  = ey  ; col[7].loc  = LEFT;     col[7].c  = 0; valA[7]  = -      etaLeft  / (hx*hy); /* up left x edge */
          col[8].i = ex  ; col[8].j  = ey  ; col[8].loc  = RIGHT;    col[8].c  = 0; valA[8]  =        etaRight / (hx*hy); /* up right x edge */
          col[9].i = ex  ; col[9].j  = ey-1; col[9].loc  = ELEMENT;  col[9].c  = 0; valA[9]  =  Kcont / hy;
          col[10].i = ex ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0; valA[10] = -Kcont / hy;
        }

        /* Insert Y-momentum entries */
        ierr = DMStagMatSetValuesStencil(dm,A  ,1,&row,nEntries,col, valA  , INSERT_VALUES); CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm,rhs,1,&row,             &valRhs, INSERT_VALUES); CHKERRQ(ierr);
      }

      /* Right Boundary velocity Dirichlet */
      if (ex == Nx-1) {
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = Kbound;

        row.i = ex; row.j = ey; row.loc = RIGHT; row.c = 0;
        valRhs = 0.0;

        ierr   = DMStagMatSetValuesStencil(dm,A  ,1,&row,1,&row,&valA  ,INSERT_VALUES); CHKERRQ(ierr);
        ierr   = DMStagVecSetValuesStencil(dm,rhs,1,&row,       &valRhs,INSERT_VALUES); CHKERRQ(ierr);
      }

      /* Left velocity Dirichlet */
      if (ex == 0) {
        DMStagStencil row;
        PetscScalar   valRhs;
        const PetscScalar valA = Kbound;

        row.i = ex; row.j = ey; row.loc = LEFT; row.c = 0;
        valRhs = 0.0;

        ierr = DMStagMatSetValuesStencil(dm,A  ,1,&row,1,&row,&valA  ,INSERT_VALUES); CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm,rhs,1,&row,       &valRhs,INSERT_VALUES); CHKERRQ(ierr);
      }

      /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
      else {
        PetscInt nEntries;
        DMStagStencil row,col[11];
        PetscScalar   valRhs,valA[11];
        PetscScalar   etaLeft,etaRight,etaUp,etaDown;

        /* Get eta values */
        etaLeft  = ETA(x_center-hx,y_center);
        etaRight = ETA(x_center,y_center);
        etaUp    = ETA(x_left,y_up);
        etaDown  = ETA(x_left,y_down);

        /* Bottom boundary x velocity stencil (with zero vel deriv) */
        if (ey == 0) {
          nEntries = 10;
          row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0; valA[0]  = -2.0 * (etaLeft + etaRight) / (hx*hx) -(etaUp) / (hy*hy);
          /* Missing element below */
          col[1].i  = ex  ; col[1].j  = ey+1; col[1].loc  = LEFT;    col[1].c   = 0; valA[1]  =        etaUp    / (hy*hy);
          col[2].i  = ex-1; col[2].j  = ey  ; col[2].loc  = LEFT;    col[2].c   = 0; valA[2]  =  2.0 * etaLeft  / (hx*hx);
          col[3].i  = ex+1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0; valA[3]  =  2.0 * etaRight / (hx*hx);
          col[4].i  = ex-1; col[4].j  = ey  ; col[4].loc  = DOWN;    col[4].c   = 0; valA[4]  =        etaDown  / (hx*hy); /* down left */
          col[5].i  = ex  ; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0; valA[5]  = -      etaDown  / (hx*hy); /* down right */
          col[6].i  = ex-1; col[6].j  = ey  ; col[6].loc  = UP;      col[6].c   = 0; valA[6]  = -      etaUp    / (hx*hy); /* up left */
          col[7].i  = ex  ; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  =        etaUp    / (hx*hy); /* up right */
          col[8].i  = ex-1; col[8].j  = ey  ; col[8].loc  = ELEMENT; col[8].c   = 0; valA[8]  =  Kcont / hx;
          col[9].i = ex   ; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0; valA[9]  = -Kcont / hx;
          valRhs = 0.0;
        }

        /* Top boundary x velocity stencil */
        else if (ey == Ny-1) {
          nEntries = 10;
          row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0; valA[0]  = -2.0 * (etaLeft + etaRight) / (hx*hx) -(etaDown) / (hy*hy);
          col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc  = LEFT;    col[1].c   = 0; valA[1]  =        etaDown  / (hy*hy);
          /* Missing element above */
          col[2].i  = ex-1; col[2].j  = ey  ; col[2].loc  = LEFT;    col[2].c   = 0; valA[2]  =  2.0 * etaLeft  / (hx*hx);
          col[3].i  = ex+1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0; valA[3]  =  2.0 * etaRight / (hx*hx);
          col[4].i  = ex-1; col[4].j  = ey  ; col[4].loc  = DOWN;    col[4].c   = 0; valA[4]  =        etaDown  / (hx*hy); /* down left */
          col[5].i  = ex  ; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0; valA[5]  = -      etaDown  / (hx*hy); /* down right */
          col[6].i  = ex-1; col[6].j  = ey  ; col[6].loc  = UP;      col[6].c   = 0; valA[6]  = -      etaUp    / (hx*hy); /* up left */
          col[7].i  = ex  ; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  =        etaUp    / (hx*hy); /* up right */
          col[8].i  = ex-1; col[8].j  = ey  ; col[8].loc  = ELEMENT; col[8].c   = 0; valA[8]  =  Kcont / hx;
          col[9].i = ex   ; col[9].j  = ey   ; col[9].loc = ELEMENT;  col[9].c  = 0; valA[9]  = -Kcont / hx;
          valRhs = 0.0;
        }

        /* U_x interior equation */
        else {
          nEntries = 11;
          row.i     = ex  ; row.j     = ey  ; row.loc     = LEFT;    row.c      = 0;
          col[0].i  = ex  ; col[0].j  = ey  ; col[0].loc  = LEFT;    col[0].c   = 0; valA[0]  = -2.0 * (etaLeft + etaRight) / (hx*hx) -(etaUp + etaDown) / (hy*hy);
          col[1].i  = ex  ; col[1].j  = ey-1; col[1].loc  = LEFT;    col[1].c   = 0; valA[1]  =        etaDown  / (hy*hy);
          col[2].i  = ex  ; col[2].j  = ey+1; col[2].loc  = LEFT;    col[2].c   = 0; valA[2]  =        etaUp    / (hy*hy);
          col[3].i  = ex-1; col[3].j  = ey  ; col[3].loc  = LEFT;    col[3].c   = 0; valA[3]  =  2.0 * etaLeft  / (hx*hx);
          col[4].i  = ex+1; col[4].j  = ey  ; col[4].loc  = LEFT;    col[4].c   = 0; valA[4]  =  2.0 * etaRight / (hx*hx);
          col[5].i  = ex-1; col[5].j  = ey  ; col[5].loc  = DOWN;    col[5].c   = 0; valA[5]  =        etaDown  / (hx*hy); /* down left */
          col[6].i  = ex  ; col[6].j  = ey  ; col[6].loc  = DOWN;    col[6].c   = 0; valA[6]  = -      etaDown  / (hx*hy); /* down right */
          col[7].i  = ex-1; col[7].j  = ey  ; col[7].loc  = UP;      col[7].c   = 0; valA[7]  = -      etaUp    / (hx*hy); /* up left */
          col[8].i  = ex  ; col[8].j  = ey  ; col[8].loc  = UP;      col[8].c   = 0; valA[8]  =        etaUp    / (hx*hy); /* up right */
          col[9].i  = ex-1; col[9].j  = ey  ; col[9].loc  = ELEMENT; col[9].c   = 0; valA[9]  =  Kcont / hx;
          col[10].i = ex  ; col[10].j = ey  ; col[10].loc = ELEMENT; col[10].c  = 0; valA[10] = -Kcont / hx;
          valRhs = 0.0;
        }

        /* Insert X-momentum entries */
        ierr = DMStagMatSetValuesStencil(dm,A  ,1,&row,nEntries,col, valA  ,INSERT_VALUES); CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm,rhs,1,&row,             &valRhs,INSERT_VALUES); CHKERRQ(ierr);
      }

      /* P equation : u_x + v_y = 0 */
      /* Pin the first pressure node to zero, if requested */
      if (pinPressure && ex == 0 && ey == 0) {
        DMStagStencil row;
        PetscScalar valA,valRhs;

        row.i = ex; row.j = ey; row.loc = ELEMENT; row.c = 0;
        valA   = Kbound;
        valRhs = 0.0;

        ierr = DMStagMatSetValuesStencil(dm,A  ,1,&row,1,&row,&valA  ,INSERT_VALUES); CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm,rhs,1,&row,       &valRhs,INSERT_VALUES); CHKERRQ(ierr);
      }
      else {
        DMStagStencil row,col[5];
        PetscScalar   valA[5],valRhs;

        row.i    = ex; row.j    = ey; row.loc    = ELEMENT; row.c    = 0;
        col[0].i = ex; col[0].j = ey; col[0].loc = LEFT;    col[0].c = 0; valA[0] = -Kcont / hx;
        col[1].i = ex; col[1].j = ey; col[1].loc = RIGHT;   col[1].c = 0; valA[1] =  Kcont / hx;
        col[2].i = ex; col[2].j = ey; col[2].loc = DOWN;    col[2].c = 0; valA[2] = -Kcont / hy;
        col[3].i = ex; col[3].j = ey; col[3].loc = UP;      col[3].c = 0; valA[3] =  Kcont / hy;
        col[4] = row;                                                     valA[4] = 0.0;
        valRhs = 0.0;

        /* Insert P-equation entries */
        ierr = DMStagMatSetValuesStencil(dm,A  ,1,&row,5,col, valA  ,INSERT_VALUES); CHKERRQ(ierr);
        ierr = DMStagVecSetValuesStencil(dm,rhs,1,&row,      &valRhs,INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }

  /* Matrix and Vector Assembly */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
  ierr = VecAssemblyEnd  (rhs); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode MatCreatePreallocator_private(Mat A,Mat *p)
{
  Mat                    preallocator;
  PetscInt               M,N,m,n,bs;
  DM                     dm;
  ISLocalToGlobalMapping l2g[] = { NULL, NULL };
  PetscErrorCode         ierr;

  ierr = MatGetSize(A,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
  ierr = MatGetBlockSize(A,&bs);CHKERRQ(ierr);
  ierr = MatGetDM(A,&dm);CHKERRQ(ierr);
  ierr = MatGetLocalToGlobalMapping(A,&l2g[0],&l2g[1]);CHKERRQ(ierr);

  ierr = MatCreate(PetscObjectComm((PetscObject)A),&preallocator);CHKERRQ(ierr);
  ierr = MatSetType(preallocator,MATPREALLOCATOR);CHKERRQ(ierr);
  ierr = MatSetSizes(preallocator,m,n,M,N);CHKERRQ(ierr);
  ierr = MatSetBlockSize(preallocator,bs);CHKERRQ(ierr);
  ierr = MatSetDM(preallocator,dm);CHKERRQ(ierr);
  if (l2g[0] && l2g[1]) { ierr = MatSetLocalToGlobalMapping(preallocator,l2g[0],l2g[1]);CHKERRQ(ierr); }
  ierr = MatSetUp(preallocator);CHKERRQ(ierr);

  ierr = PetscObjectCompose((PetscObject)A,"__mat_preallocator__",(PetscObject)preallocator);CHKERRQ(ierr);
  if (p) {
    *p = preallocator;
  }
  PetscFunctionReturn(0);
}

/* may return a NULL pointer */
static PetscErrorCode MatGetPreallocator(Mat A,Mat *preallocator)
{
  PetscErrorCode ierr;
  Mat            p = NULL;

  ierr = PetscObjectQuery((PetscObject)A,"__mat_preallocator__",(PetscObject*)&p);CHKERRQ(ierr);
  *preallocator = p;
  PetscFunctionReturn(0);
}

/*
 Returns preallocator, a matrix of type "preallocator".
 The user should not call MatDestroy() on preallocator;
*/
PetscErrorCode MatPreallocatePhaseBegin(Mat A,Mat *preallocator)
{
  PetscErrorCode ierr;
  Mat            p = NULL;
  PetscInt       bs;

  ierr = MatGetPreallocator(A,&p);CHKERRQ(ierr);
  if (p) {
    ierr= MatDestroy(&p);CHKERRQ(ierr);
    p = NULL;
    ierr = PetscObjectCompose((PetscObject)A,"__mat_preallocator__",(PetscObject)p);CHKERRQ(ierr);
  }
  ierr = MatCreatePreallocator_private(A,&p);CHKERRQ(ierr);

  /* zap existing non-zero structure in A */
  /*
   It is a good idea to remove any exisiting non-zero structure in A to
   (i) reduce memory immediately
   (ii) to facilitate raising an error if someone trys to insert values into A after
   MatPreallocatorBegin() has been called - which signals they are doing something wrong/inconsistent
   */
  ierr = MatGetBlockSize(A,&bs);CHKERRQ(ierr);
  ierr = MatXAIJSetPreallocation(A,bs,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  *preallocator = p;
  PetscFunctionReturn(0);
}

static PetscErrorCode MatPreallocatePhaseEnd(Mat A)
{
  PetscErrorCode ierr;
  Mat            p = NULL;

  ierr = MatGetPreallocator(A,&p);CHKERRQ(ierr);
  if (!p) SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_USER,"Must call MatPreallocatorBegin() first");
  ierr = MatAssemblyBegin(p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* create new non-zero structure */
  ierr = MatPreallocatorPreallocate(p,PETSC_TRUE,A);CHKERRQ(ierr);

  /* clean up and remove the preallocator object from A */
  ierr= MatDestroy(&p);CHKERRQ(ierr);
  p = NULL;
  ierr = PetscObjectCompose((PetscObject)A,"__mat_preallocator__",(PetscObject)p);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
