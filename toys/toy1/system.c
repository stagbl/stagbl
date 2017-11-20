#include "system.h"
#include "output.h"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>

static PetscErrorCode AssembleStokesOperator(Mat,Ctx);
static PetscErrorCode AssembleRHS(Vec,Ctx);

/* =========================================================================== */
PetscErrorCode CreateSystem(Mat *pA,Vec *pb, Ctx ctx)
{
  PetscErrorCode ierr;
  Mat            A;
  Vec            b;

  PetscFunctionBeginUser;

  // Scaling constants
  {
    PetscReal hxAvgInv = 2.0/(ctx->hxCharacteristic + ctx->hyCharacteristic);
    ctx->Kcont         = ctx->etaCharacteristic*hxAvgInv;
    ctx->Kbound        = ctx->etaCharacteristic*hxAvgInv*hxAvgInv;
  }

  // RHS
  ierr = DMCreateGlobalVector(ctx->dmStokes,pb);CHKERRQ(ierr);
  b = *pb;
  ierr = PetscObjectSetName((PetscObject)b,"RHS");CHKERRQ(ierr);
  ierr = AssembleRHS(b,ctx);CHKERRQ(ierr);

  // Operator
  ierr = DMCreateMatrix(ctx->dmStokes,pA);CHKERRQ(ierr); // Or use known nonzero pattern
  A = *pA;
  ierr = PetscObjectSetName((PetscObject)A,"Stokes");CHKERRQ(ierr);
  ierr = AssembleStokesOperator(A,ctx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* =========================================================================== */
// Assemble a monolithic Stokes operator
// This routine shamefully does index arithmetic instead of trying to use the PETSc API for everything
// it does not work in parallel because the equations are numbered differently.
static PetscErrorCode AssembleStokesOperator(Mat A,Ctx ctx)
{
  PetscErrorCode         ierr;
  DM                     dmStokes,daC,daX,daY;
  PetscInt               MC,NC,MX,NX,MY,NY,offsetC,offsetX,offsetY;
  PetscReal              Kbound,Kcont;

  PetscFunctionBeginUser;
  {
    PetscMPIInt size;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Assembled operator not implemented in parallel");
  }
  dmStokes = ctx->dmStokes;
  Kbound = ctx->Kbound;
  Kcont = ctx->Kcont;
  offsetC = ctx->offsetC; offsetX = ctx->offsetX; offsetY = ctx->offsetY;
  MC = ctx->MC; NC=ctx->NC; MX=ctx->MX; NX=ctx->NX; MY=ctx->MY; NY=ctx->NY;

  ierr = DMCompositeGetEntries(ctx->dmStokes,&daC,&daX,&daY);CHKERRQ(ierr);

  // Preallocation
  ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
  // TODO  preallocation (and get rid of getting matrix from DM)

  // Mass
  {
    PetscInt i,j,si,sj,ni,nj;
    PetscMPIInt rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = DMDAGetCorners(daC,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
    for (i=si;i<si+ni;++i) {
      for (j=sj;j<sj+nj;++j) {
        PetscInt row,cols[4];
        PetscScalar vals[4];

        row = i + MC*j + offsetC;
        /* We impose 4 pressure symmetry BCs, at the corners, 
           and a fixed node, as implied by Exercise 7.1.
           Using an implied ghost velocity node (as also described in the text)
           seems cleaner, though */
        if (i==ctx->pFixX && j==ctx->pFixY){
          /* Fix pressure node */
          cols[0] = row;
          vals[0] = Kbound;
          ierr = MatSetValues(A,1,&row,1,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        } else if (i == 0 && j == 0) {
          cols[0] = row;
          cols[1] = (i+1) + MC*j + offsetC;
          vals[0] = Kbound;
          vals[1] = -Kbound;
          ierr = MatSetValues(A,1,&row,2,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        } else if (i == 0 && j == NC-1) {
          cols[0] = row;
          cols[1] = (i+1) + MC*j + offsetC;
          vals[0] = Kbound;
          vals[1] = -Kbound;
          ierr = MatSetValues(A,1,&row,2,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        } else if (i == MC-1 && j == 0) {
          cols[0] = (i-1) + MC*j + offsetC;
          cols[1] = row;
          vals[0] = -Kbound;
          vals[1] = Kbound;
          ierr = MatSetValues(A,1,&row,2,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        } else if (i == MC-1 && j == NC-1) {
          cols[0] = (i-1) + MC*j + offsetC;
          cols[1] = row;
          vals[0] = -Kbound;
          vals[1] = Kbound;
          ierr = MatSetValues(A,1,&row,2,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          // nat_A(i,j) is the global eqn number for element i,j on grid A
          // row nat_C(i,j) + pressure grid offset
          // cols:
          //   nat_X(i,j)   + x grid offset
          const PetscReal hx = ctx->hxCharacteristic; // assume constant
          const PetscReal hy = ctx->hyCharacteristic; // assume constant
          cols[0] = i + MX*j + offsetX;
          vals[0] = -Kcont/hx;
          //   nat_X(i+1,j) + x grid offset
          cols[1] = i+1 + MX*j + offsetX;
          vals[1] = Kcont/hx;
          //   nat_Y(i,j)   + y grid offset
          cols[2] = i + MY*j + offsetY;
          vals[2] = -Kcont/hy;
          //   nat_Y(i,j+1) + y grid offset
          cols[3] = i + MY*(j+1) + offsetY;
          vals[3] = Kcont/hy;
          ierr = MatSetValues(A,1,&row,4,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
  }

  // X Momentum
  {
    PetscInt    i,j,si,sj,ni,nj;
    PetscMPIInt rank;
    PetscScalar **arrEtaC,**arrEtaCor;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = DMDAGetCorners(daX,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(ctx->daC,ctx->etaC,&arrEtaC);CHKERRQ(ierr); // Will not work in parallel
    ierr = DMDAVecGetArrayRead(ctx->daCor,ctx->etaCor,&arrEtaCor);CHKERRQ(ierr); // Will not work in parallel
    for (i=PetscMax(si,1);i<PetscMin(si+ni,MX-1);++i) {
      for (j=PetscMax(sj,1);j<PetscMin(sj+nj,NX-1);++j) {
        PetscInt row,cols[11]; 
        PetscScalar vals[11];

        const PetscReal hx       = ctx->hxCharacteristic; // assume constant
        const PetscReal hy       = ctx->hyCharacteristic; // assume constant
        const PetscReal hxi2     = 1.0/(hx*hx);
        const PetscReal hyi2     = 1.0/(hy*hy);
        const PetscReal hxhyi    = 1.0/(hx*hy);
        const PetscReal etaRight = arrEtaC[j][i];
        const PetscReal etaLeft  = arrEtaC[j][i-1];
        const PetscReal etaUp    = arrEtaCor[j+1][i];
        const PetscReal etaDown  = arrEtaCor[j][i];
        row = i + MX*j + offsetX;
        // Couples to (see Fig 7.19) :
        //   pressures on each side : (i-1,j) and (i,j) in C grid 
        cols[0] = (i-1) + MC*j     + offsetC;
        vals[0] = Kcont/hx;
        cols[1] = i     + MC*j     + offsetC;
        vals[1] = -Kcont/hx;
        // 5-point stencil in X-velocity grid (i,j),(i,j\pm 1),(i\pm 1,j) in X grid
        // Note: my i's and j's are backwards from Figure 7.19
        cols[2] = i     + MX*(j-1) + offsetX;
        vals[2] = etaDown*hyi2;
        cols[3] = (i-1) + MX*j     + offsetX;
        vals[3] = 2.0*etaLeft*hxi2;
        cols[4] = i     + MX*j     + offsetX;
        vals[4] = -2.0*(etaLeft+etaRight)*hxi2 -(etaUp+etaDown)*hyi2;
        cols[5] = (i+1) + MX*j     + offsetX;
        vals[5] = 2.0*etaRight*hxi2;
        cols[6] = i     + MX*(j+1) + offsetX;
        vals[6] = etaUp*hyi2;
        // Neighboring 4 points in Y-velocity grid (i-1,j) (i,j) (i-1,j+1) (i,j+1) in Y grid 
        // Note: the indices of these points differ from those in Figure 7.19,
        // because we use i for x and j for y
        cols[7] = (i-1) + MY*j     + offsetY;
        vals[7] = etaDown*hxhyi;
        cols[8] = i     + MY*j     + offsetY;
        vals[8] = -etaDown*hxhyi;
        cols[9] = (i-1) + MY*(j+1) + offsetY;
        vals[9] = -etaUp*hxhyi;
        cols[10]= i     + MY*(j+1) + offsetY;
        vals[10]= etaUp*hxhyi;
        ierr = MatSetValues(A,1,&row,11,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    // Boundary Conditions : X velocity
    // Zero boundary conditions on left and right boundaries - insert Kbound on diagonal
    // TODO - this is doing redundant work if you have several procs..
    for (j=sj;j<sj+nj;++j) {
      PetscScalar val = Kbound;
      PetscInt row;
      i = 0;
      row = i + MX*j + offsetX;
      ierr = MatSetValues(A,1,&row,1,&row,&val,INSERT_VALUES);CHKERRQ(ierr);
      i = MX-1;
      row = i + MX*j + offsetX;
      ierr = MatSetValues(A,1,&row,1,&row,&val,INSERT_VALUES);CHKERRQ(ierr);
    }
    // Zero gradient functions on the top and bottom
    // TODO - this is doing redundant work if you have several procs..
    for (i=PetscMax(si,1);i<PetscMin(si+ni,MX-1);++i) {
        PetscInt row,cols[2];
        PetscScalar vals[2];
        j = 0;
        row = i + MX*j + offsetX;
        cols[0] = row;
        vals[0] = Kbound;
        cols[1] = i + MX*(j+1) + offsetX;
        vals[1] = -Kbound;
        ierr = MatSetValues(A,1,&row,2,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        j = (NX-1);
        row = i + MX*j + offsetX;
        cols[0] = i + MX*(j-1) + offsetX;
        vals[0] = -Kbound;
        cols[1] = row;
        vals[1] = Kbound;
        ierr = MatSetValues(A,1,&row,2,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = DMDAVecRestoreArrayRead(ctx->daC,ctx->etaC,&arrEtaC);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(ctx->daCor,ctx->etaCor,&arrEtaCor);CHKERRQ(ierr);
  }

  // Y Momentum
  {
    PetscInt    i,j,si,sj,ni,nj;
    PetscMPIInt rank;
    PetscScalar **arrEtaC,**arrEtaCor;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = DMDAVecGetArrayRead(ctx->daC,ctx->etaC,&arrEtaC);CHKERRQ(ierr); // Note this won't work in parallel (need global->local)
    ierr = DMDAVecGetArrayRead(ctx->daCor,ctx->etaCor,&arrEtaCor);CHKERRQ(ierr); // Note this won't work in parallel!
    ierr = DMDAGetCorners(daY,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
    for (i=PetscMax(si,1);i<PetscMin(si+ni,MY-1);++i) {
      for (j=PetscMax(sj,1);j<PetscMin(sj+nj,NY-1);++j) {
        PetscInt row,cols[11]; 
        PetscScalar vals[11];
        const PetscReal hx = ctx->hxCharacteristic; // assume constant
        const PetscReal hy = ctx->hyCharacteristic; // assume constant
        const PetscReal hxi2 = 1.0/(hx*hx);
        const PetscReal hyi2 = 1.0/(hy*hy);
        const PetscReal hxhyi = 1.0/(hx*hy);
        const PetscReal etaUp = arrEtaC[j][i];
        const PetscReal etaDown = arrEtaC[j-1][i];
        const PetscReal etaLeft = arrEtaCor[j][i];
        const PetscReal etaRight = arrEtaCor[j][i+1];
        // TODO : check use of these etas, as our test case wouldn't reveal all issues
        row = i + MY*j + offsetY;
        // Couples to (see Fig 7.19) :
        //  pressure above and below (i,j-1) (i,j) in C grid
        cols[0] = i     + MC*(j-1) + offsetC;
        vals[0] = Kcont/hy;
        cols[1] = i     + MC*j     + offsetC;
        vals[1] = -Kcont/hy;
        // 4 neighboring points in X-velocity grid
        //  (i,j-1) (i+1,j-1) (i,j) (i+1,j)
        cols[2] = i     + MX*(j-1)  + offsetX;
        vals[2] = etaLeft*hxhyi;
        cols[3] = (i+1) + MX*(j-1) + offsetX;
        vals[3] = -etaRight*hxhyi;
        cols[4] = i     + MX*j     + offsetX;
        vals[4] = -etaLeft*hxhyi;
        cols[5] = (i+1) + MX*j     + offsetX;
        vals[5] = etaRight*hxhyi;
        // 5-point stencil of Y-velocities (i,j) (i ,j\pm 1) (i \pm 1,j) in Y grid
        cols[6] = i     + MY*(j-1) + offsetY;
        vals[6] = 2.0*etaDown*hyi2;
        cols[7] = (i-1) + MY*j     + offsetY;
        vals[7] = etaLeft*hxi2;
        cols[8] = i     + MY*j     + offsetY;
        vals[8] = -2.0*(etaUp+etaDown)*hyi2 -(etaLeft+etaRight)*hxi2;
        cols[9] = (i+1) + MY*j     + offsetY;
        vals[9] = etaRight*hxi2;
        cols[10]= i     + MY*(j+1) + offsetY;
        vals[10]= 2.0*etaUp*hyi2;
        ierr = MatSetValues(A,1,&row,11,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    // Boundary Conditions: Y velocity
    // TODO - this is doing redundant work on multiple procs
    // Zero boundary conditions on top and bottom
    for (i=si;i<si+ni;++i) {
      PetscScalar val = Kbound;
      PetscInt row;
      j = 0;
      row = i + MY*j + offsetY;
      ierr = MatSetValues(A,1,&row,1,&row,&val,INSERT_VALUES);CHKERRQ(ierr);
      j = NY-1;
      row = i + MY*j + offsetY;
      ierr = MatSetValues(A,1,&row,1,&row,&val,INSERT_VALUES);CHKERRQ(ierr);
    }
    // Zero-gradient boundary conditions on left and right
    // TODO - this is doing redundant work
    for (j=PetscMax(sj,1);j<PetscMin(sj+nj,NY-1);++j) {
        PetscInt row,cols[2];
        PetscScalar vals[2];
        i = 0;
        row = i + MY*j + offsetY;
        cols[0] = row;
        vals[0] = Kbound;
        cols[1] = i+1 + MY*j + offsetY;
        vals[1] = -Kbound;
        ierr = MatSetValues(A,1,&row,2,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        i = MY-1;
        row = i + MY*j + offsetY;
        cols[0] = (i-1) + MY*j + offsetY;
        vals[0] = -Kbound;
        cols[1] = row;
        vals[1] = Kbound;
        ierr = MatSetValues(A,1,&row,2,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = DMDAVecRestoreArrayRead(ctx->daC,ctx->etaC,&arrEtaC);CHKERRQ(ierr); // Note this won't work in parallel (need global->local)
    ierr = DMDAVecRestoreArrayRead(ctx->daCor,ctx->etaCor,&arrEtaCor);CHKERRQ(ierr); // Note this won't work in parallel!
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* =========================================================================== */
PetscErrorCode AssembleRHS(Vec b,Ctx ctx)
{
  PetscErrorCode ierr;
  Vec             bY;
  PetscReal       Kbound,Kcont;

  PetscFunctionBegin;
  Kbound = ctx->Kbound;
  Kcont = ctx->Kcont;
  ierr = VecZeroEntries(b);CHKERRQ(ierr); //slightly inefficient
  ierr = DMCompositeGetAccess(ctx->dmStokes,b,NULL,NULL,&bY);CHKERRQ(ierr);

  // Mass (all zero)

  // X mom (all zero with no gx component)
  if (ctx->gx != 0.0 ) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Not implemented for non-zero gx");

  // Y mom
  {
    PetscScalar **vArr,**rhoArr;
    PetscInt    si,sj,ni,nj,i,j;
    ierr = DMDAVecGetArray(ctx->daY,bY,&vArr);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(ctx->daCor,ctx->rho,&rhoArr);CHKERRQ(ierr);
    ierr = DMDAGetCorners(ctx->daY,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
    for (i=PetscMax(si,1); i<PetscMin(si+ni,ctx->MY-1);++i){
      for (j=PetscMax(sj,1); j<PetscMin(sj+nj,ctx->NY-1);++j){
        vArr[j][i] = -ctx->gy * 0.5 * (rhoArr[j][i] + rhoArr[j][i+1]);
      }
    }
    ierr = DMDAVecRestoreArrayRead(ctx->daCor,ctx->rho,&rhoArr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(ctx->daY,bY,&vArr);CHKERRQ(ierr);
  }

  ierr = DMCompositeRestoreAccess(ctx->dmStokes,b,NULL,NULL,&bY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
