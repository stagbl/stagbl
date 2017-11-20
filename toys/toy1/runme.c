static char help[] = "Naive 3-DMDA-based implementation of Gerya Ch. 7 Exercises\n\
                      -nx,-ny : number of cells in x and y direction\n\
                      -isoviscous : don't use variable viscosity\n";

#include <petsc.h>
#include "ctx.h"
#include "system.h"
#include "output.h"

static PetscErrorCode CreateGrids(Ctx);
static PetscErrorCode CreateMaterialProperties(Ctx);
static PetscErrorCode DumpOutput(Ctx,Vec);

/* =========================================================================== */
int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  Ctx            ctx;
  Mat            A;
  Vec            x,b;
  KSP            ksp;
  PC             pc;
  PetscMPIInt    size;

  /* --- Initialize and Create Context --------------------------------------- */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = CreateCtx(&ctx);CHKERRQ(ierr);

  /* --- Create DMDA objects and Composite DM -------------------------------- */
  ierr = CreateGrids(ctx);CHKERRQ(ierr);
  ierr = DMView(ctx->dmStokes,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); // print DM info

  /* --- Populate Material Properties ---------------------------------------- */
  ierr = CreateMaterialProperties(ctx);CHKERRQ(ierr);

  /* --- Define the Linear System -------------------------------------------- */
  ierr = CreateSystem(&A,&b,ctx);CHKERRQ(ierr);

  /* --- Solve the Linear System --------------------------------------------- */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERUMFPACK);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)x,"Solution");CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /* --- Simple Diagnostics -------------------------------------------------- */
  {
    PetscReal pMax,pMin,vxMax,vxMin,vyMax,vyMin;
    Vec       p,vx,vy;
    ierr = DMCompositeGetAccess(ctx->dmStokes,x,&p,&vx,&vy);CHKERRQ(ierr);
    ierr = VecMin(p,NULL,&pMin);CHKERRQ(ierr);
    ierr = VecMax(p,NULL,&pMax);CHKERRQ(ierr);
    ierr = VecMin(vx,NULL,&vxMin);CHKERRQ(ierr);
    ierr = VecMax(vx,NULL,&vxMax);CHKERRQ(ierr);
    ierr = VecMin(vy,NULL,&vyMin);CHKERRQ(ierr);
    ierr = VecMax(vy,NULL,&vyMax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"-- Info -----\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Kbound = %g\n",ctx->Kbound);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Kcont = %g\n",ctx->Kcont);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"hx = %g\n",ctx->hxCharacteristic);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"hy = %g\n",ctx->hyCharacteristic);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"-- Min / Max --\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"p' %g / %g \n",pMin,pMax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"p  %g / %g\n",ctx->Kcont*pMin,ctx->Kcont*pMax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"vx %g / %g\n",vxMin,vxMax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"vy %g / %g\n",vyMin,vyMax);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"-------------\n");CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccess(ctx->dmStokes,x,&p,&vx,&vy);CHKERRQ(ierr);
  }

  /* --- Output -------------------------------------------------------------- */

  // Dump raw matrix and vectors for debugging and scoping purposes
  ierr = OutputMatBinary(A,"A.petscbin");CHKERRQ(ierr);
  ierr = OutputVecBinary(x,"x.petscbin",PETSC_FALSE);CHKERRQ(ierr);
  ierr = OutputVecBinary(b,"b.petscbin",PETSC_FALSE);CHKERRQ(ierr);
  
  // Create xml/binary XDMFs
  { 
    Vec pprime,vX,vY;
    ierr = DMCompositeGetAccess(ctx->dmStokes,x,&pprime,&vX,&vY);CHKERRQ(ierr);
    {
      PetscViewer viewer;
      ierr = OutputDMCoordsBinary(ctx->daC,"coordsC.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFStart(ctx->daC,"C.xmf","coordsC.bin",&viewer);CHKERRQ(ierr);
      ierr = OutputVecBinary(pprime,"pprime.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFAddAttribute(ctx->daC,"pprime.bin","pprime",viewer);CHKERRQ(ierr);
      ierr = OutputVecBinary(ctx->etaC,"etaC.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFAddAttribute(ctx->daC,"etaC.bin","eta",viewer);CHKERRQ(ierr);
      ierr = DMDA2dXDMFFinish(&viewer);CHKERRQ(ierr);
    }
    {
      PetscViewer viewer;
      ierr = OutputDMCoordsBinary(ctx->daCor,"coordsCor.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFStart(ctx->daCor,"Cor.xmf","coordsCor.bin",&viewer);CHKERRQ(ierr);
      ierr = OutputVecBinary(ctx->etaC,"etaCor.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFAddAttribute(ctx->daCor,"etaCor.bin","eta",viewer);CHKERRQ(ierr);
      ierr = OutputVecBinary(ctx->rho,"rho.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFAddAttribute(ctx->daCor,"rho.bin","rho",viewer);CHKERRQ(ierr);
      ierr = DMDA2dXDMFFinish(&viewer);CHKERRQ(ierr);
    }
    {
      PetscViewer viewer;
      ierr = OutputDMCoordsBinary(ctx->daX,"coordsX.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFStart(ctx->daX,"X.xmf","coordsX.bin",&viewer);CHKERRQ(ierr);
      ierr = OutputVecBinary(vX,"vX.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFAddAttribute(ctx->daX,"vX.bin","vX",viewer);CHKERRQ(ierr);
      ierr = DMDA2dXDMFFinish(&viewer);CHKERRQ(ierr);
    }
    {
      PetscViewer viewer;
      ierr = OutputDMCoordsBinary(ctx->daY,"coordsY.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFStart(ctx->daY,"Y.xmf","coordsY.bin",&viewer);CHKERRQ(ierr);
      ierr = OutputVecBinary(vY,"vY.bin",PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMDA2dXDMFAddAttribute(ctx->daY,"vY.bin","vY",viewer);CHKERRQ(ierr);
      ierr = DMDA2dXDMFFinish(&viewer);CHKERRQ(ierr);
    }
    ierr = DMCompositeRestoreAccess(ctx->dmStokes,x,&pprime,&vX,&vY);CHKERRQ(ierr);
  }

  // Dump our main output XDMF file
  ierr = DumpOutput(ctx,x);CHKERRQ(ierr);

  /* --- Clean up and Finalize ----------------------------------------------- */
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = DestroyCtx(&ctx);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/* =========================================================================== */
PetscErrorCode CreateGrids(Ctx ctx)
{
  PetscErrorCode ierr;
  PetscReal hx,hy,xminCell,xmaxCell,yminCell,ymaxCell; 

  PetscFunctionBeginUser;

  // Default grid size (cells in the physical domain)
  ctx->MC = 30;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&ctx->MC,NULL);CHKERRQ(ierr);
  ctx->NC = 20;
  ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&ctx->NC,NULL);CHKERRQ(ierr);

  // Default grid dimensions and helpers
  ctx->xmin = 0.0;
  ctx->xmax = 1e6;
  ierr = PetscOptionsGetReal(NULL,NULL,"-xmax",&ctx->xmax,NULL);CHKERRQ(ierr);
  ctx->ymin = 0.0;
  ctx->ymax = 1.5e6;
  ierr = PetscOptionsGetReal(NULL,NULL,"-ymax",&ctx->ymax,NULL);CHKERRQ(ierr);
  hx = (ctx->xmax - ctx->xmin)/ctx->MC;
  ctx->hxCharacteristic = hx;
  hy = (ctx->ymax - ctx->ymin)/ctx->NC;
  ctx->hyCharacteristic = hy;
  ctx->MCor = ctx->MC+1;
  ctx->NCor = ctx->NC+1;
  ctx->MX   = ctx->MC+1;
  ctx->NX   = ctx->NC;
  ctx->MY   = ctx->MC;
  ctx->NY   = ctx->NC+1;

  if (ctx->MC < 3 || ctx->NC < 3) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Grids with less than three cells in any direction are not supported"); 

  // Which pressure node to fix
  if (ctx->MC > 2 && ctx->NC > 3) {
    ctx->pFixX = 1; ctx->pFixY = 2; /* To correspond with Chapter 7 solutions */ 
  } else {
    ctx->pFixX = 0; ctx->pFixY = 1;
  }

  // (UGLY. depends on ordering of DMs in Composite and the indexing of a DMDA)
  ctx->offsetC = 0;
  ctx->offsetX = ctx->offsetC + ctx->MC*ctx->NC;
  ctx->offsetY = ctx->offsetX + ctx->MX*ctx->NX;

  xminCell = ctx->xmin + hx/2.0; xmaxCell= ctx->xmax - hx/2.0;
  yminCell = ctx->ymin + hy/2.0; ymaxCell= ctx->ymax - hy/2.0;

  ierr = DMCompositeCreate(PETSC_COMM_WORLD,&ctx->dmStokes);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)ctx->dmStokes,"Composite DM");CHKERRQ(ierr);

  // Cell-centered grid
  ierr = DMDACreate2d(
      PETSC_COMM_WORLD,
      DM_BOUNDARY_GHOSTED,
      DM_BOUNDARY_GHOSTED,
      DMDA_STENCIL_BOX,
      ctx->MC,            // number of cells (x)
      ctx->NC,            // number of cells (y)
      PETSC_DECIDE,       // Local size (x)
      PETSC_DECIDE,       // Local size (y)
      1,                  // One dof (pressure)
      1,                  // Stencil Width of 1
      NULL,
      NULL,
      &ctx->daC);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)ctx->daC,"Center DMDA");CHKERRQ(ierr);
  ierr = DMSetUp(ctx->daC);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(ctx->daC,xminCell,xmaxCell,yminCell,ymaxCell,0,0);CHKERRQ(ierr); // After SetUp 
  ierr = DMCompositeAddDM(ctx->dmStokes,ctx->daC);CHKERRQ(ierr);

  // Corners grid
  ierr = DMDACreate2d(
      PETSC_COMM_WORLD,
      DM_BOUNDARY_NONE,
      DM_BOUNDARY_NONE,
      DMDA_STENCIL_BOX,
      ctx->MCor,          // number of nodes (x)
      ctx->NCor,          // number of nodes (y)
      PETSC_DECIDE,       // Local size (x)
      PETSC_DECIDE,       // Local size (y)
      1,                  // One dof 
      1,                  // Stencil Width of 1
      NULL,
      NULL,
      &ctx->daCor);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)ctx->daCor,"Corner DMDA");CHKERRQ(ierr);
  ierr = DMSetUp(ctx->daCor);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(ctx->daCor,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0,0);CHKERRQ(ierr); // After SetUp 
  // (not part of the DMComposite, as we don't solve for quanitites on the corners)

  // Edges grids
  ierr = DMDACreate2d(
      PETSC_COMM_WORLD,
      DM_BOUNDARY_GHOSTED,
      DM_BOUNDARY_GHOSTED,
      DMDA_STENCIL_BOX,
      ctx->MX,            // number of nodes (x)
      ctx->NX,            // number of cells (y)
      PETSC_DECIDE,       // Local size (x)
      PETSC_DECIDE,       // Local size (y)
      1,                  // One dof 
      1,                  // Stencil Width of 1
      NULL,
      NULL,
      &ctx->daX);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)ctx->daC,"X-edges DMDA");CHKERRQ(ierr);
  ierr = DMSetUp(ctx->daX);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(ctx->daX,ctx->xmin,ctx->xmax,yminCell,ymaxCell,0,0);CHKERRQ(ierr); // After SetUp 
  ierr = DMCompositeAddDM(ctx->dmStokes,ctx->daX);CHKERRQ(ierr);

  ierr = DMDACreate2d(
      PETSC_COMM_WORLD,
      DM_BOUNDARY_GHOSTED,
      DM_BOUNDARY_GHOSTED,
      DMDA_STENCIL_BOX,
      ctx->MY,            // number of nodes (x)
      ctx->NY,            // number of cells (y)
      PETSC_DECIDE,       // Local size (x)
      PETSC_DECIDE,       // Local size (y)
      1,                  // One dof 
      1,                  // Stencil Width of 1
      NULL,
      NULL,
      &ctx->daY);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)ctx->daY,"Y-edges DMDA");CHKERRQ(ierr);
  ierr = DMSetUp(ctx->daY);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(ctx->daY,xminCell,xmaxCell,ctx->ymin,ctx->ymax,0,0);CHKERRQ(ierr); // After SetUp 
  ierr = DMCompositeAddDM(ctx->dmStokes,ctx->daY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* =========================================================================== */
PetscErrorCode DumpOutput(Ctx ctx,Vec x)
{
  PetscErrorCode ierr;
  DM             dmStokes,daC,daC2,daX,daY;
  Vec            vInterp,vX,vY,vC;
  PetscScalar    **vXArr,**vYArr,***vInterpArr;
  PetscInt       i,j,si,sj,ni,nj;

  PetscFunctionBeginUser;

  dmStokes = ctx->dmStokes;
  ierr = DMCompositeGetEntries(dmStokes,&daC,&daX,&daY);CHKERRQ(ierr);

  // Create a new DMDA, like daC, but with 2 dof
  ierr = DMDAGetReducedDMDA(daC,2,&daC2);CHKERRQ(ierr);
  ierr = DMView(daC2,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMSetUp(daC2);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)daC2,"Cell-centered Interpolation DMDA");CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(daC2,&vInterp);CHKERRQ(ierr);
  {
    Vec vXLocal,vYLocal;

    ierr = DMCompositeGetAccess(dmStokes,x,&vC,&vX,&vY);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(daX,&vXLocal);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(daY,&vYLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(daX,vX,INSERT_VALUES,vXLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(daY,vY,INSERT_VALUES,vYLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(daX,vX,INSERT_VALUES,vXLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(daY,vY,INSERT_VALUES,vYLocal);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(daX,vXLocal,&vXArr);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(daY,vYLocal,&vYArr);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(daC2,vInterp,&vInterpArr);CHKERRQ(ierr);
    ierr = DMDAGetCorners(daC,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
    // Note : this has NOT been tested in parallel
    for (i=si;i<si+ni;++i) {
      for (j=sj;j<sj+nj;++j) {
        vInterpArr[j][i][0] = (vXArr[j][i] + vXArr[j][i+1])*0.5;
        vInterpArr[j][i][1] = (vYArr[j][i] + vYArr[j+1][i])*0.5;
      }
    }
    ierr = DMDAVecRestoreArrayRead(daX,vXLocal,&vXArr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(daY,vYLocal,&vYArr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(daC2,vInterp,&vInterpArr);CHKERRQ(ierr);
    ierr = VecDestroy(&vXLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&vYLocal);CHKERRQ(ierr);
  }

  // Dump headerless binary of cell-centered 2-dof field
  {
    Vec         vInterp_natural;
    PetscViewer viewer;

    ierr = DMDACreateNaturalVector(daC2,&vInterp_natural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalBegin(daC2,vInterp,INSERT_VALUES,vInterp_natural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalEnd(daC2,vInterp,INSERT_VALUES,vInterp_natural);CHKERRQ(ierr);
    ierr = PetscViewerCreate(PetscObjectComm((PetscObject)daC2),&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
    ierr = PetscViewerBinarySetSkipHeader(viewer,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer,"vInterp.bin");CHKERRQ(ierr);
    ierr = VecView(vInterp_natural,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = VecDestroy(&vInterp_natural);CHKERRQ(ierr);
  }
  ierr = DMCompositeRestoreAccess(ctx->dmStokes,x,&vInterp,NULL,NULL);CHKERRQ(ierr);

  // Scale Pressure field to dimensional units and dump as headerless binary
  {
    Vec vCNatural;

    ierr = DMDACreateNaturalVector(daC,&vCNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalBegin(daC,vC,INSERT_VALUES,vCNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalEnd(daC,vC,INSERT_VALUES,vCNatural);CHKERRQ(ierr);
    ierr = VecScale(vCNatural,ctx->Kcont);CHKERRQ(ierr);
    ierr = OutputVecBinary(vCNatural,"pressure.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecDestroy(&vCNatural);CHKERRQ(ierr);
  }

  // Dump xdmf of cell-centered fields
  // Assumes that the required .bin files have already been generated
  {
    PetscViewer viewer;
    ierr = DMDA2dXDMFStart(ctx->daC,"paraviewMe.xmf","coordsC.bin",&viewer);CHKERRQ(ierr);
    ierr = DMDA2dXDMFAddAttribute(daC2,"vinterp.bin","vinterp",viewer);CHKERRQ(ierr);
    ierr = DMDA2dXDMFAddAttribute(ctx->daC,"pressure.bin","p",viewer);CHKERRQ(ierr);
    ierr = DMDA2dXDMFAddAttribute(ctx->daC,"etaC.bin","eta",viewer);CHKERRQ(ierr);
    ierr = DMDA2dXDMFFinish(&viewer);CHKERRQ(ierr);
  }

  // Clean Up
  ierr = DMCompositeRestoreAccess(dmStokes,x,&vC,&vX,&vY);CHKERRQ(ierr);
  ierr = VecDestroy(&vInterp);CHKERRQ(ierr);
  ierr = DMDestroy(&daC2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ========================================================================== */
static PetscReal getRho(Ctx ctx,PetscReal x,PetscReal y)
{
  if (x < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->rho1;
  } else {
    return ctx->rho2;
  }
}

/* ========================================================================== */
static PetscReal getEta(Ctx ctx,PetscReal x,PetscReal y)
{
  if (ctx->isoviscous ) return ctx->etaCharacteristic;
  if (x < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->eta1;
  } else {
    return ctx->eta2;
  }
}

/* =========================================================================== */
PetscErrorCode CreateMaterialProperties(Ctx ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  {
    PetscScalar **arr;
    PetscInt    si,sj,ni,nj,i,j;
    DM          dmc;
    Vec         coords;
    DMDACoor2d  **arrCoords;

    ierr = DMGetCoordinateDM(ctx->daCor,&dmc);CHKERRQ(ierr);
    ierr = DMGetCoordinates(ctx->daCor,&coords);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dmc,coords,&arrCoords);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(ctx->daCor,&ctx->rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(ctx->daCor,ctx->rho,&arr);CHKERRQ(ierr);
    ierr = DMDAGetCorners(ctx->daCor,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
    for (i=si; i<si+ni;++i){
      for (j=sj; j<sj+nj;++j){
        arr[j][i] = getRho(ctx,arrCoords[j][i].x,arrCoords[j][i].y);
      }
    }
    ierr = DMDAVecRestoreArray(ctx->daCor,ctx->rho,&arr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(dmc,coords,&arrCoords);CHKERRQ(ierr);
  }
  {
    PetscScalar **arr;
    PetscInt    si,sj,ni,nj,i,j;
    DM          dmc;
    Vec         coords;
    DMDACoor2d  **arrCoords;

    ierr = DMGetCoordinateDM(ctx->daC,&dmc);CHKERRQ(ierr);
    ierr = DMGetCoordinates(ctx->daC,&coords);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dmc,coords,&arrCoords);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(ctx->daC,&ctx->etaC);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(ctx->daC,ctx->etaC,&arr);CHKERRQ(ierr);
    ierr = DMDAGetCorners(ctx->daC,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
    for (i=si; i<si+ni;++i){
      for (j=sj; j<sj+nj;++j){
        arr[j][i] = getEta(ctx,arrCoords[j][i].x,arrCoords[j][i].y);
      }
    }
    ierr = DMDAVecRestoreArray(ctx->daC,ctx->etaC,&arr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(dmc,coords,&arrCoords);CHKERRQ(ierr);
  }
  {
    PetscScalar **arr;
    PetscInt    si,sj,ni,nj,i,j;
    DM          dmc;
    Vec         coords;
    DMDACoor2d  **arrCoords;

    ierr = DMGetCoordinateDM(ctx->daCor,&dmc);CHKERRQ(ierr);
    ierr = DMGetCoordinates(ctx->daCor,&coords);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dmc,coords,&arrCoords);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(ctx->daCor,&ctx->etaCor);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(ctx->daCor,ctx->etaCor,&arr);CHKERRQ(ierr);
    ierr = DMDAGetCorners(ctx->daCor,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
    for (i=si; i<si+ni;++i){
      for (j=sj; j<sj+nj;++j){
        arr[j][i] = getEta(ctx,arrCoords[j][i].x,arrCoords[j][i].y);
      }
    }
    ierr = DMDAVecRestoreArray(ctx->daCor,ctx->etaCor,&arr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(dmc,coords,&arrCoords);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
