static char help[] = "Prototype DMStag implementation of Gerya Ch. 7 Exercises\n\
                      -nx,-ny : number of cells in x and y direction\n\
                      -isoviscous : don't use variable viscosity\n";

#include <petsc.h>
#include "ctx.h"
#include "system.h"
#include "output.h"
#include "dmstag.h"
#include "system.h"

/* =========================================================================== */
int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  Ctx            ctx;
  DM             stokesGrid,paramGrid,elementOnlyGrid,vertexOnlyGrid;
  Vec            paramLocal;
  Mat            A;
  Vec            b;
  Vec            x;
  KSP            ksp;
  PC             pc;
  PetscMPIInt    size;
  PetscScalar    pMin,pMax,pPrimeMin,pPrimeMax,vxInterpMin,vxInterpMax,vyInterpMin,vyInterpMax;

  /* --- Initialize and Create Context --------------------------------------- */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = CreateCtx(&ctx);CHKERRQ(ierr);

  /* --- Create DMStag Objects ----------------------------------------------- */
  ierr = DMRegister(DMSTAG,DMCreate_Stag);CHKERRQ(ierr);

  // A DMStag to hold the unknowns to solve the Stokes problem
  ierr = DMStagCreate2d(PETSC_COMM_WORLD, 
      DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,       // no special boundary support yet
      ctx->M,ctx->N,PETSC_DECIDE,PETSC_DECIDE, // sizes provided as DMDA
      0,1,1,                                   // no dof per vertex: 0 per vertex, 1 per edge, 1 per face/element
      DMSTAG_GHOST_STENCIL_BOX,                // elementwise stencil pattern
      1,                                       // elementwise stencil width
      &ctx->stokesGrid);
  stokesGrid = ctx->stokesGrid;
  ierr = DMSetUp(stokesGrid);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinates(stokesGrid,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr); 

  // A DMStag to hold the coefficient fields for the Stokes problem: 1 dof per vertex and two per element/face
  ierr = DMStagCreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,ctx->M,ctx->N,PETSC_DECIDE,PETSC_DECIDE,1,0,2,DMSTAG_GHOST_STENCIL_BOX,1,&ctx->paramGrid);CHKERRQ(ierr);
  paramGrid = ctx->paramGrid;
  ierr = DMSetUp(paramGrid);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinates(paramGrid,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr); 

  /* --- Set up Problem ------------------------------------------------------ */
  {
    typedef struct {PetscReal xCorner,yCorner,xElement,yElement;} CoordinateData; // Definitely need to provide these types somehow to the user.. (specifying them wrong leads to very annoying bugs)
    typedef struct {PetscScalar etaCorner,etaElement,rho;} ElementData; // ? What's the best way to have users to do this? Provide these before hand? 
    DM             dmc;
    Vec            coordsGlobal,coordsLocal;
    CoordinateData **arrCoords;
    ElementData    **arrParam; 
    PetscInt       start[2],n[2],i,j;

    // We do this with the help of a local-global mapping (even in serial), which even with no overlap should make things easier
    // Note - in this implementation we have conflated two types of ghost nodes. It would be conceivable to have up to 4 vector types: with/without interior ghosts, with/without exterior ghosts (used just to get an even number of blocks of the same size)
    ierr = DMGetCoordinateDM(paramGrid,&dmc);CHKERRQ(ierr);
    ierr = DMGetCoordinates(paramGrid,&coordsGlobal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(dmc,&coordsLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dmc,coordsGlobal,INSERT_VALUES,coordsLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dmc,coordsGlobal,INSERT_VALUES,coordsLocal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(paramGrid,&paramLocal);CHKERRQ(ierr);
    ierr = DMStagGetGhostCorners(paramGrid,&start[0],&start[1],NULL,&n[0],&n[1],NULL);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(dmc,coordsLocal,&arrCoords);CHKERRQ(ierr); // should be a "Read" version 

    // Iterate over all elements in the ghosted region. This is redundant, but means we already have the local information for later (thus doing more computation but less communication).This is also encouraged by our conflation of the two types of ghost points, as mentioned above.
    for (j=start[1];j<start[1]+n[1];++j) {
      for(i=start[0]; i<start[0]+n[0]; ++i) {
        arrParam[j][i].etaCorner  = getEta(ctx,arrCoords[j][i].xCorner, arrCoords[j][i].yCorner);
        arrParam[j][i].etaElement = getEta(ctx,arrCoords[j][i].xElement,arrCoords[j][i].yElement);
        arrParam[j][i].rho        = getRho(ctx,arrCoords[j][i].xElement,arrCoords[j][i].yElement);
      }
    }

    ierr = DMStagVecRestoreArray(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(dmc,coordsLocal,&arrCoords);CHKERRQ(ierr); // should be a "Read" version 
    ierr = VecDestroy(&coordsLocal);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(paramGrid,&ctx->param);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(paramGrid,paramLocal,INSERT_VALUES,ctx->param);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(paramGrid,paramLocal,INSERT_VALUES,ctx->param);CHKERRQ(ierr);
  }

  /* --- Create and Solve Linear System -------------------------------------- */
  ierr = CreateSystem(ctx,&A,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  if (size == 1) {
    ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERUMFPACK);CHKERRQ(ierr);
  } else {
    ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);CHKERRQ(ierr); 
  }
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /* --- Output -------------------------------------------------------------- */
  // Output the system as binary with headers, to see in Octave/MATLAB (see loadData.m)
  ierr = OutputMatBinary(A,"A.petscbin");
  ierr = OutputVecBinary(b,"b.petscbin",PETSC_FALSE);
  ierr = OutputVecBinary(x,"x.petscbin",PETSC_FALSE);

  /* Create auxiliary DMStag objects to help with output */
  // DMStag with a single dof per element/face
  ierr = DMStagCreate2d( PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,ctx->M,ctx->N,PETSC_DECIDE,PETSC_DECIDE,0,0,1,DMSTAG_GHOST_STENCIL_BOX, 1,&elementOnlyGrid);CHKERRQ(ierr);
  ierr = DMSetUp(elementOnlyGrid);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinates(elementOnlyGrid,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr); 

  // DMStag with a single dof per vertex
  ierr = DMStagCreate2d( PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, ctx->M,ctx->N, PETSC_DECIDE,PETSC_DECIDE, 1,0,0,DMSTAG_GHOST_STENCIL_BOX,1,&vertexOnlyGrid);CHKERRQ(ierr);
  ierr = DMSetUp(vertexOnlyGrid);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinates(vertexOnlyGrid,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr); 

  // Create single dof element- or vertex-only vectors to dump with our xdmf
  // Also do some simple diagonostics
  {
    typedef struct {PetscScalar etaCorner,etaElement,rho;} ParamData; // Design question: What's the best way to have users to do this? Provide these beforehand somehow?
    typedef struct {PetscScalar vy,vx,p;} StokesData; 
    Vec         etaElNatural,etaElLocal,rhoNatural,rhoLocal,etaCornerNatural,etaCornerLocal,bLocal,pLocal,pNatural,xLocal,vxInterpLocal,vxInterpNatural,vyInterpLocal,vyInterpNatural;
    PetscScalar **arrRho,**arrEtaEl,**arrEtaCorner,**arrp,**arrvxinterp,**arrvyinterp;
    ParamData   **arrParam; 
    StokesData  **arrb,**arrx;
    PetscInt    start[2],n[2],N[2],i,j;
    PetscMPIInt rank;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // TODO add checks that all these DMs are "compatible" to assure that simultaneous iteration is safe

    ierr = DMCreateLocalVector(stokesGrid,&bLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(stokesGrid,b,INSERT_VALUES,bLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(stokesGrid,b,INSERT_VALUES,bLocal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(stokesGrid,&xLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(stokesGrid,x,INSERT_VALUES,xLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(stokesGrid,x,INSERT_VALUES,xLocal);CHKERRQ(ierr);

    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&etaElNatural);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(elementOnlyGrid,&etaElLocal);CHKERRQ(ierr);

    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&rhoNatural);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(elementOnlyGrid,&rhoLocal);CHKERRQ(ierr);

    ierr = DMStagCreateNaturalVector(vertexOnlyGrid,&etaCornerNatural);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(vertexOnlyGrid,&etaCornerLocal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(      elementOnlyGrid,&pLocal          );CHKERRQ(ierr);
    ierr = DMCreateLocalVector(      elementOnlyGrid,&vxInterpLocal   );CHKERRQ(ierr);
    ierr = DMCreateLocalVector(      elementOnlyGrid,&vyInterpLocal   );CHKERRQ(ierr);
    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&pNatural        );CHKERRQ(ierr);
    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&vxInterpNatural);CHKERRQ(ierr);
    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&vyInterpNatural);CHKERRQ(ierr);

    ierr = DMStagVecGetArray(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr); // should be Read
    ierr = DMStagVecGetArray(elementOnlyGrid,etaElLocal,&arrEtaEl);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(elementOnlyGrid,rhoLocal,&arrRho);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(vertexOnlyGrid,etaCornerLocal,&arrEtaCorner);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(stokesGrid,bLocal,&arrb);CHKERRQ(ierr); // should be Read
    ierr = DMStagVecGetArray(elementOnlyGrid,pLocal,&arrp);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(elementOnlyGrid,vxInterpLocal,&arrvxinterp);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(elementOnlyGrid,vyInterpLocal,&arrvyinterp);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(stokesGrid,xLocal,&arrx);CHKERRQ(ierr); // should be Read

    ierr = DMStagGetCorners(paramGrid,&start[0],&start[1],NULL,&n[0],&n[1],NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMStagGetGlobalSizes(paramGrid,&N[0],&N[1],NULL);CHKERRQ(ierr);
    for (j=start[1]; j<start[1]+n[1]; ++j) {
      for(i=start[0]; i<start[0]+n[0]; ++i) {
        arrEtaEl[j][i]     = arrParam[j][i].etaElement;
        arrRho[j][i]       = arrParam[j][i].rho;
        arrEtaCorner[j][i] = arrParam[j][i].etaCorner;
        arrvyinterp[j][i]  = 0.5*(arrx[j][i].vy + arrx[j+1][i].vy);
        arrvxinterp[j][i]  = 0.5*(arrx[j][i].vx + arrx[j][i+1].vx);
        arrp[j][i]         = arrx[j][i].p;
      }
    }

    ierr = DMStagVecRestoreArray(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr); // should be Read
    ierr = DMStagVecRestoreArray(elementOnlyGrid,etaElLocal,&arrEtaEl);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(elementOnlyGrid,rhoLocal,&arrRho);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(vertexOnlyGrid,etaCornerLocal,&arrEtaCorner);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(stokesGrid,bLocal,&arrb);CHKERRQ(ierr); // should be Read
    ierr = DMStagVecRestoreArray(elementOnlyGrid,pLocal,&arrp);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(elementOnlyGrid,vxInterpLocal,&arrvxinterp);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(elementOnlyGrid,vyInterpLocal,&arrvyinterp);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(stokesGrid,xLocal,&arrx);CHKERRQ(ierr); // should be Read

    ierr = DMStagLocalToNaturalBegin(elementOnlyGrid,etaElLocal,INSERT_VALUES,etaElNatural);CHKERRQ(ierr);
    ierr = DMStagLocalToNaturalEnd(elementOnlyGrid,etaElLocal,INSERT_VALUES,etaElNatural);CHKERRQ(ierr);
    ierr = OutputVecBinary(etaElNatural,"etaElNatural.bin",PETSC_TRUE);CHKERRQ(ierr);

    ierr = DMStagLocalToNaturalBegin(elementOnlyGrid,rhoLocal,INSERT_VALUES,rhoNatural);CHKERRQ(ierr);
    ierr = DMStagLocalToNaturalEnd(elementOnlyGrid,rhoLocal,INSERT_VALUES,rhoNatural);CHKERRQ(ierr);
    ierr = OutputVecBinary(rhoNatural,"rhoNatural.bin",PETSC_TRUE);CHKERRQ(ierr);

    ierr = DMStagLocalToNaturalBegin(vertexOnlyGrid,etaCornerLocal,INSERT_VALUES,etaCornerNatural);CHKERRQ(ierr);
    ierr = DMStagLocalToNaturalEnd(vertexOnlyGrid,etaCornerLocal,INSERT_VALUES,etaCornerNatural);CHKERRQ(ierr);
    ierr = OutputVecBinary(etaCornerNatural,"etaCornerNatural.bin",PETSC_TRUE);CHKERRQ(ierr);

    ierr = DMStagLocalToNaturalBegin(elementOnlyGrid,pLocal,INSERT_VALUES,pNatural);CHKERRQ(ierr);
    ierr = DMStagLocalToNaturalEnd(elementOnlyGrid,pLocal,INSERT_VALUES,pNatural);CHKERRQ(ierr);
    ierr = VecMin(pNatural,NULL,&pPrimeMin);CHKERRQ(ierr);
    ierr = VecMax(pNatural,NULL,&pPrimeMax);CHKERRQ(ierr);
    ierr = OutputVecBinary(pNatural,"pPrimeNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecScale(pNatural,ctx->Kcont);CHKERRQ(ierr);
    ierr = VecMin(pNatural,NULL,&pMin);CHKERRQ(ierr);
    ierr = VecMax(pNatural,NULL,&pMax);CHKERRQ(ierr);
    ierr = OutputVecBinary(pNatural,"pNatural.bin",PETSC_TRUE);CHKERRQ(ierr);

    ierr = DMStagLocalToNaturalBegin(elementOnlyGrid,vxInterpLocal,INSERT_VALUES,vxInterpNatural);CHKERRQ(ierr);
    ierr = DMStagLocalToNaturalEnd(elementOnlyGrid,vxInterpLocal,INSERT_VALUES,vxInterpNatural);CHKERRQ(ierr);
    ierr = VecMin(vxInterpNatural,NULL,&vxInterpMin);CHKERRQ(ierr);
    ierr = VecMax(vxInterpNatural,NULL,&vxInterpMax);CHKERRQ(ierr);
    ierr = OutputVecBinary(vxInterpNatural,"vxInterpNatural.bin",PETSC_TRUE);CHKERRQ(ierr);

    ierr = DMStagLocalToNaturalBegin(elementOnlyGrid,vyInterpLocal,INSERT_VALUES,vyInterpNatural);CHKERRQ(ierr);
    ierr = DMStagLocalToNaturalEnd(elementOnlyGrid,vyInterpLocal,INSERT_VALUES,vyInterpNatural);CHKERRQ(ierr);
    ierr = VecMin(vyInterpNatural,NULL,&vyInterpMin);CHKERRQ(ierr);
    ierr = VecMax(vyInterpNatural,NULL,&vyInterpMax);CHKERRQ(ierr);
    ierr = OutputVecBinary(vyInterpNatural,"vyInterpNatural.bin",PETSC_TRUE);CHKERRQ(ierr);

    ierr = VecDestroy(&bLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&etaElNatural);CHKERRQ(ierr);
    ierr = VecDestroy(&etaElLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&etaCornerNatural);CHKERRQ(ierr);
    ierr = VecDestroy(&etaCornerLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&rhoNatural);CHKERRQ(ierr);
    ierr = VecDestroy(&rhoLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&pLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&vxInterpLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&vyInterpLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&pNatural);CHKERRQ(ierr);
    ierr = VecDestroy(&vxInterpNatural);CHKERRQ(ierr);
    ierr = VecDestroy(&vyInterpNatural);CHKERRQ(ierr);
    ierr = VecDestroy(&xLocal);CHKERRQ(ierr);
  }

  // Elements-only Xdmf
  {
    PetscViewer viewer;
    ierr = OutputDMCoordsNaturalBinary(elementOnlyGrid,"coordsElNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = DMStag2dXDMFStart(elementOnlyGrid,"elements_ParaView_Me.xmf","coordsElNatural.bin",&viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"etaElNatural.bin","eta",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"rhoNatural.bin","rho",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"pNatural.bin","p",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"pPrimeNatural.bin","pPrime",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"vxInterpNatural.bin","vxInterp",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"vyInterpNatural.bin","vyInterp",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFFinish(&viewer);CHKERRQ(ierr);
  }

  // Vertices-only Xdmf
  {
    PetscViewer viewer;
    ierr = OutputDMCoordsNaturalBinary(vertexOnlyGrid,"coordsVerNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = DMStag2dXDMFStart(vertexOnlyGrid,"vertices.xmf","coordsVerNatural.bin",&viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(vertexOnlyGrid,"etaCornerNatural.bin","eta",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFFinish(&viewer);CHKERRQ(ierr);
  }

  /* --- Simple Diagnostics -------------------------------------------------- */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"-- Info -----\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Kbound = %g\n",ctx->Kbound);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Kcont = %g\n",ctx->Kcont);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"hx = %g\n",ctx->hxCharacteristic);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"hy = %g\n",ctx->hyCharacteristic);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"-- Min / Max --\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"p' %g / %g \n",pPrimeMin,pPrimeMax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"p  %g / %g\n",pMin,pMax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"vxInterp %g / %g\n",vxInterpMin,vxInterpMax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"vyInterp %g / %g\n",vyInterpMin,vyInterpMax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"-------------\n");CHKERRQ(ierr);

  /* --- Clean up and Finalize ----------------------------------------------- */
  ierr = DMDestroy(&elementOnlyGrid);CHKERRQ(ierr);
  ierr = DMDestroy(&vertexOnlyGrid);CHKERRQ(ierr);
  ierr = VecDestroy(&paramLocal);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = DestroyCtx(&ctx);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
