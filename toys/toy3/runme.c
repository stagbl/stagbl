static char help[] = " Extension of toy2, which relies on a branch of PETSc including DMStag,\n\
                      and adds very simple particle advection using DMSwarm\n\
                      -nx,-ny : number of cells in x and y direction\n\
                      -xmax,-ymax : domain sizes (m) \n\
                      -p : particles per cell per dimension\n\
                      -dt : timestep for simple advection\n\
                      -isoviscous : don't use variable viscosity\n";

#include <petsc.h>
#include "ctx.h"
#include "system.h"
#include "output.h"
#include "system.h"

// TODO: introduce interface to fix this (or decide that the approach here
//       is in general hacky and move more into the DMStag impl)
#include <petsc/private/dmstagimpl.h>

PetscErrorCode MaterialPoint_AdvectRK1(DM dm_vp,Vec vp,PetscReal dt,DM dm_mpoint);

/* A deprecated function, removed from DMStag but still used here */
PetscErrorCode DMStagGetVertexDMDA(DM,PetscInt,DM*);

/* =========================================================================== */
int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  Ctx            ctx;
  DM             stokesGrid,paramGrid,elementOnlyGrid,vertexOnlyGrid,swarm,daVertex,daVertex2;
  Vec            paramLocal;
  Mat            A;
  Vec            b;
  Vec            x;
  Vec            vVertex=NULL,vVertexLocal=NULL;
  KSP            ksp;
  PC             pc;
  PetscMPIInt    size;
  PetscScalar    pMin,pMax,pPrimeMin,pPrimeMax,vxInterpMin,vxInterpMax,vyInterpMin,vyInterpMax;
  PetscInt       nel,npe;
  const PetscInt *element_list;

  /* --- Initialize and Create Context --------------------------------------- */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = CreateCtx(&ctx);CHKERRQ(ierr);

  /* --- Create Main DMStag Objects ------------------------------------------ */
  // A DMStag to hold the unknowns to solve the Stokes problem
  ierr = DMStagCreate2d(PETSC_COMM_WORLD,
      DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,       // no special boundary support yet
      ctx->M,ctx->N,PETSC_DECIDE,PETSC_DECIDE, // sizes provided as DMDA (global x, global y, local x, local y)
      0,1,1,                                   // no dof per vertex: 0 per vertex, 1 per edge, 1 per face/element
      DMSTAG_GHOST_STENCIL_BOX,                // elementwise stencil pattern
      1,                                       // elementwise stencil width
      NULL,NULL,
      &ctx->stokesGrid);
  stokesGrid = ctx->stokesGrid;
  ierr = DMSetUp(stokesGrid);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesExplicit(stokesGrid,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr);  // TODO : here (and everywhere) use DMStagSetUniformCoordinates (once implemented)

  // A DMStag to hold the coefficient fields for the Stokes problem: 1 dof per vertex and two per element/face
  ierr = DMStagCreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,ctx->M,ctx->N,PETSC_DECIDE,PETSC_DECIDE,1,0,2,DMSTAG_GHOST_STENCIL_BOX,1,NULL,NULL,&ctx->paramGrid);CHKERRQ(ierr);
  paramGrid = ctx->paramGrid;
  ierr = DMSetUp(paramGrid);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesExplicit(paramGrid,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"paramGrid:\n");
  ierr = DMView(paramGrid,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /* --- Create a DMDAs representing the vertices of our grid (since the logic is already there for DMSwarm) --- */
    ierr = DMStagGetVertexDMDA(paramGrid,1,&daVertex);
    ierr = DMSetUp(daVertex);CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates(daVertex,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr); // TODO put this into the function to create the DMDA?
    ierr = DMStagGetVertexDMDA(paramGrid,2,&daVertex2);
    ierr = DMSetUp(daVertex2);CHKERRQ(ierr);

    ierr = DMDASetUniformCoordinates(daVertex2,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr);
    ierr = DMView(daVertex2,PETSC_VIEWER_STDOUT_WORLD);

    ierr = DMDAGetElements(daVertex2,&nel,&npe,&element_list);CHKERRQ(ierr);
    ierr = DMDARestoreElements(daVertex2,&nel,&npe,&element_list);CHKERRQ(ierr);

  /* --- Create a DMSwarm (particle system) object --------------------------- */
  {
    PetscInt particlesPerElementPerDim = 10;

    ierr = PetscOptionsGetInt(NULL,NULL,"-p",&particlesPerElementPerDim,NULL);CHKERRQ(ierr);
    ierr = DMCreate(PETSC_COMM_WORLD,&swarm);CHKERRQ(ierr);
    ierr = DMSetType(swarm,DMSWARM);CHKERRQ(ierr);
    ierr = DMSetDimension(swarm,2);CHKERRQ(ierr);
    ierr = DMSwarmSetType(swarm,DMSWARM_PIC);CHKERRQ(ierr);
    ierr = DMSwarmSetCellDM(swarm,daVertex);CHKERRQ(ierr);
    ierr = DMSwarmRegisterPetscDatatypeField(swarm,"rho",1,PETSC_REAL);CHKERRQ(ierr);
    ierr = DMSwarmRegisterPetscDatatypeField(swarm,"eta",1,PETSC_REAL);CHKERRQ(ierr);
    ierr = DMSwarmFinalizeFieldRegister(swarm);CHKERRQ(ierr);
    ierr = DMSwarmSetLocalSizes(swarm,nel*particlesPerElementPerDim,100);CHKERRQ(ierr);
    ierr = DMSwarmInsertPointsUsingCellDM(swarm,DMSWARMPIC_LAYOUT_REGULAR,particlesPerElementPerDim);CHKERRQ(ierr);

    // Set properties for particles.  Each particle has a viscosity and density (would probably be more efficient
    //  in this example to use material ids)
    {
      PetscReal   *array_x,*array_e,*array_r;
      PetscInt    npoints,p;

      ierr = DMSwarmGetField(swarm,DMSwarmPICField_coor,NULL,NULL,(void**)&array_x);CHKERRQ(ierr);
      ierr = DMSwarmGetField(swarm,"eta",               NULL,NULL,(void**)&array_e);CHKERRQ(ierr);
      ierr = DMSwarmGetField(swarm,"rho",               NULL,NULL,(void**)&array_r);CHKERRQ(ierr);
      ierr = DMSwarmGetLocalSize(swarm,&npoints);CHKERRQ(ierr);
      for (p = 0; p < npoints; p++) {
        PetscReal x_p[2];

        /* One could apply random noise to the particles here (and should!)
         See KSP tutorial ex70 for an example; note that one should call
         DMSwarmMigrate() since particles may move across cell boundaries */

        // Get the coordinates of point p
        x_p[0] = array_x[2*p + 0];
        x_p[1] = array_x[2*p + 1];

        // Call functions to compute eta and rho at that location
        array_e[p] = getEta(ctx,x_p[0],x_p[1]);
        array_r[p] = getRho(ctx,x_p[0],x_p[1]);

      }
      ierr = DMSwarmRestoreField(swarm,"rho",NULL,NULL,(void**)&array_r);CHKERRQ(ierr);
      ierr = DMSwarmRestoreField(swarm,"eta",NULL,NULL,(void**)&array_e);CHKERRQ(ierr);
      ierr = DMSwarmRestoreField(swarm,DMSwarmPICField_coor,NULL,NULL,(void**)&array_x);CHKERRQ(ierr);
    }

    // View DMSwarm object
    ierr = DMView(swarm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  /* --- Set up Problem ------------------------------------------------------ */

  // Project from particles, using the existing DMDA implementation for DMSwarm and then making an ugly transfer.
  // TODO do this each time step (but reuse these vectors and move decls up, and make sure that the vecs and the array of vecs are freed, as both are created here)
  {
    Vec               *pfields;
    Vec               etaVertex,rhoVertex;
    Vec               etaVertexLocal,rhoVertexLocal;
    const char        *fieldnames[]={"eta","rho"};
    const PetscBool   reuse = PETSC_FALSE;
    PetscScalar       *arrParamRaw;
    const PetscScalar *vArrEta,*vArrRho;

    PetscInt          nel,npe,eidx,nelxDA,xs,ys,xe,ye,Xs,Ys,Xe,Ye,entriesPerElementRowGhost;
    const PetscInt    *element_list;
    DM_Stag           *stag; // TODO: obviously this isn't good. We shouldn't have to access this data once we have a proper API

    PetscMPIInt       rank; // debug

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    stag = (DM_Stag*)stokesGrid->data;

    /* Note here that it would be better to implement DMStag as a base
       DM for DMSwarm, as we could then project nicely to both the corners
       and the cell centers */

    ierr = DMSwarmProjectFields(swarm,2,fieldnames,&pfields,reuse);CHKERRQ(ierr);
    etaVertex = pfields[0];
    rhoVertex = pfields[1];

    // Global->Local
    ierr = DMGetLocalVector(daVertex,&etaVertexLocal);CHKERRQ(ierr);
    ierr = DMGetLocalVector(daVertex,&rhoVertexLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(daVertex,etaVertex,INSERT_VALUES,etaVertexLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(daVertex,etaVertex,INSERT_VALUES,etaVertexLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(daVertex,rhoVertex,INSERT_VALUES,rhoVertexLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(daVertex,rhoVertex,INSERT_VALUES,rhoVertexLocal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(paramGrid,&paramLocal);CHKERRQ(ierr);
    ierr = VecGetArray(paramLocal,&arrParamRaw);CHKERRQ(ierr);

    ierr = VecGetArrayRead(etaVertex,&vArrEta);CHKERRQ(ierr);
    ierr = VecGetArrayRead(rhoVertex,&vArrRho);CHKERRQ(ierr);

    // Logic copied from DMDAGetElements_2D()
    // TODO : check element type is correct
    ierr   = DMDAGetCorners(daVertex,&xs,&ys,0,&xe,&ye,0);CHKERRQ(ierr);
    ierr   = DMDAGetGhostCorners(daVertex,&Xs,&Ys,0,&Xe,&Ye,0);CHKERRQ(ierr);
    xe    += xs; Xe += Xs; if (xs != Xs) xs -= 1;
    ye    += ys; Ye += Ys; if (ys != Ys) ys -= 1;
    nelxDA = xe-xs-1;
    ierr = DMDAGetElements(daVertex,&nel,&npe,&element_list);CHKERRQ(ierr);

    entriesPerElementRowGhost =  stag->nGhost[0] * stag->entriesPerElement;

    for (eidx = 0; eidx < nel; eidx++) {
      PetscInt localRowDA,localColDA;

      // This lists the local element numbers
      const PetscInt *element = &element_list[npe*eidx];
      PetscInt indTo,indFrom;

      localRowDA = eidx / nelxDA; // Integer div
      localColDA = eidx % nelxDA;

      // Lower left
      indFrom = element[0];
      indTo = (localRowDA  ) * entriesPerElementRowGhost + (localColDA  ) * stag->entriesPerElement+ 0; // etaCorner
      arrParamRaw[indTo] = vArrEta[indFrom];
      indTo = (localRowDA  ) * entriesPerElementRowGhost + (localColDA  ) * stag->entriesPerElement+ 1; // etaElement
      arrParamRaw[indTo] = vArrEta[indFrom];
      indTo = (localRowDA  ) * entriesPerElementRowGhost + (localColDA  ) * stag->entriesPerElement+ 2; // rho (element)
      arrParamRaw[indTo] = vArrRho[indFrom];

      // Lower Right
      indFrom = element[1];
      indTo = (localRowDA  ) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 0; // etaCorner
      arrParamRaw[indTo] = vArrEta[indFrom];
      indTo = (localRowDA  ) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 1; // etaElement
      arrParamRaw[indTo] = vArrEta[indFrom];
      indTo = (localRowDA  ) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 2; // rho (element)
      arrParamRaw[indTo] = vArrRho[indFrom];

      // Upper Right
      indFrom = element[2];
      indTo = (localRowDA+1) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 0; // etaCorner
      arrParamRaw[indTo] = vArrEta[indFrom];
      indTo = (localRowDA+1) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 1; // etaElement
      arrParamRaw[indTo] = vArrEta[indFrom];
      indTo = (localRowDA+1) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 2; // rho (element)
      arrParamRaw[indTo] = vArrRho[indFrom];

      // Upper Left
      indFrom = element[3];
      indTo = (localRowDA+1) * entriesPerElementRowGhost + (localColDA ) * stag->entriesPerElement+ 0; // etaCorner
      arrParamRaw[indTo] = vArrEta[indFrom];
      indTo = (localRowDA+1) * entriesPerElementRowGhost + (localColDA ) * stag->entriesPerElement+ 1; // etaElement
      arrParamRaw[indTo] = vArrEta[indFrom];
      indTo = (localRowDA+1) * entriesPerElementRowGhost + (localColDA ) * stag->entriesPerElement+ 2; // rho (element)
      arrParamRaw[indTo] = vArrRho[indFrom];
    }

    ierr = VecRestoreArray(paramLocal,&arrParamRaw);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(etaVertex,&vArrEta);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(rhoVertex,&vArrRho);CHKERRQ(ierr);

    // Local->global
    ierr = DMCreateGlobalVector(paramGrid,&ctx->param);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(paramGrid,paramLocal,INSERT_VALUES,ctx->param);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(paramGrid,paramLocal,INSERT_VALUES,ctx->param);CHKERRQ(ierr);

    // Inefficiently, fill in any ghost/dummy points that we missed by scattering back again
    ierr = DMGlobalToLocalBegin(paramGrid,ctx->param,INSERT_VALUES,paramLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(paramGrid,ctx->param,INSERT_VALUES,paramLocal);CHKERRQ(ierr);

    ierr = VecDestroy(&etaVertexLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&rhoVertexLocal);CHKERRQ(ierr);
  }

#if 0
  // This is how one might directly populate a DMStag array (keep this code around, as it'll be useful when implementing the proper usage of DMStag as a base DM)
  {
    typedef struct {PetscReal xCorner,yCorner,xElement,yElement;} CoordinateData; // Definitely need to provide these types somehow to the user.. (specifying them wrong leads to very annoying bugs)
    typedef struct {PetscScalar etaCorner,etaElement,rho;} ElementData; // ? What's the best way to have users to do this? Provide these before hand?
    DM             dmc;
    Vec            coordsGlobal,coordsLocal;
    CoordinateData **arrCoords;
    ElementData    **arrParam;
    PetscInt       start[2],n[2],i,j;

    ierr = DMGetCoordinateDM(paramGrid,&dmc);CHKERRQ(ierr);
    ierr = DMGetCoordinates(paramGrid,&coordsGlobal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(dmc,&coordsLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dmc,coordsGlobal,INSERT_VALUES,coordsLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dmc,coordsGlobal,INSERT_VALUES,coordsLocal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(paramGrid,&paramLocal);CHKERRQ(ierr);
    ierr = DMStagGetGhostCorners(paramGrid,&start[0],&start[1],NULL,&n[0],&n[1],NULL);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr);
    ierr = DMStagVecGetArrayRead(dmc,coordsLocal,&arrCoords);CHKERRQ(ierr);

    // Iterate over all elements in the ghosted region. This is redundant, but means we already have the local information for later (thus doing more computation but less communication).This is also encouraged by our conflation of the two types of ghost points, as mentioned above.
    // TODO change this with getting stuff from etaVertex and rhoVertex (and interpolating)
    for (j=start[1];j<start[1]+n[1];++j) {
      for(i=start[0]; i<start[0]+n[0]; ++i) {
        arrParam[j][i].etaCorner  = getEta(ctx,arrCoords[j][i].xCorner, arrCoords[j][i].yCorner);
        arrParam[j][i].etaElement = getEta(ctx,arrCoords[j][i].xElement,arrCoords[j][i].yElement);
        arrParam[j][i].rho        = getRho(ctx,arrCoords[j][i].xElement,arrCoords[j][i].yElement);
      }
    }

    ierr = DMStagVecRestoreArray(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArrayRead(dmc,coordsLocal,&arrCoords);CHKERRQ(ierr);
    ierr = VecDestroy(&coordsLocal);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(paramGrid,&ctx->param);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(paramGrid,paramLocal,INSERT_VALUES,ctx->param);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(paramGrid,paramLocal,INSERT_VALUES,ctx->param);CHKERRQ(ierr);

  }
#endif

  /* --- Create and Solve Linear System -------------------------------------- */

  // TODO put the solve, and recomputation of rho/eta from tracers, inside the time loop
  ierr = CreateSystem(ctx,&A,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  if (size == 1) {
    ierr = PCFactorSetMatSolverType(pc,MATSOLVERUMFPACK);CHKERRQ(ierr);
  } else {
    ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);CHKERRQ(ierr);
  }
  ierr = PetscOptionsSetValue(NULL,"-ksp_converged_reason","");CHKERRQ(ierr); // To get info on direct solve success
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  {
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
    if (reason < 0) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,"Linear solve failed");CHKERRQ(ierr);
  }

  /* Transfer velocities to arrays living on the vertex-only DMDA.
  TODO: this is hacky and it should be possible to do this through the API */
  // TODO : introduce compatibility check for this (or perhaps make it so that we can do this with a vertex-only DMStag)
  ierr = DMCreateGlobalVector(daVertex2,&vVertex);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(daVertex2,&vVertexLocal);CHKERRQ(ierr);
  {
    PetscScalar       *vArr;
    Vec               xLocal;
    const PetscScalar *arrxLocalRaw;
    PetscInt          nel,npe,eidx,nelxDA,xs,ys,xe,ye,Xs,Ys,Xe,Ye,entriesPerElementRowGhost;
    const PetscInt    *element_list;
    DM_Stag           *stag; // TODO: obviously this isn't good. We shouldn't have to access this data once we have a proper API

    PetscMPIInt       rank; // debug

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    stag = (DM_Stag*)stokesGrid->data;

    ierr = DMCreateLocalVector(stokesGrid,&xLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(stokesGrid,x,INSERT_VALUES,xLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(stokesGrid,x,INSERT_VALUES,xLocal);CHKERRQ(ierr);
    ierr = VecGetArrayRead(xLocal,&arrxLocalRaw);CHKERRQ(ierr);

    // Logic copied from DMDAGetElements_2D()
    // TODO : check element type is correct
    ierr   = DMDAGetCorners(daVertex2,&xs,&ys,0,&xe,&ye,0);CHKERRQ(ierr);
    ierr   = DMDAGetGhostCorners(daVertex2,&Xs,&Ys,0,&Xe,&Ye,0);CHKERRQ(ierr);
    xe    += xs; Xe += Xs; if (xs != Xs) xs -= 1;
    ye    += ys; Ye += Ys; if (ys != Ys) ys -= 1;
    nelxDA = xe-xs-1;
    ierr = DMDAGetElements(daVertex2,&nel,&npe,&element_list);CHKERRQ(ierr);

    ierr = VecGetArray(vVertexLocal,&vArr);CHKERRQ(ierr);

    entriesPerElementRowGhost =  stag->nGhost[0] * stag->entriesPerElement;

    for (eidx = 0; eidx < nel; eidx++) {
      PetscInt localRowDA,localColDA;

      // This lists the local element numbers
      const PetscInt *element = &element_list[npe*eidx];
      const PetscInt NSD = 2;
      PetscInt indTo,indFrom;

      localRowDA = eidx / nelxDA; // Integer div
      localColDA = eidx % nelxDA;
      //PetscPrintf(PETSC_COMM_SELF,"[%d] DMDA el %d (%d,%d) : %d %d %d %d\n",rank,eidx,localRowDA,localColDA,element[0],element[1],element[2],element[3]);

      // TODO wrap this in a function in the API (thus obviating the need for access to the DM_Stag object)
      indTo = element[0]*NSD+0;
      indFrom = (localRowDA  ) * entriesPerElementRowGhost + (localColDA  ) * stag->entriesPerElement+ 1; // vx is the first entry
      vArr[indTo] = arrxLocalRaw[indFrom];
      //ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d --> %d (%g) \n",rank,indFrom,indTo,vArr[indTo]);CHKERRQ(ierr);

      indTo = element[0]*NSD+1;
      indFrom = (localRowDA  ) * entriesPerElementRowGhost + (localColDA  ) * stag->entriesPerElement+ 0;
      vArr[indTo] = arrxLocalRaw[indFrom];
      //ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d --> %d (%g) \n",rank,indFrom,indTo,vArr[indTo]);CHKERRQ(ierr);

      indTo = element[1]*NSD+0;
      indFrom = (localRowDA  ) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 1;
      vArr[indTo] = arrxLocalRaw[indFrom];
      //ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d --> %d (%g) \n",rank,indFrom,indTo,vArr[indTo]);CHKERRQ(ierr);

      indTo = element[1]*NSD+1;
      indFrom = (localRowDA  ) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 0;
      vArr[indTo] = arrxLocalRaw[indFrom];
      //ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d --> %d (%g) \n",rank,indFrom,indTo,vArr[indTo]);CHKERRQ(ierr);

      // This is top-right in DMDA ordering
      indTo = element[2]*NSD+0;
      indFrom = (localRowDA+1) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 1;
      vArr[indTo] = arrxLocalRaw[indFrom];
      //ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d --> %d (%g) \n",rank,indFrom,indTo,vArr[indTo]);CHKERRQ(ierr);

      indTo = element[2]*NSD+1;
      indFrom = (localRowDA+1) * entriesPerElementRowGhost + (localColDA+1) * stag->entriesPerElement+ 0;
      vArr[indTo] = arrxLocalRaw[indFrom];
      //ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d --> %d (%g) \n",rank,indFrom,indTo,vArr[indTo]);CHKERRQ(ierr);

      // This is top-left in DMDA ordering
      indTo = element[3]*NSD+0;
      indFrom = (localRowDA+1) * entriesPerElementRowGhost + (localColDA  ) * stag->entriesPerElement+ 1;
      vArr[indTo] = arrxLocalRaw[indFrom];
      //ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d --> %d (%g) \n",rank,indFrom,indTo,vArr[indTo]);CHKERRQ(ierr);

      indTo = element[3]*NSD+1;
      indFrom = (localRowDA+1) * entriesPerElementRowGhost + (localColDA  ) * stag->entriesPerElement+ 0;
      vArr[indTo] = arrxLocalRaw[indFrom];
      //ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d --> %d (%g) \n",rank,indFrom,indTo,vArr[indTo]);CHKERRQ(ierr);
    }
    // TODO: improve the above to use averaged values, not just grabbing from one neighboring edge

    ierr = VecRestoreArray(vVertexLocal,&vArr);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(xLocal,&arrxLocalRaw);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(daVertex2,vVertexLocal,INSERT_VALUES,vVertex);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(daVertex2,vVertexLocal,INSERT_VALUES,vVertex);CHKERRQ(ierr);
    ierr = VecDestroy(&xLocal);CHKERRQ(ierr);
  }

  /* --- Push Particles and Dump Xdmf ---------------------------------------- */
  // A naive forward Euler step
  // Note: open these in Paraview 5.3.0 on OS X with "Xdmf Reader", not "Xdmf 3.0 reader"
  {
    PetscInt  step;
    char      filename[PETSC_MAX_PATH_LEN];
    PetscInt  nSteps;
    PetscReal dt;

    step = 0;
    nSteps = 10;
    dt = 1e11;
      // TODO pick a timestep based on max grid velocity (see ex70)
    ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"-- Advection -----\n",step);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-nsteps",&nSteps,NULL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Step %D of %D\n",step,nSteps);CHKERRQ(ierr);
    ierr = DMSwarmViewXDMF(swarm,"swarm_0000.xmf");CHKERRQ(ierr);
    for (step=1; step<=nSteps; ++step) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Step %D of %D\n",step,nSteps);CHKERRQ(ierr);

      // Advect and Migrate
      ierr = MaterialPoint_AdvectRK1(daVertex2,vVertex,dt,swarm);CHKERRQ(ierr);
      ierr = DMSwarmMigrate(swarm,PETSC_TRUE);CHKERRQ(ierr);

      // TODO allow for points outside of the domain as in 8.4 p. 121?
      // TODO update grid quanitities from tracers (once other stuff works)

      // Output
      ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"swarm_%.4D.xmf",step);CHKERRQ(ierr);
      ierr = DMSwarmViewXDMF(swarm,filename);CHKERRQ(ierr);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\rDone (%d steps).\n",nSteps);CHKERRQ(ierr);
  }

  /* --- Additional Output --------------------------------------------------- */
  // Output the system as binary with headers, to see in Octave/MATLAB (see loadData.m)
  ierr = OutputMatBinary(A,"A.pbin");
  ierr = OutputVecBinary(b,"b.pbin",PETSC_FALSE);
  ierr = OutputVecBinary(x,"x.pbin",PETSC_FALSE);

  // Create auxiliary DMStag object, with one dof per element, to help with output
  ierr = DMStagCreate2d( PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,ctx->M,ctx->N,PETSC_DECIDE,PETSC_DECIDE,0,0,1,DMSTAG_GHOST_STENCIL_BOX,1,NULL,NULL,&elementOnlyGrid);CHKERRQ(ierr);
  ierr = DMSetUp(elementOnlyGrid);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesExplicit(elementOnlyGrid,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr);

  // Another auxiliary DMStag to help dumping vertex/corner - only data
  ierr = DMStagCreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,ctx->M,ctx->N,PETSC_DECIDE,PETSC_DECIDE,1,0,0,DMSTAG_GHOST_STENCIL_BOX,1,NULL,NULL,&vertexOnlyGrid);CHKERRQ(ierr);
  ierr = DMSetUp(vertexOnlyGrid);CHKERRQ(ierr);
  ierr = DMStagSetUniformCoordinatesExplicit(vertexOnlyGrid,ctx->xmin,ctx->xmax,ctx->ymin,ctx->ymax,0.0,0.0);CHKERRQ(ierr);

  // Create single dof element- or vertex-only vectors to dump with our Xdmf
  {
    typedef struct {PetscScalar etaCorner,etaElement,rho;} ParamData; // Design question: What's the best way to have users to do this? Provide these beforehand somehow?
    typedef struct {PetscScalar vy,vx,p;} StokesData;
    Vec         xLocal;
    Vec         etaElGlobal,etaElNatural,etaElLocal;
    Vec         rhoGlobal,rhoNatural,rhoLocal;
    Vec         pGlobal,pLocal,pNatural;
    Vec         vxInterpGlobal,vxInterpLocal,vxInterpNatural,vyInterpGlobal,vyInterpLocal,vyInterpNatural;
    PetscScalar **arrRho,**arrEtaEl,**arrp,**arrvxinterp,**arrvyinterp;
    ParamData   **arrParam;
    StokesData  **arrx;
    PetscInt    start[2],n[2],N[2],i,j;
    PetscMPIInt rank;

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // TODO add checks that all these DMs are "compatible" to assure that simultaneous iteration is safe

    ierr = DMCreateLocalVector(stokesGrid,&xLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(stokesGrid,x,INSERT_VALUES,xLocal);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(stokesGrid,x,INSERT_VALUES,xLocal);CHKERRQ(ierr);

    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&etaElNatural);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(elementOnlyGrid,&etaElLocal);CHKERRQ(ierr);

    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&rhoNatural);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(elementOnlyGrid,&rhoLocal);CHKERRQ(ierr);

    ierr = DMCreateLocalVector(      elementOnlyGrid,&pLocal          );CHKERRQ(ierr);
    ierr = DMCreateLocalVector(      elementOnlyGrid,&vxInterpLocal   );CHKERRQ(ierr);
    ierr = DMCreateLocalVector(      elementOnlyGrid,&vyInterpLocal   );CHKERRQ(ierr);
    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&pNatural        );CHKERRQ(ierr);
    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&vxInterpNatural);CHKERRQ(ierr);
    ierr = DMStagCreateNaturalVector(elementOnlyGrid,&vyInterpNatural);CHKERRQ(ierr);

    ierr = DMStagVecGetArrayRead(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(elementOnlyGrid,etaElLocal,&arrEtaEl);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(elementOnlyGrid,rhoLocal,&arrRho);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(elementOnlyGrid,pLocal,&arrp);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(elementOnlyGrid,vxInterpLocal,&arrvxinterp);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(elementOnlyGrid,vyInterpLocal,&arrvyinterp);CHKERRQ(ierr);
    ierr = DMStagVecGetArrayRead(stokesGrid,xLocal,&arrx);CHKERRQ(ierr);

    ierr = DMStagGetCorners(paramGrid,&start[0],&start[1],NULL,&n[0],&n[1],NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMStagGetGlobalSizes(paramGrid,&N[0],&N[1],NULL);CHKERRQ(ierr);
    for (j=start[1]; j<start[1]+n[1]; ++j) {
      for(i=start[0]; i<start[0]+n[0]; ++i) {
        arrEtaEl[j][i]     = arrParam[j][i].etaElement;
        arrRho[j][i]       = arrParam[j][i].rho;
        arrvyinterp[j][i]  = 0.5*(arrx[j][i].vy + arrx[j+1][i].vy);
        arrvxinterp[j][i]  = 0.5*(arrx[j][i].vx + arrx[j][i+1].vx);
        arrp[j][i]         = arrx[j][i].p;
      }
    }

    ierr = DMStagVecRestoreArrayRead(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(elementOnlyGrid,etaElLocal,&arrEtaEl);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(elementOnlyGrid,rhoLocal,&arrRho);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(elementOnlyGrid,pLocal,&arrp);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(elementOnlyGrid,vxInterpLocal,&arrvxinterp);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArray(elementOnlyGrid,vyInterpLocal,&arrvyinterp);CHKERRQ(ierr);
    ierr = DMStagVecRestoreArrayRead(stokesGrid,xLocal,&arrx);CHKERRQ(ierr);

    ierr = VecDestroy(&xLocal);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(elementOnlyGrid,&etaElGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(elementOnlyGrid,etaElLocal,INSERT_VALUES,etaElGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(elementOnlyGrid,etaElLocal,INSERT_VALUES,etaElGlobal);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalBegin(elementOnlyGrid,etaElGlobal,INSERT_VALUES,etaElNatural);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalEnd(elementOnlyGrid,etaElGlobal,INSERT_VALUES,etaElNatural);CHKERRQ(ierr);
    ierr = OutputVecBinary(etaElNatural,"etaElNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecDestroy(&etaElGlobal);CHKERRQ(ierr);
    ierr = VecDestroy(&etaElNatural);CHKERRQ(ierr);
    ierr = VecDestroy(&etaElLocal);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(elementOnlyGrid,&rhoGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(elementOnlyGrid,rhoLocal,INSERT_VALUES,rhoGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(elementOnlyGrid,rhoLocal,INSERT_VALUES,rhoGlobal);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalBegin(elementOnlyGrid,rhoGlobal,INSERT_VALUES,rhoNatural);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalEnd(elementOnlyGrid,rhoGlobal,INSERT_VALUES,rhoNatural);CHKERRQ(ierr);
    ierr = OutputVecBinary(rhoNatural,"rhoNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecDestroy(&rhoGlobal);CHKERRQ(ierr);
    ierr = VecDestroy(&rhoNatural);CHKERRQ(ierr);
    ierr = VecDestroy(&rhoLocal);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(elementOnlyGrid,&pGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(elementOnlyGrid,pLocal,INSERT_VALUES,pGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(elementOnlyGrid,pLocal,INSERT_VALUES,pGlobal);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalBegin(elementOnlyGrid,pGlobal,INSERT_VALUES,pNatural);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalEnd(elementOnlyGrid,pGlobal,INSERT_VALUES,pNatural);CHKERRQ(ierr);
    ierr = VecMin(pNatural,NULL,&pPrimeMin);CHKERRQ(ierr);
    ierr = VecMax(pNatural,NULL,&pPrimeMax);CHKERRQ(ierr);
    ierr = OutputVecBinary(pNatural,"pPrimeNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecScale(pNatural,ctx->Kcont);CHKERRQ(ierr);
    ierr = VecMin(pNatural,NULL,&pMin);CHKERRQ(ierr);
    ierr = VecMax(pNatural,NULL,&pMax);CHKERRQ(ierr);
    ierr = OutputVecBinary(pNatural,"pNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecDestroy(&pLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&pNatural);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(elementOnlyGrid,&vxInterpGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(elementOnlyGrid,vxInterpLocal,INSERT_VALUES,vxInterpGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(elementOnlyGrid,vxInterpLocal,INSERT_VALUES,vxInterpGlobal);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalBegin(elementOnlyGrid,vxInterpGlobal,INSERT_VALUES,vxInterpNatural);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalEnd(elementOnlyGrid,vxInterpGlobal,INSERT_VALUES,vxInterpNatural);CHKERRQ(ierr);
    ierr = VecMin(vxInterpNatural,NULL,&vxInterpMin);CHKERRQ(ierr);
    ierr = VecMax(vxInterpNatural,NULL,&vxInterpMax);CHKERRQ(ierr);
    ierr = OutputVecBinary(vxInterpNatural,"vxInterpNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecDestroy(&vxInterpGlobal);CHKERRQ(ierr);
    ierr = VecDestroy(&vxInterpLocal);CHKERRQ(ierr);
    ierr = VecDestroy(&vxInterpNatural);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(elementOnlyGrid,&vyInterpGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(elementOnlyGrid,vyInterpLocal,INSERT_VALUES,vyInterpGlobal);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(elementOnlyGrid,vyInterpLocal,INSERT_VALUES,vyInterpGlobal);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalBegin(elementOnlyGrid,vyInterpGlobal,INSERT_VALUES,vyInterpNatural);CHKERRQ(ierr);
    ierr = DMStagGlobalToNaturalEnd(elementOnlyGrid,vyInterpGlobal,INSERT_VALUES,vyInterpNatural);CHKERRQ(ierr);
    ierr = VecMin(vyInterpNatural,NULL,&vyInterpMin);CHKERRQ(ierr);
    ierr = VecMax(vyInterpNatural,NULL,&vyInterpMax);CHKERRQ(ierr);
    ierr = OutputVecBinary(vyInterpNatural,"vyInterpNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecDestroy(&vyInterpGlobal);CHKERRQ(ierr);
    ierr = VecDestroy(&vyInterpNatural);CHKERRQ(ierr);
    ierr = VecDestroy(&vyInterpLocal);CHKERRQ(ierr);
  }

  // Elements-only Xdmf
  {
    PetscViewer viewer;
    ierr = OutputDMCoordsNaturalBinary(elementOnlyGrid,"coordsElNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = DMStag2dXDMFStart(elementOnlyGrid,"ParaViewMe.xmf","coordsElNatural.bin",&viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"etaElNatural.bin","eta",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"rhoNatural.bin","rho",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"pNatural.bin","p",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"pPrimeNatural.bin","pPrime",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"vxInterpNatural.bin","vxInterp",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFAddAttribute(elementOnlyGrid,"vyInterpNatural.bin","vyInterp",viewer);CHKERRQ(ierr);
    ierr = DMStag2dXDMFFinish(&viewer);CHKERRQ(ierr);
  }

  // Vertex-only Xdmf
  {
    PetscViewer viewer;
    ierr = OutputDMCoordsNaturalBinary(vertexOnlyGrid,"coordsVertexNatural.bin",PETSC_TRUE);CHKERRQ(ierr);
    ierr = DMStag2dXDMFStart(vertexOnlyGrid,"ParaViewMe_Too.xmf","coordsVertexNatural.bin",&viewer);CHKERRQ(ierr);
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
  ierr = DMDestroy(&daVertex);CHKERRQ(ierr);
  ierr = DMDestroy(&daVertex2);CHKERRQ(ierr);
  ierr = DMDestroy(&elementOnlyGrid);CHKERRQ(ierr);
  ierr = DMDestroy(&vertexOnlyGrid);CHKERRQ(ierr);
  ierr = VecDestroy(&paramLocal);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  if (vVertex)      {ierr = VecDestroy(&vVertex);CHKERRQ(ierr);}
  if (vVertexLocal) {ierr = VecDestroy(&vVertexLocal);CHKERRQ(ierr);}
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = DestroyCtx(&ctx);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/* Routine taken from KSP tutorial ex70 */

#define NSD            2 /* number of spatial dimensions */
#define NODES_PER_EL   4 /* nodes per element */

static void EvaluateBasis_Q1(PetscScalar _xi[],PetscScalar N[])
{
  PetscScalar xi  = _xi[0];
  PetscScalar eta = _xi[1];

  N[0] = 0.25*(1.0-xi)*(1.0-eta);
  N[1] = 0.25*(1.0+xi)*(1.0-eta);
  N[2] = 0.25*(1.0+xi)*(1.0+eta);
  N[3] = 0.25*(1.0-xi)*(1.0+eta);
}

PetscErrorCode MaterialPoint_AdvectRK1(DM dm_vp,Vec vp,PetscReal dt,DM dm_mpoint)
{
  PetscErrorCode    ierr;
  Vec               vp_l,coor_l;
  const PetscScalar *LA_vp;
  PetscInt          i,p,e,npoints,nel,npe;
  PetscInt          *mpfield_cell;
  PetscReal         *mpfield_coor;
  const PetscInt    *element_list;
  const PetscInt    *element;
  PetscScalar       xi_p[NSD],Ni[NODES_PER_EL];
  const PetscScalar *LA_coor;
  PetscScalar       dx[NSD];

  PetscFunctionBeginUser;
  ierr = DMGetCoordinatesLocal(dm_vp,&coor_l);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor_l,&LA_coor);CHKERRQ(ierr);

  ierr = DMGetLocalVector(dm_vp,&vp_l);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm_vp,vp,INSERT_VALUES,vp_l);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm_vp,vp,INSERT_VALUES,vp_l);CHKERRQ(ierr);
  ierr = VecGetArrayRead(vp_l,&LA_vp);CHKERRQ(ierr);

  ierr = DMDAGetElements(dm_vp,&nel,&npe,&element_list);CHKERRQ(ierr);
  ierr = DMSwarmGetLocalSize(dm_mpoint,&npoints);CHKERRQ(ierr);
  ierr = DMSwarmGetField(dm_mpoint,DMSwarmPICField_coor,NULL,NULL,(void**)&mpfield_coor);CHKERRQ(ierr);
  ierr = DMSwarmGetField(dm_mpoint,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);
  for (p=0; p<npoints; p++) {
    PetscReal         *coor_p;
    PetscScalar       vel_n[NSD*NODES_PER_EL],vel_p[NSD];
    const PetscScalar *x0;
    const PetscScalar *x2;

    e       = mpfield_cell[p];
    coor_p  = &mpfield_coor[NSD*p];
    element = &element_list[NODES_PER_EL*e];

    /* compute local coordinates: (xp-x0)/dx = (xip+1)/2 */
    x0 = &LA_coor[NSD*element[0]];
    x2 = &LA_coor[NSD*element[2]];

    dx[0] = x2[0] - x0[0];
    dx[1] = x2[1] - x0[1];

    xi_p[0] = 2.0 * (coor_p[0] - x0[0])/dx[0] - 1.0;
    xi_p[1] = 2.0 * (coor_p[1] - x0[1])/dx[1] - 1.0;
    if (PetscRealPart(xi_p[0]) < -1.0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"value (xi) too small %1.4e [e=%D]\n",(double)PetscRealPart(xi_p[0]),e);
    if (PetscRealPart(xi_p[0]) >  1.0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"value (xi) too large %1.4e [e=%D]\n",(double)PetscRealPart(xi_p[0]),e);
    if (PetscRealPart(xi_p[1]) < -1.0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"value (eta) too small %1.4e [e=%D]\n",(double)PetscRealPart(xi_p[1]),e);
    if (PetscRealPart(xi_p[1]) >  1.0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"value (eta) too large %1.4e [e=%D]\n",(double)PetscRealPart(xi_p[1]),e);

    /* evaluate basis functions */
    EvaluateBasis_Q1(xi_p,Ni);

    /* get cell nodal velocities */
    for (i=0; i<NODES_PER_EL; i++) {
      PetscInt nid;

      nid = element[i];
      // Note: this differs from ex70, since we only have 2 dof (not 3 - they have pressure)
      vel_n[NSD*i+0] = LA_vp[(NSD)*nid+0];
      vel_n[NSD*i+1] = LA_vp[(NSD)*nid+1];
    }

    /* interpolate velocity */
    vel_p[0] = vel_p[1] = 0.0;
    for (i=0; i<NODES_PER_EL; i++) {
      vel_p[0] += Ni[i] * vel_n[NSD*i+0];
      vel_p[1] += Ni[i] * vel_n[NSD*i+1];
    }

    coor_p[0] += dt * PetscRealPart(vel_p[0]);
    coor_p[1] += dt * PetscRealPart(vel_p[1]);
  }

  ierr = DMSwarmRestoreField(dm_mpoint,DMSwarmPICField_cellid,NULL,NULL,(void**)&mpfield_cell);CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(dm_mpoint,DMSwarmPICField_coor,NULL,NULL,(void**)&mpfield_coor);CHKERRQ(ierr);
  ierr = DMDARestoreElements(dm_vp,&nel,&npe,&element_list);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(vp_l,&LA_vp);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm_vp,&vp_l);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(coor_l,&LA_coor);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* A deprecated function, removed from DMStag but still used here */
PetscErrorCode DMStagGetVertexDMDA(DM dmstag,PetscInt dof,DM *dmda)
{
  PetscErrorCode  ierr;
  PetscInt        dim,i,j,stencilWidth;
  DM_Stag * const stag = (DM_Stag*)dmstag->data;
  DMDAStencilType stencilType;
  PetscInt        *l[DMSTAG_MAX_DIM];

  PetscFunctionBegin;
  ierr = DMGetDimension(dmstag,&dim);CHKERRQ(ierr);

  switch (stag->ghostStencilType) {
    case DMSTAG_GHOST_STENCIL_NONE: stencilType = DMDA_STENCIL_STAR; stencilWidth= 0                       ; break;
    case DMSTAG_GHOST_STENCIL_STAR: stencilType = DMDA_STENCIL_STAR; stencilWidth = stag->ghostStencilWidth; break;
    case DMSTAG_GHOST_STENCIL_BOX : stencilType = DMDA_STENCIL_BOX ; stencilWidth = stag->ghostStencilWidth; break;
    default: SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported Stencil Type %d",stag->ghostStencilType);
  }

  for (i=0;i<dim;++i) {
    ierr = PetscMalloc1(stag->nRanks[i],&l[i]);CHKERRQ(ierr);
    for (j=0; j<stag->nRanks[i]; ++j) l[i][j] = stag->l[i][j];
    l[i][stag->nRanks[i]-1] += 1; /* extra vertex on last rank in each dimension */
  }

  switch (dim) {
    case 2 :
      ierr = DMDACreate2d(PetscObjectComm((PetscObject)dmstag),
            stag->boundaryType[0],
            stag->boundaryType[1],
            stencilType,
            stag->N[0] + 1,       /* Note + 1 */
            stag->N[1] + 1,       /* Note + 1 */
            stag->nRanks[0],
            stag->nRanks[1],
            dof,
            stencilWidth,
            l[0],
            l[1],
            dmda);CHKERRQ(ierr);
      break;
    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"DMStagGetVertexDMDA not implemented for dim %d (but simple to add)",dim);
  }

  for (i=0;i<dim;++i) {
    ierr = PetscFree(l[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}