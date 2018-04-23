#include <petsc/private/dmstagimpl.h>
#include "system.h"
#include "petscdmstag.h"

PetscReal getRho(Ctx ctx,PetscReal x,PetscReal y)
{
  if (x < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->rho1;
  } else {
    return ctx->rho2;
  }
}

PetscReal getEta(Ctx ctx,PetscReal x,PetscReal y)
{
  if (ctx->isoviscous ) return ctx->etaCharacteristic;
  if (x < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->eta1;
  } else {
    return ctx->eta2;
  }
}

PetscErrorCode CreateSystem(Ctx ctx,Mat *pA,Vec *pb)
{
  typedef struct {PetscScalar etaCorner,etaElement,rho;} ParamData;
  typedef struct {PetscScalar vy,vx,p;} StokesData; 
  PetscErrorCode  ierr;
  DM_Stag         *stag;
  ParamData       **arrParam;
  StokesData      **arrb;
  PetscInt        startGhost[2],nGhost[2],n[2],start[2],N[2],nDummy[2],i,j;
  Vec             b,bLocal,paramLocal;
  Mat             A;
  DM              stokesGrid = ctx->stokesGrid,paramGrid = ctx->paramGrid;
  const PetscReal hx = ctx->hxCharacteristic; // assume constant
  const PetscReal hy = ctx->hyCharacteristic; // assume constant
  const PetscReal hxi2 = 1.0/(hx*hx);
  const PetscReal hyi2 = 1.0/(hy*hy);
  const PetscReal hxhyi = 1.0/(hx*hy);

  PetscFunctionBeginUser;
  stag = (DM_Stag*)stokesGrid->data;

  // Create RHS vector and local version
  ierr = DMCreateGlobalVector(stokesGrid,pb);CHKERRQ(ierr);
  b = *pb;

  // Create Matrix
  // Here, we drastically overestimate the number of nonzeros
  {
    PetscInt sizeLocal;
    ISLocalToGlobalMapping ltogmap;
    const PetscInt width = 11; // largest stencil (momentum)
    ierr = VecGetLocalSize(b,&sizeLocal);CHKERRQ(ierr);
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,sizeLocal,sizeLocal,PETSC_DETERMINE,PETSC_DETERMINE,width,NULL,width,NULL,pA);CHKERRQ(ierr);
    A = *pA;
    ierr = DMGetLocalToGlobalMapping(stokesGrid,&ltogmap);CHKERRQ(ierr);
    ierr = MatSetLocalToGlobalMapping(A,ltogmap,ltogmap);CHKERRQ(ierr);
  }

  // Get access to parameters array
  ierr = DMCreateLocalVector(paramGrid,&paramLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(paramGrid,ctx->param,INSERT_VALUES,paramLocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(paramGrid,ctx->param,INSERT_VALUES,paramLocal);CHKERRQ(ierr);
  ierr = DMStagVecGetArray(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr); // should be Read

  ierr = DMStagGetCorners(paramGrid,&start[0],&start[1],NULL,&n[0],&n[1],NULL,&nDummy[0],&nDummy[1],NULL);CHKERRQ(ierr);
  ierr = DMStagGetGhostCorners(paramGrid,&startGhost[0],&startGhost[1],NULL,&nGhost[0],&nGhost[1],NULL);CHKERRQ(ierr);
  ierr = DMStagGetGlobalSizes(stokesGrid,&N[0],&N[1],NULL);CHKERRQ(ierr);

  // Populate RHS
  ierr = DMCreateLocalVector(stokesGrid,&bLocal);CHKERRQ(ierr);
  ierr = VecZeroEntries(bLocal);CHKERRQ(ierr);
  ierr = DMStagVecGetArray(stokesGrid,bLocal,&arrb);CHKERRQ(ierr);
  if (ctx->gx != 0.0 ) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Not implemented for non-zero gx");
  // Note: we waste some time here computing things on the ghost nodes. 
  for (j=PetscMax(1,start[1]); j<start[1]+n[1]; ++j) {
    for(i=PetscMax(1,start[0]); i<PetscMin(N[0]-1,start[0]+n[0]); ++i) {
      arrb[j][i].vy = -ctx->gy * 0.5 * (arrParam[j][i].rho + arrParam[j-1][i].rho);
      arrb[j][i].vx = 0; // assume gx = 0
      arrb[j][i].p = 0;
    }
  }
  ierr = DMStagVecRestoreArray(stokesGrid,bLocal,&arrb);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(stokesGrid,bLocal,INSERT_VALUES,b);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(stokesGrid,bLocal,INSERT_VALUES,b);CHKERRQ(ierr);
  ierr = VecDestroy(&bLocal);CHKERRQ(ierr);

  // Populate matrix
  // Note: this is not tested for extremely small domains! If you have less than 3 elements in any direction, this might not work
  // We loop over the entire local domain at once
  // TODO Introduce a way to set with "stencil" objects
  for (j=startGhost[1]; j<startGhost[1]+nGhost[1]; ++j) { // i and j refer to global (natural) element numbers
    const PetscInt jLocal = j-startGhost[1];
    for (i=startGhost[0]; i<startGhost[0]+nGhost[0]; ++i) {
      const PetscInt iLocal = i-startGhost[0];
      PetscScalar    vals[11]; // 11 is most entries we will add in one row
      PetscInt       row,cols[11],epr,poffset,vyoffset,vxoffset,e,ncols;
      PetscReal      Kbound=ctx->Kbound,Kcont=ctx->Kcont;

      e = stag->entriesPerElement;
      epr = stag->entriesPerElementRowGhost;

      vyoffset = stag->dof[0]; // 2d only
      vxoffset = stag->dof[0] + stag->dof[1]; // 2d only
      poffset = stag->dof[0] + 2*stag->dof[1]; // 2d only
      
      // vy
      // For each element, we update the edge below. Thus, we ignore all ghost/dummy cells except the dummies on the top
      if (i >= start[0] && i < start[0] + n[0] && j>= start[1] && j < start[1] + n[1] +nDummy[1]) {
        row = e*iLocal + epr*jLocal + vyoffset;
        if (j == N[1] || j == 0) { // if we are at the top boundary dummy element, or at the bottom
          ncols = 1;
          cols[0] = row;
          vals[0] = Kbound;
        } else if (i == N[0]-1) { // right boundary (not the dummy element) (and not top or bottom)
          ncols = 2;
          cols[0] = e*(iLocal-1) + epr*jLocal + vyoffset;
          cols[1] = row; 
          vals[0] = -Kbound;
          vals[1] = Kbound;
        } else if (i == 0) { // left boundary (and not top or bottom)
          ncols = 2;
          cols[0] = row; 
          cols[1] = e*(iLocal+1) + epr*jLocal + vyoffset;
          vals[0] = -Kbound;
          vals[1] = Kbound;
        } else { // interior
          ncols = 11;
          const PetscReal etaUp    = arrParam[j  ][i  ].etaElement;
          const PetscReal etaDown  = arrParam[j-1][i  ].etaElement;
          const PetscReal etaLeft  = arrParam[j  ][i  ].etaCorner;
          const PetscReal etaRight = arrParam[j  ][i+1].etaCorner;

          cols[0] = e*iLocal     + epr*(jLocal-1) + vyoffset;
          cols[1] = e*(iLocal-1) + epr*jLocal     + vyoffset;
          cols[2] = e*iLocal     + epr*jLocal     + vyoffset;
          cols[3] = e*(iLocal+1) + epr*jLocal     + vyoffset;
          cols[4] = e*iLocal     + epr*(jLocal+1) + vyoffset;
          vals[0] =  2.0 * etaDown * hyi2;
          vals[1] =        etaLeft * hxi2;
          vals[2] = -2.0 * (etaDown + etaUp) * hyi2 - (etaLeft + etaRight) * hxi2;
          vals[3] =        etaRight * hxi2;
          vals[4] =  2.0 * etaUp    * hyi2;

          cols[5] = e*iLocal     + epr*(jLocal-1) + vxoffset;
          cols[6] = e*(iLocal+1) + epr*(jLocal-1) + vxoffset;
          cols[7] = e*iLocal     + epr*jLocal     + vxoffset;
          cols[8] = e*(iLocal+1) + epr*jLocal     + vxoffset;
          vals[5] =  etaLeft  * hxhyi;
          vals[6] = -etaRight * hxhyi;
          vals[7] = -etaLeft  * hxhyi;
          vals[8] =  etaRight * hxhyi;

          cols[9]  = e*iLocal + epr*(jLocal-1) + poffset;
          cols[10] = e*iLocal + epr*jLocal     + poffset;
          vals[9]  = Kcont/hy;
          vals[10] = -Kcont/hy;
        }
        ierr = MatSetValuesLocal(A,1,&row,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
      }
      
      // vx 
      if (i >= start[0] && i < start[0] + n[0] + nDummy[0]  && j>= start[1] && j < start[1] + n[1]) {
        row = e*iLocal + epr*jLocal + vxoffset;
        if (i == 0 || i == N[0] ){
          ncols = 1;
          cols[0] = row;
          vals[0] = Kbound;
        } else if (j == N[1]-1) {
          ncols = 2;
          cols[0] = e*iLocal + epr*(jLocal-1) + vxoffset;
          cols[1] = row; 
          vals[0] = -Kbound;
          vals[1] = Kbound;
        } else if (j == 0) {
          ncols = 2;
          cols[0] = row; 
          cols[1] = e*iLocal + epr*(jLocal+1) + vxoffset;
          vals[0] = -Kbound;
          vals[1] = Kbound;
        } else {
          ncols = 11;
          const PetscReal etaUp    = arrParam[j+1][i ].etaCorner;
          const PetscReal etaDown  = arrParam[j ][i  ].etaCorner;
          const PetscReal etaLeft  = arrParam[j ][i-1].etaElement;
          const PetscReal etaRight = arrParam[j ][i  ].etaElement;

          cols[0] = e*(iLocal-1) + epr*jLocal     + vyoffset;
          cols[1] = e*iLocal     + epr*jLocal     + vyoffset;
          cols[2] = e*(iLocal-1) + epr*(jLocal+1) + vyoffset;
          cols[3] = e*iLocal     + epr*(jLocal+1) + vyoffset;
          vals[0] =  etaDown * hxhyi;
          vals[1] = -etaDown * hxhyi;
          vals[2] = -etaUp   * hxhyi;
          vals[3] =  etaUp   * hxhyi;

          cols[4] =  e*iLocal     + epr*(jLocal-1) + vxoffset;
          cols[5] =  e*(iLocal-1) + epr*jLocal     + vxoffset;
          cols[6] =  e*iLocal     + epr*jLocal     + vxoffset;
          cols[7] =  e*(iLocal+1) + epr*jLocal     + vxoffset;
          cols[8] =  e*iLocal     + epr*(jLocal+1) + vxoffset;
          vals[4] =        etaDown  * hyi2;
          vals[5] =  2.0 * etaLeft  * hxi2;
          vals[6] = -2.0 * (etaLeft + etaRight) * hxi2 - (etaUp + etaDown) * hyi2;
          vals[7] =  2.0 * etaRight * hxi2;
          vals[8] =        etaUp    * hyi2;

          cols[9]  = e*(iLocal-1) + epr*jLocal + poffset;
          cols[10] = e*iLocal     + epr*jLocal + poffset;
          vals[9]  = Kcont/hx;
          vals[10] = -Kcont/hx;
        }
        ierr = MatSetValuesLocal(A,1,&row,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
      }

      // Pressure/mass
      if (i >= start[0] && i < start[0] + n[0] && j>= start[1] && j < start[1] + n[1]) { // don't consider any ghost/dummy nodes
      row = e*iLocal + epr*jLocal + poffset;
        /* We impose 4 pressure BCs, at the corners, since our boundary conditions
           leave the corner nodes out of any momentum equation,
           and a fixed node, as implied by Exercise 7.1.
           Using an implied ghost velocity node (as also described in the text)
           seems cleaner, though */
        if (i == ctx->pFixx && j == ctx->pFixy){
          cols[0] = row;
          vals[0] = Kbound;
          ncols = 1;
        } else if (i == 0 && j == 0) {
          cols[0] = row;
          cols[1] = e*(iLocal+1) + epr*jLocal + poffset;
          vals[0] = Kbound;
          vals[1] = -Kbound;
          ncols = 2;
        } else if (i == 0 && j == N[1]-1) {
          cols[0] = row;
          cols[1] = e*(iLocal+1) + epr*jLocal + poffset;
          vals[0] = Kbound;
          vals[1] = -Kbound;
          ncols = 2;
        } else if (i == N[0]-1 && j == 0) {
          cols[0] = e*(iLocal-1) + epr*jLocal + poffset;
          cols[1] = row;
          vals[0] = -Kbound;
          vals[1] = Kbound;
          ncols = 2;
        } else if (i == N[0]-1 && j == N[1]-1) {
          cols[0] = e*(iLocal-1) + epr*jLocal + poffset;
          cols[1] = row;
          vals[0] = -Kbound;
          vals[1] = Kbound;
          ncols = 2;
        } else {
          const PetscReal hx = ctx->hxCharacteristic; // assume constant
          const PetscReal hy = ctx->hyCharacteristic; // assume constant
          cols[0] = e*iLocal + epr*jLocal + vxoffset;
          vals[0] = -Kcont/hx;
          cols[1] = e*(iLocal+1) + epr*jLocal + vxoffset;
          vals[1] = Kcont/hx;
          cols[2] = e*iLocal + epr*jLocal + vyoffset;
          vals[2] = -Kcont/hy;
          cols[3] = e*iLocal + epr*(jLocal+1) + vyoffset;
          vals[3] = Kcont/hy;
          ncols = 4;
        }
        ierr = MatSetValuesLocal(A,1,&row,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = DMStagVecRestoreArray(paramGrid,paramLocal,&arrParam);CHKERRQ(ierr); // should be Read
  ierr = VecDestroy(&paramLocal);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
