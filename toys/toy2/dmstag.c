#include <petsc/private/dmimpl.h>
#include "dmstag.h"

PetscErrorCode DMStagCreate2d( MPI_Comm comm, DMBoundaryType bndx ,DMBoundaryType bndy, PetscInt M,PetscInt N, PetscInt m,PetscInt n, PetscInt dof0,PetscInt dof1,PetscInt dof2, DMStagGhostStencilType ghostStencil, PetscInt ghostStencilWidth, DM* dm)
{
  PetscErrorCode ierr;
  DM_Stag        *stag;
  PetscInt       i;

  PetscFunctionBegin;

  ierr = DMCreate(comm,dm);CHKERRQ(ierr);
  ierr = DMSetType(*dm,DMSTAG);CHKERRQ(ierr);
  stag = (DM_Stag*)(*dm)->data;

  /* Global sizes and flags (derivative quantities set in DMSetUp_Stag) */
  stag->dim = 2; // We don't wrap everything with setters/getters in this prototype
  stag->boundaryType[0] = bndx; stag->boundaryType[1] = bndy;
  stag->N[0] = M; stag->N[1] = N;
  stag->nprocs[0] = m; stag->nprocs[1] = n; // could be PETSC_DECIDE
  stag->ghostStencil = ghostStencil;
  stag->ghostStencilWidth = ghostStencilWidth;
  stag->dof[0] = dof0; stag->dof[1] = dof1; stag->dof[2] = dof2; stag->dof[3] = 0;
  for (i=0; i<DMSTAG_MAX_STRATA; ++i) stag->stratumActive[i] = stag->dof[i] > 0; // for now, we always set these based on whether dof is positive
  stag->nprocs[0] = m; stag->nprocs[1] = n; // These are adjusted later in DMSetUp_Stag
  PetscFunctionReturn(0);
}

static PetscErrorCode DMDestroy_Stag(DM dm)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  if (stag->gton) {
    ierr = VecScatterDestroy(&stag->gton);CHKERRQ(ierr);
  }
  if (stag->ntol) {
    ierr = VecScatterDestroy(&stag->ntol);CHKERRQ(ierr);
  }
  if (stag->lton) {
    ierr = VecScatterDestroy(&stag->lton);CHKERRQ(ierr);
  }
  if (stag->gtol) {
    ierr = VecScatterDestroy(&stag->gtol);CHKERRQ(ierr);
  }

  ierr = PetscFree(stag);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMCreateGlobalVector_Stag(DM dm,Vec *vec)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  ierr = VecCreateMPI(PetscObjectComm((PetscObject)dm),stag->entries,PETSC_DECIDE,vec);CHKERRQ(ierr);
  // Could set some ops, as DMDA does
  // Note that there they set an ltogmap in the vector as well
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreateLocalVector_Stag(DM dm,Vec *vec)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidPointer(vec,2);
  ierr = VecCreateSeq(PETSC_COMM_SELF,stag->entriesGhost,vec);CHKERRQ(ierr);
  ierr = VecSetDM(*vec,dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
Note there are several orderings in play here.
In all cases, non-element dof are associated with the element that they are below/left/behind, and the order in 2D proceeds vertex/bottom edge/left edge/element (with all dof on each together).
Also in all cases, only subdomains which are the last in their dimension have partial elements.

1) "Natural" Ordering. Number adding each full or partial (on the right or top) element, starting at the bottom left (i=0,j=0) and proceeding across the entire domain, row by row to get a global numbering.
2) Global ("PETSc") ordering. The same as natural, but restricted to each domain. So, traverse all elements (again starting at the bottom left and going row-by-row) on proc 0, then continue numbering with proc 1, and so on.
3) Local ordering. Including ghost elements (both interior and on the right/top/front to complete partial elements), use the same convention to create a local numbering.
*/

PetscErrorCode DMLocalToGlobalBegin_Stag(DM dm,Vec l,InsertMode mode,Vec g)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(l,VEC_CLASSID,2);
  PetscValidHeaderSpecific(g,VEC_CLASSID,4);
  if (mode == ADD_VALUES) {
    SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Unsupported InsertMode");
    // Would be easy ( see DMDA)
  } else if (mode == INSERT_VALUES) {
    ierr = VecScatterBegin(stag->gtol,l,g,mode,SCATTER_REVERSE_LOCAL);CHKERRQ(ierr);
  } else SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Unsupported InsertMode");
  PetscFunctionReturn(0);
}

PetscErrorCode DMLocalToGlobalEnd_Stag(DM dm,Vec l,InsertMode mode,Vec g)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(l,VEC_CLASSID,2);
  PetscValidHeaderSpecific(g,VEC_CLASSID,4);
  if (mode == ADD_VALUES) {
    SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Unsupported InsertMode");
    // Would be easy ( see DMDA)
  } else if (mode == INSERT_VALUES) {
    ierr = VecScatterEnd(stag->gtol,l,g,mode,SCATTER_REVERSE_LOCAL);CHKERRQ(ierr);
  } else SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Unsupported InsertMode");
  PetscFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalBegin_Stag(DM dm,Vec g,InsertMode mode,Vec l)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(l,VEC_CLASSID,2);
  PetscValidHeaderSpecific(g,VEC_CLASSID,4);
  ierr = VecScatterBegin(stag->gtol,g,l,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMGlobalToLocalEnd_Stag(DM dm,Vec g,InsertMode mode,Vec l)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(l,VEC_CLASSID,2);
  PetscValidHeaderSpecific(g,VEC_CLASSID,4);
  ierr = VecScatterEnd(stag->gtol,g,l,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMCreateCoordinateDM_Stag(DM dm,DM *dmc)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBegin;
  if (stag->dim == 2) {
    ierr = DMStagCreate2d(
         PetscObjectComm((PetscObject)dm),
         stag->boundaryType[0],stag->boundaryType[1],
         stag->N[0],stag->N[1],
         stag->nprocs[0],stag->nprocs[1],
         stag->stratumActive[0] ? stag->dim : 0, // dim (2) dof on each active stratum
         stag->stratumActive[1] ? stag->dim : 0,
         stag->stratumActive[2] ? stag->dim : 0,
         stag->ghostStencil,
         stag->ghostStencilWidth,
         dmc);
  } else {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Set up not implemented for dim!=2");
  }
  ierr = DMSetUp(*dmc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagSetUniformCoordinates(DM dm,PetscReal xmin,PetscReal xmax,PetscReal ymin,PetscReal ymax,PetscReal zmin,PetscReal zmax)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data, *stagc;
  DM             dmc;
  PetscReal      hx,hy,hz;
  PetscInt       entriesPerElement,entriesPerEdge,entriesPerCorner,entriesPerElementRow,i,j,d;
  Vec            coords;
  PetscScalar    *arr;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDM(dm, &dmc);CHKERRQ(ierr);
  stagc = (DM_Stag*)dmc->data; // Ugly!! For next impl, protect everything..

  if (stag->dim != 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"DMSTagSetUniformCoordinates only implemented for dim=2");

  entriesPerElement    = stagc->entriesPerElement;
  entriesPerEdge       = stagc->entriesPerEdge;
  entriesPerCorner     = stagc->entriesPerCorner;
  entriesPerElementRow = stagc->entriesPerElementRow;

  if (stag->boundaryType[0] != DM_BOUNDARY_NONE || stag->boundaryType[1] != DM_BOUNDARY_NONE || stag->boundaryType[2] != DM_BOUNDARY_NONE) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Uniform coordinates only implemented for boundary type DM_BOUNDARY_NONE");
  hx = (xmax-xmin) / stag->N[0];
  hy = (ymax-ymin) / stag->N[1];
  hz = (zmax-zmin) / stag->N[2];

  ierr = DMCreateGlobalVector(dmc,&coords);CHKERRQ(ierr);
  ierr = VecGetArray(coords,&arr);CHKERRQ(ierr);

  // Note that this assumes 2 dof per point
  if ((stagc->stratumActive[0] && stagc->dof[0]!=2) || (stagc->stratumActive[1] && stagc->dof[1]!=2) || (stagc->stratumActive[2] && stagc->dof[2]!=2)) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Uniform coordinate routine assumes 2 dof per stratum");
  for (j=0; j<stag->n[1]; ++j){ // local element index
    const PetscInt jg = j + stag->start[1]; // global element index
    for (i=0; i<stag->n[0]; ++i){
      const PetscInt ig = i + stag->start[0];
      d = 0;
      if (stagc->stratumActive[0]) { // Note: all the branches are ugly, but would be perfectly predicted (and could/should be replaced bby a marked-constant value)
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hx*(ig);      ++d;
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hy*(jg);      ++d;
      }
      if (stagc->stratumActive[1]) {
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hx*(ig+0.5);  ++d;
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hy*(jg);      ++d;
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hx*(ig);      ++d;
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hy*(jg+0.5);  ++d;
      }
      if (stagc->stratumActive[2]) {
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hx*(ig+0.5);  ++d;
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hy*(jg+0.5);  ++d;
      }
    }
  }
  // Additional points on right boundary
  if (stag->lastproc[0]){
    i = stag->n[0];
    for (j=0; j<stag->n[1]; ++j) {
      const PetscInt ig = i + stag->start[0];
      const PetscInt jg = j + stag->start[1];
      d = 0;
      if (stagc->stratumActive[0]) {
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hx*(ig);      ++d;
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hy*(jg);      ++d;
      }
      if (stagc->stratumActive[1]) {
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hx*(ig);      ++d;
        arr[j*entriesPerElementRow + entriesPerElement*i+d]  = hy*(jg+0.5);  ++d;
      }
    }
  }
  // Additional points on top boundary
  if (stag->lastproc[1]) {
    j = stag->n[1];
    for (i=0; i<stag->n[0]; ++i){
      const PetscInt ig = i + stag->start[0];
      const PetscInt jg = j + stag->start[1];
      d = 0;
      if (stagc->stratumActive[0]) {
        arr[j*entriesPerElementRow + entriesPerEdge*i+d]  = hx*(ig);         ++d;
        arr[j*entriesPerElementRow + entriesPerEdge*i+d]  = hy*(jg);         ++d;
      }
      if (stagc->stratumActive[1]) {
        arr[j*entriesPerElementRow + entriesPerEdge*i+d]  = hx*(ig+0.5);     ++d;
        arr[j*entriesPerElementRow + entriesPerEdge*i+d]  = hy*(jg);         ++d;
      }
    }
  }
  // Additional top right points
  if (stag->lastproc[0] && stag->lastproc[1]){
    i = stag->n[0];
    j = stag->n[1];
    const PetscInt ig = i + stag->start[0];
    const PetscInt jg = j + stag->start[1];
    d = 0;
    if (stagc->stratumActive[0]) {
      arr[j*entriesPerElementRow + entriesPerEdge*i+d]  = hx*(ig);           ++d;
      arr[j*entriesPerElementRow + entriesPerEdge*i+d]  = hy*(jg);           ++d;
    }
  }
  ierr = VecRestoreArray(coords,&arr);CHKERRQ(ierr);
  ierr = DMSetCoordinates(dm,coords);CHKERRQ(ierr);
  ierr = PetscLogObjectParent((PetscObject)dm,(PetscObject)coords);CHKERRQ(ierr);
  ierr = VecDestroy(&coords);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMSetUp_Stag(DM dm)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag *) dm->data;
  PetscMPIInt    size,rank;
  PetscInt       i,j,m,n,M,N;
  PetscInt       *l[DMSTAG_MAX_DIM]; /* Arrays of number of elements/proc in each dimension (could accept as args like DMDA creation routines do)*/
  MPI_Comm       comm;
  PetscMPIInt    *neighbors;
  PetscInt       *globalOffsets;

  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

  if (stag->dim != 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Set up not implemented for dim!=2");

  m = stag->nprocs[0]; n = stag->nprocs[1]; // this is confusing, but we do it to as to borrow DMDA code below
  M = stag->N[0]; N = stag->N[1];

  /* Processor grid sizes  (from da2.c) */
  {
    if (m != PETSC_DECIDE) {
      if (m < 1) SETERRQ1(comm,PETSC_ERR_ARG_OUTOFRANGE,"Non-positive number of processors in X direction: %D",m);
      else if (m > size) SETERRQ2(comm,PETSC_ERR_ARG_OUTOFRANGE,"Too many processors in X direction: %D %d",m,size);
    }
    if (n != PETSC_DECIDE) {
      if (n < 1) SETERRQ1(comm,PETSC_ERR_ARG_OUTOFRANGE,"Non-positive number of processors in Y direction: %D",n);
      else if (n > size) SETERRQ2(comm,PETSC_ERR_ARG_OUTOFRANGE,"Too many processors in Y direction: %D %d",n,size);
    }

    if (m == PETSC_DECIDE || n == PETSC_DECIDE) {
      if (n != PETSC_DECIDE) {
        m = size/n;
      } else if (m != PETSC_DECIDE) {
        n = size/m;
      } else {
        /* try for squarish distribution */
        m = (PetscInt)(0.5 + PetscSqrtReal(((PetscReal)M)*((PetscReal)size)/((PetscReal)N)));
        if (!m) m = 1;
        while (m > 0) {
          n = size/m;
          if (m*n == size) break;
          m--;
        }
        if (M > N && m < n) {PetscInt _m = m; m = n; n = _m;}
      }
      if (m*n != size) SETERRQ(comm,PETSC_ERR_PLIB,"Unable to create partition, check the size of the communicator and input m and n ");
    } else if (m*n != size) SETERRQ(comm,PETSC_ERR_ARG_OUTOFRANGE,"Given Bad partition");

    if (M < m) SETERRQ2(comm,PETSC_ERR_ARG_OUTOFRANGE,"Partition in x direction is too fine! %D %D",M,m);
    if (N < n) SETERRQ2(comm,PETSC_ERR_ARG_OUTOFRANGE,"Partition in y direction is too fine! %D %D",N,n);

    stag->nprocs[0] = m; stag->nprocs[1] = n;
  }

  // TODO : make a proper View function at some point and/or supply some PetscInfo
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Partition: %d x %d\n",stag->nprocs[0],stag->nprocs[1]);CHKERRQ(ierr);

  /* Determine location of process in grid (these get extra boundary points on the last element)
     Order is x-fast, as usual */
  {
    stag->proc[0] = rank % stag->nprocs[0];
    stag->proc[1] = (rank - stag->proc[0])/stag->nprocs[0];
    for(i=0;i<2;++i) stag->firstproc[i] = !stag->proc[i];
    // These don't depend on the previous two
    stag->lastproc[0] =  !((rank+1)%stag->nprocs[0]);
    stag->lastproc[1] = rank >= size - stag->nprocs[0];
    ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] I am processor (%d,%d), last (%d,%d)\n",rank,stag->proc[0],stag->proc[1],(int)stag->lastproc[0],(int)stag->lastproc[1]);CHKERRQ(ierr);
  }

  /* Determine Locally owned region */
  // Divide equally, giving lower procs in each dimension more elements if needbe.
  // We don't accept stag->l from the user here but easily could.
  for (i=0;i<stag->dim;++i) {
    PetscInt j;
    const PetscInt Ni = stag->N[i], nprocsi = stag->nprocs[i];
    ierr = PetscMalloc1(stag->nprocs[i],&l[i]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Dim %d decomp: ",i);CHKERRQ(ierr);
    for (j=0;j<stag->nprocs[i];++j){
      l[i][j] = Ni/nprocsi + ((Ni % nprocsi) > j);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"%d ",l[i][j]);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }

  // Retrieve local size in stag->n
  for (i=0;i<stag->dim;++i) stag->n[i] = l[i][stag->proc[i]];
#if defined(PETSC_USE_DEBUG)
  for (i=0;i<stag->dim;++i) {
    PetscInt Ncheck,j;
    Ncheck = 0;
    for (j=0;j<stag->nprocs[i];++j) Ncheck += l[i][j];
    if (Ncheck != stag->N[i]) SETERRQ3(comm,PETSC_ERR_ARG_SIZ,"Local sizes in dimension %d don't add up. %d != %d\n",i,Ncheck,stag->N[i]);
  }
#endif

  // Compute starting elements (used in lton.c, but not otherwise essential)
  for (i=0; i<stag->dim; ++i) {
    stag->start[i] = 0;
    for(j=0;j<stag->proc[i];++j){
      stag->start[i] += l[i][j];
    }
  }

  /* Determine ranks of neighboring processes, using DMDA's convention

     n6 n7 n8
     n3    n5
     n0 n1 n2                                               */
  for (i=0; i<stag->dim; ++i) {
    if (stag->boundaryType[i] != DM_BOUNDARY_NONE) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Neighbor determination assumes boundary_none");
  }
  ierr = PetscMalloc1(9,&neighbors);CHKERRQ(ierr);
  neighbors[0] = stag->firstproc[0] || stag->firstproc[1] ? -1 : rank - 1 - stag->nprocs[0]; // down left
  neighbors[1] =                       stag->firstproc[1] ? -1 : rank     - stag->nprocs[0]; // down
  neighbors[2] = stag->lastproc[0]  || stag->firstproc[1] ? -1 : rank + 1 - stag->nprocs[0]; // down right
  neighbors[3] = stag->firstproc[0]                       ? -1 : rank - 1                  ; // left
  neighbors[4] =                                                 rank                      ; // here
  neighbors[5] = stag->lastproc[0]                        ? -1 : rank + 1                  ; // right
  neighbors[6] = stag->firstproc[0] || stag->lastproc[1]  ? -1 : rank - 1 + stag->nprocs[0]; // up left 
  neighbors[7] =                       stag->lastproc[1]  ? -1 : rank     + stag->nprocs[0]; // up
  neighbors[8] = stag->lastproc[0]  || stag->lastproc[1]  ? -1 : rank + 1 + stag->nprocs[0]; // up right

  /* Define useful sizes */
  if (stag->dim == 2){
    stag->entriesPerElement          = stag->dof[0] + 2*stag->dof[1] + stag->dof[2];
    stag->entriesPerEdge             = stag->dof[0] + stag->dof[1];
    stag->entriesPerCorner           = stag->dof[0];
    stag->entriesPerElementRow       = stag->n[0]*stag->entriesPerElement + (stag->lastproc[0] ? stag->entriesPerEdge : 0);
    stag->entries                    = stag->n[1]*stag->entriesPerElementRow +  (stag->lastproc[1] ? stag->n[0]*stag->entriesPerEdge : 0 ) + (stag->lastproc[0] && stag->lastproc[1] ? stag->entriesPerCorner: 0);
    stag->entriesPerElementRowGlobal = stag->N[0]*stag->entriesPerElement + stag->entriesPerEdge;
  } else {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"not implemented for dim!=2");
  }   

  /* Compute offsets for each rank into global vectors */
  // Note that it get a bit silly computing som many different twiddly things, but in the next impl we will at least know ahed of time what we need,
  // and be able to assemble a smaller set of things, and during this setup organize them better.
  {
    ierr = PetscMalloc1(size,&globalOffsets);CHKERRQ(ierr);
    globalOffsets[0] = 0;
    PetscInt count = 1; // note the count is offset by 1 here. We add the size of the previous rank
    for (j=0; j<stag->nprocs[1]-1; ++j) {
      const PetscInt nnj = l[1][j];
      for (i=0; i<stag->nprocs[0]-1; ++i) {
        const PetscInt nni = l[0][i];
        globalOffsets[count] = globalOffsets[count-1] + nnj*nni*stag->entriesPerElement; // No right/top/front boundaries
        ++count;
      }
      {
        i = stag->nprocs[0]-1;
        const PetscInt nni = l[0][i];
        globalOffsets[count] = globalOffsets[count-1] + nnj*(nni*stag->entriesPerElement + stag->entriesPerEdge); // extra vert/edge per row, on the right
        ++count;
      }
    }
    {
      j = stag->nprocs[1]-1;
      const PetscInt nnj = l[1][j];
      for (i=0; i<stag->nprocs[0]-1; ++i) {
        const PetscInt nni = l[0][i];
        globalOffsets[count] = globalOffsets[count-1] + (nnj*stag->entriesPerElement + stag->entriesPerEdge)*nni; // extra vert/edge per column, on the top
        ++count;
      }
      {
        /* Don't need this last one!
        i = stag->nprocs[0]-1;
        const PetscInt nni = l[1][i];
        */
      }
    }
  }

  /* Define ghosted/local sizes */
  {
    // Note - for a elements-only DMStag, the extra elements on the edges aren't necessary (but whatever since then you might as well be using DMDA)
    switch (stag->ghostStencil) {
      case DMSTAG_GHOST_STENCIL_NONE : // only the extra one on the right/top edges
        for (i=0;i<stag->dim;++i) {
          stag->nghost[i] = stag->n[i];
          stag->startGhost[i] = stag->start[i];
          if (stag->lastproc[i]) stag->nghost[i] += 1;
        }
        break;
      case DMSTAG_GHOST_STENCIL_STAR : // allocate the corners but don't use them
      case DMSTAG_GHOST_STENCIL_BOX :
        for (i=0;i<stag->dim;++i) {
          stag->nghost[i] = stag->n[i];
          stag->startGhost[i] = stag->start[i];
          if (!stag->firstproc[i]) {
            stag->nghost[i]     += stag->ghostStencilWidth; // add interior ghost elements
            stag->startGhost[i] -= stag->ghostStencilWidth;
          }
          if (!stag->lastproc[i]) {
            stag->nghost[i] += stag->ghostStencilWidth; // add interior ghost elements
          } else {
            stag->nghost[i] += 1; // one element on the boundary to complete blocking
          }
        }
        break;
      default :
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unrecognized ghost stencil type %d",stag->ghostStencil);
    }
    stag->entriesGhost = stag->nghost[0]*stag->nghost[1]*stag->entriesPerElement; // ghost size is easy, as all elements are the same
    stag->entriesPerElementRowGhost = stag->nghost[0]*stag->entriesPerElement; // local/ghost rep is missing the top/right/front edges/verts
  }

  ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] Local dimensions %d x %d (%d x %d ghosted), start %d,%d (%d,%d ghosted), global offset %d\n",rank,stag->n[0],stag->n[1],stag->nghost[0],stag->nghost[1],stag->start[0],stag->start[1],stag->startGhost[0],stag->startGhost[1],globalOffsets[rank]);CHKERRQ(ierr);

  /* Create the default section (note this uses local sizes) 
     We don't use this for now, but keep this for future reference. 
     */
#if 0 
  {
    // Note that these are all local sizes, not taking any ghost points into account
    // Note that we don't use this to do local->global scatters. Maybe we should, but for this prototype we will just do a more direct vec scatter approach.
    PetscInt npoints[DMSTAG_MAX_STRATA],nPointsTotal,pointsPerElement,pointsPerEdge,pointsPerRow;
    PetscInt extrax = stag->lastproc[0] ? 1 : 0, extray = stag->lastproc[1] ? 1 : 0;
    PetscSection   section;

    // Numbers of points per element and per edge on the boundary
    if (stag->dim == 2){
      pointsPerElement = ((PetscInt) stag->stratumActive[0]) + 2*((PetscInt) stag->stratumActive[1]) + ((PetscInt) stag->stratumActive[2]);
      pointsPerEdge = ((PetscInt) stag->stratumActive[0]) + ((PetscInt) stag->stratumActive[1]);
      pointsPerRow = pointsPerElement * stag->n[0] + (stag->lastproc[0] ? pointsPerEdge : 0);
    } else {
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"not implemented for dim!=2");
    }   

    // Compute some extra things just for output/debugging
    npoints[0] = stag->stratumActive[0] ? (stag->n[0] + extrax) * (stag->n[1] + extray)                                     : 0; // number of vertices
    npoints[1] = stag->stratumActive[1] ? (stag->n[0]    ) * (stag->n[1] + extray)  + (stag->n[0] + extrax) * ( stag->n[1]) : 0; // number of edges
    npoints[2] = stag->stratumActive[2] ? (stag->n[0]    ) * (stag->n[1]    )                                               : 0; // number of elements/faces
    npoints[3] = 0;
    nPointsTotal = 0; for (i=0;i<DMSTAG_MAX_STRATA;++i) nPointsTotal += npoints[i];
    {
      PetscMPIInt rank;
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d vertices, %d edges, %d faces %d total\n",rank,npoints[0],npoints[1],npoints[2],nPointsTotal);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] %d entries, %d / elrow, %d / +bnd edge %d ++bnd corner\n",rank,stag->entries,stag->entriesPerElementRow,stag->entriesPerEdge,stag->entriesPerCorner);CHKERRQ(ierr);
    }

    ierr = PetscSectionCreate(PetscObjectComm((PetscObject)(dm)),&section);CHKERRQ(ierr);
    ierr = PetscSectionSetChart(section,0,nPointsTotal);CHKERRQ(ierr);

    if (stag->dim !=2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"not implemented for dim!=2");
    if (stag->stratumActive[0] && stag->stratumActive[1] && stag->stratumActive[2]) {
      // Points per row of cells, except the last/top one
      for (j=0; j<stag->n[1]; ++j) {
        for (i=0; i<stag->n[0]; ++i) {
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 0, stag->dof[0]);CHKERRQ(ierr); // vertex
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 1, stag->dof[1]);CHKERRQ(ierr); // edge (y normal)
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 2, stag->dof[1]);CHKERRQ(ierr); // edge (x normal)
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 3, stag->dof[2]);CHKERRQ(ierr); // element
        }
        if (stag->lastproc[0]) {
          i = stag->n[0];
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 0, stag->dof[0]);CHKERRQ(ierr); // vertex
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 1, stag->dof[1]);CHKERRQ(ierr); // edge (x normal)
        }
      }
      if (stag->lastproc[1]) {
        j = stag->n[1];
        for (i=0; i<stag->n[0]; ++i) {
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerEdge*i + 0, stag->dof[0]);CHKERRQ(ierr); // vertex
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerEdge*i + 1, stag->dof[1]);CHKERRQ(ierr); // edge (y normal)
        }
        if (stag->lastproc[0]) {
          i = stag->n[0];
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerEdge*i + 0, stag->dof[0]);CHKERRQ(ierr); // vertex
        }
      }
    } else if (!stag->stratumActive[0] && stag->stratumActive[1] && stag->stratumActive[2]) { // "Stokes"
      for (j=0; j<stag->n[1]; ++j) {
        for (i=0; i<stag->n[0]; ++i) {
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 0, stag->dof[1]);CHKERRQ(ierr); // edge (y normal)
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 1, stag->dof[1]);CHKERRQ(ierr); // edge (x normal)
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 2, stag->dof[2]);CHKERRQ(ierr); // element
        }
        if (stag->lastproc[0]) {
          i = stag->n[0];
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 0, stag->dof[1]);CHKERRQ(ierr); // edge (x normal)
        }
      }
      if (stag->lastproc[1]) {
        j = stag->n[1];
        for (i=0; i<stag->n[0]; ++i) {
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerEdge*i + 0, stag->dof[1]);CHKERRQ(ierr); // edge (y normal)
        }
      }
    } else if (stag->stratumActive[0] && !stag->stratumActive[1] && stag->stratumActive[2]) { // Vertices and elements/faces
      // Points per row of cells, except the last/top one
      for (j=0; j<stag->n[1]; ++j) {
        for (i=0; i<stag->n[0]; ++i) {
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 0, stag->dof[0]);CHKERRQ(ierr); // vertex
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 1, stag->dof[2]);CHKERRQ(ierr); // element
        }
        if (stag->lastproc[0]) {
          i = stag->n[0];
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 0, stag->dof[0]);CHKERRQ(ierr); // vertex
        }
      }
      if (stag->lastproc[1]) {
        j = stag->n[1];
        for (i=0; i<stag->n[0]; ++i) {
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerEdge*i + 0, stag->dof[0]);CHKERRQ(ierr); // vertex
        }
        if (stag->lastproc[0]) {
          i = stag->n[0];
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerEdge*i + 0, stag->dof[0]);CHKERRQ(ierr); // vertex
        }
      }
    } else if (stag->stratumActive[0] && !stag->stratumActive[1] && !stag->stratumActive[2]) { // Vertices only 
      // This one so simple we don't follow the pattern (but we should, and condense all this even more)
      PetscInt rowsize = stag->n[0] + (stag->lastproc[0] ? 1 : 0); // 1 points per cell, 1 extra on the right edge
      PetscInt colsize = stag->n[1] + (stag->lastproc[1] ? 1 : 0); // 1 points per cell, 1 extra on the top
      for (j=0; j<colsize; ++j) {
        for (i=0; i<rowsize; ++i) {
          ierr = PetscSectionSetDof(section,j*rowsize+i,stag->dof[0]);CHKERRQ(ierr); // vertex
        }
      }
    } else if (!stag->stratumActive[0] && stag->stratumActive[1] && !stag->stratumActive[2]) { // Edges only
      for (j=0; j<stag->n[1]; ++j) {
        for (i=0; i<stag->n[0]; ++i) {
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 0, stag->dof[1]);CHKERRQ(ierr); // edge (y normal)
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 1, stag->dof[1]);CHKERRQ(ierr); // edge (x normal)
        }
        if (stag->lastproc[0]) {
          i = stag->n[0];
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerElement*i + 0, stag->dof[1]);CHKERRQ(ierr); // edge (x normal)
        }
      }
      if (stag->lastproc[1]) {
        j = stag->n[1];
        for (i=0; i<stag->n[0]; ++i) {
          ierr = PetscSectionSetDof(section,j*pointsPerRow + pointsPerEdge*i + 0, stag->dof[1]);CHKERRQ(ierr); // edge (y normal)
        }
      }
    } else if (!stag->stratumActive[0] && !stag->stratumActive[1] && stag->stratumActive[2]) { // elements/faces only
      for (j=0; j<stag->n[1]; ++j) {
        for (i=0; i<stag->n[0]; ++i) {
          // pointsPerElement = 1 here, of course
          ierr = PetscSectionSetDof(section,j*stag->n[0] + pointsPerElement*i + 0, stag->dof[2]);CHKERRQ(ierr); // face/element
        }
      }
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Not implemented for all options of deactivated strata yet!");
    }

    ierr = PetscSectionSetUp(section);CHKERRQ(ierr);
    ierr = DMSetDefaultSection(dm,section);CHKERRQ(ierr);
    ierr = PetscSectionDestroy(&section);CHKERRQ(ierr); // reference remains in DM
  }
#endif

  /* Create global-->local VecScatter and local->global ISLocalToGlobalMapping

     Note that this could (should?) be done by defining a PetscSF to complement
     the section above (which should include the ghost points?), but this 
     seems confusing and at first glance seems to require one to know the local
     indices of points on other processors.

     Note further here that the local vectors:
     - Are always an integral number of elements-worth of points. That is,
     they are a point in a cell, and all boundary points that are down/left/behind.
     This allows easy iteration over the vector, as its has blocks of the same size,
     at the cost of having unused points on the right/top/front. Thus, even
     with a stencil width of zero, on 1 processor, local and global vectors
     are different!
     - Include both these extra points to make complete cells, and ghost
     cells on the interior boundaries (fictitious processor boundaries)
     Thus, on the top/right/front boundaries, there are local entries that don't 
     participate in the local-global maps

     This has NOT been tested for a stencil width greater than one. 
     We assume that the size on each rank is greater than or equal to the 
     stencil width.
     */
  {
    PetscInt *idxLocal,*idxGlobal,*idxGlobalAll;
    PetscInt count,entriesToTransferTotal,i,j,d,ighostoffset,jghostoffset,ighostoffsetright,jghostoffsettop;
    IS isLocal,isGlobal;
    PetscInt jghost,ighost;
    PetscInt nNeighbor[2];

    // These offsets should always be non-negative, and describe how many
    // ghost elements precede the first non-ghost element in each dimensions
    // For am (elementwise) stencil width of 1, these would be 1 on interior boundaries
    ighostoffset = stag->start[0] - stag->startGhost[0]; // not the greatest names
    jghostoffset = stag->start[1] - stag->startGhost[1];
    ighostoffsetright = stag->startGhost[0]+stag->nghost[0] - (stag->start[0]+stag->n[0]); // not the greatest names
    jghostoffsettop = stag->startGhost[1]+stag->nghost[1] - (stag->start[1]+stag->n[1]);

    entriesToTransferTotal  = stag->nghost[1]*stag->nghost[0]*stag->entriesPerElement // overestimate
      + (stag->lastproc[0]                      ? stag->nghost[1] * (-stag->entriesPerElement + stag->entriesPerEdge  ) : 0)
      + (stag->lastproc[1]                      ? stag->nghost[0] * (-stag->entriesPerElement + stag->entriesPerEdge  ) : 0)
      + (stag->lastproc[0] && stag->lastproc[1] ? stag->dof[2] : 0); // replace double-counted element on corner (edges were each removed once)

    // allocate idxTo and idxFrom of the correct size
    ierr = PetscMalloc1(entriesToTransferTotal,&idxLocal);CHKERRQ(ierr);
    ierr = PetscMalloc1(entriesToTransferTotal,&idxGlobal);CHKERRQ(ierr);

    // Populate index arrays, one rank at a time
    // Note that boundary ghost elements are processed with this rank (neighbor 4).
    count = 0;

    // Here and below, we work with (i,j) describing element numbers within a neighboring rank's global ordering,
    // and ighost,jghost referring to element numbers within this ranks local (ghosted) ordering

    // Neighbor 0 (down left)
    if (!stag->firstproc[0] && !stag->firstproc[1]) {
    // We may be a ghosted boundary in x, but the neighbor never is
      const PetscInt neighborRank = neighbors[0];
      nNeighbor[0] = l[0][neighborRank % stag->nprocs[0]];
      nNeighbor[1] = l[1][(neighborRank/stag->nprocs[0])];
      const PetscInt entriesPerElementRowNeighbor = stag->entriesPerElement * l[0][stag->proc[0]-1]; // no possibility of a down left neighbor being a right/top boundary!
      for (jghost = 0; jghost<jghostoffset; ++jghost) {
        const PetscInt j = nNeighbor[1] - jghostoffset + jghost;
        for (ighost = 0; ighost<ighostoffset; ++ighost) {
          const PetscInt i = nNeighbor[0] - ighostoffset + ighost;
          for (d=0; d<stag->entriesPerElement; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
      }
    }

    // Neighbor 1 (down)
    if (!stag->firstproc[1]) {
    // We may be a ghosted boundary in x, in which case the nieghbor is also
      const PetscInt neighborRank = neighbors[1];
      nNeighbor[0] = l[0][neighborRank % stag->nprocs[0]];
      nNeighbor[1] = l[1][(neighborRank/stag->nprocs[0])];
      const PetscInt entriesPerElementRowNeighbor = stag->entriesPerElementRow; // rank below has same width and boundary status!
      for (jghost = 0; jghost<jghostoffset; ++jghost) {
        const PetscInt j = nNeighbor[1] - jghostoffset + jghost;
        for (ighost = ighostoffset; ighost<stag->nghost[0]-ighostoffsetright; ++ighost) {
          const PetscInt i = ighost - ighostoffset;
          for (d=0; d<stag->entriesPerElement; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
        if (stag->lastproc[0]) {
          // assumes a single ghost element on the boundaries
          const PetscInt ighost = stag->nghost[0]-1; // careful with these!!
          const PetscInt i = stag->n[0]; // or nNeighbor[0]
          // vertex
          for (d=0; d<stag->dof[0]; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
          // edge
          for (d=0; d<stag->dof[1]; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + stag->dof[0]                + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + stag->dof[0] + stag->dof[1] + d;
          }
        }
      }
    }

    // Neighbor 2 (down right)
    {
    // We are never a ghosted boundary in x, but the neighbor may be
    // Note that here we assume that none of the boudnary points on the neighboring rank are involved (only normal elements are in the ghost region)
    if (!stag->lastproc[0] && !stag->firstproc[1]) {
      const PetscInt neighborRank = neighbors[2];
      nNeighbor[0] = l[0][neighborRank % stag->nprocs[0]];
      nNeighbor[1] = l[1][(neighborRank/stag->nprocs[0])];
      const PetscInt entriesPerElementRowNeighbor = stag->proc[0] == stag->nprocs[0]-2 ?  nNeighbor[0]*stag->entriesPerElement + stag->entriesPerEdge : nNeighbor[0]*stag->entriesPerElement; // are we second-to-last in x?
      for (jghost = 0; jghost<jghostoffset; ++jghost) {
        const PetscInt j = nNeighbor[1] - jghostoffset + jghost;
        for (i=0; i<ighostoffsetright; ++i) {
          const PetscInt ighost = stag->nghost[0] - ighostoffsetright + i;
          for (d=0; d<stag->entriesPerElement; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
      }
    }
    }

    // Neighbor 3 (left)
    if (!stag->firstproc[0]) {
    // Our neighbor is never a ghosted boundary in x, but we may be
    // Here, we may be a ghosted boundary in y and thus so will our neighbor be
      const PetscInt neighborRank = neighbors[3];
      nNeighbor[0] = l[0][neighborRank % stag->nprocs[0]];
      nNeighbor[1] = l[1][(neighborRank/stag->nprocs[0])];
      const PetscInt entriesPerElementRowNeighbor = nNeighbor[0]*stag->entriesPerElement;
      for (jghost = jghostoffset;jghost<stag->nghost[1]-jghostoffsettop; ++jghost) {
        const PetscInt j = jghost-jghostoffset;
        for (ighost = 0; ighost<ighostoffset; ++ighost) {
          const PetscInt i = nNeighbor[0] - ighostoffset + ighost;
          for (d=0; d<stag->entriesPerElement; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
      }
      if (stag->lastproc[1]) {
        // assumes a single ghost element on the boundary in ecah column
        const PetscInt jghost = stag->nghost[1]-1; // careful with these
        const PetscInt j = stag->n[1]; // or nNeighbor[1];
        for (ighost = 0; ighost<ighostoffset; ++ighost) {
          const PetscInt i = nNeighbor[1] - ighostoffset + ighost;
          for (d=0; d<stag->entriesPerEdge; ++d) { // only vertices and horizontal edge (which are the first dof)
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerEdge    + d; // i moves by edge here
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
      }
    }

    // Interior/Resident-here-in-global elements ("Neighbor 4" - same rank)
    // *including* boundary ghosts
    for (j=0; j<stag->n[1]; ++j) {
      const PetscInt jghost = j + jghostoffset;
      for (i=0; i<stag->n[0]; ++i) {
        const PetscInt ighost = i + ighostoffset;
        for (d=0; d<stag->entriesPerElement; ++d) {
          idxGlobal[count ] = globalOffsets[neighbors[4]] + j     *stag->entriesPerElementRow       + i     *stag->entriesPerElement + d;
          idxLocal[count++] =                                     jghost*stag->entriesPerElementRowGhost  + ighost*stag->entriesPerElement + d;
        }
      }
      if (stag->lastproc[0]) {
        i = stag->n[0];
        const PetscInt ighost = i + ighostoffset;
        for (d=0; d<stag->dof[0]; ++d) { // vertex first
          idxGlobal[count ] = globalOffsets[neighbors[4]] // here
            + j     *stag->entriesPerElementRow       + i     *stag->entriesPerElement + d;
          idxLocal[count++] = jghost*stag->entriesPerElementRowGhost  + ighost*stag->entriesPerElement + d;
        }
        for (d=0; d<stag->dof[1]; ++d) { // then left ege (skipping bottom edge)
          idxGlobal[count ] = globalOffsets[neighbors[4]] + j     *stag->entriesPerElementRow       + i     *stag->entriesPerElement + stag->dof[0]                + d; // global doesn't have bottom edge
          idxLocal[count++] =                                     jghost*stag->entriesPerElementRowGhost  + ighost*stag->entriesPerElement + stag->dof[0] + stag->dof[1] + d; // local does have bottom edge (unused)
        }
      }
    }
    if (stag->lastproc[1]) {
      j = stag->n[1];
      const PetscInt jghost = j + jghostoffset;
      for (i=0; i<stag->n[0]; ++i) {
        const PetscInt ighost = i + ighostoffset;
        for (d=0; d<stag->entriesPerEdge; ++d) { // vertex and bottom edge (which are the first entries)
          idxGlobal[count ] = globalOffsets[neighbors[4]] + j     *stag->entriesPerElementRow       + i     *stag->entriesPerEdge    + d; // note i increment by entriesPerEdge
          idxLocal[count++] =                                     jghost*stag->entriesPerElementRowGhost  + ighost*stag->entriesPerElement + d; // note i increment by entriesPerElement, still
        }

      }
      if (stag->lastproc[0]) {
        i = stag->n[0];
        const PetscInt ighost = i + ighostoffset;
        for (d=0; d<stag->entriesPerCorner; ++d) { // vertex only
          idxGlobal[count ] = globalOffsets[neighbors[4]] + j     *stag->entriesPerElementRow       + i     *stag->entriesPerEdge    + d; // note i increment by entriesPerEdge
          idxLocal[count++] =                                           jghost*stag->entriesPerElementRowGhost  + ighost*stag->entriesPerElement + d; // note i increment by entriesPerElement, still
        }
      }
    }

    // Neighbor 5 (right)
    if (!stag->lastproc[0]) {
    // We can never be right boundary, but the right neighbor may be
    // we may be a top boundary, along with the right neighbor
      const PetscInt neighborRank = neighbors[5];
      nNeighbor[0] = l[0][neighborRank % stag->nprocs[0]];
      nNeighbor[1] = l[1][(neighborRank/stag->nprocs[0])];
      const PetscInt entriesPerElementRowNeighbor = stag->proc[0] == stag->nprocs[0]-2 ?  nNeighbor[0]*stag->entriesPerElement + stag->entriesPerEdge : nNeighbor[0]*stag->entriesPerElement; // are we second-to-last in x?
      for (jghost = jghostoffset;jghost<stag->nghost[1]-jghostoffsettop; ++jghost) {
        const PetscInt j = jghost-jghostoffset;
        for (i=0; i<ighostoffsetright; ++i) {
          const PetscInt ighost = stag->nghost[0] - ighostoffsetright + i;
          for (d=0; d<stag->entriesPerElement; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
      }
      if (stag->lastproc[1]) {
        // assumes a single ghost element on the boundary in ecah column
        const PetscInt jghost = stag->nghost[1]-1; // careful with these
        const PetscInt j = stag->n[1]; // or nNeighbor[1];
        for (i=0; i<ighostoffsetright; ++i) { // as everywhere around here, we assume that this width is small enough that we don't hit the boundary (hence stencil width greater than 1 not supported, though it might work if you have at least that many elements per rank)
          const PetscInt ighost = stag->nghost[0] - ighostoffsetright + i;
          for (d=0; d<stag->entriesPerEdge; ++d) { // only vertices and horizontal edge (which are the first dof)
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerEdge    + d; // i moves by edge here
            idxLocal[count++] =                                     jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
      }
    }

    // Neighbor 6 (up left)
    if (!stag->firstproc[0] && !stag->lastproc[1]) {
      // We can never be a top boundary, but our neighbor may be
      // We may be a right boundary, but our neighbor cannot be
      const PetscInt neighborRank = neighbors[6];
      nNeighbor[0] = l[0][neighborRank % stag->nprocs[0]];
      nNeighbor[1] = l[1][(neighborRank/stag->nprocs[0])];
      const PetscInt entriesPerElementRowNeighbor =  nNeighbor[0]*stag->entriesPerElement; // the neighbor cannot be a right boundary
      for (j=0; j<jghostoffsettop; ++j) {
        const PetscInt jghost = stag->nghost[1] - jghostoffsettop + j;
        for (ighost = 0; ighost<ighostoffset; ++ighost) {
          const PetscInt i = nNeighbor[0] - ighostoffset + ighost;
          for (d=0; d<stag->entriesPerElement; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
      }
    }

    // Neighbor 7 (up) 
    if (!stag->lastproc[1]) {
      // We cannot be the last proc in y, though our neighbor may be
      // We may be the last proc in x, in which case our neighbor is also
      const PetscInt neighborRank = neighbors[7];
      nNeighbor[0] = l[0][neighborRank % stag->nprocs[0]];
      nNeighbor[1] = l[1][(neighborRank/stag->nprocs[0])];
      const PetscInt entriesPerElementRowNeighbor = stag->entriesPerElementRow; // rank above has same width and boundary status!
      for (j=0; j<jghostoffsettop; ++j) {
        const PetscInt jghost = stag->nghost[1] - jghostoffsettop + j;
        for (ighost = ighostoffset; ighost<stag->nghost[0]-ighostoffsetright; ++ighost) {
          const PetscInt i = ighost - ighostoffset;
          for (d=0; d<stag->entriesPerElement; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                                     jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
        if (stag->lastproc[0]) {
          // assumes a single ghost element on the boundaries
          const PetscInt ighost = stag->nghost[0]-1; // careful with these!!
          const PetscInt i = stag->n[0]; // or nNeighbor[0]
          // vertex
          for (d=0; d<stag->dof[0]; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
          // edge
          for (d=0; d<stag->dof[1]; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + stag->dof[0]                + d;
            idxLocal[count++] =                                     jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + stag->dof[0] + stag->dof[1] + d;
          }
        }
      }
    }

    // Neighbor 8 (up right)
    if (!stag->lastproc[0] && !stag->lastproc[1]) {
      // We can never be a ghosted boundary
      // Our neighbor may be a top boundary, a right boundary, or both, but we assume here that the stencil width isn't greater than the rank size, so this doesn't complicate things too much.
      const PetscInt neighborRank = neighbors[8];
      nNeighbor[0] = l[0][neighborRank % stag->nprocs[0]];
      nNeighbor[1] = l[1][(neighborRank/stag->nprocs[0])];
      const PetscInt entriesPerElementRowNeighbor = stag->proc[0] == stag->nprocs[0]-2 ?  nNeighbor[0]*stag->entriesPerElement + stag->entriesPerEdge : nNeighbor[0]*stag->entriesPerElement; // are we second-to-last in x?
      for (j=0; j<jghostoffsettop; ++j) {
        const PetscInt jghost = stag->nghost[1] - jghostoffsettop + j;
        for (i=0; i<ighostoffsetright; ++i) {
          const PetscInt ighost = stag->nghost[0] - ighostoffsetright + i;
          for (d=0; d<stag->entriesPerElement; ++d) {
            idxGlobal[count ] = globalOffsets[neighborRank] + j     *entriesPerElementRowNeighbor    + i     *stag->entriesPerElement + d;
            idxLocal[count++] =                               jghost*stag->entriesPerElementRowGhost + ighost*stag->entriesPerElement + d;
          }
        }
      }
    }

    // Create Local IS (transferring pointer ownership)
    ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm),entriesToTransferTotal,idxLocal,PETSC_OWN_POINTER,&isLocal);CHKERRQ(ierr);

    // ltog requires a mapping from [0,n) to global degrees of freedom. This means we need to be in local order, and include the dummy nodes (which we will map to global -1).
    // This is horrendously inefficient, TODO replace with an in-order local sweep of all local points (including dummies) as in da2.c
    // (we do this in betwen creating the ISs so we can use isLocal to look up, and still access idxGlobal before losing ownership)
    ierr = PetscMalloc1(stag->entriesGhost,&idxGlobalAll);CHKERRQ(ierr);
    for (i=0;i<stag->entriesGhost;++i) {
      PetscInt loc;
      ierr = ISLocate(isLocal,i,&loc);CHKERRQ(ierr); // !! probably slow! (could at least use a hash of some kind before implementing the more efficient option mentioned above)
      if (loc >= 0) {
        idxGlobalAll[i] = idxGlobal[loc];
      } else {
        idxGlobalAll[i] =  -1;
      }
    }

    // Create Global IS (transferring pointer ownership)
    ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm),entriesToTransferTotal,idxGlobal,PETSC_OWN_POINTER,&isGlobal);CHKERRQ(ierr);

    // Create stag->gtol (order arbitary (computed as PETSc ordering), doesn't include dummy entries))
    {
      Vec local,global;
      ierr = VecCreateMPIWithArray(PetscObjectComm((PetscObject)dm),1,stag->entries,PETSC_DECIDE,NULL,&global);CHKERRQ(ierr); // Here and everywhere else we use a blocksize of 1, but for local/ghosted vecs, maybe should have a block size of stag->entriesPerElement
      ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,stag->entriesGhost,NULL,&local);CHKERRQ(ierr);
      ierr = VecScatterCreate(global,isGlobal,local,isLocal,&stag->gtol);CHKERRQ(ierr);
      ierr = VecDestroy(&global);CHKERRQ(ierr);
      ierr = VecDestroy(&local);CHKERRQ(ierr);
    }

    // Destroy ISs
    ierr = ISDestroy(&isLocal);CHKERRQ(ierr);
    ierr = ISDestroy(&isGlobal);CHKERRQ(ierr);

    // Create local-to-global map (in local ordering, includes maps to -1 for dummy points)
    ierr = ISLocalToGlobalMappingCreate(comm,1,stag->entriesGhost,idxGlobalAll,PETSC_OWN_POINTER,&dm->ltogmap);CHKERRQ(ierr); // 1 dof
    ierr = PetscLogObjectParent((PetscObject)dm,(PetscObject)dm->ltogmap);CHKERRQ(ierr);
  }

  // Get rid of the local sizes array
  for (i=0;i<stag->dim;++i) {
    ierr = PetscFree(l[i]);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DMStagVecGetArray(DM dm,Vec vec,void *array)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;
  PetscInt       nLocal;

  PetscFunctionBegin;

  /* This function assumes that you want an n-d array where each entry is an elements-worth of data.
     Thus, it only works for "local" vectors, which include both types of ghost point. */

  ierr = VecGetLocalSize(vec,&nLocal);CHKERRQ(ierr);
  if (nLocal != stag->entriesGhost) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Vector local size %D is not compatible with DMStag local size %D\n",nLocal,stag->entriesGhost);

  switch (stag->dim) {
    case 2:
      ierr = VecGetArray2d(vec,stag->nghost[1],stag->nghost[0]*stag->entriesPerElement,stag->startGhost[1],stag->startGhost[0]*stag->entriesPerElement,(PetscScalar***)array);CHKERRQ(ierr);
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Not implemented for dim!=2");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagVecRestoreArray(DM dm,Vec vec,void *array)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;
  PetscInt       nLocal;

  PetscFunctionBegin;

  /* This function assumes that you want an n-d array where each entry is an elements-worth of data.
     Thus, it only works for "local" vectors, which include both types of ghost point. */

  ierr = VecGetLocalSize(vec,&nLocal);CHKERRQ(ierr);
  if (nLocal != stag->entriesGhost) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Vector local size %D is not compatible with DMStag local size %D\n",nLocal,stag->entriesGhost);

  switch (stag->dim) {
    case 2:
      ierr = VecRestoreArray2d(vec,stag->nghost[1],stag->nghost[0]*stag->entriesPerElement,stag->startGhost[1],stag->startGhost[0]*stag->entriesPerElement,(PetscScalar***)array);CHKERRQ(ierr);
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Not implemented for dim!=2");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagGetCorners(DM dm,PetscInt *x,PetscInt *y,PetscInt *z,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt *nDummyx,PetscInt *nDummyy,PetscInt *nDummyz)
{
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBegin;
  if (x) *x = stag->start[0];
  if (y) *y = stag->start[1];
  if (z) *z = stag->start[2];
  if (m) *m = stag->n[0];
  if (n) *n = stag->n[1];
  if (p) *p = stag->n[2];
  if (nDummyx) *nDummyx = stag->lastproc[0] ? 1 : 0;
  if (nDummyy) *nDummyy = stag->lastproc[1] ? 1 : 0;
  if (nDummyz) *nDummyz = stag->lastproc[2] ? 1 : 0;
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagGetGhostCorners(DM dm,PetscInt *x,PetscInt *y,PetscInt *z,PetscInt *m,PetscInt *n,PetscInt *p)
{
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBegin;
  if (x) *x = stag->startGhost[0];
  if (y) *y = stag->startGhost[1];
  if (z) *z = stag->startGhost[2];
  if (m) *m = stag->nghost[0];
  if (n) *n = stag->nghost[1];
  if (p) *p = stag->nghost[2];
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagGetGlobalSizes(DM dm,PetscInt* M,PetscInt* N,PetscInt* P)
{
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBegin;
  if (M) *M = stag->N[0];
  if (N) *N = stag->N[1];
  if (P) *P = stag->N[2];
  PetscFunctionReturn(0);
}

PetscErrorCode DMCreate_Stag(DM dm)
{
  PetscErrorCode ierr;
  DM_Stag        *stag;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidPointer(dm,1);
  ierr = PetscNewLog(dm,&stag);CHKERRQ(ierr);
  dm->data = stag;

  ierr = PetscObjectChangeTypeName((PetscObject)dm,DMSTAG);CHKERRQ(ierr);

  stag->gton = NULL;
  stag->gtol = NULL;
  stag->lton = NULL;
  stag->ntol = NULL;

  for (i=0;i<DMSTAG_MAX_STRATA;++i) stag->dof[i] = 0;

  dm->ops->destroy            = DMDestroy_Stag;
  dm->ops->createglobalvector = DMCreateGlobalVector_Stag;
  dm->ops->createlocalvector  = DMCreateLocalVector_Stag;
  dm->ops->creatematrix       = NULL;
  dm->ops->view               = NULL;
  dm->ops->load               = NULL;
  dm->ops->globaltolocalbegin = DMGlobalToLocalBegin_Stag;
  dm->ops->globaltolocalend   = DMGlobalToLocalEnd_Stag;
  dm->ops->localtoglobalbegin = DMLocalToGlobalBegin_Stag;
  dm->ops->localtoglobalend   = DMLocalToGlobalEnd_Stag;
  dm->ops->localtolocalbegin  = NULL;
  dm->ops->localtolocalend    = NULL;
  dm->ops->createsubdm        = NULL;
  dm->ops->createcoordinatedm = DMCreateCoordinateDM_Stag;
  dm->ops->clone              = NULL;
  dm->ops->setup              = DMSetUp_Stag;
  PetscFunctionReturn(0);
}
