#include <petsc/private/dmimpl.h>
#include "dmstag.h"

PetscErrorCode DMStagCreateNaturalVector(DM dm,Vec *n)
{
  PetscErrorCode ierr;
  PetscInt       totalLocalEntries,entriesPerEdge,entriesPerElement,entriesPerCorner;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;

  ierr = VecCreate(PetscObjectComm((PetscObject)dm),n);CHKERRQ(ierr);
  if (stag->dim != 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"vector creation not implemented for dim!=2");
    entriesPerElement = stag->dof[0] + 2*stag->dof[1] + stag->dof[2];
    entriesPerEdge    = stag->dof[0] + stag->dof[1];
    entriesPerCorner  = stag->dof[0];
    totalLocalEntries = entriesPerElement * stag->n[0] * stag->n[1] 
      + (stag->lastproc[0]                      ? entriesPerEdge * stag->n[1] : 0) 
      + (stag->lastproc[1]                      ? entriesPerEdge * stag->n[0] : 0) 
      + (stag->lastproc[0] && stag->lastproc[1] ? entriesPerCorner            : 0);
    ierr = VecSetSizes(*n,totalLocalEntries,PETSC_DETERMINE);CHKERRQ(ierr); 
    // don't set block size, unlike DMDA (should we? can we?)

    ierr = VecSetType(*n,VECMPI);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagGlobalToNatural_Create(DM dm)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;
  PetscInt       *lidx;
  PetscInt       lict,i,j,Nlocal,sizeLocal,start;
  IS             from,to;
  PetscInt entriesPerElementRowNatural,entriesPerElement,entriesPerEdge,entriesPerCorner;
  Vec            globalDummy,naturalDummy;

  PetscFunctionBegin;
  if (stag->dim != 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"global to natural not implemented for dim != 2");

  entriesPerElement           = (stag->dof[0] + 2*stag->dof[1] + stag->dof[2]); // assumes that if the stratum is deactivaed, dof is 0
  entriesPerEdge              = stag->dof[0] + stag->dof[1];
  entriesPerCorner            = stag->dof[0];
  entriesPerElementRowNatural = stag->N[0] * entriesPerElement + stag->dof[0] + stag->dof[1];

  Nlocal = entriesPerElement * stag->n[0] * stag->n[1] 
    + (stag->lastproc[0]                      ? entriesPerEdge * stag->n[1] : 0) 
    + (stag->lastproc[1]                      ? entriesPerEdge * stag->n[0] : 0) 
    + (stag->lastproc[0] && stag->lastproc[1] ? entriesPerCorner            : 0);

  ierr = PetscMalloc1(Nlocal,&lidx);CHKERRQ(ierr);

  // Note that this is only for dimension 2!
  lict = 0;
  for (j=stag->start[1];j<stag->start[1]+stag->n[1];++j){
    for (i=stag->start[0];i<stag->start[0]+stag->n[0];++i){
      PetscInt d,offsetLocal;
      offsetLocal = 0;
      for(d=0; d<stag->dof[0]; ++d) { // vertices
        lidx[lict++] = j * entriesPerElementRowNatural + i * entriesPerElement + offsetLocal + d;
      }
      offsetLocal += stag->dof[0];
      for(d=0; d<2*stag->dof[1]; ++d) { // edges
        lidx[lict++] = j * entriesPerElementRowNatural + i * entriesPerElement + offsetLocal + d;
      }
      offsetLocal += 2*stag->dof[1];
      for(d=0; d<stag->dof[2]; ++d) { // elements
        lidx[lict++] = j * entriesPerElementRowNatural + i * entriesPerElement + offsetLocal + d;
      }
    }
    if (stag->lastproc[0]) {
      PetscInt d,offsetLocal;
      i = stag->N[0];
      offsetLocal = 0;
      for(d=0; d<stag->dof[0]; ++d) { // vertices
        lidx[lict++] = j * entriesPerElementRowNatural + i * entriesPerElement + offsetLocal + d;
      }
      offsetLocal += stag->dof[0];
      for(d=0; d<stag->dof[1]; ++d) { // edge
        lidx[lict++] = j * entriesPerElementRowNatural + i * entriesPerElement + offsetLocal + d;
      }
    }
  }
  if (stag->lastproc[1]) {
    j = stag->N[1];
    for (i=stag->start[0];i<stag->start[0]+stag->n[0];++i) {
      PetscInt d,offsetLocal;
      offsetLocal = 0;
      for(d=0; d<stag->dof[0]; ++d) { // vertices
        lidx[lict++] = j * entriesPerElementRowNatural + i * entriesPerEdge + offsetLocal + d;
      }
      offsetLocal += stag->dof[0];
      for(d=0; d<stag->dof[1]; ++d) { // edge
        lidx[lict++] = j * entriesPerElementRowNatural + i * entriesPerEdge + offsetLocal + d;
      }
    }
    if (stag->lastproc[0]) {
      i = stag->N[0];
      PetscInt d,offsetLocal;
      offsetLocal=0;
      for(d=0; d<stag->dof[0]; ++d) { // vertex
        lidx[lict++] = j * entriesPerElementRowNatural + i * entriesPerEdge + offsetLocal + d;
      }
    }
  }
  ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm),Nlocal,lidx,PETSC_OWN_POINTER,&to);CHKERRQ(ierr);

  // Create the from IS and scatter with the help of two dummy vectors (as everywhere here, bs=1)
  ierr = VecCreateMPIWithArray(PetscObjectComm((PetscObject)dm),1,stag->entries,PETSC_DETERMINE,NULL,&globalDummy);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(PetscObjectComm((PetscObject)dm),1,stag->entries,PETSC_DETERMINE,NULL,&naturalDummy);CHKERRQ(ierr);
  ierr = VecGetLocalSize(globalDummy,&sizeLocal);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(globalDummy,&start,NULL);CHKERRQ(ierr);
  ierr = ISCreateStride(PetscObjectComm((PetscObject)dm),sizeLocal,start,1,&from);CHKERRQ(ierr);

  // debug
#if 0 
  ierr = ISView(to,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = ISView(from,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

  ierr = VecScatterCreate(globalDummy,from,naturalDummy,to,&stag->gton);CHKERRQ(ierr); 
  ierr = VecDestroy(&globalDummy);CHKERRQ(ierr);
  ierr = VecDestroy(&naturalDummy);CHKERRQ(ierr);
  ierr = ISDestroy(&to);CHKERRQ(ierr);
  ierr = ISDestroy(&from);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMStagGlobalToNaturalBegin(DM dm,Vec globalVec,InsertMode mode, Vec naturalVec)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(globalVec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(naturalVec,VEC_CLASSID,4);
  if (!stag->gton) {
    ierr = DMStagGlobalToNatural_Create(dm);CHKERRQ(ierr);
  }
  ierr = VecScatterBegin(stag->gton,globalVec,naturalVec,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode  DMStagGlobalToNaturalEnd(DM dm,Vec globalVec,InsertMode mode,Vec naturalVec)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(naturalVec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(globalVec,VEC_CLASSID,4);
  ierr = VecScatterEnd(stag->gton,globalVec,naturalVec,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagNaturalToGlobalBegin(DM dm,Vec naturalVec,InsertMode mode, Vec globalVec)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(naturalVec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(globalVec,VEC_CLASSID,4);
  if (!stag->gton) {
    ierr = DMStagGlobalToNatural_Create(dm);CHKERRQ(ierr);
  }
  ierr = VecScatterBegin(stag->gton,naturalVec,globalVec,mode,SCATTER_REVERSE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagNaturalToGlobalEnd(DM dm,Vec naturalVec,InsertMode mode, Vec globalVec)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(naturalVec,VEC_CLASSID,2);
  PetscValidHeaderSpecific(globalVec,VEC_CLASSID,4);
  if (!stag->gton) {
    ierr = DMStagGlobalToNatural_Create(dm);CHKERRQ(ierr);
  }
  ierr = VecScatterEnd(stag->gton,naturalVec,globalVec,mode,SCATTER_REVERSE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagLocalToNatural_Create(DM dm)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;
  IS             to,from;
  PetscInt       *idxTo,*idxFrom;
  PetscInt       ighostoffset,jghostoffset,count,i,j;

  ighostoffset = stag->start[0] - stag->startGhost[0];
  jghostoffset = stag->start[1] - stag->startGhost[1];

  // allocate idxTo and idxFrom of the correct size
  ierr = PetscMalloc1(stag->entries,&idxTo);CHKERRQ(ierr);
  ierr = PetscMalloc1(stag->entries,&idxFrom);CHKERRQ(ierr);

  // iterate over interior points to collect global indices to map to in "to", and local indices to map "from"
  count = 0;
  for (j=0;j<stag->n[1];++j) {
    const PetscInt jghost  = j + jghostoffset;
    const PetscInt jglobal = j + stag->start[1];
    for (i=0;i<stag->n[0];++i) {
      const PetscInt ighost  = i + ighostoffset;
      const PetscInt iglobal = i + stag->start[0];
      PetscInt d;
      for (d=0;d<stag->entriesPerElement;++d) {
        idxFrom[count  ] = jghost *stag->entriesPerElementRowGhost  + ighost *stag->entriesPerElement + d;
        idxTo  [count++] = jglobal*stag->entriesPerElementRowGlobal + iglobal*stag->entriesPerElement + d;
      }
    }
    if (stag->lastproc[0]) {
      i = stag->n[0];
      const PetscInt ighost  = i + ighostoffset;
      const PetscInt iglobal = i + stag->start[0];
      PetscInt d,dGhost;
      for (d=0; d<stag->dof[0]; ++d) {
        idxFrom[count  ] = jghost *stag->entriesPerElementRowGhost  + ighost *stag->entriesPerElement + d;
        idxTo  [count++] = jglobal*stag->entriesPerElementRowGlobal + iglobal*stag->entriesPerElement + d;
      }
      // skip bottom edge
      for(dGhost = stag->dof[0]+stag->dof[1];dGhost < stag->dof[0] + 2*stag->dof[1]; ++dGhost) {
        const PetscInt dGlobal = dGhost-stag->dof[1];
        idxFrom[count  ] = jghost *stag->entriesPerElementRowGhost  + ighost *stag->entriesPerElement + dGhost;
        idxTo  [count++] = jglobal*stag->entriesPerElementRowGlobal + iglobal*stag->entriesPerElement + dGlobal;
      } 
    }
  }
  if (stag->lastproc[1]) {
    j = stag->n[1];
    const PetscInt jghost  = j + jghostoffset;
    const PetscInt jglobal = j + stag->start[1];
    for (i=0;i<stag->n[0];++i) {
      const PetscInt ighost  = i + ighostoffset;
      const PetscInt iglobal = i + stag->start[0];
      PetscInt d;
      for (d=0;d<stag->dof[0]+stag->dof[1];++d) {
        idxFrom[count  ] = jghost *stag->entriesPerElementRowGhost + ighost *stag->entriesPerElement + d; // note top ghost row has full elements
        idxTo  [count++] = jglobal*stag->entriesPerElementRowGlobal+ iglobal*stag->entriesPerEdge    + d; // note top global row has only edges
      }
    }
    if (stag->lastproc[0]) {
      i = stag->n[0];
      const PetscInt ighost  = i + ighostoffset;
      const PetscInt iglobal = i + stag->start[0];
      PetscInt d;
      for (d=0;d<stag->dof[0];++d) {
        idxFrom[count  ] = jghost *stag->entriesPerElementRowGhost + ighost *stag->entriesPerElement + d; // note top ghost row has full elements
        idxTo  [count++] = jglobal*stag->entriesPerElementRowGlobal+ iglobal*stag->entriesPerEdge    + d; // note top global row has only edges
      }
    }
  }

  if (count != stag->entries) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_INCOMP,"Count is wrong generating local->global scatter");

  // Create GeneralISs for to and from (passing ownership of pointers)
  ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm),stag->entries,idxFrom,PETSC_OWN_POINTER,&from);CHKERRQ(ierr);
  ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm),stag->entries,idxTo  ,PETSC_OWN_POINTER,&to  );CHKERRQ(ierr);

#if 0 
  ierr = PetscPrintf(PETSC_COMM_WORLD," >> Scatter From\n");
  ierr = ISView(from,PETSC_VIEWER_STDOUT_WORLD);
  ierr = PetscPrintf(PETSC_COMM_WORLD," >> Scatter To\n");
  ierr = ISView(to,PETSC_VIEWER_STDOUT_WORLD);
#endif

  // Create Vec scatter and set in stag->lton
  {
    Vec local,natural;
    ierr = VecCreateMPIWithArray(PetscObjectComm((PetscObject)dm),1,stag->entries,PETSC_DECIDE,NULL,&natural);CHKERRQ(ierr); // Here and everywhere else we use a blocksize of 1, but for local/ghosted vecs, maybe should have a block size of stag->entriesPerElement
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,stag->entriesGhost,NULL,&local);CHKERRQ(ierr);
    ierr = VecScatterCreate(local,from,natural,to,&stag->lton);CHKERRQ(ierr);
    ierr = VecDestroy(&natural);CHKERRQ(ierr);
    ierr = VecDestroy(&local);CHKERRQ(ierr);
  }

  // Destroy ISs
  ierr = ISDestroy(&to);CHKERRQ(ierr);
  ierr = ISDestroy(&from);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMStagLocalToNaturalBegin(DM dm,Vec l,InsertMode mode,Vec n)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(l,VEC_CLASSID,2);
  PetscValidHeaderSpecific(n,VEC_CLASSID,4);
  if (!stag->lton) {
    ierr = DMStagLocalToNatural_Create(dm);CHKERRQ(ierr);
  }
  ierr = VecScatterBegin(stag->lton,l,n,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagLocalToNaturalEnd(DM dm,Vec l,InsertMode mode,Vec n)
{
  
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(l,VEC_CLASSID,2);
  PetscValidHeaderSpecific(n,VEC_CLASSID,4);
  ierr = VecScatterEnd(stag->lton,l,n,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagNaturalToLocal_Create(DM dm)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;
  IS       to,from;
  PetscInt *idxTo,*idxFrom;
  PetscInt count,entriesToTransferTotal,i,j;

  PetscFunctionBegin;
  entriesToTransferTotal  = stag->nghost[1]*stag->nghost[0]*stag->entriesPerElement // overestimate
    + (stag->lastproc[0]                      ? stag->nghost[1] * (-stag->entriesPerElement + stag->entriesPerEdge  ) : 0)
    + (stag->lastproc[1]                      ? stag->nghost[0] * (-stag->entriesPerElement + stag->entriesPerEdge  ) : 0)
    + (stag->lastproc[0] && stag->lastproc[1] ? stag->dof[2]  : 0); // replace double-counted element on corner (edges were each removed once)

  // allocate idxTo and idxFrom of the correct size
  ierr = PetscMalloc1(entriesToTransferTotal,&idxTo);CHKERRQ(ierr);
  ierr = PetscMalloc1(entriesToTransferTotal,&idxFrom);CHKERRQ(ierr);

  // Populate index arrays
  count = 0;
  {
    const PetscInt jmax = stag->nghost[1] + (stag->lastproc[1] ? -1 : 0);
    const PetscInt imax = stag->nghost[0] + (stag->lastproc[0] ? -1 : 0);
    PetscInt d;
    for (j=0;j<jmax;++j) {
      const PetscInt jglobal = j + stag->startGhost[1];
      for (i=0;i<imax;++i) {
        const PetscInt iglobal = i + stag->startGhost[0];
        for (d=0;d<stag->entriesPerElement;++d) {
          idxFrom[count  ] = jglobal*stag->entriesPerElementRowGlobal + iglobal*stag->entriesPerElement + d;
          idxTo  [count++] = j      *stag->entriesPerElementRowGhost  + i      *stag->entriesPerElement + d;
        } 
      }
      if (stag->lastproc[0]) {
        i = stag->nghost[0]-1;
        const PetscInt iglobal = i + stag->startGhost[0];
        for (d=0;d<stag->dof[0];++d) {
          idxFrom[count  ] = jglobal*stag->entriesPerElementRowGlobal + iglobal*stag->entriesPerElement + d;
          idxTo  [count++] = j      *stag->entriesPerElementRowGhost  + i      *stag->entriesPerElement + d;
        }
        for (d=stag->dof[0];d<stag->dof[0]+stag->dof[1];++d) {
          const PetscInt dGhost = d + stag->dof[1]; // skip bottom edge
          idxFrom[count  ] = jglobal*stag->entriesPerElementRowGlobal + iglobal*stag->entriesPerElement + d;
          idxTo  [count++] = j      *stag->entriesPerElementRowGhost  + i      *stag->entriesPerElement + dGhost;
        }
      }
    }
    if (stag->lastproc[1]) {
      j = stag->nghost[1]-1;
      const PetscInt jglobal = j + stag->startGhost[1];
      for (i=0;i<imax;++i) {
        const PetscInt iglobal = i + stag->startGhost[0];
        for (d=0;d<stag->dof[0]+stag->dof[1];++d) { // vertex and bottom edge
          idxFrom[count  ] = jglobal*stag->entriesPerElementRowGlobal + iglobal*stag->entriesPerEdge    + d; // note i moves by edge in global
          idxTo  [count++] = j      *stag->entriesPerElementRowGhost  + i      *stag->entriesPerElement + d; // note i still moves by element in ghost
        }
      }
      if (stag->lastproc[0]) {
        i = stag->nghost[0]-1;
        const PetscInt iglobal = i + stag->startGhost[0];
        for (d=0;d<stag->dof[0];++d) { // vertex
          idxFrom[count  ] = jglobal*stag->entriesPerElementRowGlobal + iglobal*stag->entriesPerEdge    + d; // note i moves by edge in global
          idxTo  [count++] = j      *stag->entriesPerElementRowGhost  + i      *stag->entriesPerElement + d; // note i still moves by element in ghost
        }
      }
    }
  }

  // Create GeneralISs for to and from (passing ownership of pointers)
  ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm),entriesToTransferTotal,idxFrom,PETSC_OWN_POINTER,&from);CHKERRQ(ierr);
  ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm),entriesToTransferTotal,idxTo  ,PETSC_OWN_POINTER,&to  );CHKERRQ(ierr);

#if 0
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Natural-->Local ISs\n");
  ierr = PetscPrintf(PETSC_COMM_WORLD," >> Scatter From\n");
  ierr = ISView(from,PETSC_VIEWER_STDOUT_WORLD);
  ierr = PetscPrintf(PETSC_COMM_WORLD," >> Scatter To\n");
  ierr = ISView(to,PETSC_VIEWER_STDOUT_WORLD);
#endif

  // Create vec scatter and set in stag->ntol
  {
    Vec local,global;
    ierr = VecCreateMPIWithArray(PetscObjectComm((PetscObject)dm),1,stag->entries,PETSC_DECIDE,NULL,&global);CHKERRQ(ierr); // Here and everywhere else we use a blocksize of 1, but for local/ghosted vecs, maybe should have a block size of stag->entriesPerElement

    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,stag->entriesGhost,NULL,&local);CHKERRQ(ierr);
    ierr = VecScatterCreate(global,from,local,to,&stag->ntol);CHKERRQ(ierr);

    ierr = VecDestroy(&global);CHKERRQ(ierr);
    ierr = VecDestroy(&local);CHKERRQ(ierr);

#if 0
    ierr = PetscPrintf(comm,"Natural-to-Local Scatter:\n");CHKERRQ(ierr);
    ierr = VecScatterView(stag->ntol,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
  }

  // Destroy ISs
  ierr = ISDestroy(&to);CHKERRQ(ierr);
  ierr = ISDestroy(&from);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DMStagNaturalToLocalBegin(DM dm,Vec n,InsertMode mode,Vec l)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(n,VEC_CLASSID,2);
  PetscValidHeaderSpecific(l,VEC_CLASSID,4);
  if (!stag->ntol) {
    ierr = DMStagNaturalToLocal_Create(dm);CHKERRQ(ierr);
  }
  ierr = VecScatterBegin(stag->ntol,n,l,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  // Note: it's redundant to have both ntol and lton. Could just do SCATTER_REVERSE..
  PetscFunctionReturn(0);
}

PetscErrorCode DMStagNaturalToLocalEnd(DM dm,Vec n,InsertMode mode,Vec l)
{
  PetscErrorCode ierr;
  DM_Stag        *stag = (DM_Stag*)dm->data;

  PetscFunctionBeginUser;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(n,VEC_CLASSID,2);
  PetscValidHeaderSpecific(l,VEC_CLASSID,4);
  ierr = VecScatterEnd(stag->ntol,n,l,mode,SCATTER_FORWARD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
