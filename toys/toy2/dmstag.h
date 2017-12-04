#if !defined(DMSTAG_H__)
#define DMSTAG_H__

#include <petscdm.h>

/* Note: we define things internally based on this, but the interfaces assume it's 
         equal to 3 */
#define DMSTAG_MAX_DIM 3
#define DMSTAG_MAX_STRATA DMSTAG_MAX_DIM + 1

#define DMSTAG "stag"

typedef enum{DMSTAG_GHOST_STENCIL_NONE=0,DMSTAG_GHOST_STENCIL_STAR,DMSTAG_GHOST_STENCIL_BOX} DMStagGhostStencilType;

typedef struct {
  PetscInt               dim;
  PetscInt               N[DMSTAG_MAX_DIM];                /* Global dimensions (elements)    */
  PetscInt               n[DMSTAG_MAX_DIM];                /* Local dimensions (elements)     */
  PetscInt               nghost[DMSTAG_MAX_DIM];           /* Local dimensions (with ghosts)  */
  PetscInt               start[DMSTAG_MAX_DIM];            /* First element number            */
  PetscInt               startGhost[DMSTAG_MAX_DIM];       /* First element number (ghosted)  */
  PetscMPIInt            proc[DMSTAG_MAX_DIM];             /* Location in processor grid      */
  PetscBool              lastproc[DMSTAG_MAX_DIM];         /* Last proc in this dim?          */ 
  PetscBool              firstproc[DMSTAG_MAX_DIM];        /* First proc in this dim?         */ 
  PetscMPIInt            nprocs[DMSTAG_MAX_DIM];           /* Procs in each direction         */
  PetscInt               dof[DMSTAG_MAX_STRATA];           /* dof per point for each stratum  */
  PetscBool              stratumActive[DMSTAG_MAX_STRATA]; 
  DMStagGhostStencilType ghostStencil;                     /* element-wise ghost stencil      */
  PetscInt               ghostStencilWidth;                /* elementwise ghost width         */
  DMBoundaryType         boundaryType[DMSTAG_MAX_DIM];
  VecScatter             gton;                             /* Global  --> Natural             */
  VecScatter             gtol;                             /* Global  --> Local               */
  VecScatter             ntol;                             /* Natural --> Local               */
  VecScatter             lton;                             /* Local   --> Natural (redundant) */
  PetscInt entriesPerElement,entriesPerEdge,entriesPerCorner,entriesPerElementRow,
           entries,entriesGhost,entriesPerElementRowGhost,entriesPerElementRowGlobal;
} DM_Stag;

PetscErrorCode DMStagCreate2d(MPI_Comm,DMBoundaryType,DMBoundaryType,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,DMStagGhostStencilType,PetscInt,DM*);
PetscErrorCode DMCreate_Stag(DM);
PetscErrorCode DMStagCreateNaturalVector(DM,Vec*);
PetscErrorCode DMStagGlobalToNaturalBegin(DM,Vec,InsertMode,Vec);
PetscErrorCode DMStagGlobalToNaturalEnd(DM,Vec,InsertMode,Vec);
PetscErrorCode DMStagNaturalToGlobalBegin(DM,Vec,InsertMode,Vec);
PetscErrorCode DMStagNaturalToGlobalEnd(DM,Vec,InsertMode,Vec);
PetscErrorCode DMStagSetUniformCoordinates(DM,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal);
PetscErrorCode DMStagVecGetArray(DM,Vec,void*);
PetscErrorCode DMStagVecRestoreArray(DM,Vec,void*);
PetscErrorCode DMStagGetCorners(DM,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PetscErrorCode DMStagGetGhostCorners(DM,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PetscErrorCode DMStagGetGlobalSizes(DM,PetscInt*,PetscInt*,PetscInt*);
PetscErrorCode DMStagNaturalToLocal_Create(DM);
PetscErrorCode DMStagNaturalToLocalBegin(DM,Vec,InsertMode,Vec);
PetscErrorCode DMStagNaturalToLocalEnd(DM,Vec,InsertMode,Vec);
PetscErrorCode DMStagLocalToNatural_Create(DM);
PetscErrorCode DMStagLocalToNaturalBegin(DM,Vec,InsertMode,Vec);
PetscErrorCode DMStagLocalToNaturalEnd(DM,Vec,InsertMode,Vec);

#endif
