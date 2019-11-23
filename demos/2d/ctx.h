#ifndef CTX_H_
#define CTX_H_

#include "stagbl.h"

typedef struct {
  MPI_Comm    comm;
  char        mode[1024];
  StagBLGrid  stokesGrid,coeffGrid;
  StagBLArray coeffArray;
  PetscScalar  xmax,ymax,xmin,ymin,hxCharacteristic,hyCharacteristic;
  PetscScalar  eta1,eta2,rho1,rho2,gy,Kbound,Kcont,etaCharacteristic;
  PetscInt   pinx,piny;
  PetscScalar (*getEta)(void*,PetscScalar,PetscScalar);
  PetscScalar (*getRho)(void*,PetscScalar,PetscScalar);
} CtxData;
typedef CtxData* Ctx;

PetscErrorCode CtxCreate(MPI_Comm,const char* mode,Ctx*);
PetscErrorCode CtxSetupFromGrid(Ctx);

#endif
