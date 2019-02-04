#ifndef CTX_H_
#define CTX_H_

#include "stagbl.h"

typedef struct {
  MPI_Comm    comm;
  StagBLGrid  stokesGrid,coeffGrid;
  StagBLArray coeffArray;
  StagBLReal  xmax,ymax,xmin,ymin,hxCharacteristic,hyCharacteristic;
  StagBLReal  eta1,eta2,rho1,rho2,gy,Kbound,Kcont,etaCharacteristic;
  StagBLInt   pinx,piny;
  StagBLReal (*getEta)(void*,StagBLReal,StagBLReal);
  StagBLReal (*getRho)(void*,StagBLReal,StagBLReal);
} CtxData;
typedef CtxData* Ctx;

StagBLErrorCode CtxCreate(MPI_Comm,Ctx*);
StagBLErrorCode CtxSetupFromGrid(Ctx);

#endif
