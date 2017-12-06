#ifndef CTX_H__
#define CTX_H__

#include <petscdm.h>

// An Application Context
typedef struct CtxData {
  DM stokesGrid,paramGrid;
  Vec param;
  PetscBool isoviscous;
  PetscInt  M,N; // global numbers of cells
  PetscReal xmin,ymin,xmax,ymax; // domain size
  PetscReal rho1,rho2,eta1,eta2; // for simple coefficients
  PetscReal etaCharacteristic;                 // For scaling
  PetscReal hxCharacteristic,hyCharacteristic; // For scaling
  PetscReal gx,gy;
  PetscInt pFixx,pFixy;
  PetscReal Kbound,Kcont;
} CtxData;
typedef CtxData* Ctx;

PetscErrorCode CreateCtx(Ctx *pctx);
PetscErrorCode DestroyCtx(Ctx *pctx);

#endif
