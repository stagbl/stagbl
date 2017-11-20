#ifndef CTX_H__
#define CTX_H__

#include <petscdm.h>

// An Application Context
typedef struct CtxData {
  PetscBool isoviscous;
  PetscReal xmin,xmax,ymin,ymax;               // Domain Bounds
  DM        dmStokes;                          // Packer
  DM        daC,daCor,daX,daY;                 // Cells, Corners, X-vel edges, Y-vel edges
  Vec       rho,etaC,etaCor;                   // Coefficient Fields
  PetscReal rho1,rho2,eta1,eta2;               // For simple coefficient fields
  PetscReal gx,gy;                             // Gravity
  PetscReal etaCharacteristic;                 // For scaling
  PetscReal hxCharacteristic,hyCharacteristic; // For scaling
  PetscReal Kcont,Kbound;                      // Scaling Constants
  PetscInt  MC,NC,MX,NX,MY,NY,MCor,NCor;       // Sizes (redundant)
  PetscInt  offsetC,offsetX,offsetY;           // Offsets for naive assembly 
  PetscInt  pFixX,pFixY;                       // Fixes pressure node
} CtxData;
typedef CtxData* Ctx;

PetscErrorCode CreateCtx(Ctx *pctx);
PetscErrorCode DestroyCtx(Ctx *pctx);

#endif
