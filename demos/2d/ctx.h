#ifndef CTX_H_
#define CTX_H_

#include <petscdm.h>

typedef struct {
  MPI_Comm    comm;
  DM          dmStokes,dmCoeff;
  Vec         coeff;
  PetscReal   xmax,ymax,xmin,ymin,hxCharacteristic,hyCharacteristic;
  PetscScalar eta1,eta2,rho1,rho2,gy,Kbound,Kcont,etaCharacteristic;
  PetscInt    pinx,piny;
} CtxData;
typedef CtxData* Ctx;

#endif
