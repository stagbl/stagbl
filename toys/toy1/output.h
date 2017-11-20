#ifndef OUTPUT_H__
#define OUTPUT_H__

#include <petscmat.h>
#include <petscdm.h>

PetscErrorCode OutputVecBinary(Vec,const char *,PetscBool);
PetscErrorCode OutputMatBinary(Mat,const char *);
PetscErrorCode OutputDMCoordsBinary(DM,const char*,PetscBool);
PetscErrorCode DMDA2dXDMFStart(DM,const char*,const char*,PetscViewer *);
PetscErrorCode DMDA2dXDMFAddAttribute(DM,const char*,const char*,PetscViewer);
PetscErrorCode DMDA2dXDMFFinish(PetscViewer *);

#endif
