#ifndef OUTPUT_H__
#define OUTPUT_H__

#include <petscmat.h>
#include <petscdm.h>

PetscErrorCode OutputVecBinary(Vec,const char *,PetscBool);
PetscErrorCode OutputMatBinary(Mat,const char *);
PetscErrorCode OutputDMCoordsNaturalBinary(DM,const char*,PetscBool);
PetscErrorCode DMStag2dXDMFStart(DM,const char*,const char*,PetscViewer *);
PetscErrorCode DMStag2dXDMFAddAttribute(DM,const char*,const char*,PetscViewer);
PetscErrorCode DMStag2dXDMFFinish(PetscViewer *);

#endif
