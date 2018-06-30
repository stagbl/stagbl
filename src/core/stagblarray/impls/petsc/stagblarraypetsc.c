#include "../../stagblarrayimpl.h"
#include "stagblarraypetscimpl.h"
#include <stdlib.h>

void StagBLArrayDestroy_PETSc(StagBLArray stagblarray)
{
  StagBLArray_PETSc *data= (StagBLArray_PETSc*) stagblarray->data;
  if (data->vec) {
    VecDestroy(&data->vec);
  }
  free(stagblarray->data);
  stagblarray->data = NULL;
}

void StagBLArrayCreate_PETSc(StagBLArray stagblarray)
{
  StagBLArray_PETSc *data;
  stagblarray->data = (void*) malloc(sizeof(StagBLArray_PETSc));
  data = (StagBLArray_PETSc*) stagblarray->data;
  data->vec = NULL;
  stagblarray->ops->destroy = StagBLArrayDestroy_PETSc;
}

void StagBLArrayPETScGetVecPointer(StagBLArray stagblarray,Vec **vec)
{
  StagBLArray_PETSc * const data = (StagBLArray_PETSc*) stagblarray->data;
  *vec = &(data->vec);
}
