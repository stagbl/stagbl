#include "stagbl/private/stagblarrayimpl.h"
#include "stagblarraysimpleimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLArrayDestroy_Simple(StagBLArray stagblarray)
{
  StagBLArray_Simple *data= (StagBLArray_Simple*) stagblarray->data;

  if (data->local) {
    free(data->local);
  }
  if (data->global) {
    free(data->global);
  }
  free(stagblarray->data);
  stagblarray->data = NULL;
  return 0;
}

PetscErrorCode StagBLArrayCreate_Simple(StagBLArray stagblarray)
{
  StagBLArray_Simple *data;

  stagblarray->data = (void*) malloc(sizeof(StagBLArray_Simple));
  data = (StagBLArray_Simple*) stagblarray->data;
  data->local = NULL;
  data->global = NULL;
  stagblarray->ops->destroy = StagBLArrayDestroy_Simple;
  return 0;
}
