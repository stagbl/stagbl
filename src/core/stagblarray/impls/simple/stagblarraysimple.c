#include "stagbl/private/stagblarrayimpl.h"
#include "stagblarraysimpleimpl.h"
#include <stdlib.h>

static PetscErrorCode StagBLArraySimpleCreateGlobalVector_Private(StagBLArray);
static PetscErrorCode StagBLArraySimpleCreateLocalVector_Private(StagBLArray);

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

static PetscErrorCode StagBLArraySimpleCreateLocalVector_Private(StagBLArray array)
{
  PetscErrorCode    ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;

  PetscFunctionBegin;
  if (!data->local) {
    DM       dm;
    PetscInt local_size;

    ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
    ierr = DMStagGetEntriesLocal(dm,&local_size);CHKERRQ(ierr);
    data->local = (PetscScalar*) malloc(local_size * sizeof(PetscScalar));
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLArraySimpleCreateGlobalVector_Private(StagBLArray array)
{
  PetscErrorCode    ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;

  PetscFunctionBegin;
  if (!data->global) {
    DM       dm;
    PetscInt native_size;

    ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
    ierr = DMStagGetEntries(dm,&native_size);CHKERRQ(ierr);
    data->global = (PetscScalar*) malloc(native_size * sizeof(PetscScalar));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayPrint_Simple(StagBLArray array)
{
  PetscErrorCode    ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;

  PetscFunctionBegin;
  if (data->local) {
    DM          dm;
    PetscInt    local_size;
    PetscMPIInt rank;

    ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
    ierr = DMStagGetEntriesLocal(dm,&local_size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm),&rank);CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PetscObjectComm((PetscObject)dm),"[%d] local portion of global entries:\n",rank);CHKERRQ(ierr);
    for (PetscInt i=0; i<local_size; ++i) {
      ierr = PetscSynchronizedPrintf(PetscObjectComm((PetscObject)dm),"%g\n", (double)data->local[i]);CHKERRQ(ierr);
    }
    ierr = PetscSynchronizedFlush(PetscObjectComm((PetscObject)dm),PETSC_STDOUT);CHKERRQ(ierr);
  }
  if (data->global) {
    DM          dm;
    PetscInt    native_size;
    PetscMPIInt rank;

    ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
    ierr = DMStagGetEntries(dm,&native_size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm),&rank);CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PetscObjectComm((PetscObject)dm),"[%d] local portion of global entries:\n",rank);CHKERRQ(ierr);
    for (PetscInt i=0; i<native_size; ++i) {
      ierr = PetscSynchronizedPrintf(PetscObjectComm((PetscObject)dm),"%g\n", (double)data->global[i]);CHKERRQ(ierr);
    }
    ierr = PetscSynchronizedFlush(PetscObjectComm((PetscObject)dm),PETSC_STDOUT);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArraySetLocalConstant_Simple(StagBLArray array, PetscScalar value)
{
  PetscErrorCode    ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;

  PetscFunctionBegin;
  if (!data->local) {
    ierr = StagBLArraySimpleCreateLocalVector_Private(array);CHKERRQ(ierr);
  }
  {
    DM       dm;
    PetscInt local_size;

    ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
    ierr = DMStagGetEntriesLocal(dm,&local_size);CHKERRQ(ierr);
    for (PetscInt i=0; i<local_size; ++i) data->local[i] = value;
  }
  array->current_local = PETSC_TRUE;
  array->current_global = PETSC_FALSE;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArrayCreate_Simple(StagBLArray stagblarray)
{
  StagBLArray_Simple *data;

  stagblarray->data = (void*) malloc(sizeof(StagBLArray_Simple));
  data = (StagBLArray_Simple*) stagblarray->data;
  data->local = NULL;
  data->global = NULL;
  stagblarray->ops->destroy = StagBLArrayDestroy_Simple;
  stagblarray->ops->print = StagBLArrayPrint_Simple;
  stagblarray->ops->setlocalconstant = StagBLArraySetLocalConstant_Simple;
  return 0;
}
