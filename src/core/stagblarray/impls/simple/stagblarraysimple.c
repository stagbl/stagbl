#include "stagbl/private/stagblarrayimpl.h"
#include "stagblarraysimpleimpl.h"
#include <stdlib.h>

/* For DMStagStencilToIndexLocal(), which should be made public */
#include <petsc/private/dmstagimpl.h>

static PetscErrorCode StagBLArraySimpleCreateGlobalVector_Private(StagBLArray);
static PetscErrorCode StagBLArraySimpleCreateLocalVector_Private(StagBLArray);

static PetscErrorCode StagBLArrayDestroy_Simple(StagBLArray stagblarray)
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

static PetscErrorCode StagBLArrayGetLocalValuesStencil_Simple(StagBLArray array,PetscInt n,const DMStagStencil *pos,PetscScalar *values)
{
  PetscErrorCode    ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;
  DM                dm;
  PetscInt          *indices_local;
  PetscInt           dim;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  if (!data->local) StagBLError(PetscObjectComm((PetscObject)dm),"Local array data not defined");
  if (!array->current_local) StagBLError(PetscObjectComm((PetscObject)dm),"Local array data not current");
  indices_local = (PetscInt*) malloc(n * sizeof(PetscInt));
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  ierr = DMStagStencilToIndexLocal(dm,dim,n,pos,indices_local);CHKERRQ(ierr);
  for (PetscInt i=0; i<n; ++i) values[i] = data->local[indices_local[i]];
  free(indices_local);
  PetscFunctionReturn(0);
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

static PetscErrorCode StagBLArraySetLocalValuesStencil_Simple(StagBLArray array,PetscInt n,const DMStagStencil *pos,const PetscScalar *values)
{
  PetscErrorCode    ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;
  DM                dm;
  PetscInt          *indices_local;
  PetscInt           dim;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  if (!data->local) StagBLError(PetscObjectComm((PetscObject)dm),"Local array data not defined");
  if (!array->current_local) StagBLError(PetscObjectComm((PetscObject)dm),"Local array data not current");
  indices_local = (PetscInt*) malloc(n * sizeof(PetscInt));
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  ierr = DMStagStencilToIndexLocal(dm,dim,n,pos,indices_local);CHKERRQ(ierr);
  for (PetscInt i=0; i<n; ++i) data->local[indices_local[i]] = values[i];
  free(indices_local);
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

static PetscErrorCode StagBLArraySimpleCreateShellGlobalVec_Private(StagBLArray array, Vec *p_vec)
{
  PetscErrorCode     ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;
  DM                 dm;
  PetscInt           native_size;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  ierr = DMStagGetEntries(dm,&native_size);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(PetscObjectComm((PetscObject)dm),1,native_size,PETSC_DECIDE,data->global,p_vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLArraySimpleCreateShellLocalVec_Private(StagBLArray array, Vec *p_vec)
{
  PetscErrorCode     ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;
  DM                 dm;
  PetscInt           local_size;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  ierr = DMStagGetEntriesLocal(dm,&local_size);CHKERRQ(ierr);
  ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,local_size,data->local,p_vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Note that this function proceeds by wrapping the raw pointers in PETSc Vec
   objects. If one wanted to do things "manually", one can get access to the
   actual local-to-global map with functions like
      - DMGetLocalToGlobalMapping()
      - ISLocalToGlobalMappingGetIndices()
 */
static PetscErrorCode StagBLArrayGlobalToLocal_Simple(StagBLArray array)
{
  PetscErrorCode     ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;
  DM                 dm;
  Vec                shell_local, shell_global;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  if (!data->global) StagBLError(PetscObjectComm((PetscObject)dm),"Cannot perform global-to-local before global data is defined");
  if (!array->current_global) StagBLError(PetscObjectComm((PetscObject)dm),"Refusing to perform global-to-local with stale global data");
  if (!data->local) {
    ierr = StagBLArraySimpleCreateLocalVector_Private(array);CHKERRQ(ierr);
  }
  ierr = StagBLArraySimpleCreateShellGlobalVec_Private(array, &shell_global);CHKERRQ(ierr);
  ierr = StagBLArraySimpleCreateShellLocalVec_Private(array, &shell_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocal(dm,shell_global,INSERT_VALUES,shell_local);CHKERRQ(ierr);
  ierr = VecDestroy(&shell_global);CHKERRQ(ierr);
  ierr = VecDestroy(&shell_local);CHKERRQ(ierr);
  array->current_local = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/* Note that this function proceeds by wrapping the raw pointers in PETSc Vec
   objects. If one wanted to do things "manually", one can get access to the
   actual local-to-global map with functions like
      - DMGetLocalToGlobalMapping()
      - ISLocalToGlobalMappingGetIndices()
 */
static PetscErrorCode StagBLArrayLocalToGlobal_Simple(StagBLArray array)
{
  PetscErrorCode     ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;
  DM                 dm;
  Vec                shell_local, shell_global;

  PetscFunctionBegin;
  ierr = StagBLGridPETScGetDM(array->grid,&dm);CHKERRQ(ierr);
  if (!data->local) StagBLError(PetscObjectComm((PetscObject)dm),"Cannot perform local-to-global before local data is defined");
  if (!array->current_local) StagBLError(PetscObjectComm((PetscObject)dm),"Refusing to perform global-to-local with stale global data");
  if (!data->global) {
    ierr = StagBLArraySimpleCreateGlobalVector_Private(array);CHKERRQ(ierr);
  }
  ierr = StagBLArraySimpleCreateShellGlobalVec_Private(array, &shell_global);CHKERRQ(ierr);
  ierr = StagBLArraySimpleCreateShellLocalVec_Private(array, &shell_local);CHKERRQ(ierr);
  ierr = DMLocalToGlobal(dm,shell_local,INSERT_VALUES,shell_global);CHKERRQ(ierr);
  ierr = VecDestroy(&shell_global);CHKERRQ(ierr);
  ierr = VecDestroy(&shell_local);CHKERRQ(ierr);
  array->current_global = PETSC_TRUE;
  PetscFunctionReturn(0);
}

static PetscErrorCode StagBLArrayPrint_Simple(StagBLArray array)
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

static PetscErrorCode StagBLArraySetLocalConstant_Simple(StagBLArray array, PetscScalar value)
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
  stagblarray->ops->getlocalvaluesstencil = StagBLArrayGetLocalValuesStencil_Simple;
  stagblarray->ops->globaltolocal = StagBLArrayGlobalToLocal_Simple;
  stagblarray->ops->localtoglobal = StagBLArrayLocalToGlobal_Simple;
  stagblarray->ops->print = StagBLArrayPrint_Simple;
  stagblarray->ops->setlocalvaluesstencil = StagBLArraySetLocalValuesStencil_Simple;
  stagblarray->ops->setlocalconstant = StagBLArraySetLocalConstant_Simple;
  return 0;
}

PetscErrorCode StagBLArraySimpleGetGlobalRaw(StagBLArray array,PetscScalar **raw)
{
  PetscErrorCode     ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;

  PetscFunctionBegin;
  // FIXME definitely need type check here!
  if (!data->global) {
    ierr = StagBLArraySimpleCreateGlobalVector_Private(array);CHKERRQ(ierr);
  }
  *raw = data->global;
  array->current_local = PETSC_FALSE;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLArraySimpleGetLocalRaw(StagBLArray array,PetscScalar **raw)
{
  PetscErrorCode     ierr;
  StagBLArray_Simple *data = (StagBLArray_Simple*) array->data;

  PetscFunctionBegin;
  // FIXME definitely need type check here!
  if (!data->local) {
    ierr = StagBLArraySimpleCreateLocalVector_Private(array);CHKERRQ(ierr);
  }
  *raw = data->local;
  array->current_global = PETSC_FALSE;
  PetscFunctionReturn(0);
}
