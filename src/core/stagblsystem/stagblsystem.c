#include "stagbl/private/stagblsystemimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLSystemCreate(StagBLGrid grid,StagBLSystem *system,StagBLSystemType type)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,system);CHKERRQ(ierr);
  ierr = PetscCalloc1(1,&(*system)->ops);CHKERRQ(ierr);

  (*system)->type = type;
  (*system)->grid = grid;

  /* Set the creation function and call it, which sets other ops */
  if (StagBLCheckType(type,STAGBLSYSTEMPETSC)) {
      (*system)->ops->create = StagBLSystemCreate_PETSc;
  } else if (StagBLCheckType(type,STAGBLSYSTEMSIMPLE)) {
      (*system)->ops->create = StagBLSystemCreate_Simple;
  } else StagBLError1(PETSC_COMM_WORLD,"System creation not implemented for type %s",type);
  ierr = ((*system)->ops->create)(*system);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemCreateStagBLSolver(StagBLSystem system,StagBLSolver *p_solver)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = StagBLSolverCreate(system,p_solver,system->solver_type);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemDestroy(StagBLSystem *system)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*system) PetscFunctionReturn(0);
  if ((*system)->ops->destroy) {
    ierr = ((*system)->ops->destroy)(*system);CHKERRQ(ierr);
  }
  ierr = PetscFree((*system)->ops);CHKERRQ(ierr);
  ierr = PetscFree(*system);CHKERRQ(ierr);
  *system = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemGetGrid(StagBLSystem system,StagBLGrid *grid)
{
  PetscFunctionBegin;
  *grid = system->grid;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemOperatorSetValuesStencil(StagBLSystem system,PetscInt nrows,const DMStagStencil *rows,PetscInt ncols,const DMStagStencil *cols, const PetscScalar *values)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (system->ops->operatorsetvaluesstencil) {
    ierr = (system->ops->operatorsetvaluesstencil)(system,nrows,rows,ncols,cols,values);CHKERRQ(ierr);
  } else StagBLError2(PETSC_COMM_WORLD,"%s not implemented for StagBLSystem object of type %s",__func__,system->type);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemRHSSetConstant(StagBLSystem system,PetscScalar value)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (system->ops->rhssetconstant) {
    ierr = (system->ops->rhssetconstant)(system,value);CHKERRQ(ierr);
  } else StagBLError2(PETSC_COMM_WORLD,"%s not implemented for StagBLSystem object of type %s",__func__,system->type);
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSystemRHSSetValuesStencil(StagBLSystem system,PetscInt nrows,const DMStagStencil *rows,const PetscScalar *values)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (system->ops->operatorsetvaluesstencil) {
    ierr = (system->ops->rhssetvaluesstencil)(system,nrows,rows,values);CHKERRQ(ierr);
  } else StagBLError2(PETSC_COMM_WORLD,"%s not implemented for StagBLSystem object of type %s",__func__,system->type);
  PetscFunctionReturn(0);
}
