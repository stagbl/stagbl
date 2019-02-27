#include "stagbl/private/stagblsolverimpl.h"
#include "stagblsolverpetscimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLSolverDestroy_PETSc(StagBLSolver solver)
{
  StagBLSolver_PETSc *data = (StagBLSolver_PETSc*) solver->data;
  if (data->ksp) {
    KSPDestroy(&data->ksp);
  }
  free(solver->data);
  solver->data = NULL;
  return 0;
}

StagBLErrorCode StagBLSolverPETScGetKSPPointer(StagBLSolver stagblsolver,KSP **ksp)
{
  StagBLSolver_PETSc * const data = (StagBLSolver_PETSc*) stagblsolver->data;
  *ksp = &(data->ksp);
  return 0;
}

StagBLErrorCode StagBLSolverSolve_PETSc(StagBLSolver solver,StagBLArray sol)
{
  StagBLSolver_PETSc * const data = (StagBLSolver_PETSc*) solver->data;
  StagBLGrid                 grid;
  PetscErrorCode ierr;
  Vec vecRHS,vecSol;
  DM dm;

  ierr = StagBLArrayGetStagBLGrid(sol,&grid);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(grid,&dm);CHKERRQ(ierr);

  ierr = StagBLArrayPETScGetGlobalVec(sol,&vecSol);CHKERRQ(ierr);

  /* Create the solution Vec, if needbe */
  if (!vecSol)
  {
    Vec *pvecSol;
    ierr = StagBLArrayPETScGetGlobalVecPointer(sol,&pvecSol);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm,pvecSol);CHKERRQ(ierr);
    vecSol = *pvecSol;
  }

  /* Create the KSP object from the system, if needbe */
  if (!data->ksp) {
    Mat mat;
    ierr = StagBLSystemPETScGetMat(solver->system,&mat);CHKERRQ(ierr);
    ierr = KSPCreate(PetscObjectComm((PetscObject)dm),&data->ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(data->ksp,mat,mat);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(data->ksp);CHKERRQ(ierr); // TODO this might become problematic
  }

  ierr = StagBLSystemPETScGetVec(solver->system,&vecRHS);CHKERRQ(ierr);

  ierr = KSPSolve(data->ksp,vecRHS,vecSol);CHKERRQ(ierr);
  {
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(data->ksp,&reason);CHKERRQ(ierr);
    if (reason < 0) SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_CONV_FAILED,"Linear solve failed: %s",KSPConvergedReasons[reason]);CHKERRQ(ierr);
  }

  return 0;
}

StagBLErrorCode StagBLSolverCreate_PETSc(StagBLSolver stagblsolver)
{
  StagBLSolver_PETSc *data;
  stagblsolver->data = (void*) malloc(sizeof(StagBLSolver_PETSc));
  data = (StagBLSolver_PETSc*) stagblsolver->data;
  data->ksp = NULL;
  stagblsolver->ops->destroy = StagBLSolverDestroy_PETSc;
  stagblsolver->ops->solve   = StagBLSolverSolve_PETSc;
  return 0;
}
