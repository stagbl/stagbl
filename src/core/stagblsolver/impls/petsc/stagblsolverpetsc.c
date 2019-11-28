#include "stagbl/private/stagblsolverimpl.h"
#include "stagblsolverpetscimpl.h"
#include <stdlib.h>

PetscErrorCode StagBLSolverDestroy_PETSc(StagBLSolver solver)
{
  StagBLSolver_PETSc *data = (StagBLSolver_PETSc*) solver->data;

  PetscFunctionBegin;
  if (data->snes) {
    SNESDestroy(&data->snes);
  }
  free(solver->data);
  solver->data = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSolverSolve_PETSc(StagBLSolver solver,StagBLArray sol)
{
  StagBLSolver_PETSc * const data = (StagBLSolver_PETSc*) solver->data;
  StagBLGrid                 grid;
  PetscErrorCode             ierr;
  Vec                        vecRHS,vecSol;
  DM                         dm;

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
    /* Set default solve as direct, if possible */
    {
      PetscMPIInt size;
      PC          pc;


      ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm),&size);CHKERRQ(ierr);
      if (size == 1) {
#ifdef PETSC_HAVE_SUITESPARSE
        ierr = KSPSetType(data->ksp,KSPFGMRES);CHKERRQ(ierr);
        ierr = KSPGetPC(data->ksp,&pc);CHKERRQ(ierr);
        ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERUMFPACK);CHKERRQ(ierr);
#endif
      } else {
#ifdef PETSC_HAVE_SUPERLU_DIST
        ierr = KSPSetType(data->ksp,KSPFGMRES);CHKERRQ(ierr);
        ierr = KSPGetPC(data->ksp,&pc);CHKERRQ(ierr);
        ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU_DIST);CHKERRQ(ierr);
#endif
      }
    }
    ierr = KSPSetFromOptions(data->ksp);CHKERRQ(ierr); // TODO this might become problematic - need to figure out prefixes
  }

  ierr = StagBLSystemPETScGetVec(solver->system,&vecRHS);CHKERRQ(ierr);

  ierr = KSPSolve(data->ksp,vecRHS,vecSol);CHKERRQ(ierr);
  {
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(data->ksp,&reason);CHKERRQ(ierr);
    if (reason < 0) SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_CONV_FAILED,"Linear solve failed: %s",KSPConvergedReasons[reason]);
  }

  return 0;
}

PetscErrorCode StagBLSolverCreate_PETSc(StagBLSolver stagblsolver)
{
  StagBLSolver_PETSc *data;
  stagblsolver->data = (void*) malloc(sizeof(StagBLSolver_PETSc));
  data = (StagBLSolver_PETSc*) stagblsolver->data;
  data->ksp = NULL;
  stagblsolver->ops->destroy = StagBLSolverDestroy_PETSc;
  stagblsolver->ops->solve   = StagBLSolverSolve_PETSc;
  return 0;
}
