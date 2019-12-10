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
  Vec                        vec_sol;
  DM                         dm;

  PetscFunctionBegin;
  ierr = StagBLArrayGetStagBLGrid(sol,&grid);CHKERRQ(ierr);
  ierr = StagBLGridPETScGetDM(grid,&dm);CHKERRQ(ierr);

  ierr = StagBLArrayPETScGetGlobalVec(sol,&vec_sol);CHKERRQ(ierr);

  /* Create the solution Vec, if needbe */
  if (!vec_sol)
  {
    Vec *p_vec_sol;

    ierr = StagBLArrayPETScGetGlobalVecPointer(sol,&p_vec_sol);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm,p_vec_sol);CHKERRQ(ierr);
    vec_sol = *p_vec_sol;
  }

  /* Create the SNES object from the system, if needbe */
  if (!data->snes) {
    Mat mat;
    PetscErrorCode (*residual_function)(SNES,Vec,Vec,void*);
    PetscErrorCode (*jacobian_function)(SNES,Vec,Mat,Mat,void*);

    ierr = StagBLSystemPETScGetResidualFunction(solver->system,&residual_function);CHKERRQ(ierr);
    ierr = StagBLSystemPETScGetJacobianFunction(solver->system,&jacobian_function);CHKERRQ(ierr);

    ierr = StagBLSystemPETScGetMat(solver->system,&mat);CHKERRQ(ierr);
    ierr = SNESCreate(PetscObjectComm((PetscObject)dm),&data->snes);CHKERRQ(ierr);
    ierr = SNESSetFunction(data->snes,NULL,residual_function,(void*)solver->system);CHKERRQ(ierr);
    ierr = SNESSetJacobian(data->snes,NULL,NULL,jacobian_function,(void*)solver->system);CHKERRQ(ierr);
    {
      PetscMPIInt size;
      KSP         ksp;
      PC          pc;

      ierr = SNESSetType(data->snes,SNESKSPONLY);CHKERRQ(ierr);
      ierr = SNESGetKSP(data->snes,&ksp);CHKERRQ(ierr);
      ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);
      ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

      ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm),&size);CHKERRQ(ierr);
      if (size == 1) {
#ifdef PETSC_HAVE_SUITESPARSE
        ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERUMFPACK);CHKERRQ(ierr);
#endif
      } else {
#ifdef PETSC_HAVE_SUPERLU_DIST
        ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU_DIST);CHKERRQ(ierr);
#endif
      }
    }
    ierr = SNESSetFromOptions(data->snes);CHKERRQ(ierr);
  }

  ierr = SNESSolve(data->snes,NULL,vec_sol);CHKERRQ(ierr);
  {
    SNESConvergedReason reason;

    ierr = SNESGetConvergedReason(data->snes,&reason);CHKERRQ(ierr);
    if (reason < 0) SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_CONV_FAILED,"Solve failed: %s",SNESConvergedReasons[reason]);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode StagBLSolverCreate_PETSc(StagBLSolver stagblsolver)
{
  StagBLSolver_PETSc *data;

  PetscFunctionBegin;
  stagblsolver->data = (void*) malloc(sizeof(StagBLSolver_PETSc));
  data = (StagBLSolver_PETSc*) stagblsolver->data;
  data->snes = NULL;
  stagblsolver->ops->destroy = StagBLSolverDestroy_PETSc;
  stagblsolver->ops->solve   = StagBLSolverSolve_PETSc;
  PetscFunctionReturn(0);
}
