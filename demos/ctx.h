#ifndef CTX_H_
#define CTX_H_

#include "stagbl.h"

/**
  * A single object which represents settings and state for a 2D application
  *  code.
  *
  * For simplicity, it combines both settings and "data", mainly pointer to
  * objects representing the state of the simulation.
  */
typedef struct {
  /* A string which controls creation of the Ctx */
  char         mode[1024];

  /* Domain and boundary conditions */
  PetscBool    uniform_grid;
  PetscScalar  xmax,ymax,xmin,ymin,hxCharacteristic,hyCharacteristic;

  /* Settings */
  PetscBool    boussinesq_forcing;
  PetscBool    compute_nusselt_number;
  PetscScalar (*getEta)(void*,PetscScalar,PetscScalar);
  PetscScalar (*getRho)(void*,PetscScalar,PetscScalar);
  PetscReal    temperature_top,temperature_bottom,kappa,alpha;
  PetscScalar  eta1,eta2,rho1,rho2,gy;

  /* Timestepping */
  PetscInt     totalTimesteps; /* can be 0, only compute initial temperature field and solve Stokes once */
  PetscReal    dt;

  /* Derived/Computed Parameters */
  PetscScalar  Kbound,Kcont,etaCharacteristic,KTemp;
  PetscInt     pinx,piny;

  /* Data */
  MPI_Comm     comm;
  StagBLGrid   stokes_grid,coefficient_grid,temperature_grid;
  StagBLArray  coefficient_array,stokes_array,temperature_array;
  StagBLSystem temperature_system,stokes_system;
  StagBLSolver temperature_solver,stokes_solver;
  DM           dm_particles;
} data_Ctx;
typedef data_Ctx* Ctx;

PetscErrorCode CtxCreate(MPI_Comm,const char* mode,Ctx*);
PetscErrorCode CtxDestroy(Ctx*);
PetscErrorCode CtxSetupFromGrid(Ctx);

#endif
