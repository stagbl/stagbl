#ifndef CTX_H_
#define CTX_H_

#include "stagbl.h"

/**
  * A single object which represents settings and state for a 2D application
  *  code.
  *
  * For simplicity, it combines both settings and pointers to objects
  * representing the state of the simulation.
  */
typedef struct {
  MPI_Comm     comm;
  char         mode[1024];
  PetscInt     totalTimesteps; /* can be 0, only compute initial temperature field and solve Stokes once */
  PetscReal    dt;
  StagBLGrid   stokes_grid,coefficient_grid,temperature_grid;
  StagBLArray  coefficient_array,stokes_array,temperature_array;
  StagBLSystem temperature_system,stokes_system;
  StagBLSolver temperature_solver,stokes_solver;
  PetscBool    uniform_grid,boussinesq_forcing;
  PetscScalar  xmax,ymax,xmin,ymin,hxCharacteristic,hyCharacteristic;
  PetscScalar  eta1,eta2,rho1,rho2,gy,Kbound,Kcont,etaCharacteristic,KTemp;
  PetscReal    temperature_top,temperature_bottom,kappa,alpha;
  PetscInt     pinx,piny;
  PetscScalar (*getEta)(void*,PetscScalar,PetscScalar);
  PetscScalar (*getRho)(void*,PetscScalar,PetscScalar);
  DM           dm_particles;
} data_Ctx;
typedef data_Ctx* Ctx;

PetscErrorCode CtxCreate(MPI_Comm,const char* mode,Ctx*);
PetscErrorCode CtxDestroy(Ctx*);
PetscErrorCode CtxSetupFromGrid(Ctx);

#endif
