#if !defined(STAGBL_H_)
#define STAGBL_H_

#include <petsc.h>

PetscErrorCode StagBLInitialize(int,char**,const char*,MPI_Comm);
PetscErrorCode StagBLFinalize();

#define StagBLError(comm,str) SETERRQ(comm,PETSC_ERR_LIB,"StagBL Error: "str);
#define StagBLError1(comm,str,arg) SETERRQ1(comm,PETSC_ERR_LIB,"StagBL Error: "str,arg);

#define STAGBL_UNUSED(x) (void) x

/* StagBLGrid Data */
struct data_StagBLGrid;
typedef struct data_StagBLGrid *StagBLGrid;

/* StagBLArray Data */
struct data_StagBLArray;
typedef struct data_StagBLArray *StagBLArray;

/* StagBLSystem Data */
struct data_StagBLSystem;
typedef struct data_StagBLSystem *StagBLSystem;

/* StagBLSolver Data */
struct data_StagBLSolver;
typedef struct data_StagBLSolver *StagBLSolver;

/* StagBLGrid Functions */
PetscErrorCode StagBLGridCreate(StagBLGrid*);
PetscErrorCode StagBLGridCreateCompatibleStagBLGrid(StagBLGrid,PetscInt,PetscInt,PetscInt,PetscInt,StagBLGrid*);
PetscErrorCode StagBLGridCreateStagBLArray(StagBLGrid,StagBLArray*);
PetscErrorCode StagBLGridCreateStagBLSystem(StagBLGrid,StagBLSystem*);
PetscErrorCode StagBLGridDestroy(StagBLGrid*);

PetscErrorCode StagBLGridPETScGetDM(StagBLGrid,DM*);
PetscErrorCode StagBLGridPETScGetDMPointer(StagBLGrid,DM**);

/* StagBLArray Functions */
PetscErrorCode StagBLArrayCreate(StagBLGrid,StagBLArray*); // TODO this needn't be public (or even exist?) because we can only create StagBLArrays from StagBLGrids
PetscErrorCode StagBLArrayDestroy(StagBLArray*);
PetscErrorCode StagBLArrayGetStagBLGrid(StagBLArray,StagBLGrid*);

/* StagBLArray impls */
#define STAGBLARRAYPETSC "petsc"
PetscErrorCode StagBLArrayCreate_PETSc(StagBLArray);
PetscErrorCode StagBLArrayPETScGetLocalVec(StagBLArray,Vec*);
PetscErrorCode StagBLArrayPETScGetGlobalVec(StagBLArray,Vec*);
PetscErrorCode StagBLArrayPETScGetLocalVecPointer(StagBLArray,Vec**);
PetscErrorCode StagBLArrayPETScGetGlobalVecPointer(StagBLArray,Vec**);

/* StagBLSystem Functions */
PetscErrorCode StagBLSystemCreate(StagBLGrid,StagBLSystem*);
PetscErrorCode StagBLSystemCreateStagBLSolver(StagBLSystem,StagBLSolver*);
PetscErrorCode StagBLSystemDestroy(StagBLSystem*);

/* StagBLSystem impls */
#define STAGBLSYSTEMPETSC "petsc"
PetscErrorCode StagBLSystemCreate_PETSc(StagBLSystem);
PetscErrorCode StagBLSystemGetGrid(StagBLSystem,StagBLGrid*);
PetscErrorCode StagBLSystemPETScGetMat(StagBLSystem,Mat*);
PetscErrorCode StagBLSystemPETScGetMatPointer(StagBLSystem,Mat**);
PetscErrorCode StagBLSystemPETScGetVec(StagBLSystem,Vec*);
PetscErrorCode StagBLSystemPETScGetVecPointer(StagBLSystem,Vec**);
PetscErrorCode StagBLSystemPETScGetResidualFunction(StagBLSystem,PetscErrorCode (**f)(SNES,Vec,Vec,void*));
PetscErrorCode StagBLSystemPETScGetJacobianFunction(StagBLSystem,PetscErrorCode (**f)(SNES,Vec,Mat,Mat,void*));

/* StagBLSolver Functions */
PetscErrorCode StagBLSolverCreate(StagBLSystem,StagBLSolver*);
PetscErrorCode StagBLSolverDestroy(StagBLSolver*);
PetscErrorCode StagBLSolverGetSystem(StagBLSolver,StagBLSystem*);
PetscErrorCode StagBLSolverSolve(StagBLSolver,StagBLArray);

/* StagBLSolver impls */
#define STAGBLSOLVERPETSC "petsc"
PetscErrorCode StagBLSolverCreate_PETSc(StagBLSolver);

/* Stokes */
struct data_StagBLStokesParameters {
  StagBLGrid  stokes_grid,temperature_grid;
  StagBLArray coefficient_array,temperature_array;
  PetscBool   uniform_grid, boussinesq_forcing;
  PetscReal   xmin,xmax,ymin,ymax,zmin,zmax;
  PetscReal   gy,eta_characteristic,alpha;
};
typedef struct data_StagBLStokesParameters *StagBLStokesParameters;

PetscErrorCode StagBLStokesParametersCreate(StagBLStokesParameters*);
PetscErrorCode StagBLStokesParametersDestroy(StagBLStokesParameters*);

PetscErrorCode StagBLGridCreateStokes2DBox(MPI_Comm,PetscInt,PetscInt,PetscScalar,PetscScalar,PetscScalar,PetscScalar,StagBLGrid*);
PetscErrorCode StagBLGridCreateStokes3DBox(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar,StagBLGrid*);
PetscErrorCode StagBLCreateStokesSystem(StagBLStokesParameters,StagBLSystem*);

PetscErrorCode StagBLDumpStokes(StagBLStokesParameters,StagBLArray,PetscInt);

#endif
