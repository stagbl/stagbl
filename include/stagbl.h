#if !defined(STAGBL_H_)
#define STAGBL_H_

#include <petsc.h>

PetscErrorCode StagBLInitialize(int,char**,const char*,MPI_Comm);
PetscErrorCode StagBLFinalize();

#define StagBLError(comm,msg) SETERRQ1(comm,PETSC_ERR_LIB,"StagBL Error: %s",msg);

#define STAGBL_UNUSED(x) (void) x

// StagBLGrid Data
struct _p_StagBLGrid;
typedef struct _p_StagBLGrid *StagBLGrid;

// StagBLArray Data
struct _p_StagBLArray;
typedef struct _p_StagBLArray *StagBLArray;

// StagBLSystem Data
struct _p_StagBLSystem;
typedef struct _p_StagBLSystem *StagBLSystem;

// StagBLSolver Data
struct _p_StagBLSolver;
typedef struct _p_StagBLSolver *StagBLSolver;

// StagBLGrid Functions
PetscErrorCode StagBLGridCreate(StagBLGrid*);
PetscErrorCode StagBLGridCreateCompatibleStagBLGrid(StagBLGrid,PetscInt,PetscInt,PetscInt,PetscInt,StagBLGrid*);
PetscErrorCode StagBLGridCreateStagBLArray(StagBLGrid,StagBLArray*);
PetscErrorCode StagBLGridCreateStagBLSystem(StagBLGrid,StagBLSystem*);
PetscErrorCode StagBLGridDestroy(StagBLGrid*);

PetscErrorCode StagBLGridPETScGetDM(StagBLGrid,DM*);
PetscErrorCode StagBLGridPETScGetDMPointer(StagBLGrid,DM**);

// StagBLArray Functions
PetscErrorCode StagBLArrayCreate(StagBLGrid,StagBLArray*); // TODO this needn't be public (or even exist?) because we can only create StagBLArrays from StagBLGrids
PetscErrorCode StagBLArrayDestroy(StagBLArray*);
PetscErrorCode StagBLArrayGetStagBLGrid(StagBLArray,StagBLGrid*);

// StagBLArray impls
#define STAGBLARRAYPETSC "petsc"
PetscErrorCode StagBLArrayCreate_PETSc(StagBLArray);
PetscErrorCode StagBLArrayPETScGetLocalVec(StagBLArray,Vec*);
PetscErrorCode StagBLArrayPETScGetGlobalVec(StagBLArray,Vec*);
PetscErrorCode StagBLArrayPETScGetLocalVecPointer(StagBLArray,Vec**);
PetscErrorCode StagBLArrayPETScGetGlobalVecPointer(StagBLArray,Vec**);

// StagBLSystem Functions
PetscErrorCode StagBLSystemCreate(StagBLGrid,StagBLSystem*);
PetscErrorCode StagBLSystemCreateStagBLSolver(StagBLSystem,StagBLSolver*);
PetscErrorCode StagBLSystemDestroy(StagBLSystem*);

// StagBLSystem impls
#define STAGBLSYSTEMPETSC "petsc"
PetscErrorCode StagBLSystemCreate_PETSc(StagBLSystem);
PetscErrorCode StagBLSystemGetGrid(StagBLSystem,StagBLGrid*);
PetscErrorCode StagBLSystemPETScGetMat(StagBLSystem,Mat*);
PetscErrorCode StagBLSystemPETScGetMatPointer(StagBLSystem,Mat**);
PetscErrorCode StagBLSystemPETScGetVec(StagBLSystem,Vec*);
PetscErrorCode StagBLSystemPETScGetVecPointer(StagBLSystem,Vec**);
PetscErrorCode StagBLSystemPETScGetResidualFunction(StagBLSystem,PetscErrorCode (**f)(SNES,Vec,Vec,void*));
PetscErrorCode StagBLSystemPETScGetJacobianFunction(StagBLSystem,PetscErrorCode (**f)(SNES,Vec,Mat,Mat,void*));

// StagBLSolver Functions
PetscErrorCode StagBLSolverCreate(StagBLSystem,StagBLSolver*);
PetscErrorCode StagBLSolverDestroy(StagBLSolver*);
PetscErrorCode StagBLSolverGetSystem(StagBLSolver,StagBLSystem*);
PetscErrorCode StagBLSolverSolve(StagBLSolver,StagBLArray);

// StagBLSolver impls
#define STAGBLSOLVERPETSC "petsc"
PetscErrorCode StagBLSolverCreate_PETSc(StagBLSolver);

// Stokes
PetscErrorCode StagBLGridCreateStokes2DBox(MPI_Comm,PetscInt,PetscInt,PetscScalar,PetscScalar,PetscScalar,PetscScalar,StagBLGrid*);

#endif
