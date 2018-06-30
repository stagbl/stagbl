#if !defined(STAGBL_H_)
#define STAGBL_H_

// Hard-coding for now (later put in generated header)
#define STAGBL_HAVE_PETSC

#if defined(STAGBL_HAVE_PETSC)
#include <petsc.h>
#endif
#include <mpi.h>

void StagBLInitialize(int,char**,MPI_Comm);
void StagBLFinalize();

// StagBLGrid
struct _p_StagBLGrid;
typedef struct _p_StagBLGrid *StagBLGrid;
void StagBLGridCreate(StagBLGrid*);
void StagBLGridDestroy(StagBLGrid*);

// StagBLGrid impls
#if defined (STAGBL_HAVE_PETSC)
#define STAGBLGRIDPETSC "petsc"
void StagBLGridCreate_PETSc(StagBLGrid);
void StagBLGridPETScGetDMPointer(StagBLGrid,DM**);
#endif

// StagBLArray
struct _p_StagBLArray;
typedef struct _p_StagBLArray *StagBLArray;
void StagBLArrayCreate(StagBLArray*);
void StagBLArrayDestroy(StagBLArray*);

// StagBLArray impls
#if defined (STAGBL_HAVE_PETSC)
#define STAGBLARRAYPETSC "petsc"
void StagBLArrayCreate_PETSc(StagBLArray);
void StagBLArrayPETScGetVecPointer(StagBLArray,Vec**);
#endif

// StagBLOperator
struct _p_StagBLOperator;
typedef struct _p_StagBLOperator *StagBLOperator;
void StagBLOperatorCreate(StagBLOperator*);
void StagBLOperatorDestroy(StagBLOperator*);

// StagBLOperator impls
#if defined (STAGBL_HAVE_PETSC)
#define STAGBLOPERATORPETSC "petsc"
void StagBLOperatorCreate_PETSc(StagBLOperator);
void StagBLOperatorPETScGetMatPointer(StagBLOperator,Mat**);
#endif

// StagBLLinearSolver
struct _p_StagBLLinearSolver;
typedef struct _p_StagBLLinearSolver *StagBLLinearSolver;
void StagBLLinearSolverCreate(StagBLLinearSolver*);
void StagBLLinearSolverDestroy(StagBLLinearSolver*);

// StagBLLinearSolver impls
#if defined (STAGBL_HAVE_PETSC)
#define STAGBLLINEARSOLVERPETSC "petsc"
void StagBLLinearSolverCreate_PETSc(StagBLLinearSolver);
void StagBLLinearSolverPETScGetKSPPointer(StagBLLinearSolver,KSP**);
#endif

#endif
