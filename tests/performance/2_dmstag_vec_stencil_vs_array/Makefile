all : runme

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

runme: runme.o
	-${CLINKER} -o $@ $^ ${PETSC_LIB}
