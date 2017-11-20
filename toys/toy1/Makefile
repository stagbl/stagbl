EXNAME=runme

SRC_C = ${EXNAME}.c system.c ctx.c output.c
SRC_O = ${SRC_C:%.c=%.o}
SRC_D = ${SRC_C:%.c=%.d}

CFLAGS+=${C_DEPFLAGS}

all :: ${EXNAME}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

${EXNAME}: ${SRC_O}
	-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB}

# Aggressive clean
clean :: 
	${RM} ${EXNAME} ${SRC_O} ${SRC_D} *.xmf *.petscbin *.petscbin.info *.h5 *.bin *.bin.info *.vts octave-workspace

# Indicate that SRC_D is all up to date. This apparently prevents the include
#   from having quadratic complexity. 
$(SRC_D) : ; 
-include $(SRC_D)

.DELETE_ON_ERROR:
