# Convenience Makefile for StagBL
#
# See the README for information on how to configure and build StagBL.
#
# This assumes that you have PETSC_DIR and PETSC_ARCH set,
# and want to use STAGBL_ARCH=PETSC_ARCH

all:
	+${MAKE} -C ${PETSC_ARCH} $@

%:
	+${MAKE} -C ${PETSC_ARCH} $@

.PHONY: all
