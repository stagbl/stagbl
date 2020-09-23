# This file is meant to be included from $(STAGBL_ARCH)/Makefile, after
# $(STAGBL_ARCH)/variables.mk

BINDIR  ?= bin
INCDIR  ?= include
LIBDIR  ?= lib
OBJDIR  ?= obj
TESTDIR ?= test

all : library tests

# function to prefix directory that contains most recently-parsed
# makefile (current) if that directory is not ./
thisdir = $(addprefix $(dir $(lastword $(MAKEFILE_LIST))),$(1))
incsubdirs = $(addsuffix /local.mk,$(call thisdir,$(1)))
srctoobj = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(filter-out $(OBJDIR)/%,$(1)))

# Rules
include $(STAGBL_DIR)/rules.mk

# Initialize set of targets and recursively include files for all targets
libstagbl-y.c :=
include $(SRCDIR)/local.mk

# Build libstagbl from sources here
libstagbl = $(LIBDIR)/libstagbl.$(AR_LIB_SUFFIX)
libstagbl : $(libstagbl)
$(libstagbl) : $(call srctoobj,$(libstagbl-y.c))

library : $(libstagbl) includes

.PHONY : library

$(OBJDIR)/%.o: $(OBJDIR)/%.c
	$(STAGBL_COMPILE.c) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $$(@D)/.DIR
	$(STAGBL_COMPILE.c) $< -o $@

# Configuration-specific includes
includes : $(INCDIR)/.DIR

.PHONY: includes

clean:
	rm -rf $(BINDIR) $(INCDIR) $(LIBDIR) $(OBJDIR) $(TESTDIR)

.PHONY: all clean print

srcs.c := $(libstagbl-y.c)
srcs.o := $(call srctoobj,$(srcs.c))
srcs.d := $(srcs.o:%.o=%.d)
# Tell make that srcs.d are all up to date.  Without this, the include
# below has quadratic complexity
$(srcs.d) : ;

-include $(srcs.d)

# Demo applications for testing and convenience. Since these are mainly
# intended as standalone example of using the library, we use a phony target
# to clean and rebuild them with their own makefiles, clumsily moving the resulting
# binary.
demos : library $(BINDIR)/.DIR
	$(MAKE) -C $(STAGBL_DIR)/demos clean
	$(MAKE) -C $(STAGBL_DIR)/demos all
	mv $(STAGBL_DIR)/demos/stagbldemo2d $(BINDIR)
	mv $(STAGBL_DIR)/demos/stagbldemo3d $(BINDIR)

.PHONY: demos

# Additional Test Executables
# This is currently done manually, and would be more maintainable with an automated process
tests : \
	$(BINDIR)/test_dmstag_vs_dmda \
	$(BINDIR)/test_dmstag_vs_dmda_mf_op \
	$(BINDIR)/test_dmstag_vs_dmda_matstencil \
	$(BINDIR)/test_dmstag_vec_stencil_vs_array \
	$(BINDIR)/test_dmstag_preallocate \
	$(BINDIR)/test_stagbl_array \
	$(BINDIR)/test_stokes_operator \
	$(BINDIR)/test_stokes_assembly_and_solve \

$(BINDIR)/test_dmstag_vs_dmda : $(OBJDIR)/src/tests/performance/test_dmstag_vs_dmda.o library | $$(@D)/.DIR
	$(STAGBL_LINK) $< $(STAGBL_LIB)
$(BINDIR)/test_dmstag_vs_dmda_mf_op : $(OBJDIR)/src/tests/performance/test_dmstag_vs_dmda_mf_op.o library | $$(@D)/.DIR
	$(STAGBL_LINK) $< $(STAGBL_LIB)
$(BINDIR)/test_dmstag_vs_dmda_matstencil : $(OBJDIR)/src/tests/performance/test_dmstag_vs_dmda_matstencil.o library | $$(@D)/.DIR
	$(STAGBL_LINK) $< $(STAGBL_LIB)
$(BINDIR)/test_dmstag_vec_stencil_vs_array : $(OBJDIR)/src/tests/performance/test_dmstag_vec_stencil_vs_array.o library | $$(@D)/.DIR
	$(STAGBL_LINK) $< $(STAGBL_LIB)
$(BINDIR)/test_dmstag_preallocate : $(OBJDIR)/src/tests/performance/test_dmstag_preallocate.o library | $$(@D)/.DIR
	$(STAGBL_LINK) $< $(STAGBL_LIB)
$(BINDIR)/test_stagbl_array : $(OBJDIR)/src/tests/unit/test_stagbl_array.o library | $$(@D)/.DIR
	$(STAGBL_LINK) $< $(STAGBL_LIB)
$(BINDIR)/test_stokes_operator : $(OBJDIR)/src/tests/unit/test_stokes_operator.o library | $$(@D)/.DIR
	$(STAGBL_LINK) $< $(STAGBL_LIB)
$(BINDIR)/test_stokes_assembly_and_solve : $(OBJDIR)/src/tests/integration/test_stokes_assembly_and_solve.o library | $$(@D)/.DIR
	$(STAGBL_LINK) $< $(STAGBL_LIB)

.PHONY : tests

# Run tests
STAGBL_SCIATH_COMMAND= cd ${TEST_DIR} && STAGBL_DIR=$(STAGBL_DIR) STAGBL_ARCH=$(STAGBL_ARCH) python -m sciath $(STAGBL_DIR)/tests/tests.yml -w sciath.conf

test : tests demos $(TESTDIR)/.DIR
	${STAGBL_SCIATH_COMMAND}

test_check : tests demos $(TESTDIR)/.DIR
	${STAGBL_SCIATH_COMMAND} -v

test_clean : $(TESTDIR)/.DIR
	${STAGBL_SCIATH_COMMAND} -p

.PHONY: test test_clean
