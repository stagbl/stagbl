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
stagbltests-y.c :=
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

srcs.c := $(libstagbl-y.c) $(stagbltests-y.c)
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
$(BINDIR)/test_% : $(OBJDIR)/src/tests/test_%.o library | $$(@D)/.DIR
	$(STAGBL_LINK) $< $(STAGBL_LIB)

tests : $(patsubst $(SRCDIR)/src/tests/%.c,$(BINDIR)/%,$(stagbltests-y.c))

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
