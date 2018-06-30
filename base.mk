# This file is meant to be included from $(STAGBL_ARCH)/Makefile

.SECONDEXPANSION:		# to expand $$(@D)/.DIR
.SUFFIXES:	        # Clear .SUFFIXES because we don't use implicit rules
.DELETE_ON_ERROR:   # Delete likely-corrupt target file if rule fails

OBJDIR ?= obj
LIBDIR ?= lib

# function to prefix directory that contains most recently-parsed
# makefile (current) if that directory is not ./
thisdir = $(addprefix $(dir $(lastword $(MAKEFILE_LIST))),$(1))
incsubdirs = $(addsuffix /local.mk,$(call thisdir,$(1)))
srctoobj = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(filter-out $(OBJDIR)/%,$(1)))

libstagbl-y.c :=

all : library

.PHONY : library

# Recursively include files for all targets
include $(SRCDIR)/local.mk

#### Rules ####
ifeq ($(V),)
  quiet_HELP := "Use \"$(MAKE) V=1\" to see the verbose compile lines.\n"
  quiet = @printf $(quiet_HELP)$(eval quiet_HELP:=)"  %10s %s\n" "$1$2" "$@"; $($1)
else ifeq ($(V),0)		# Same, but do not print any help
  quiet = @printf "  %10s %s\n" "$1$2" "$@"; $($1)
else				# Show the full command line
  quiet = $($1)
endif

libstagbl = $(LIBDIR)/libstagbl.$(AR_LIB_SUFFIX)
libstagbl : $(libstagbl)
$(libstagbl) : $(call srctoobj,$(libstagbl-y.c))

library : $(libstagbl)


# TODO This doesn't actually work without PETSc. Need to define these AR-style things to build the library.

%.$(AR_LIB_SUFFIX) : | $$(@D)/.DIR
	$(call quiet,AR) $(AR_FLAGS) $@ $^
	$(call quiet,RANLIB) $@

# gcc/gfortran style dependency flags; these are set in petscvariables starting with petsc-3.5
C_DEPFLAGS ?= $(if $(CONFIG_XLCOMPILER),-qmakedep=gcc,-MMD -MP)

# GCC-style syntax for C99.  Use "make C99FLAGS=-qlanglvl=extc99" or similar
# on systems that use different syntax to specify C99.
C99FLAGS := $(if $(findstring c99,$(PCC_FLAGS) $(STAGBL_CFLAGS) $(CFLAGS)),,$(if $(CONFIG_XLCOMPILER),-qlanglvl=extc99,-std=c99))

# TODO change the name of SRCDIR perhaps to STAGBL_ROOT_DIR or something?
STAGBL_INCLUDE = -I$(SRCDIR)/include 

STAGBL_COMPILE.c = $(call quiet,CC) -c $(C99FLAGS) $(STAGBL_CPPFLAGS) $(CPPFLAGS) $(STAGBL_CFLAGS) $(CFLAGS) $(C_DEPFLAGS) $(STAGBL_INCLUDE)
STAGBL_LINK = $(call quiet,CCLD) $(STAGBL_CFLAGS) $(CFLAGS) $(STAGBL_LDFLAGS) $(LDFLAGS) -o $@
CC = $(STAGBL_CC)
CCLD = $(if $(PCC_LINKER),$(PCC_LINKER),$(STAGBL_CC))

$(OBJDIR)/%.o: $(OBJDIR)/%.c
	$(STAGBL_COMPILE.c) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $$(@D)/.DIR
	$(STAGBL_COMPILE.c) $< -o $@

%/.DIR :
	@mkdir -p $(@D)
	@touch $@

.PRECIOUS: %/.DIR

.PHONY: all clean print

clean:
	rm -rf $(OBJDIR) $(LIBDIR)

# make print VAR=the-variable
print:
	@echo $($(VAR))

srcs.c := $(libstagbl-y.c)
srcs.o := $(call srctoobj,$(srcs.c))
srcs.d := $(srcs.o:%.o=%.d)
# Tell make that srcs.d are all up to date.  Without this, the include
# below has quadratic complexity, taking more than one second for a
# do-nothing build of PETSc (much worse for larger projects)
$(srcs.d) : ;

-include $(srcs.d)
