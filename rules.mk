# This file is meant to be included from stagbl.mk, which is turn included
# by $(STAGBL_ARCH)/Makefile, along with $(STAGBL_ARCH)/variables.mk
# This pattern can be copied to more easily define makefiles for compatible applications.

.SECONDEXPANSION:		# to expand $$(@D)/.DIR
.SUFFIXES:	        # Clear .SUFFIXES because we don't use implicit rules
.DELETE_ON_ERROR:   # Delete likely-corrupt target file if rule fails

ifeq ($(V),)
  quiet_HELP := "Use \"$(MAKE) V=1\" to see the verbose compile lines.\n"
  quiet = @printf $(quiet_HELP)$(eval quiet_HELP:=)"  %10s %s\n" "$1$2" "$@"; $($1)
else ifeq ($(V),0)		# Same, but do not print any help
  quiet = @printf "  %10s %s\n" "$1$2" "$@"; $($1)
else				# Show the full command line
  quiet = $($1)
endif

# gcc/gfortran style dependency flags; these are set in petscvariables starting with petsc-3.5
C_DEPFLAGS ?= $(if $(CONFIG_XLCOMPILER),-qmakedep=gcc,-MMD -MP)

# GCC-style syntax for C99.  Use "make C99FLAGS=-qlanglvl=extc99" or similar
# on systems that use different syntax to specify C99.
C99FLAGS := $(if $(findstring c99,$(PCC_FLAGS) $(STAGBL_CFLAGS) $(CFLAGS)),,$(if $(CONFIG_XLCOMPILER),-qlanglvl=extc99,-std=c99))

# Public headers for the library
STAGBL_INCLUDE = -I$(STAGBL_DIR)/include

# FIXME: consider putting info about libstagbl itself in here, so that it can be more easily linked to. stagbl.mk picks up this info and uses it to actually build the library.

# Compile and Link
STAGBL_COMPILE.c = $(call quiet,CC) -c $(C99FLAGS) $(STAGBL_CPPFLAGS) $(CPPFLAGS) $(STAGBL_CFLAGS) $(CFLAGS) $(C_DEPFLAGS) $(STAGBL_INCLUDE)
STAGBL_LINK = $(call quiet,CCLD) $(STAGBL_CFLAGS) $(CFLAGS) $(STAGBL_LDFLAGS) $(LDFLAGS) -o $@
CC = $(STAGBL_CC)
CCLD = $(if $(PCC_LINKER),$(PCC_LINKER),$(STAGBL_CC))

# Build libraries
%.$(AR_LIB_SUFFIX) : | $$(@D)/.DIR
	$(call quiet,AR) $(AR_FLAGS) $@ $^
	$(call quiet,RANLIB) $@

%/.DIR :
	@mkdir -p $(@D)
	@touch $@

.PRECIOUS: %/.DIR

# make print VAR=the-variable
print:
	@echo $($(VAR))
