# This file is meant to be included from $(STAGBL_ARCH)/Makefile, after
# $(STAGBL_ARCH)/variables.mk

OBJDIR ?= obj
LIBDIR ?= lib

all : library

.PHONY : library

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

library : $(libstagbl)

$(OBJDIR)/%.o: $(OBJDIR)/%.c
	$(STAGBL_COMPILE.c) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $$(@D)/.DIR
	$(STAGBL_COMPILE.c) $< -o $@

.PHONY: all clean print

clean:
	rm -rf $(OBJDIR) $(LIBDIR)

srcs.c := $(libstagbl-y.c)
srcs.o := $(call srctoobj,$(srcs.c))
srcs.d := $(srcs.o:%.o=%.d)
# Tell make that srcs.d are all up to date.  Without this, the include
# below has quadratic complexity, taking more than one second for a
# do-nothing build of PETSc (much worse for larger projects)
$(srcs.d) : ;

-include $(srcs.d)
