# Expects STAGBL_DIR and STAGBL_ARCH to be set
STAGBL_DIR ?= ../..
STAGBL_ARCH ?= $(PETSC_ARCH)
ifeq ($(STAGBL_ARCH),)
STAGBL_ARCH:=build
endif

EXNAME = stagbldemo3d

all : $(EXNAME)

# For ease of compatibility, use the same makefiles that libstagbl uses
include $(STAGBL_DIR)/$(STAGBL_ARCH)/variables.mk
include $(STAGBL_DIR)/rules.mk

%.o: %.c
	$(STAGBL_COMPILE.c) $< -o $@

$(EXNAME) : main.o
	$(STAGBL_LINK) $< $(STAGBL_LIB)

clean :
	rm -f  $(EXNAME)
	rm -f  main.o main.d
	rm -rf $(EXNAME).dSYM

.PHONY : all clean
