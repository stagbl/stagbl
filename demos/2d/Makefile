# Expects STAGBL_DIR and STAGBL_ARCH to be set
STAGBL_DIR ?= ../..
STAGBL_ARCH ?= $(PETSC_ARCH)
ifeq ($(STAGBL_ARCH),)
STAGBL_ARCH:=build
endif

EXNAME = stagbldemo2d

OBJ = args.o ctx.o dump.o main.o system.o system2.o

all : $(EXNAME)

# For ease of compatibility, use the same makefiles that libstagbl uses
include $(STAGBL_DIR)/$(STAGBL_ARCH)/variables.mk
include $(STAGBL_DIR)/rules.mk

%.o: %.c
	$(STAGBL_COMPILE.c) $< -o $@

$(EXNAME) : $(OBJ)
	$(STAGBL_LINK) $^ $(STAGBL_LIB)

clean :
	rm -f $(EXNAME)
	rm -f $(OBJ)
	rm -rf $(EXNAME).dSYM

.PHONY : all clean
