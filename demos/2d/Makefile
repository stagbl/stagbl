STAGBL_DIR ?= ../..
STAGBL_ARCH ?= $(PETSC_ARCH)
ifeq ($(STAGBL_ARCH),)
$(error STAGBL_ARCH or PETSC_ARCH must be set)
endif

EXNAME = stagbldemo2d

OBJ = args.o coeff.o ctx.o dump.o main.o nusselt.o particles.o temperature.o

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
