STAGBL_DIR ?= ..
STAGBL_ARCH ?= $(PETSC_ARCH)
ifeq ($(STAGBL_ARCH),)
$(error STAGBL_ARCH or PETSC_ARCH must be set)
endif

all: 2d 3d

EXNAME2D = stagbldemo2d
EXNAME3D = stagbldemo3d

2d: $(EXNAME2D)
3d: $(EXNAME3D)

OBJ2D = args.o coeff.o ctx.o dump.o main_2d.o nusselt.o particles.o temperature.o
OBJ3D = args.o coeff.o ctx.o dump.o main_3d.o

all : $(EXNAME)

# For ease of compatibility, use the same makefiles that libstagbl uses
include $(STAGBL_DIR)/$(STAGBL_ARCH)/variables.mk
include $(STAGBL_DIR)/rules.mk

%.o: %.c
	$(STAGBL_COMPILE.c) $< -o $@

$(EXNAME2D) : $(OBJ2D)
	$(STAGBL_LINK) $^ $(STAGBL_LIB)

$(EXNAME3D) : $(OBJ3D)
	$(STAGBL_LINK) $^ $(STAGBL_LIB)

clean_output:
	rm -f *.vtr *.pbin *.xmf

clean : clean_output
	rm -f $(EXNAME2D) $(EXNAME3D)
	rm -f $(OBJ2D) $(OBJ2D:.o=.d) $(OBJ3D) $(OBJ3D:.o=.d)
	rm -rf $(EXNAME2D).dSYM $(EXNAME3D).dSYM

.PHONY : all 2d 3d clean clean_output
