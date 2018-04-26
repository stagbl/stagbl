StagBL Layout Proposal (WIP)
----------------------

A proposal for the layout of the StagBL project
Guidelines: keep it lean, keep it self-documenting as much as possible

config/

demos/

developer/                            # Everything that only the developer would care about (toys, util scripts, etc.)
  toys/

documentation/

LICENSE.txt

Makefile                              # or some other root configure script

README.md

src/                                  # Each subdirectory contains a collection of classes
  core/
    stagblgrid/
        stagblgrid.h
        stagblgrid.c
        impls/
          stagblgridpetsc.c
          stagblgridbasic.c
    stagblarray/
        stagblarray.h
        stagblarray.c
        impls/
          stagblarraypetsc.c
          stagblarraybasic.c
  io/
  stokes/
  particles/

tests/
