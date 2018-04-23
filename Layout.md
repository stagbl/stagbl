StagBL Layout Proposal (WIP)
----------------------

A proposal for the layout of the StagBL project
Guidelines: keep it lean, keep it self-documenting as much as possible

Root :
  README.md                          # Contains quickstart w/ pictures
  Makefile [or some other config..]
  documentation
  examples
  src
  include
  developer


developer/                            # Everything that only the developer would care about (toys, util scripts, etc.)
  toys/

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
