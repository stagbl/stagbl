Toy 3
-----

Similar to toy 2, but 
 - Relies on a development branch of PETSc, which includes DMStag : psanan/dmstag at bitbucket.org/psanan/petsc-private/
 - Relying on this same branch, interacts with DMSwarm to passively advect particles

This can and probably should be extended to solve problem 8.4 in Gerya's textbook.

Note: Expect a little bit of frustration loading XDMF files with paraview. Loading with the XDMF reader (not the XDMF 3 reader) seems to work on Ubuntu 16.04 with Paraview 5.4.1. (Note further that installing paraview from APT didn't work - I had to use a download from the Paraview website)
