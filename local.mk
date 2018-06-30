# This file defines the top of the source tree.
# Other directories deeper in the tree define files also
# include local.mk files to define which source files to add to which targets,
# and which subdirectories to descend into

include $(call incsubdirs,src)
