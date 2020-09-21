libstagbl-y.c += $(call thisdir,\
  stagbl.c\
)

include $(call incsubdirs,control core stokes)
