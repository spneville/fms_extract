########################################################################
# Compiler flags
########################################################################
FLAGS=-llapack \
      -lblas

########################################################################
# Modules
########################################################################
MODPATH=./lib
MODS=$(MODPATH)/iomod.f90 \
     $(MODPATH)/errormod.f90 \
     $(MODPATH)/parsemod.f90 \
     $(MODPATH)/sysdef.f90 \
     $(MODPATH)/trajdef.f90 \
     $(MODPATH)/expec.f90 \
     $(MODPATH)/intcoo.f90 \
     $(MODPATH)/gausstools.f90 \
     $(MODPATH)/density.f90 \
     $(MODPATH)/adpop.f90 \
     $(MODPATH)/dysonmod.f90 \
     $(MODPATH)/dysonprep.f90 \
     $(MODPATH)/cirmsd.f90 \
     $(MODPATH)/trpes.f90

########################################################################
# Compilation
########################################################################
fms_extract:	fms_extract.F90 $(MODS)
	gfortran -o fms_extract $(MODS) fms_extract.F90 $(FLAGS)

clean: 
	-rm *.mod fms_extract
