########################################################################
# Compiler flags
########################################################################
#F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -fopenmp -O2 -fbacktrace -fcheck=all

F90OPTS = -g

FLAGS=-llapack \
      -lblas \

#FLAGS=/usr/lib64/libblas.so.3 /usr/lib64/liblapack.so.3

########################################################################
# Modules
########################################################################
MODPATH=./lib
MODS=$(MODPATH)/iomod.f90 \
     $(MODPATH)/errormod.f90 \
     $(MODPATH)/parsemod.f90 \
     $(MODPATH)/mathlib.f90 \
     $(MODPATH)/sysdef.f90 \
     $(MODPATH)/kabschmod.f90 \
     $(MODPATH)/permutemod.f90 \
     $(MODPATH)/trajdef.f90 \
     $(MODPATH)/centdef.f90 \
     $(MODPATH)/expec.f90 \
     $(MODPATH)/intcoo.f90 \
     $(MODPATH)/gausstools.f90 \
     $(MODPATH)/density.f90 \
     $(MODPATH)/density_2d.f90 \
     $(MODPATH)/density_mom.f90 \
     $(MODPATH)/adpop.f90 \
     $(MODPATH)/dysonmod.f90 \
     $(MODPATH)/postprep.f90 \
     $(MODPATH)/dysonprep.f90 \
     $(MODPATH)/adcprep.f90 \
     $(MODPATH)/cirmsd.f90 \
     $(MODPATH)/projmod.f90 \
     $(MODPATH)/trpes.f90 \
     $(MODPATH)/mcspline.f90 \
     $(MODPATH)/adcspec.f90 \
     $(MODPATH)/tspsgmod.f90 \
     $(MODPATH)/gamessprep.f90 \
     $(MODPATH)/gamesstrxas.f90 \
     $(MODPATH)/dipoleprep.f90 \
     $(MODPATH)/dipole.f90

########################################################################
# Compilation
########################################################################
fms_extract:	fms_extract.f90 $(MODS)
	gfortran $(F90OPTS) -o fms_extract $(MODS) fms_extract.f90 $(FLAGS)
	-rm *.mod

clean: 
	-rm *.mod fms_extract
