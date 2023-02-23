# Compiler and flags

.SUFFIXES: .o .f .f90 .mod


VPATH = src

#compiler
F90 = f95
F77 = $(F90)

#flags
#F90FLAGS =  -xW -O4 #-check bounds 
#F77FLAGS = $(F90FLAGS)
###GRIMESBY###
#LINKFLAGS = -static -lsvml  #-lgfortran  -lguide -lpthread #-check bounds
###HELIOS###
#LINKFLAGS = -static -lsvml #-lacml -L/local/acml3.6.0/ifort64/lib
###KARI###
#LINKFLAGS =  -lsvml  -L/opt/acml3.6.0/ifort64/lib -lacml -lacml_mv -i_dynamic -static
TARGET = atom

#make objects from the source 
%.o:%.f90
	$(F90) -c  $(F90FLAGS) $<
%.o:%.F90
	$(F90) -c  $(F90FLAGS) $<

default: $(TARGET)

#main object
PROGRAM = atom.o

MODULES = constants.o globals.o initmod.o output.o \
           pxcfuncs3d.o dfsubs.o solvers.o 


# module dependencies
globals.o: constants.o
pxcfuncs3d.o: constants.o 
output.o: globals.o constants.o pxcfuncs.o
initmod.o: constants.o globals.o output.o
dfsubs.o: constants.o globals.o pxcfuncs3d.o output.o
solvers.o: constants.o globals.o dfsubs.o output.o

#program dependencies
$(PROGRAM): $(MODULES)


$(TARGET): $(MODULES) $(PROGRAM)
	 $(F90) $(MODULES) $(PROGRAM) -o  $(TARGET) $(LINKFLAGS)


clean:
	rm *.o *.mod

realclean:
	rm atom

