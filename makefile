#F90 = mpif90 -i4 -real-size 32 -O2 -mcmodel=medium
#F90 = mpifort -mcmodel=medium -i4 -real-size 32 -O2
#Fone = ifort -i4 -real-size 32 -O2 -mcmodel=medium -openmp

#F90 = mpif90 -ffree-line-length-none -O2 -mcmodel=medium -g -pedantic -fbacktrace -fbounds-check -ffree-form -frecursive #-i4 #-real-size 32 -O2 -mmacosx-version-min=10.6

F90 =  mpif90 -i4 -real-size 32 -O2 -mcmodel=medium -f90=ifort# -O2 -mmacosx-version-min=10.6
#Fone = ifort -i4 -real-size 32 -O2 -mcmodel=medium -openmp
#F90 = mpiifort -march=alderlake -O3 -mcmodel=medium -ffree-line-length-none

FILES = dimensions.f90 mult_proc.f90 grid.f90 boundary.f90 var_arrays.f90 inputs.f90 misc.f90  grid_interp.f90 gutsf.f90 gutsp.f90  initial.f90 part_init.f90 maind.f90
DEBUG = -g -pedantic -fbacktrace -fbounds-check -fstack-check
FILESO = dimensions.f90 boundary.f90 grid_interp.f90 mult_proc.f90 var_arrays.f90 inputs.f90 grid.f90 initial.f90 gutsf.f90 misc.f90 part_init.f90 gutsp.f90 maind.f90

INCLUDE = dimensions.o inputs.o grid.o mult_proc.o var_arrays.o
INCLUDE2 = dimensions.o inputs.o grid.o mult_proc.o boundary.o var_arrays.o
INCLUDE3 = dimensions.o inputs.o grid.o mult_proc.o boundary.o misc.o grid_interp.o gutsp.o gutsf.o var_arrays.o
INCLUDE4 = dimensions.o inputs.o grid.o mult_proc.o boundary.o misc.o grid_interp.o gutsp.o gutsf.o initial.o part_init.o var_arrays.o
OBJECTS = dimensions.o inputs.o grid.o mult_proc.o boundary.o misc.o grid_interp.o gutsp.o gutsf.o initial.o part_init.o var_arrays.o maind.o
MODS = dimensions.mod mult_proc.mod var_array.mod inputs.mod grid.mod initial.mod gutsf.mod misc.mod boundary.mod part_init.mod grid_interp.mod gutsp.mod

hybrid: $(OBJECTS) 
	$(FCFLAGS) $(F90) -o hybrid $(OBJECTS) 

debug: $(OBJECTS) 
	$(F90) -o hybrid $(OBJECTS) $(DEBUG)

clean:
	rm *.o hybrid *.out *.mod disp hybrid_d

dimensions.o:dimensions.f90;$(F90) -c dimensions.f90
mult_proc.o:mult_proc.f90;$(F90) -c mult_proc.f90
grid.o:grid.f90;$(F90) -c grid.f90
inputs.o:inputs.f90 dimensions.o mult_proc.o grid.o var_arrays.o;$(F90) -c inputs.f90
boundary.o:boundary.f90 $(INCLUDE);$(F90) -c boundary.f90
misc.o:misc.f90 $(INCLUDE) boundary.o grid_interp.o;$(F90) -c misc.f90
grid_interp.o:grid_interp.f90 $(INCLUDE2);$(F90) -c grid_interp.f90
gutsf.o:gutsf.f90 $(INCLUDE2) grid_interp.o;$(F90) -c gutsf.f90
gutsp.o:gutsp.f90 $(INCLUDE2) grid_interp.o;$(F90) -c gutsp.f90
initial.o:initial.f90 inputs.o grid.o mult_proc.o;$(F90) -c initial.f90
part_init.o:part_init.f90 $(INCLUDE3);$(F90) -c part_init.f90
var_arrays.o:var_arrays.f90;$(F90) -c var_arrays.f90
maind.o:maind.f90 $(INCLUDE4);$(F90) -c maind.f90

dimensions.mod:dimensions.f90 $(INCLUDE);$(F90) -c dimensions.f90
mult_proc.mod:mult_proc.f90 $(INCLUDE);$(F90) -c mult_proc.f90
grid.mod:grid.f90 $(INCLUDE);$(F90) -c grid.f90
grid_interp.mod:grid_interp.f90 $(INCLUDE);$(F90) -c grid_interp.f90
gutsf.mod:gutsf.f90 $(INCLUDE);$(F90) -c gutsf.f90
gutsp.mod:gutsp.f90 $(INCLUDE);$(F90) -c gutsp.f90
misc.mod:misc.f90 $(INCLUDE);$(F90) -c misc.f90
boundary.mod:boundary.f90 $(INCLUDE);$(F90) -c boundary.f90
part_init.mod:part_init.f90 $(INCLUDE);$(F90) -c part_init.f90
initial.mod:initial.f90 $(INCLUDE);$(F90) -c initial.f90
inputs.mod:inputs.f90 $(INCLUDE);$(F90) -c inputs.f90
var_arrays.mod:var_arrays.f90 $(INCLUDE);$(F90) -c var_arrays.f90
