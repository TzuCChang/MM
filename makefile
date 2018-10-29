FC90 = ifort -traceback -O3 -i8 -qopenmp -parallel -threads -fpp
FCFLAGS=   -I${MKLROOT}/include -mkl=parallel 
MKLROOT = /home/pec1065/intel/compilers_and_libraries_2016.1.150/linux/mkl
LIBS =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -ldl 
TARGETS= clean ShearCell
OBJSOC = track.o segment_prop_calc.o motion2.o M_SOLVERS.o ini_positioning.o EXCL_VOL_WALLS.o EXCL_VOL.o data_output.o bending.o break.o tipos.o set_up_hinges_and_grid.o
OBJMOD = $(wildcard *.mod)

ShearCell: $(OBJSOC) main.f90 
		$(FC90) $(FCFLAGS) -o ShearCell main.f90  $(LIBS) $(OBJSOC)
track.o: tipos.o track.f90
		$(FC90) -c track.f90

set_up_hinges_and_grid.o: data_structures.mod excl_volume.mod set_up_hinges_and_grid.f90
		$(FC90) -c set_up_hinges_and_grid.f90
		
segment_prop_calc.o: data_structures.mod segment_prop_calc.f90
		$(FC90) -c segment_prop_calc.f90
		
motion2.o: data_structures.mod MATRIX_SOLVERS.mod motion2.f90
		$(FC90) -c motion2.f90

MATRIX_SOLVERS.mod: M_SOLVERS.f90
		$(FC90) $(FCFLAGS)  -c  M_SOLVERS.f90 $(FCFLAGS) $(LIBS) 

ini_positioning.o: data_structures.mod ini_positioning.f90
		$(FC90) -c ini_positioning.f90 

EXCL_VOL_WALLS.o: data_structures.mod EXCL_VOL.o bending.o EXCL_VOL_WALLS.f90
		$(FC90) -c EXCL_VOL_WALLS.f90 

EXCL_VOL.o: data_structures.mod bending.o EXCL_VOL.f90
		$(FC90) -c EXCL_VOL.f90 

data_output.o: data_structures.mod data_output.f90
		$(FC90) -c data_output.f90 

bending.o: data_structures.mod bending.f90
		$(FC90) -c bending.f90 
		
break.o: data_structures.mod break.f90
		$(FC90) -c break.f90 
		
tipos.o: tipos.f90
		$(FC90) -c tipos.f90
		
clean: 
	rm -f *.o *.mod

