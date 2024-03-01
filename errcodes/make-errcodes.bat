gfortran -c program_constants0.f90
gfortran -c error_handler0.f90
gfortran -c dynalloc0.f90
gfortran -c tabulate0.f90
gfortran errcodes.f90 program_constants0.o error_handler0.o dynalloc0.o tabulate0.o
rm *.o *.mod
