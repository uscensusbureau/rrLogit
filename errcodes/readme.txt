The errorcodes program is a development tool that enables the
Fortran procedures called in an R package to pass text highly informative text
messages to the R user. The messages not only tell the user what error
occurred, but also the name of the procedure where it occurred.

How error handling works

The error_handler Fortran module for the R package, whose source code is in
    rrLogit/src/error_handler.f90
    rrLogit/src/icodes_1.h
    rrLogit/src/icodes_2.h
    rrLogit/src/icodes_3.h
manages error messages for all Fortran procedures. When a fatal error
occurs, the program does not crash, but gracefully exit and provide
detailed information. 

However, .Fortran function in R does not like character arguments. So
the text messages are converted to integer codes, which are passed to
R as integer arguments to .Fortran.

The text corresponding to the integer codes is stored in the three icodes_*.h
files. An R function defined in 
    rrLogit/R/msg.R
does the conversion back to text, using lookup tables stored in an
internal package data object named
    icodesError
Manually updating these header files and the R lookup tables would be
extremely tedious and difficult, so the errorcodes program does it
automatically. 

What errorcodes does

When executed, the errorcodes program scans the source code of every
Fortran module in the src directory, updating error_handler.f90 and
icodes_*.h. It creates a file named
   rrLogit/src/icodes.R
which contains R code for creating the lookup tables, which are then
saved as R objects in
   rrLogit/src/sysdata.rda

Scripts

This directory has two shell scripts. Please examine these scripts and
modify their paths if necessary before running. The script
    make-errcodes
compiles the errorcodes program from source. The script
    run-errcodes
runs the program to update the rrLogit Fortran source, and then creates
the lookup tables in sysdata.rda
