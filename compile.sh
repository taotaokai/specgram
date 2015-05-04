#!/bin/bash

wkdir=$(pwd)

#incdir=include
#libdir=lib
#subdir=sub

## add search path into C_INCLUDE_PATH, LIBRARY_PATH
##NOTE: your need to add netcdf and fftw here, if thery are not 
##      in the search path already.
#export C_INCLUDE_PATH=$wkdir/include:${C_INCLUDE_PATH}
#export LIBRARY_PATH=$wkdir/lib:${LIBRARY_PATH}
#
## compile libraries in sub/
#for libname in par sac txt
#do
#  cd $wkdir/$subdir/lib${libname}/src
#  make clean
#  make
#  cp ../lib/*.a $wkdir/$libdir
#  cp ../include/*.h $wkdir/$incdir 
#done

# compile main program
cd $wkdir/src
make clean
make

# compile auxiliary program
cd $wkdir/aux
make clean
make
