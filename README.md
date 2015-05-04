# specgram

Dependent libaries:

  netcdf, fftw


Programs:

  1. src/specgram.c: make S-transform of a given sac file.

  2. aux/specgram_stack.c: auxiliary program to stack several output nc files from
      specgram


How to use:

  execute binaries in bin/ without any argument will print out help document.


How to compile:

  excute ./compile.sh 

#  You need to modify compile.sh by adding netcdf and fftw paths 
#    in to C_INCLUDE_PATH and LIBRARY_PATH, if they are not in the path already.
