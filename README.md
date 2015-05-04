# specgram

make spectrogram from a sac file

__Dependent libaries:__

  netcdf, fftw


__Programs:__

  1. src/specgram.c: make S-transform of a given sac file.

  2. aux/specgram\_stack.c: auxiliary program to stack several output nc files from
      specgram


__How to use:__

  execute binaries in bin/ without any argument will print out help document.


__How to compile:__

  excute ./compile.sh 

  You need to modify compile.sh by adding netcdf and fftw paths 
    in to C_INCLUDE_PATH and LIBRARY_PATH, if they are not in the path already.
