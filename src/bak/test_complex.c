#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

int main()
{

  fftw_complex a, c;
  double b;

  b = 3;
  a = 1 + I*2;

  c = b*a;

  fprintf(stderr, "b=%e \n",b);
  fprintf(stderr, "a=%e+I%e \n",creal(a),cimag(a));
  fprintf(stderr, "c=%e+I%e \n",creal(c),cimag(c));

  return 0;
}