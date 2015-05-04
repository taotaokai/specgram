#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "par.h"
#include "dbg.h"
#include "sac.h"
#include "sacutil.h"

#define PI M_PI

/* self documentation */
char *sdoc[] = {
"NAME",
"    specgram_st - make spectrogram by S-transform",
"",
"SYNOPSIS",
"    specgram_st in= out= dt=0.1 a=2*PI/log(10) freqs=1,100,1",
"",
"DESCRIPTION",
"    in=  input sac file name ",
"    out=  output file name",
"    dt=0.1  output time sampling interval(second)",
"    a=2*PI/log(10)  time window parameter T(w,a) (see below S-transfrom)",
"    freqs=1,100,1  output frequency samples(Hz) begin,end,interval",
"",
"S-transform",
"    The Stockwell transform of f(t) is defined as ",
"      Sf(t,w) = Int[f(s) * g(t-s;T(w,a)) * exp(-I*w*s), {s,-inf,inf}], ",
"      where g(t;T(w,a)) = exp(-(t/T(w,a))^2), and ",
"      , and T is charateristic time window length depending on frequency.",
"    In frequency domian S-transform is: ",
"      Sf(w',w) = f(w'+w) * G(w';T(w,a), ",
"      where G(w';T) = |T| * exp(-(w'*T)^2)",
"",
"Time window function",
"    In frequency domain: G(w';T(fc,a)) = T*exp(-(w'*T)^2), ",
"      where fc is the central frequency, and T(fc,a) = a/fc. ",
"    Then Dt*Dw ~ 1 => fc[i] = exp(2*PI/a*i) would ",
"      give the optimal time-frequency distribution. ",
"    If choose a=2*PI/log(10)~2.73, then fc[i] = 10^i, ",
"      i.e. a logarithm frequency axis. ", 
"",
"About output sampling interval",
"    Time window length varies from a/fc to dt with increasing frequency, ",
"    use dt instead of a*T(wc) to account for output sampling interval, ",
"    which limits the minimum time window length.",
"",
NULL};

/**** shift fft result F to center zero frequency ****/
int fftshift(fftw_complex *F, int nf)
{
  int i, n_nyq;
  fftw_complex *dump=NULL;

  check_mem(dump = malloc(sizeof(fftw_complex)*nf));
  for (i=0;i<nf;i++) dump[i] = F[i];

  n_nyq = (nf+1)/2;

  for (i=0; i<n_nyq; i++) F[i+nf-n_nyq] = dump[i];
  for (i=n_nyq; i<nf; i++) F[i-n_nyq] = dump[i];

  fftw_free(dump);
  return 0;

error:
  fftw_free(dump);
  return -1;
}

/**** shift F back to fft order ****/
int ifftshift(fftw_complex *F, int nf)
{
  int i, n_nyq;
  fftw_complex *dump=NULL;

  check_mem(dump = malloc(sizeof(fftw_complex)*nf));
  for (i=0;i<nf;i++) dump[i] = F[i];

  n_nyq = (nf+1)/2;

  for (i=0; i<n_nyq; i++) F[i] = dump[i+nf-n_nyq];
  for (i=n_nyq; i<nf; i++) F[i] = dump[i-n_nyq];

  fftw_free(dump);
  return 0;

error:
  fftw_free(dump);
  return -1;
}

/**** main program ****/
int main(int argc, char *argv[])
{
  /**** define variables ****/

  /* command line arguments */
  char *sacin=NULL;
  char *outfile=NULL;
  double dt1, gaussa;
  //double freqs[3];
  int nfreq;

  /* local variables */
  //sac trace
  sac *tr=NULL;
  //input
  int n0, n0_nyq;
  double dt0, dw0;
  double *w0=NULL;
  fftw_complex *x0=NULL;
  fftw_complex *F0=NULL;
  //output
  int n1, n1_nyq;
  double dw1;
  double *w1=NULL;
  fftw_complex *x1=NULL;
  fftw_complex *F1=NULL;
  //fft
  fftw_plan p_fft=NULL, p_bfft=NULL;
  //other
  int i, idx;
  double fc, wc, Twin;
  double base, exp1, exp0, dexp;
  int iexp;
  FILE * fp_out=NULL;

  /**** get command line arguments ****/
  //debug("get command line arguments.");

  initargs(argc, argv);
  requestdoc(1);
  //IO files
  if (!getparstring("in",&sacin)) sentinel("must specify in=");
  if (!getparstring("out",&outfile)) sentinel("must specify out=");
  //S-transform control
  if (!getpardouble("a",&gaussa)) gaussa = 2*PI/log(10);
  //output
  if (!getpardouble("dt",&dt1)) dt1 = 0.1;
  if (!getparint("nfreq",&nfreq)) nfreq = 100;
  //if (countparval("freqs")==3) {
  //  getpardouble("freqs",freqs);
  //} else {
  //  freqs[0]=1; freqs[1]=100; freqs[2]=1;
  //}
  //open output file
  fp_out = fopen(outfile,"w");

  /**** read input sac file ****/
  //debug("read input sac file.");

  check_mem(tr = sacnewn(1));
  sacread(tr, sacin);

  /**** fft of input data ****/
  dt0 = tr->delta;
  n0 = tr->npts;
  dw0 = 2.0*PI/(n0*dt0);

  //fprintf(fp_out,"# dt0=%f n0=%d dw0=%f\n",dt0,n0,dw0);

  check_mem(x0 = fftw_malloc(sizeof(fftw_complex) * n0));
  check_mem(F0 = fftw_malloc(sizeof(fftw_complex) * n0));
  for (i=0; i<n0; i++) x0[i] = tr->data[i];
  p_fft = fftw_plan_dft_1d(n0, x0, F0, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p_fft);

  //fprintf(fp_out,"# x0\n");
  //for (i=0; i<n0; i++)
  //  fprintf(fp_out,"%e+%ei ", creal(x0[i]), cimag(x0[i]));
  //fprintf(fp_out,"\n");

  //fprintf(fp_out,"# F0\n");
  //for (i=0; i<n0; i++)
  //  fprintf(fp_out,"%e %e \n", creal(F0[i]), cimag(F0[i]));
  //fprintf(fp_out,"\n");

  //w0: frequency samples in natural order 
  check_mem(w0 = malloc(sizeof(double) * n0));
  n0_nyq = (n0+1)/2;
  for (i=0; i<n0; i++) w0[i] = (i+n0_nyq-n0)*dw0;

  //fprintf(fp_out,"# w0\n");
  //for (i=0; i<n0; i++)
  //  fprintf(fp_out,"%e \n", w0[i]);
  //fprintf(fp_out,"\n");

  //fftshift: shift F0 to center zero frequency
  fftshift(F0, n0);
  //for (i=0; i<n0_nyq; i++) {
  //  dumpc = F0[i];
  //  F0[i] = F0[i+n0_nyq];
  //  F0[i+n0_nyq] = dumpc;
  //}

  //for (i=0; i<n0; i++)
  //  fprintf(fp_out,"%e %e \n", creal(F0[i]), cimag(F0[i]));

  /**** prepare backward fft of output data ****/
  n1 = n0*dt0/dt1;
  dw1 = 2.0*PI/(n1*dt1);

  //w1: frequency samples in fft order
  check_mem(w1 = malloc(sizeof(double) * n1));
  n1_nyq = (n1+1)/2;
  for (i=0; i<n1; i++) w1[i] = (i+n1_nyq-n1)*dw1;

  //for (i=0; i<n1; i++)
  //  fprintf(fp_out,"%e \n", w1[i]);

  check_mem(x1 = fftw_malloc(sizeof(fftw_complex) * n1));
  check_mem(F1 = fftw_malloc(sizeof(fftw_complex) * n1));
  p_bfft = fftw_plan_dft_1d(n1, F1, x1, FFTW_BACKWARD, FFTW_ESTIMATE);

  /**** Stockwell transform at each frequency ****/

  //output frequency range
  base = 2*PI/gaussa;
  fprintf(stderr,"log frequency axis base = %f\n", exp(base));
  exp0 = log(1/dt0/n0)/base;
  exp1 = log(1/dt0/2)/base;
  dexp = (exp1-exp0)/nfreq;

  //for (fc=freqs[0]; fc<=freqs[1]; fc+=freqs[2]) {
  //for (fc=2.0; fc<=2.1; fc+=1.0) {
  for (iexp=0; iexp<nfreq; iexp++) {
    //central frequency
    fc = exp(base*(exp0+iexp*dexp));
    //debug("fc = %f", fc);
    wc = 2*PI*fc;

    //interpolate F0(w0) -> F1(w1) = F0(w1+wc)
    for (i=0; i<n1; i++) {
      //locate w0[idx] <= w1[i] <= w0[idx+1]
      idx = (w1[i] + wc - w0[0])/dw0;
      if ((idx < 0) || (idx >=n0-1)) {
        F1[i] = 0.0; //out of data frequency range
      } else {
        //linear interpolation
        F1[i] = (F0[idx]*(w0[idx+1] - w1[i] - wc) 
               + F0[idx+1]*(w1[i] + wc - w0[idx]))/dw0;
      }
      //fprintf(fp_out,"%d %e %e \n", idx, creal(F1[i]), cimag(F1[i]));
    }

    //Multiply window function G(w';T(wc,a)) = T*exp(-(w'*T)^2).
    if (fc == 0.0) {
      for (i=1; i<n1; i++) F1[i] = 0.0;
    } else {
      //time window length varies from dt1 to a*T(wc) with decreasing frequency
      //use dt1 instead of a*T(wc) to account for output sampling interval. 
      Twin = gaussa/fc;
      //fprintf(stderr,"Twin=%f\n", Twin);
      Twin = sqrt(1.0 + pow(Twin/dt1,2))*dt1;
      //fprintf(stderr,"Twin=%f\n", Twin);
      for (i=0; i<n1; i++) {
        F1[i] *= Twin * exp(-pow(w1[i]*Twin,2));
        //fprintf(fp_out,"%e %e \n", creal(F1[i]), cimag(F1[i]));
      }
    }

    //shift F1 back to put zero frequency back at first element
    ifftshift(F1,n1);

    //for (i=0; i<n1; i++) fprintf(fp_out,"%e %e \n", creal(F1[i]), cimag(F1[i]));

    //ifft(F1)
    fftw_execute(p_bfft);

    //output abs(x1) = abs(Sf(t,fc))
    fprintf(fp_out, "%f ", fc);
    for (i=0; i<n1; i++) fprintf(fp_out, "%e ", cabs(x1[i]));
    fprintf(fp_out, "\n");
  }

  /**** clean up ****/
  //sacfreen(tr,1);
  //fftw_destroy_plan(p_fft); fftw_destroy_plan(p_bfft);
  //free(w0); fftw_free(x0); fftw_free(F0);
  //free(w1); fftw_free(x1); fftw_free(F1);
  //fclose(fp_out);

error:
  sacfreen(tr,1);
  fftw_destroy_plan(p_fft); fftw_destroy_plan(p_bfft);
  free(w0); fftw_free(x0); fftw_free(F0);
  free(w1); fftw_free(x1); fftw_free(F1);
  if (fp_out) fclose(fp_out);
}
