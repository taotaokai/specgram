#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <netcdf.h>

#include "par.h"
#include "dbg.h"
#include "sac.h"
#include "sacutil.h"

/***** self documentation *****/
char *sdoc[] = {
"NAME",
"    specgram - make spectrogram using S-transform",
"",
"SYNOPSIS",
"    specgram in= out=out.nc dt=0.1 fmin=0.1 fmax=10 nf=10 twina=1 sinca=20",
"",
"DESCRIPTION",
"    in=  input sac file name ",
"    out=out.nc  output file name (format: netcdf)",
"    dt=0.1  output time sampling interval (second)",
"    fmin=0.1  minimum frequency (Hz)",
"    fmax=10  maximum frequency (Hz)",
"    nf=10  number of frequency samples",
"    twina=1  scaling factor of time window width T(f,a)=a/f",
"    sinca=20  sinc window width for lanczos interpolation",
"",
"OUTPUT FILE CONTENT",
"    the output netcdf file consists of ",
"    - axis: time(nt), log10_frequency(nf)",
"    - data: ",
"      stf_abs(nf,nt) = |X(t,f)|, absolute value of S-transform",
"",
"COMMETS",
" 1) S-transform, moving time window T(f,a):",
"    X(t,f) = Int[x(s) * g(t-s;T(f,a)) * exp(-I*2*PI*s), {s,-inf,inf}], ",
"      where g(t;T(f,a)) = exp(-2*(t/T(f,a))^2)*sqrt(2)/sqrt(PI)/abs(T), and ",
"      , and T is approximatly the window width(|g|>0.6).",
"    The time window width scales with frequency as T(f,a) = a/f",
"    In frequency domian X(w,f) = Fx(w+2*PI*f) * G(w;T(f,a)), ",
"      where G(w;T) = exp(-(w*T)^2/2), ",
"      Fx(w+2*PI*f) is futher interpolated from discrete frequency samples of ",
"      fft(x(t)). see {Lanczos sinc interpolation} ",
"",
" 2) About further tuning the time window width T(f,a):",
"    Time window width varies from T = a/f to dt with increasing frequency, ",
"    use dt instead of a/f at high frequency end to account for output time sampling interval, ",
"    which limits the minimum time resolution.",
"",
" 3) Output time-frequency samples: ",
"    time: same begin time as input sac file, sampling interval dt, ",
"          the end time is closest to the end time as the input sac file.",
"    frequency: eqidistant log10 axis from fmin to fmax, with nf samples including",
"               fmin and fmax",
"",
" 4) Lanczos sinc interpolation(reconstruction):",
"    windowed sinc: L(x) = sinc(x)*sinc(x/a), |x|<=a ",
"                        = 0, |x|>a ",
"    S(x) = sum(L(x-i) * S[i]) ",
"",
" 5) Relation between S-transfrom and bandpass filtered signal:",
"    Gaussian bandpass filter of central frequency f is defined as G(w;T(f)), which",
"    is also Fourier transform of the time window used in S-transform. ",
"    Let Gx(t,f) be the bandpass filtered x(t), and X(t,f) the S-transform of x(t), then ",
"    Gx(t,f) = Re[X(t,f) * exp(I*2*PI*f*t)] = |X(t,f)|*cos(phi(t)+2*PI*f*t), ",
"    so |X(t,f)| and phi(t) correspond to the instantaneous amplitude and phase of ",
"    the bandpass filtered signal, respectively.",
NULL};


/***** subroutines *****/

int fftshift(fftw_complex *F, int nf)
/* shift fft result F to center zero frequency */
{
  int i, n_nyq;
  fftw_complex *dump=NULL;

  check_mem(dump = fftw_malloc(sizeof(fftw_complex) * nf));
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


int ifftshift(fftw_complex *F, int nf)
/* shift F back to fft order */
{
  int i, n_nyq;
  fftw_complex *dump=NULL;

  check_mem(dump = fftw_malloc(sizeof(fftw_complex) * nf));
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

/* macro: double dsinc(x) */
#define PI M_PI 
/* note the usage of dsincarg=PI*(a) in macro DSINC(a) */
static double dsincarg;
#define DSINC(a) ((dsincarg=PI*(a)) == 0.0 ? 1 : sin(dsincarg)/dsincarg )

int lanczos_interp(
  fftw_complex *F,  double w0,  double dw,  int nw, 
  fftw_complex *Fi, double wi0, double dwi, int nwi, int na)
/* windowed sinc interpolation of spectrum.
 * interpolate F[w0+i*dw],i=0..nw-1 onto wi0+i*dwi, i=0..nwi-1
 * na: window half width in points */
{
  fftw_complex sum;
  double d0, nf, resp;
  int i, n0, n1, nn;

  d0 = wi0 - w0;

  /* loop each interpolating point */
  for (i=0; i<nwi; i++) {

    /* get window index range: n0..n1-1 */
    nf = (d0 + i*dwi)/dw;
    n0 = nf - na;
    n1 = nf + na + 1;
    if (n0 < 0) n0 = 0;
    if (n1 > nw) n1 = nw;

    /* loop each point in sinc window */
    sum = 0.0;
    for (nn=n0; nn<n1; nn++) {
      /* resp: band-limited response at position nf from sample nn
       * DSINC is executed with a sequence point in between
       * to avoid undefined operation on global variable dsincarg */
      resp = DSINC(nn-nf);
      sum += resp * DSINC((nn-nf)/na) * F[nn]; 
    }

    Fi[i] = sum;
  }

  return 0;
}


/***** main program  *****/

/* Handle errors for netcdf API */
static int nc_retval;
#define nc_check(A) if((nc_retval=(A))!= NC_NOERR) {fprintf(stderr,  \
    "[ERROR_NC] (%s:%d: errno: %s) \n", __FILE__, __LINE__, \
    nc_strerror(nc_retval)); goto error;}

/* netcdf constants */
#define NDIMS 2
#define UNITS "units"

int main(int argc, char *argv[])
{
  /*** declaration of variables ***/

  /* shared */
  int i;
  fftw_plan p_fft=NULL, p_bfft=NULL;

  /* command line arguments */
  char *sacin=NULL;
  char *outfile=NULL;
  double deltat, twina, fmin, fmax;
  int nfreq, sinca;
  
  /* input data */
  sac *tr=NULL;
  struct tm sacgmt;
  int n0, n0_nyq;
  double dt0, dw0, f0_nyq;
  double *w0=NULL;
  fftw_complex *x0=NULL;
  fftw_complex *F0=NULL;

  /* S-transform */
  int n1, n1_nyq;
  double dt1, dw1;
  double *w1=NULL;
  fftw_complex *x1=NULL;
  fftw_complex *F1=NULL;

  /* freuqency samples */
  int ifreq;
  double fc, wc, winlen;
  double log10_fmin, log10_fmax, dlog10_freq; 
  double *log10_freqs=NULL;

  /* output netcdf */
  int ncid=-1;
  int dimids[NDIMS], time_dimid, freq_dimid;
  int time_varid, freq_varid; /* coordinate variables */
  int stf_abs_varid; 
  /* int stf_arg_varid; // data variables */
  int ntime; double *times=NULL;
  float *stf_abs=NULL; /* magnitude of X(t,f)*/
  /* float *stf_arg=NULL; // phase angle of X(t,f)*/
  const char time_name[]="time";
  /* time_units format: "seconds since yyyy-mm-dd hh:mm:ss.sss" */
  char time_units[38]="\0"; /* to be determined from sac file */
  const char freq_name[]="log10_frequency";
  const char freq_units[]="Hz";
  const char stf_abs_name[]="stf_abs"; /* magnitude */
  const char stf_abs_units[]="unknown";
  /*
  const char stf_arg_name[]="stf_arg"; // phase angle //
  const char stf_arg_units[]="rad";
  */

  /*** get command line arguments ***/

  initargs(argc, argv);
  requestdoc(1);
  /* IO files */
  if (!getparstring("in",&sacin)) sentinel("must specify in=");
  if (!getparstring("out",&outfile)) outfile = "out.nc";
  /* output time-frequency samples */
  if (!getpardouble("dt",&deltat)) deltat = 0.1;
  if (!getpardouble("fmin",&fmin)) fmin = 1;
  if (!getpardouble("fmax",&fmax)) fmax = 10;
  if (!getparint("nf",&nfreq)) nfreq = 10;
  /* S-transform time window length  */
  if (!getpardouble("twina",&twina)) twina = 1;
  /* sinc window length for interpolation */
  if (!getparint("sinca",&sinca)) sinca = 20;

  /**** check command line arguments ****/

  check(deltat > 0, "dt(%f) must be greater than 0", deltat);
  if (nfreq < 1) {
    log_warn("nf(%d) should .ge. 1; use 1 instead", nfreq);
    nfreq = 1;
  }
  check(fmin > 0, "fmin(%f) must be greater than 0", fmin);
  check(twina > 0, "twina(%f) must be greater than 0", twina);
  check(sinca >= 5, "sinca(%d) should .ge. 5", sinca);

  fprintf(stderr, "# parameters used in this run:\n");
  fprintf(stderr, "# in=%s out=%s dt=%f fmin=%f fmax=%f nf=%d twina=%f " \
    "sinca=%d \n", sacin, outfile, deltat, fmin, fmax, nfreq, twina, sinca);

  /**** read input sac file ****/

  check_mem(tr = sacnewn(1));
  check(sacread(tr, sacin)==0, "ERROR: failed to read sac file: %s", sacin);

  /*** fft input data (suffix: 0) ***/

  dt0 = tr->delta;
  n0 = tr->npts;
  dw0 = 2.0*PI/(n0*dt0);
  f0_nyq = 0.5/dt0;

  if (fmax > f0_nyq) 
    log_warn("fmax(%f) is greater than nyquist frequency(%f) of input signal.",
        fmax, f0_nyq);

  fprintf(stderr,"# input sac: dt=%f npts=%d \n",dt0,n0);

  check_mem(x0 = fftw_malloc(sizeof(fftw_complex) * n0));
  check_mem(F0 = fftw_malloc(sizeof(fftw_complex) * n0));
  for (i=0; i<n0; i++) x0[i] = tr->data[i];
  p_fft = fftw_plan_dft_1d(n0, x0, F0, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p_fft);

  /* w0: frequency samples from minus to positive nyquist frequency */
  check_mem(w0 = malloc(sizeof(double) * n0));
  n0_nyq = (n0+1)/2;
  for (i=0; i<n0; i++) w0[i] = (i+n0_nyq-n0)*dw0;

  /* fftshift: shift F0 to center zero frequency */
  fftshift(F0, n0);

  /*** prepare backward fft for S-transform ***/
  dt1 = deltat;
  n1 = n0*dt0/dt1;
  dw1 = 2.0*PI/(n1*dt1);

  fprintf(stderr,"# output: dt=%f npts=%d\n",dt1,n1);

  /* w1: frequency samples from minus to positive nyquist frequency */
  check_mem(w1 = malloc(sizeof(double) * n1));
  n1_nyq = (n1+1)/2;
  for (i=0; i<n1; i++) w1[i] = (i+n1_nyq-n1)*dw1;

  check_mem(x1 = fftw_malloc(sizeof(fftw_complex) * n1));
  check_mem(F1 = fftw_malloc(sizeof(fftw_complex) * n1));
  p_bfft = fftw_plan_dft_1d(n1, F1, x1, FFTW_BACKWARD, FFTW_ESTIMATE);

  /*** define frequency samples in logarithmic scale ***/
  log10_fmin = log10(fmin);
  log10_fmax = log10(fmax);
  if (nfreq <= 1) {
    nfreq = 1;
    dlog10_freq = 0;
  } else {
    dlog10_freq = (log10_fmax - log10_fmin)/(nfreq-1);
  }
  check_mem(log10_freqs = malloc(sizeof(double)*nfreq));
  for (ifreq=0; ifreq<nfreq; ifreq++)
    log10_freqs[ifreq] = log10_fmin+ifreq*dlog10_freq;

  /*** S-transform ***/

  check_mem(stf_abs = malloc(sizeof(float) * n1 * nfreq));
  /* check_mem(stf_arg = malloc(sizeof(float) * n1 * nfreq)); */

  /* S-transfrom for each frequency */
  for (ifreq=0; ifreq<nfreq; ifreq++) {

    /* central frequency */
    fc = pow(10, log10_freqs[ifreq]);
    wc = 2*PI*fc;

    /* interpolate F0(w0) -> F1(w1) = F0(w1+wc) */
    lanczos_interp(F0, w0[0], dw0, n0, F1, w1[0]+wc, dw1, n1, sinca);

    /* Multiply window function G(w1;T(fc,a)) = exp(-(w1*T)^2/2) */
    winlen = twina/fc;
    winlen = sqrt(1.0 + pow(winlen/dt1,2))*dt1;
    for (i=0; i<n1; i++) {
      F1[i] *= exp(-pow(w1[i]*winlen,2)/2);
    }

    /* shift F1 back to put zero frequency back at first element */
    ifftshift(F1,n1);

    /* ifft(F1) */
    fftw_execute(p_bfft);

    /* correct the amplitude (1/npts)
     * because fftw_backward is the same as fft except an opposite sign */
    for (i=0; i<n1; i++) x1[i] = x1[i]/n1;

    /* convert x1(t,fc) to stf_abs(t,fc) and stf_arg(t,fc) */
    /* stf_abs,arg[ifreq][itime] */
    for (i=0; i<n1; i++) stf_abs[i + n1*ifreq] = cabs(x1[i]);
    /*for (i=0; i<n1; i++) stf_arg[i + n1*ifreq] = carg(x1[i]); */
  }


  /*** Output netcdf file ***/

  /* create the file. */
  nc_check(nc_create(outfile, NC_CLOBBER, &ncid));

  /* define dimensions */
  ntime = n1;
  nc_check(nc_def_dim(ncid, time_name, ntime, &time_dimid));
  nc_check(nc_def_dim(ncid, freq_name, nfreq, &freq_dimid));

  /* define coordinate variables */
  nc_check(nc_def_var(ncid, time_name, NC_DOUBLE, 1, &time_dimid, &time_varid));
  nc_check(nc_def_var(ncid, freq_name, NC_DOUBLE, 1, &freq_dimid, &freq_varid));

  /* define units attributes for coordinate vars */
  sac_get_tm(tr, &sacgmt);
  snprintf(time_units, 38, "seconds since %04d-%02d-%02d %02d:%02d:%02d.%03d", 
           sacgmt.tm_year+1900, sacgmt.tm_mon+1, sacgmt.tm_mday,
           sacgmt.tm_hour, sacgmt.tm_min, sacgmt.tm_sec, tr->nzmsec);
  nc_check(nc_put_att_text(ncid, time_varid, UNITS, strlen(time_units),
        time_units));
  nc_check(nc_put_att_text(ncid, freq_varid, UNITS, strlen(freq_units),
        freq_units));

  /* define data variables */
  dimids[1] = time_dimid;
  dimids[0] = freq_dimid;

  nc_check(nc_def_var(ncid, stf_abs_name, NC_FLOAT, NDIMS, dimids,
      &stf_abs_varid));
  /*
  nc_check(nc_def_var(ncid, stf_arg_name, NC_FLOAT, NDIMS, dimids,
      &stf_arg_varid));
  */

  /* define units attributes for data variables */
  nc_check(nc_put_att_text(ncid, stf_abs_varid, UNITS, strlen(stf_abs_units),
      stf_abs_units));
  /*
  nc_check(nc_put_att_text(ncid, stf_arg_varid, UNITS, strlen(stf_arg_units),
      stf_arg_units)); */

  /* end define mode */
  nc_check(nc_enddef(ncid));

  /* write coordinate variable data */
  check_mem(times = malloc(sizeof(double) * ntime));
  for (i=0; i<ntime; i++) times[i] = i*dt1 + tr->b; /* starts at sac b time */
  nc_check(nc_put_var_double(ncid, time_varid, times));
  nc_check(nc_put_var_double(ncid, freq_varid, log10_freqs));

  /* write data variable */
  nc_check(nc_put_var_float(ncid, stf_abs_varid, stf_abs));
  /* nc_check(nc_put_var_float(ncid, stf_arg_varid, stf_arg)); */

  /* close file */
  nc_check(nc_close(ncid));

  /*** finalize ***/
  sacfreen(tr,1);
  fftw_destroy_plan(p_fft); fftw_destroy_plan(p_bfft);
  free(w0); fftw_free(x0); fftw_free(F0);
  free(w1); fftw_free(x1); fftw_free(F1);
  free(times); free(log10_freqs);
  free(stf_abs); /* free(stf_arg); */
  return 0;

  /*** error handling ***/

error:
  sacfreen(tr,1);
  fftw_destroy_plan(p_fft); fftw_destroy_plan(p_bfft);
  free(w0); fftw_free(x0); fftw_free(F0);
  free(w1); fftw_free(x1); fftw_free(F1);
  free(times); free(log10_freqs);
  free(stf_abs); /* free(stf_arg); */
  nc_close(ncid);
  return -1;
}
