#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <netcdf.h>

#include "dbg.h"
#include "par.h"
#include "txt.h"

/***** self documentation *****/
char *sdoc[] = {
"NAME",
"   specgram_stack - stack spectrograms",
"",
"SYNOPSIS",
"   specgram_stack in=nc.lst out=stack.nc nroot=1",
"",
"DESCRIPTION",
"   (string) in=  list of nc files to be stacked",
"   (string) out=stack.nc  output file name",
"   (int) nroot=2  n-th root stacking",
"",
"OUTPUT FILE CONTENT",
"   the output netcdf file consists of ",
"    - axis: time(nt), log10_frequency(nf)",
"    - data: ",
"      stf_abs: n-th root stack of stf_abs",
"      stf_dB: n-th root stack of dB(nf,nt)=10*log10(stf_abs.^2/P0), ",
"      where P0 = mean(stf_abs.^2)",
"",
"COMMENTS",
"   For all input files, only the axis names('time','log10_frequency') ",
"   and corresponding dimension size(nt,nf) are checked to be the same.",
"",
NULL};


/***** main program *****/

/* error handling for netcdf API */
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

  /* command line arguments */
  char *inlist=NULL;
  char *outfile=NULL;
  int nroot;

  /* common  */
  size_t itf, ntime, nfreq, ntf;
  float inv_nroot, P0, P0_dB;

  /* input list */
  char **filelist=NULL;
  size_t ifile, nfile=0;
  
  /* netcdf */
  int ncid=-1;

  /* netcdf: dimension */
  const char time_name[]="time";
  const char freq_name[]="log10_frequency";
  int time_dimid, freq_dimid, dimids[NDIMS];
  size_t time_dimlen, freq_dimlen;

  /* netcdf: coordinate variables */
  int time_varid, freq_varid;
  double *time_coord=NULL, *freq_coord=NULL;
  size_t time_units_len=0;
  char *time_units=NULL; /* to be determined from sac file */
  const char freq_units[]="Hz";

  /* netcdf: data variables */
  int stf_abs_varid, stf_dB_varid;
  const char stf_abs_name[]="stf_abs";
  const char stf_abs_units[]="unkown";
  float *stf_abs=NULL;

  /* netcdf: stacked data variables */
  float *stf_abs_nstk=NULL;
  const char stf_dB_name[]="stf_dB";
  const char stf_dB_units[]="dB";
  float *stf_dB_nstk=NULL;

  /*** get command line arguments ***/

  initargs(argc, argv);
  requestdoc(1);

  if (!getparstring("in",&inlist)) sentinel("must specify in=");
  if (!getparstring("out",&outfile)) outfile = "stack.nc";
  if (!getparint("nroot",&nroot)) nroot = 2;

  /* check command line arguments */
  if (nroot < 1) {
    log_warn("nroot(%d) should .ge. 1; use 2 instead", nroot);
    nroot = 2;
  }

  fprintf(stderr, "# parameters used in this run:\n");
  fprintf(stderr, "# in=%s out=%s nroot=%d\n" \
    ,inlist, outfile, nroot);

  /**** read input list ****/
  filelist = readlines(inlist, &nfile);
  check(nfile >= 1, "nothing read from the list: %s", inlist);

  /**** determine dimensions from the first file ****/
  nc_check(nc_open(filelist[0], NC_NOWRITE, &ncid));

  /* get dimension id */
  nc_check(nc_inq_dimid(ncid, time_name, &time_dimid));
  nc_check(nc_inq_dimid(ncid, freq_name, &freq_dimid));

  /* get dimension length */
  nc_check(nc_inq_dimlen(ncid, time_dimid, &ntime));
  nc_check(nc_inq_dimlen(ncid, freq_dimid, &nfreq));
    
  /* get coordinates */
  nc_check(nc_inq_varid(ncid, time_name, &time_varid));

  check_mem(time_coord = malloc(sizeof(double) * ntime));
  nc_check(nc_get_var_double(ncid, time_varid, time_coord));

  nc_check(nc_inq_attlen(ncid, time_varid, UNITS, &time_units_len));
  check_mem(time_units = calloc(sizeof(char), time_units_len));
  nc_check(nc_get_att_text(ncid, time_varid, UNITS, time_units));

  check_mem(freq_coord = malloc(sizeof(double) * nfreq));
  nc_check(nc_inq_varid(ncid, freq_name, &freq_varid));
  nc_check(nc_get_var_double(ncid, freq_varid, freq_coord));

  /* initialize data arrays */
  ntf = ntime * nfreq;
  check_mem(stf_abs = malloc(sizeof(float) * ntf));
  check_mem(stf_abs_nstk = malloc(sizeof(float) * ntf));
  check_mem(stf_dB_nstk = malloc(sizeof(float) * ntf));
  for (itf=0; itf<ntf; itf++) {
    stf_abs_nstk[itf] = 0.0;
    stf_dB_nstk[itf] = 0.0;
  }

  /*** stack all files ***/
  inv_nroot = 1.0/nroot;
  for (ifile=0; ifile<nfile; ifile++) {

    /* open input file. */
    nc_check(nc_open(filelist[ifile], NC_NOWRITE, &ncid));
  
    /* get dimension id */
    nc_check(nc_inq_dimid(ncid, time_name, &time_dimid));
    nc_check(nc_inq_dimid(ncid, freq_name, &freq_dimid));

    /* get dimension length */
    nc_check(nc_inq_dimlen(ncid, time_dimid, &time_dimlen));
    nc_check(nc_inq_dimlen(ncid, freq_dimid, &freq_dimlen));
      
    /* check dimension same as the first file */
    check(time_dimlen==ntime, "wrong size of time dimension[%zu] in %s"
        , time_dimlen, filelist[ifile]);
    check(freq_dimlen==nfreq, "wrong size of frequency dimension[%zu] in %s"
        , freq_dimlen, filelist[ifile]);

    /* get data variable */
    nc_check(nc_inq_varid(ncid, stf_abs_name, &stf_abs_varid));
    nc_check(nc_get_var_float(ncid, stf_abs_varid, stf_abs));

    /* close file */
    nc_check(nc_close(ncid));

    /* stack of stf_abs=|X(t,f)|: mean(sgn(xi)*|xi|^(1/nroot), i=0..nfile) */
    for (itf=0; itf<ntf; itf++) {
      stf_abs_nstk[itf] += copysignf( \
          powf(fabsf(stf_abs[itf]),inv_nroot)/nfile, stf_abs[itf]);
    }

    /* stack of stf_dB=10*log10(|X(t,f)|^2/P0) */
    P0 = 0.0;
    for (itf=0; itf<ntf; itf++)
        P0 += stf_abs[itf] * stf_abs[itf];
    P0 = P0/ntf;
    P0_dB = 10*log10f(P0);
    for (itf=0; itf<ntf; itf++) {
      P0 = 20*log10f(stf_abs[itf]) - P0_dB;
      stf_dB_nstk[itf] += copysign( \
          powf(fabsf(P0),inv_nroot)/nfile, P0);
    }

  }

  /* finish n-th root stack */
  for (itf=0; itf<ntf; itf++) {
    stf_abs_nstk[itf] = copysign( \
        powf(fabsf(stf_abs_nstk[itf]), nroot), stf_abs_nstk[itf]);
    stf_dB_nstk[itf] = copysign( \
        powf(fabsf(stf_dB_nstk[itf]), nroot), stf_dB_nstk[itf]);
  }

  /*** write out stacked file ***/
  nc_check(nc_create(outfile, NC_CLOBBER, &ncid));

  /* define dimensions */
  nc_check(nc_def_dim(ncid, time_name, ntime, &time_dimid));
  nc_check(nc_def_dim(ncid, freq_name, nfreq, &freq_dimid));

  /* define coordinate variables */
  nc_check(nc_def_var(ncid, time_name, NC_DOUBLE, 1, &time_dimid, &time_varid));
  nc_check(nc_def_var(ncid, freq_name, NC_DOUBLE, 1, &freq_dimid, &freq_varid));

  /* define units attributes for coordinate vars */
  nc_check(nc_put_att_text(ncid, time_varid, UNITS, time_units_len
        , time_units));
  nc_check(nc_put_att_text(ncid, freq_varid, UNITS, strlen(freq_units),
        freq_units));

  /* define data variables */
  dimids[1] = time_dimid;
  dimids[0] = freq_dimid;

  nc_check(nc_def_var(ncid, stf_abs_name, NC_FLOAT, NDIMS, dimids,
      &stf_abs_varid));
  nc_check(nc_def_var(ncid, stf_dB_name, NC_FLOAT, NDIMS, dimids,
      &stf_dB_varid));

  /* define units attributes for data variables */
  nc_check(nc_put_att_text(ncid, stf_abs_varid, UNITS, strlen(stf_abs_units),
      stf_abs_units));
  nc_check(nc_put_att_text(ncid, stf_dB_varid, UNITS, strlen(stf_dB_units),
      stf_dB_units));

  /* end define mode */
  nc_check(nc_enddef(ncid));

  /* write coordinate variable data */
  nc_check(nc_put_var_double(ncid, time_varid, time_coord));
  nc_check(nc_put_var_double(ncid, freq_varid, freq_coord));

  /* write data variable */
  nc_check(nc_put_var_float(ncid, stf_abs_varid, stf_abs_nstk));
  nc_check(nc_put_var_float(ncid, stf_dB_varid, stf_dB_nstk));

  /* close file */
  nc_check(nc_close(ncid));

  /*** finalize ***/
  free(filelist[0]); free(filelist);
  free(time_coord); free(freq_coord); free(time_units);
  free(stf_abs); free(stf_abs_nstk); free(stf_dB_nstk);
  return 0;

  /*** error handling ***/
error:
  free(filelist[0]); free(filelist);
  free(time_coord); free(freq_coord); free(time_units);
  free(stf_abs); free(stf_abs_nstk); free(stf_dB_nstk);
  nc_close(ncid);
  return -1;
}
