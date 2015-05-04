#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <netcdf.h>

#include "dbg.h"

/***** main program  *****/

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

  char filename[]="out.nc";

  /* netcdf */
  int ncid=-1;
  /* axis */
  const char time_name[]="time";
  const char freq_name[]="log10_frequency";
  size_t time_units_len=0;
  char *time_units=NULL;
  size_t freq_units_len=0;
  //char *freq_units=NULL;
  int time_dimid, freq_dimid;
  size_t time_dimlen, freq_dimlen;
  int time_varid, freq_varid; /* coordinate variables */
  double *time_var=NULL;
  //double *freq_var=NULL;

  /* open input file. */
  nc_check(nc_open(filename, NC_NOWRITE, &ncid));

  /* get dimension id */
  nc_check(nc_inq_dimid(ncid, time_name, &time_dimid));
  nc_check(nc_inq_dimid(ncid, freq_name, &freq_dimid));

  /* get dimension lengths */
  nc_check(nc_inq_dimlen(ncid, time_dimid, &time_dimlen));
  nc_check(nc_inq_dimlen(ncid, freq_dimid, &freq_dimlen));
    
  printf("time_dimlen=%zu freq_dimlen=%zu\n", time_dimlen, freq_dimlen);

  /* get coordinates (data and units): only use the first file */

  /* time */
  nc_check(nc_inq_varid(ncid, time_name, &time_varid));

  check_mem(time_var = malloc(sizeof(double) * time_dimlen));
  nc_check(nc_get_var_double(ncid, time_varid, time_var));

  nc_check(nc_inq_attlen(ncid, time_varid, UNITS, &time_units_len));
  check_mem(time_units = malloc(sizeof(char) * time_units_len));
  nc_check(nc_get_att_text(ncid, time_varid, UNITS, time_units));

  printf("time_units_len=%zu time_units=%s\n", time_units_len, time_units);

  /* clean up */
  free(time_var); free(time_units);
  //free(freq_var); free(freq_units);
  return 0;

error:
  free(time_var); free(time_units);
  //free(freq_var); free(freq_units);
  return -1;
}
