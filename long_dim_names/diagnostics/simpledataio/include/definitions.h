#include "string.h"
#include <time.h>
#include "stdio.h"
#include <stdlib.h>
#include <netcdf.h>

#ifdef PARALLEL 
#include "netcdf_par.h"
#endif
#if HAVE_MPI
#include "mpi.h"
#else
typedef int MPI_Comm;
typedef int MPI_Fint;
#endif

#define SDATIO_INT 0
#define SDATIO_FLOAT 1
#define SDATIO_DOUBLE 2
#define SDATIO_COMPLEX_DOUBLE 3
#define SDATIO_CHAR 4

#define SDATIO_UNLIMITED NC_UNLIMITED

#define QUOTE(str) #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)
/*#define VERSION 5.8*/
#define SDATIO_VERSION_STRING EXPAND_AND_QUOTE( PACKAGE_VERSION )

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION 1.0.0
#endif



struct sdatio_dimension {
	char * name;
	int size;
	int nc_id;
	int start;
};

struct sdatio_variable {
	char * name;
	int nc_id;
	int type;
	char * dimension_list;
	int * dimension_ids;
	int type_size;
	int * manual_counts;
	int * manual_starts;
	/* Only used for Fortran:*/
	int * manual_offsets;
  /* Necessary now we can have long dim names */
  int ndims;
};


struct sdatio_file {
	int nc_file_id;
	int is_parallel;
  int is_open;
	int n_dimensions;
	struct sdatio_dimension ** dimensions;
	int n_variables;
	struct sdatio_variable ** variables;
	int data_written;
  MPI_Comm * communicator;
  int mode;
  char * name;
  int has_long_dim_names;
};


int sdatio_debug;

