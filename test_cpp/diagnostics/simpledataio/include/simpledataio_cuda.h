#include "stdio.h"
#include <stdlib.h>
#include <netcdf.h>

#define SDATIO_INT 0
#define SDATIO_FLOAT 1
#define SDATIO_DOUBLE 2
#define SDATIO_COMPLEX_DOUBLE 3

#define SDATIO_UNLIMITED NC_UNLIMITED


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
};


struct sdatio_file {
	int nc_file_id;
	int is_parallel;
	int n_dimensions;
	struct sdatio_dimension ** dimensions;
	int n_variables;
	struct sdatio_variable ** variables;
	int data_written;
};


int sdatio_debug;

/* Open a new datafile for writing. fname is the name of the file 
 * The stuct sfile is used to store the state information
 * of the file.*/
extern "C" void sdatio_createfile(struct sdatio_file * sfile, char * fname);

/* Create a new dimension in the file sfile. Dimension names must
 * be a single letter. */
extern "C" void sdatio_add_dimension(struct sdatio_file * sfile, 
													 char * dimension_name, 
													 int size,
													 char * description,
													 char * units);

/* Print out a nice list of all the dimensions defined so far*/
extern "C" void sdatio_print_dimensions(struct sdatio_file * sfile);


/* Close the file and free all memory associated with sfile*/
extern "C" void sdatio_close(struct sdatio_file * sfile);

/* Ensure all variables are written to disk in case of crashes*/
extern "C" void sdatio_sync(struct sdatio_file * sfile);

/* Define a variable in the given file. Dimension list 
 * is a character string listing (in order) the dimension names
 * (which are all single characters) e.g. "xyx".*/
extern "C" void sdatio_create_variable(struct sdatio_file * sfile,
														int variable_type,
														char * variable_name,
														char * dimension_list,
														char * description,
														char * units);

/* Write to the given variable. address should be the address of the start of the array */
extern "C" void sdatio_write_variable(struct sdatio_file * sfile, char * variable_name, void * address);

/* Write to the given variable. address should be the address of the start of the array. Indexes should be an array the same size as the number of dimensions of the variable. Using the second form is quicker as the first form requires a search for the variable at every write*/
extern "C" void sdatio_write_variable_at_index(struct sdatio_file * sfile, char * variable_name, int * indexes, void * address);
extern "C" void sdatio_write_variable_at_index_fast(struct sdatio_file * sfile, struct sdatio_variable * svar, int * indexes, void * address);

/* Print out a nice list of all the variables defined so far*/
extern "C" void sdatio_print_variables(struct sdatio_file * sfile);

/* Return a pointer the struct containing all the metadata of the given variable */
extern "C" struct sdatio_variable * sdatio_find_variable(struct sdatio_file * sfile, char * variable_name);

/* Return a pointer the struct containing all the metadata of the given dimension */
extern "C" struct sdatio_dimension * sdatio_find_dimension(struct sdatio_file * sfile, char * dimension_name);

/* Increment the start of the specified infinite dimension */
extern "C" void sdatio_increment_start(struct sdatio_file * sfile, char * dimension_name);

/* Returns 1 if the given variable exists, 0 otherwise */
extern "C" int sdatio_variable_exists(struct sdatio_file * sfile, char * variable_name);
