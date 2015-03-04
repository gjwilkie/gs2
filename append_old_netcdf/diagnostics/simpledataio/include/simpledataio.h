#include "include/definitions.h"

/* Initialise a new sdatio_file object.*/
void sdatio_init(struct sdatio_file * sfile, char * fname);

/* Free memory associated with an sdatio_file object.*/
void sdatio_free(struct sdatio_file * sfile);

/* Switch the sdatio_file object to parallel mode and 
 * set the MPI communicator */
void sdatio_set_parallel(struct sdatio_file * sfile, MPI_Comm * comm);

/* Open a new datafile for writing. fname is the name of the file 
 * The stuct sfile is used to store the state information
 * of the file.*/
void sdatio_create_file(struct sdatio_file * sfile);

/* Open an existing datafile for reading and/or appending data. fname is the name of the file 
 * The struct sfile is used to store the state information
 * of the file.*/
void sdatio_open_file(struct sdatio_file * sfile);

/* Write metadata to the file (as a netcdf global attribute)*/
void sdatio_add_metadata(struct sdatio_file * sfile, const int metadata_type, const char * key, const void * value);

/* Create a new dimension in the file sfile. Dimension names must
 * be a single letter. */
void sdatio_add_dimension(struct sdatio_file * sfile, 
                           char * dimension_name, 
                           int size,
                           char * description,
                           char * units);

/* Print out a nice list of all the dimensions defined so far*/
void sdatio_print_dimensions(struct sdatio_file * sfile);


/* Close the file and free all memory associated with sfile*/
void sdatio_close(struct sdatio_file * sfile);

/* Write default metadata such as simpledataio and library versions */
void sdatio_add_standard_metadata(struct sdatio_file * sfile);

/* Ensure all variables are written to disk in case of crashes*/
void sdatio_sync(struct sdatio_file * sfile);

/* Define a variable in the given file. Dimension list 
 * is a character string listing (in order) the dimension names
 * (which are all single characters) e.g. "xyx".*/
void sdatio_create_variable(struct sdatio_file * sfile,
                            int variable_type,
                            char * variable_name,
                            char * dimension_list,
                            char * description,
                            char * units);

/* Write to the given variable. address should be the address of the start of the array */
void sdatio_write_variable(struct sdatio_file * sfile, char * variable_name, void * address);

/* Write to the given variable. address should be the address of the start of the array. Indexes should be an array the same size as the number of dimensions of the variable. Using the second form is quicker as the first form requires a search for the variable at every write*/
void sdatio_write_variable_at_index(struct sdatio_file * sfile, char * variable_name, int * indexes, void * address);
void sdatio_write_variable_at_index_fast(struct sdatio_file * sfile, struct sdatio_variable * svar, int * indexes, void * address);

/* Return a pointer the struct containing all the metadata of the given variable */
struct sdatio_variable * sdatio_find_variable(struct sdatio_file * sfile, char * variable_name);
/* Return a pointer the struct containing all the metadata of the given dimension */
struct sdatio_dimension * sdatio_find_dimension(struct sdatio_file * sfile, char * dimension_name);

/* Print out a nice list of all the variables defined so far*/
void sdatio_print_variables(struct sdatio_file * sfile);

/* Increment the start of the specified infinite dimension */
void sdatio_increment_start(struct sdatio_file * sfile, char * dimension_name);

/* Returns 1 if the given variable exists, 0 otherwise */
int sdatio_variable_exists(struct sdatio_file * sfile, char * variable_name);

