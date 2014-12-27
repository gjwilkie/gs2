/******************************************************
 * Simple Data I/O
 *
 * A wrapper that provides a simplified interface to
 * software packages such as netcdf. It makes certain
 * assumptions about the write calls, specifically that
 * you are always writing the whole of non-infinite 
 * dimensions and only one element of infinite dimensions.
 * This is often the case!
 *
 * This is free software released under GPLv2
 *
 * Authors:
 *      Edmund Highcock (edmund.highcock@users.sourceforge.net)
 *
 ********************************************************/

#include "include/simpledataio.h"


#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define DEBUG_MESS if (sdatio_debug) printf


void sdatio_init(struct sdatio_file * sfile, char * fname){
  sfile->mode = NC_NETCDF4|NC_CLOBBER;
  sfile->is_open = 0;
  sfile->is_parallel = 0;
  sfile->has_long_dim_names = 0;
  sfile->communicator = (MPI_Comm*)malloc(sizeof(MPI_Comm));
  sfile->name = (char*)malloc(sizeof(char)*(strlen(fname)+1));
  strcpy(sfile->name, fname);
}

void sdatio_free(struct sdatio_file * sfile){
  if (sfile->is_open)
    printf("WARNING: Freeing an sdatio_file object associated with an open file; this may lead to a memory leak. Suggest closing the file first using sdatio_close\n");
  free(sfile->communicator);
  free(sfile->name);
}

void sdatio_set_parallel(struct sdatio_file * sfile, MPI_Comm * comm){
  if (sfile->is_open)
    printf("WARNING: called sdatio_set_parallel on an open file... expect the unexpected\n");
  sfile->is_parallel = 1;
  sfile->mode = NC_NETCDF4|NC_CLOBBER|NC_MPIPOSIX;
  *(sfile->communicator) = *comm;
}

void sdatio_set_parallel_fortran(struct sdatio_file * sfile, MPI_Fint fcomm){
#ifdef PARALLEL
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  sdatio_set_parallel(sfile, &comm);
#else 
  printf("sdatio was built without --enable-parallel, sdatio_set_parallel will not work\n");
  abort();
#endif
}

/* Private*/
void sdatio_end_definitions(struct sdatio_file * sfile){
  int retval;
  if ((retval = nc_enddef(sfile->nc_file_id))) ERR(retval);
}

/* Private */
void sdatio_recommence_definitions(struct sdatio_file * sfile){
  int retval;
  if (sfile->data_written) 
    printf("Warning: adding more variables or dimensions after writing data can slow performance\n");
  if ((retval = nc_redef(sfile->nc_file_id))) ERR(retval);
}


void sdatio_create_file(struct sdatio_file * sfile )  {
  /*printf("called\n");*/
  int retval;
  /*if (0){}*/
  /*else {*/

  /*char * args;*/
  retval = 0;

  if (sfile->is_open){
    printf("ERROR: The supplied sdatio_file struct corresponds to an already open file.\n");
    abort();
  }

  /*MPI_Init(&retval, &args);*/
  if (sfile->is_parallel) {
#ifdef PARALLEL 
    if ((retval = nc_create_par(sfile->name, sfile->mode, *(sfile->communicator), MPI_INFO_NULL,  &(sfile->nc_file_id)))) ERR(retval);
#else
    printf("sdatio was built without --enable-parallel, sdatio_create_file will not work for parallel files\n");
    abort();
#endif
  }
  else {
    if ((retval = nc_create(sfile->name, sfile->mode, &(sfile->nc_file_id)))) ERR(retval);
  }
  sfile->is_open = 1;
  sfile->n_dimensions = 0;
  sfile->n_variables = 0;
  sfile->data_written = 0;
    /*}*/
  sdatio_end_definitions(sfile);
}



void sdatio_add_standard_metadata(struct sdatio_file * sfile){
  time_t current_time;
  char  strtime[40];
  size_t count = 40;
  current_time = time(NULL);
  strftime(strtime, count, "%c %z", localtime(&current_time)); 
  sdatio_add_metadata(sfile, SDATIO_CHAR, "creation_time", &strtime);
  sdatio_add_metadata(sfile, SDATIO_INT, "creation_time_in_seconds_after_epoch", (int *)&current_time);
  sdatio_add_metadata(sfile, SDATIO_CHAR, "simpledatio_info", 
      "This datafile was constructed using simpledatio, a simplified netCDF interface.");
  sdatio_add_metadata(sfile, SDATIO_CHAR, "simpledatio_url", 
      "http://github.com/edmundhighcock/simpledataio");
  sdatio_add_metadata(sfile, SDATIO_CHAR, "simpledatio_version", SDATIO_VERSION_STRING);
  sdatio_add_metadata(sfile, SDATIO_CHAR, "netcdf_version", nc_inq_libvers());
  sdatio_add_metadata(sfile, SDATIO_CHAR, "netcdf_url", "http://www.unidata.ucar.edu/software/netcdf");
  
}

/*void sdatio_create_file_parallel_fortran(struct sdatio_file * sfile, char * fname, MPI_Fint  fcomm)  {*/
/*#ifdef PARALLEL*/
/*MPI_Comm comm = MPI_Comm_f2c(fcomm);*/
/*sdatio_create_file_parallel(sfile, fname, &comm);*/
/*#else */

/*#endif*/
/*}*/

/*void sdatio_create_file_with_mode(struct sdatio_file * sfile, int mode, char * fname)  {*/
/**//*printf("called\n");*/
/*int retval;*/
/**//*if (0){}*/
/**//*else {*/
/*if ((retval = nc_create(fname, NC_CLOBBER, &(sfile->nc_file_id)))) ERR(retval);*/
/**//*}*/
/*sdatio_end_definitions(sfile);*/
/*}*/
/*void sdatio_create_file(struct sdatio_file * sfile, char * fname)  {*/
/*sdatio_create_file_with_mode(sfile, NC_CLOBBER|NC_NETCDF4, fname);*/
/**//*printf("called\n");*/
/**//*int retval;*/
/**//**//*if (0){}*/
/**//**//*else {*/
/**//*if ((retval = nc_create(fname, NC_CLOBBER, &(sfile->nc_file_id)))) ERR(retval);*/
/**//**//*}*/
/**//*sfile->n_dimensions = 0;*/
/**//*sfile->n_variables = 0;*/
/**//*sfile->is_parallel = 0;*/
/**//*sfile->data_written = 0;*/
/**//*sdatio_end_definitions(sfile);*/
/*}*/

int sdatio_netcdf_variable_type(int type){
  switch (type){
  case SDATIO_INT:
    return NC_INT;
  case SDATIO_FLOAT:
    return NC_FLOAT;
  case SDATIO_DOUBLE:
    return NC_DOUBLE;
  case SDATIO_CHAR:
    return NC_CHAR;
    /*case */
    /*case SDATIO_COMPLEX_DOUBLE:*/
    /*printf("Can't do complex yet\n");*/
    /*abort();*/
  default:
    printf("Unknown data type for simple data io\n");
    abort();
  }
}

void sdatio_add_metadata(struct sdatio_file * sfile, const int metadata_type, const char * key, const void * value){
  int retval;
  sdatio_recommence_definitions(sfile);
  switch (metadata_type){
    case SDATIO_CHAR:
      if ((retval = nc_put_att_text(sfile->nc_file_id, NC_GLOBAL, key, strlen(value), value))) ERR(retval);
      break;
    default:
      if ((retval = nc_put_att(sfile->nc_file_id, NC_GLOBAL, key, 
              sdatio_netcdf_variable_type(metadata_type), 1, value))) ERR(retval);
  }
  sdatio_end_definitions(sfile);
}

/***********************************************************
 *
 * Handling Dimensions
 *
 **********************************************************/


/* Private */
void sdatio_append_dimension(struct sdatio_file * sfile, struct sdatio_dimension * sdim){
  int ndims;
  struct sdatio_dimension ** new_dimensions; 
  int i;
  ndims = sfile->n_dimensions + 1;
    

  new_dimensions = (struct sdatio_dimension **) malloc(sizeof(struct sdatio_dimension *)*ndims);

  DEBUG_MESS("Setting dimensions; %d, %d\n", ndims, sfile->n_dimensions);

  for (i=0; i < ndims-1; i++){
    DEBUG_MESS("i %d\n", i);
    new_dimensions[i] = sfile->dimensions[i];
  }

  DEBUG_MESS("Setting new\n");

  new_dimensions[ndims-1] = sdim;

  DEBUG_MESS("About to deallocate old dimensions\n");

  if (sfile->n_dimensions > 0) free(sfile->dimensions);
  sfile->n_dimensions = ndims;

  sfile->dimensions = new_dimensions;
  
}

void sdatio_add_dimension(struct sdatio_file * sfile, 
                          char * dimension_name, 
                          int size,
                          char * description,
                          char * units){

  struct sdatio_dimension  * sdim;
  int retval;
  /*printf("Dimension inputs: %s, %d, %s, %s\n", dimension_name, size, description, units);*/
  sdim = (struct sdatio_dimension *) malloc(sizeof(struct sdatio_dimension));
  sdatio_recommence_definitions(sfile);
  /*if (sfile->is_parallel){}*/
  /*else {*/
  if ((retval = nc_def_dim(sfile->nc_file_id, dimension_name, size, &(sdim->nc_id)))) ERR(retval);
  /*}*/
  sdatio_end_definitions(sfile);
  sdim->size = size;
  if (strlen(dimension_name)>1){
    sfile->has_long_dim_names = 1;
    /*printf("Dimension names can only be one character long!\n");*/
    /*abort();*/
  }
  sdim->name = (char *)malloc(sizeof(char)*(strlen(dimension_name)+1));
  strcpy(sdim->name, dimension_name);
  sdim->start = 0;
  sdatio_append_dimension(sfile, sdim);

}

void sdatio_print_dimensions(struct sdatio_file * sfile){
  int i;
  struct sdatio_dimension * sdim;
  for (i=0;i<sfile->n_dimensions;i++){
    sdim = sfile->dimensions[i];
    printf("Dimension %s, size %d, has id %d\n", sdim->name, sdim->size, sdim->nc_id);
  }
}

void sdatio_increment_start(struct sdatio_file * sfile, char * dimension_name){

  int found, j;
  struct sdatio_dimension * sdim;

  found = 0;
  for (j=0;j<sfile->n_dimensions;j++){
    sdim = sfile->dimensions[j];
    if (!strcmp(sdim->name, dimension_name)){
      if (sdim->size != SDATIO_UNLIMITED) {
        printf("Dimension %s does not have unlimited size.\n", dimension_name);
        abort();
      }   
      found = 1;
      (sdim->start)++;
    }
  }
  if (!found) {
    printf("Couldn't find dimension %s in sdatio_increment_count\n", dimension_name);
    abort();
  }
}

/* Private*/
void sdatio_free_dimension(struct sdatio_dimension * sdim){
  free(sdim->name);
  free(sdim);
}



/***************************************************
 *
 * Handling Variables
 *
 * ***********************************************/



void sdatio_get_dimension_ids(struct sdatio_file * sfile, char * dimension_list, struct sdatio_variable * svar){
  int ndims;
  int i,j;
  /*char dim_name[2];*/
  char * dim_name;
  int * dimension_ids;
  int dim_name_length;
  int counter;
  int sep_size;
  ndims  = svar->ndims;
  counter = 0;
  DEBUG_MESS("ndims %d\n", ndims);
  dimension_ids = (int *) malloc(sizeof(int)*ndims);
  for (i=0;i<ndims;i++){
    /* In the next section we set counter to 
     * the beginning of the next dimension name and
     * dim_name_length to the size of the name*/
    if (sfile->has_long_dim_names){
      sep_size = 1;
      dim_name_length = 0;
      /*The first condition checks that we haven't reached the end of the string*/
      while (dimension_list[counter] && !(dimension_list[counter]==',')){
        dim_name_length++;
        counter++;
      }
      counter++;
    }
    else {
      sep_size = 0;
      dim_name_length = 1;
      counter++;
    }
    if (dim_name_length < 1){
      printf("ERROR: zero length dimension name in dimension_list %s\n", dimension_list);
      abort();
    }

    /*dim_name[0] = dimension_list[i];*/
    /*dim_name[1] = dimension_list[ndims];*/

    /* Copy the name of the current dimension to the temporary
     * variable dim_name*/
    dim_name = (char *)malloc(sizeof(char)*(dim_name_length+1));
    strncpy(dim_name, dimension_list+(counter-dim_name_length-sep_size), dim_name_length);
    dim_name[dim_name_length] = '\0';
    DEBUG_MESS("svar %s: dim_name_length %d, dim_name: %s, counter: %d\n", svar->name, dim_name_length, dim_name, counter); 
    DEBUG_MESS("i %d\n", i);
    DEBUG_MESS("Getting id for dim %s\n", dim_name);

    /*Find the id for the dimension named dim_name*/
    dimension_ids[i] = -1;
    for (j=0;j<sfile->n_dimensions;j++){
      DEBUG_MESS("j %d\n", j);
      if (!strcmp(dim_name, sfile->dimensions[j]->name)) 
        dimension_ids[i] = sfile->dimensions[j]->nc_id;
    }
    if (dimension_ids[i]==-1){
      printf("Dimension %s is undefined!\n", dim_name);
      abort();
    }
    DEBUG_MESS("Finished loop\n");
    DEBUG_MESS("dim %s has id %d \n", dim_name, dimension_ids[i]);
    free(dim_name);
  }
  svar->dimension_ids = dimension_ids;


}
/*Private*/
void sdatio_append_variable(struct sdatio_file * sfile, struct sdatio_variable * svar){
  int nvars;
  struct sdatio_variable ** new_variables; 
  int i;
  nvars = sfile->n_variables + 1;
    

  new_variables = (struct sdatio_variable **) malloc(sizeof(struct sdatio_variable *)*nvars);

  DEBUG_MESS("Setting variable %s; %d, %d\n", svar->name, nvars, sfile->n_variables);

  for (i=0; i < nvars-1; i++){
    DEBUG_MESS("i %d\n", i);
    new_variables[i] = sfile->variables[i];
  }

  DEBUG_MESS("Setting new\n");

  new_variables[nvars-1] = svar;

  DEBUG_MESS("About to deallocate old variables\n");

  if (sfile->n_variables > 0) free(sfile->variables);

  sfile->n_variables = nvars;

  sfile->variables = new_variables;

  DEBUG_MESS("Deallocated old vars\n");
  
}

int sdatio_number_of_dimensions(struct sdatio_file * sfile, char * variable_name){
  struct sdatio_variable * svar = sdatio_find_variable(sfile, variable_name);
  return svar->ndims;
}

void sdatio_number_of_unlimited_dimensions(struct sdatio_file * sfile, char * variable_name, int * n){
  struct sdatio_variable * svar = sdatio_find_variable(sfile, variable_name);
  struct sdatio_dimension * sdim;
  int i,j;
  int found;
  *n = 0;
  for (i=0;i<svar->ndims;i++){
    found = 0;
    for (j=0;j<sfile->n_dimensions;j++){
      sdim = sfile->dimensions[j];
      if (sdim->nc_id == svar->dimension_ids[i]){
        found = 1;
        if (sdim->size == SDATIO_UNLIMITED) (*n)++; 
      }
    }
    if (!found) {
      printf("Couldn't find dimension in sdatio_get_counts_and_starts\n");
      abort();
    }
  }
  /*printf("n unlimited was %d\n", *n);*/
}

/*Private*/
int sdatio_ndims_from_string(struct sdatio_file * sfile, char * dimension_list){
  int i;
  char * s = dimension_list;
  /* This little bit of code counts the number of commas in 
   * dimension_list string. I could have made a more readable bit
   * of code, but this little fragment I copied from the web was 
   * too cute to miss*/
  for (i=0; s[i]; s[i]==',' ? i++ : *s++);

  /* An empty string means no dimensions, i.e. a scalar */
  if (strlen(dimension_list)==0) return 0;
  /* Obv the number of dimensions is one more than
   * the number of commas */
  else if (sfile->has_long_dim_names || i>0) return i + 1;
  else return strlen(dimension_list);
}

void sdatio_create_variable(struct sdatio_file * sfile,
                            int variable_type,
                            char * variable_name,
                            char * dimension_list,
                            char * description,
                            char * units){
  int ndims;
  int nunlim;
  struct sdatio_variable  * svar;
  int retval;
  /*int * dimension_ids;*/

  /*dimension_ids = (int **)malloc(sizeof(int*));*/


  /*printf("dimension_list is %s\n", dimension_list);*/
  svar = (struct sdatio_variable *) malloc(sizeof(struct sdatio_variable));

  /* Set variable name*/
  svar->name = (char *)malloc(sizeof(char)*(strlen(variable_name)+1));
  strcpy(svar->name, variable_name);

  /*ndims = strlen(dimension_list);*/
  ndims = sdatio_ndims_from_string(sfile, dimension_list);
  svar->ndims = ndims;
  DEBUG_MESS("ndims = %d for variable %s\n", ndims, variable_name);

  sdatio_get_dimension_ids(sfile, dimension_list, svar);
  /*svar->dimension_ids = dimension_ids;*/


  sdatio_recommence_definitions(sfile);
  /*if (sfile->is_parallel){}*/
  /*else {*/
    if ((retval = nc_def_var(sfile->nc_file_id, variable_name, sdatio_netcdf_variable_type(variable_type), ndims, svar->dimension_ids, &(svar->nc_id)))) ERR(retval);
    if ((retval = nc_put_att_text(sfile->nc_file_id, svar->nc_id, "description", strlen(description), description))) ERR(retval);
    if ((retval = nc_put_att_text(sfile->nc_file_id, svar->nc_id, "units", strlen(units), units))) ERR(retval);
    /*}*/
  switch (variable_type){
    case SDATIO_INT:
      svar->type_size = sizeof(int);
      break;
    case SDATIO_FLOAT:
      svar->type_size = sizeof(float);
      break;
    case SDATIO_DOUBLE:
      svar->type_size = sizeof(double);
      break;
    case SDATIO_CHAR:
      svar->type_size = sizeof(char);
      break;
    default:
      printf("Unknown type in sdatio_create_variable\n");
      abort();
  }

  sdatio_end_definitions(sfile);
  
  svar->type = variable_type;
  svar->dimension_list = (char *)malloc(sizeof(char)*(strlen(dimension_list)+1));
  strcpy(svar->dimension_list, dimension_list);

  svar->manual_starts=(int*)malloc(sizeof(int)*ndims);
  svar->manual_counts=(int*)malloc(sizeof(int)*ndims);
  svar->manual_offsets=(int*)malloc(sizeof(int)*ndims);
  int i;

  for (i=0;i<ndims;i++){
    svar->manual_starts[i]=-1;
    svar->manual_counts[i]=-1;
    svar->manual_offsets[i]=-1;
  }

  DEBUG_MESS("Starting sdatio_append_variable\n");

  sdatio_append_variable(sfile, svar);

  DEBUG_MESS("Ending sdatio_append_variable\n");

#ifdef PARALLEL
  if (sfile->is_parallel){
    sdatio_number_of_unlimited_dimensions(sfile, variable_name, &nunlim);
    if (nunlim > 0)
      if ((retval = nc_var_par_access(sfile->nc_file_id, svar->nc_id, NC_COLLECTIVE))) ERR(retval);
  }
#endif
  
}

void sdatio_print_variables(struct sdatio_file * sfile){
  int i;
  struct sdatio_variable * svar;
  for (i=0;i<sfile->n_variables;i++){
    svar = sfile->variables[i];
    printf("Variable %s, dimensions %s, has id %d\n", svar->name, svar->dimension_list, svar->nc_id);
  }
}

/* Private */
void sdatio_get_counts_and_starts(struct sdatio_file * sfile, struct sdatio_variable * svar, size_t * counts, size_t * starts){
  struct sdatio_dimension * sdim;
  int i,j;
  int found;
  for (i=0;i<svar->ndims;i++){
    found = 0;
    for (j=0;j<sfile->n_dimensions;j++){
      sdim = sfile->dimensions[j];
      if (sdim->nc_id == svar->dimension_ids[i]){
        if (svar->manual_starts[i] == -1) starts[i] = sdim->start;
        else starts[i] = svar->manual_starts[i];
        found = 1;
        if (sdim->size == SDATIO_UNLIMITED) counts[i] = 1; 
        else if (svar->manual_counts[i] == -1 ) counts[i] = sdim->size;
        else counts[i] = svar->manual_counts[i];
      }
    }
    if (!found) {
      printf("Couldn't find dimension in sdatio_get_counts_and_starts\n");
      abort();
    }
  }
}
void sdatio_get_offsets(struct sdatio_file * sfile, struct sdatio_variable * svar, size_t * starts, size_t * offsets){
  struct sdatio_dimension * sdim;
  int i,j;
  int found;
  for (i=0;i<svar->ndims;i++){
    found = 0;
    for (j=0;j<sfile->n_dimensions;j++){
      sdim = sfile->dimensions[j];
      if (sdim->nc_id == svar->dimension_ids[i]){
        if (svar->manual_offsets[i] == -1) offsets[i] = starts[i];
        else offsets[i] = svar->manual_offsets[i];
        found = 1;
      }
    }
    if (!found) {
      printf("Couldn't find dimension in sdatio_get_offsets\n");
      abort();
    }
  }
}

/*** ONLY WORKS FOR THE FORTRAN INTERFACE ***/
/* Do not set offsets when using the C interface*/
void sdatio_set_offset(struct sdatio_file * sfile, char * variable_name, char * dimension_name, int * offset){
  struct sdatio_variable * svar = sdatio_find_variable(sfile, variable_name);
  struct sdatio_dimension * sdim;
  /*struct sdatio_dimension * sdim_found;*/
  int i,j;
  int found;
  int ndim;

  

  found = 0;
  for (i=0;i<svar->ndims;i++){
    for (j=0;j<sfile->n_dimensions;j++){
      sdim = sfile->dimensions[j];
      /*printf("sdim %s, comp %d\n", sdim->name, !(strcmp(sdim->name, dimension_name)));*/
      if ((sdim->nc_id == svar->dimension_ids[i]) && !strcmp(sdim->name, dimension_name)){
        found = 1;
        /*sdim_found = sdim;*/
        ndim = i;
      }
    }
  }
    if (!found) {
      printf("Couldn't find dimension %s for variable %s in sdatio_set_offset\n", dimension_name, svar->name);
      abort();
    }
    /*printf("Start is %d\n", svar->manual_offsets[ndim]);*/
  svar->manual_offsets[ndim] = *offset;
  /*printf("Start is %d\n", svar->manual_starts[ndim]);*/


}

void sdatio_set_start(struct sdatio_file * sfile, char * variable_name, char * dimension_name, int * start){
  struct sdatio_variable * svar = sdatio_find_variable(sfile, variable_name);
  struct sdatio_dimension * sdim;
  /*struct sdatio_dimension * sdim_found;*/
  int i,j;
  int found;
  int ndim;

  

  found = 0;
  for (i=0;i<svar->ndims;i++){
    for (j=0;j<sfile->n_dimensions;j++){
      sdim = sfile->dimensions[j];
      /*printf("sdim %s, comp %d\n", sdim->name, !(strcmp(sdim->name, dimension_name)));*/
      if ((sdim->nc_id == svar->dimension_ids[i]) && !strcmp(sdim->name, dimension_name)){
        found = 1;
        /*sdim_found = sdim;*/
        ndim = i;
      }
    }
  }
    if (!found) {
      printf("Couldn't find dimension %s for variable %s in sdatio_set_start\n", dimension_name, svar->name);
      abort();
    }
    /*printf("Start is %d\n", svar->manual_starts[ndim]);*/
  svar->manual_starts[ndim] = *start;
  /*printf("Start is %d\n", svar->manual_starts[ndim]);*/


}
void sdatio_set_count(struct sdatio_file * sfile, char * variable_name, char * dimension_name, int * count){
  struct sdatio_variable * svar = sdatio_find_variable(sfile, variable_name);
  struct sdatio_dimension * sdim;
  /*struct sdatio_dimension * sdim_found;*/
  int i,j;
  int found;
  int ndim;

  

    found = 0;
  for (i=0;i<svar->ndims;i++){
    for (j=0;j<sfile->n_dimensions;j++){
      sdim = sfile->dimensions[j];
      if (sdim->nc_id == svar->dimension_ids[i] && !strcmp(sdim->name, dimension_name)){
        found = 1;
        /*sdim_found = sdim;*/
        ndim = i;
      }
    }
  }
    if (!found) {
      printf("Couldn't find dimension %s for variable %s in sdatio_set_count\n", dimension_name, variable_name);
      abort();
    }
    /*printf("count is %d\n", svar->manual_counts[ndim]);*/
  svar->manual_counts[ndim] = *count;
  /*printf("count is %d\n", svar->manual_counts[ndim]);*/


}



/* Private: used for the Fortran interface*/
void sdatio_netcdf_inputs(struct sdatio_file * sfile, char * variable_name, int * fileid, int * varid, size_t * starts, size_t * counts, size_t * offsets){
  struct sdatio_variable * svar = sdatio_find_variable(sfile, variable_name);
  sdatio_get_counts_and_starts(sfile, svar, counts, starts);
  sdatio_get_offsets(sfile, svar, starts, offsets);
  *fileid = sfile->nc_file_id;
  *varid = svar->nc_id;
  /*printf("varname %s, fileid %d, varid %d, starts[0] %d \n", variable_name, *fileid, *varid, starts[0]);  */
}


void sdatio_write_variable_private(struct sdatio_file * sfile, struct sdatio_variable * svar, size_t * counts, size_t * starts, void * address){
  int retval;
  /*if (sfile->is_parallel){}*/
  /*else {*/
  switch (svar->type){
  case (SDATIO_INT):
    DEBUG_MESS("Writing an integer\n");
    if ((retval = nc_put_vara_int(sfile->nc_file_id, svar->nc_id, starts, counts, address))) ERR(retval);
    break;
  case (SDATIO_FLOAT):
    DEBUG_MESS("Writing a float\n");
    /*if ((retval = nc_put_var_double(sfile->nc_file_id, svar->nc_id, address))) ERR(retval);*/
    if ((retval = nc_put_vara_float(sfile->nc_file_id, svar->nc_id, starts, counts, address))) ERR(retval);
    break;
  case (SDATIO_DOUBLE):
    DEBUG_MESS("Writing a double\n");
    /*if ((retval = nc_put_var_double(sfile->nc_file_id, svar->nc_id, address))) ERR(retval);*/
    if ((retval = nc_put_vara_double(sfile->nc_file_id, svar->nc_id, starts, counts, address))) ERR(retval);
    break;
  }
  
  /*}*/
  sfile->data_written = 1;
}


struct sdatio_dimension * sdatio_find_dimension(struct sdatio_file * sfile, char * dimension_name){
  int i, dimension_number;
  dimension_number = -1;

  DEBUG_MESS("Finding dimension...\n");

  for (i=0;i<sfile->n_dimensions;i++)
    if (!strcmp(sfile->dimensions[i]->name, dimension_name))
      dimension_number = i;

  if (dimension_number==-1){
    printf("Couldn't find dimension %s\n", dimension_name);
    abort();
  }
  return sfile->dimensions[dimension_number];
}

struct sdatio_variable * sdatio_find_variable(struct sdatio_file * sfile, char * variable_name){
  int i, variable_number;
  variable_number = -1;

  DEBUG_MESS("Finding variable...\n");

  for (i=0;i<sfile->n_variables;i++)
    if (!strcmp(sfile->variables[i]->name, variable_name))
      variable_number = i;

  if (variable_number==-1){
    printf("Couldn't find variable %s\n", variable_name);
    abort();
  }
  return sfile->variables[variable_number];
}


void sdatio_collective(struct sdatio_file * sfile, char * variable_name){
#ifdef PARALLEL
  int retval;
  struct sdatio_variable * svar = sdatio_find_variable(sfile, variable_name);
  if ((retval = nc_var_par_access(sfile->nc_file_id, svar->nc_id, NC_COLLECTIVE))) ERR(retval);
#endif
}
void sdatio_independent(struct sdatio_file * sfile, char * variable_name){
#ifdef PARALLEL
  int retval;
  struct sdatio_variable * svar = sdatio_find_variable(sfile, variable_name);
  if ((retval = nc_var_par_access(sfile->nc_file_id, svar->nc_id, NC_INDEPENDENT))) ERR(retval);
#endif
}

/* Returns 1 if the given variable exists, 0 otherwise */
int sdatio_variable_exists(struct sdatio_file * sfile, char * variable_name){
  int i;

  DEBUG_MESS("Finding variable in sdatio_variable_exists...\n");

  for (i=0;i<sfile->n_variables;i++)
    if (!strcmp(sfile->variables[i]->name, variable_name))
      return 1;

  return 0;

}


void sdatio_write_variable_fortran_convert(struct sdatio_file * sfile, char * variable_name, void ** address){
  /*printf("address2 is %d\n", address);*/
  printf("address3 is %d\n", *address);
  printf("address4 is %d\n", &address);
  sdatio_write_variable(sfile, variable_name, *address);
}

void sdatio_write_variable(struct sdatio_file * sfile, char * variable_name, void * address){
  int  ndims;
  struct sdatio_variable * svar;
  /*double * double_array;*/
  size_t * counts, * starts;

  /*printf("address is %d\n", address);*/
  /*printf("value is %f\n", *((float*)address));*/

  
  svar = sdatio_find_variable(sfile, variable_name);

  ndims = svar->ndims;
  counts = (size_t*)malloc(sizeof(size_t)*ndims); 
  starts = (size_t*)malloc(sizeof(size_t)*ndims); 

  sdatio_get_counts_and_starts(sfile, svar, counts, starts);

  sdatio_write_variable_private(sfile, svar, counts, starts, address);

  sdatio_sync(sfile);

  free(counts);
  free(starts);

}


void sdatio_write_variable_at_index(struct sdatio_file * sfile, char * variable_name, int * indexes, void * address){
  struct sdatio_variable * svar;
  svar = sdatio_find_variable(sfile, variable_name);
  sdatio_write_variable_at_index_fast(sfile, svar, indexes, address);
}

void sdatio_write_variable_at_index_fast(struct sdatio_file * sfile, struct sdatio_variable * svar, int * indexes, void * address){
  int i, ndims;
  /*double * double_array;*/
  size_t * counts, * starts;



  ndims = svar->ndims;
  counts = (size_t*)malloc(sizeof(size_t)*ndims); 
  starts = (size_t*)malloc(sizeof(size_t)*ndims); 

  for (i=0;i<ndims;i++){
    counts[i] = 1;
    starts[i] = indexes[i];
  }

  /*sdatio_get_counts_and_starts(sfile, svar, counts, starts);*/

  sdatio_write_variable_private(sfile, svar, counts, starts, address);

  free(counts);
  free(starts);
}
/* Private*/
void sdatio_free_variable(struct sdatio_variable * svar){
  free(svar->name);
  free(svar->dimension_list);
  free(svar->dimension_ids);
  free(svar->manual_counts);
  free(svar->manual_starts);
  free(svar->manual_offsets);
  free(svar);
}

void sdatio_close(struct sdatio_file * sfile){
  int i, retval;

  if (!sfile->is_open) {
    printf("Attempting to close a file that has not been opened in sdatio_close\n");
    abort();
  }

  /*if (sfile->is_parallel){}*/
  /*else {*/

    if ((retval = nc_close(sfile->nc_file_id))) ERR(retval);
    /*}*/

  sfile->is_open = 0;
  for (i=0;i<sfile->n_dimensions;i++){
    sdatio_free_dimension(sfile->dimensions[i]);
  }
  free(sfile->dimensions);
  for (i=0;i<sfile->n_variables;i++){
    sdatio_free_variable(sfile->variables[i]);
  }
  if (sfile->n_variables > 0) free(sfile->variables);

}

void sdatio_sync(struct sdatio_file * sfile){
  int retval;

  /*if (sfile->is_parallel){}*/
  /*else {*/
  if ((retval = nc_sync(sfile->nc_file_id))) ERR(retval);
    /*}*/
}





