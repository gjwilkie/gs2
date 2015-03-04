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
int sdatio_sdatio_variable_type(int type){
  switch (type){
  case NC_INT:
    return SDATIO_INT;
  case NC_FLOAT:
    return SDATIO_FLOAT;
  case NC_DOUBLE:
    return SDATIO_DOUBLE;
  case NC_CHAR:
    return SDATIO_CHAR;
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
void sdatio_get_dimension_list(struct sdatio_file * sfile, int * dimension_ids, struct sdatio_variable * svar){
  int ndims;
  int i,j, loop;
  /*char dim_name[2];*/
  /*char * dim_name;*/
  char * dimension_list;
  struct sdatio_dimension * sdim;
  int dim_name_length;
  int counter;
  int sep_size;
  ndims  = svar->ndims;
  DEBUG_MESS("ndims %d\n", ndims);
  /* In the first iteration, we calculate
   * the length of the dimension list string.
   * In the second iteration, we copy the names*/
  for (loop=0; loop<2; loop++){
    if (loop==1) dimension_list = (char *) malloc(sizeof(char)*(counter+1));
    counter = 0;
    for (i=0;i<ndims;i++){
      /* First we find the dimension object*/
      for (j=0;j<sfile->n_dimensions+1;j++){
        if (j==sfile->n_dimensions) {
          printf("ERROR: dimension with nc_id %d not found for variable %s\n", 
              dimension_ids[i], svar->name);
          abort();
        }
        if (sfile->dimensions[j]->nc_id==dimension_ids[i]) break;
      }
      sdim = sfile->dimensions[j];
      /* In the next section we set counter to 
       * the beginning of the dimension name and
       * dim_name_length to the size of the name*/
      if (sfile->has_long_dim_names){
        sep_size = 1;
        dim_name_length = strlen(sdim->name);
      }
      else {
        sep_size = 0;
        dim_name_length = 1;
        if (strlen(sdim->name)!=1){
          printf("Fatal logic error: dim_name_length!=1 for but not has_long_dim_names\n");
          abort();
        }
      }
      if (dim_name_length < 1){
        printf("ERROR: zero length dimension name in sdatio_open_file \n");
        abort();
      }
      DEBUG_MESS("In sdatio_get_dimension_list, dim id %d, has name %s, with counter %d\n",
          i, sdim->name, counter);

      /*dim_name[0] = dimension_list[i];*/
      /*dim_name[1] = dimension_list[ndims];*/

      /* Copy the name of the current dimension to the dimension_list*/
      if (loop==1) strncpy(dimension_list+counter, sdim->name, dim_name_length);
      counter = counter + dim_name_length;
      if (i<(ndims-1)){
        if (loop==1 && sep_size==1) strncpy(dimension_list+counter, ",", sep_size);
        counter = counter + sep_size;
      }
    }
  }
  dimension_list[counter] = '\0';

    /*dim_name[dim_name_length] = '\0';*/
  svar->dimension_list = dimension_list;


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


/* Has to be at the bottom because it uses many
 * other functions */
void sdatio_open_file(struct sdatio_file * sfile )  {
  /*printf("called\n");*/
  int retval;
  int ndims, nvars;
  int i, j;
  struct sdatio_dimension * sdim;
  struct sdatio_variable * svar;
  size_t lengthp;
  int nunlimdims;
  int *unlimdims;
  int is_unlimited;
  char * name_tmp;
  int * vardimids;
  char * dimension_list;
  nc_type vartype;
  int vartypeint;
  int dummy;
  int *nunlim;

  retval = 0;

  if (sfile->is_open){
    printf("ERROR: The supplied sdatio_file struct corresponds to an already open file.\n");
    abort();
  }

  /* We make the choice to always open the file in readwrite mode*/
  sfile->mode = sfile->mode|NC_WRITE;

  /* Open the file*/
  if (sfile->is_parallel) {
#ifdef PARALLEL 
    if ((retval = nc_open_par(sfile->name, sfile->mode, *(sfile->communicator), MPI_INFO_NULL,  &(sfile->nc_file_id)))) ERR(retval);
#else
    printf("sdatio was built without --enable-parallel, sdatio_create_file will not work for parallel files\n");
    abort();
#endif
  }
  else {
    if ((retval = nc_open(sfile->name, sfile->mode, &(sfile->nc_file_id)))) ERR(retval);
  }
  /*Initialize object file data*/
  sfile->is_open = 1;
  sfile->n_dimensions = 0;
  sfile->n_variables = 0;
  sfile->data_written = 0;
  /* Get number of dimensions in the file*/
  if ((retval = nc_inq_ndims(sfile->nc_file_id, &ndims))) ERR(retval);
  /* Allocate some temp storate*/
  name_tmp = (char*)malloc(sizeof(char*)*(NC_MAX_NAME+1));
  /* Get a list of unlimited dimensions*/
  if ((retval = nc_inq_unlimdims(sfile->nc_file_id, &nunlimdims, NULL))) ERR(retval);
  unlimdims = (int*)malloc(sizeof(int*)*nunlimdims);
  if ((retval = nc_inq_unlimdims(sfile->nc_file_id, &nunlimdims, unlimdims))) ERR(retval);
  /* Add each dimension to the sfile object*/
  for (i=0; i<ndims; i++){
    if ((retval = nc_inq_dim(sfile->nc_file_id, i, name_tmp, &lengthp))) ERR(retval);
    sdim = (struct sdatio_dimension *) malloc(sizeof(struct sdatio_dimension));
    /*if ((retval = nc_def_dim(sfile->nc_file_id, dimension_name, size, &(sdim->nc_id)))) ERR(retval);*/
    /*}*/
    /*sdatio_end_definitions(sfile);*/
    is_unlimited = 0;
    for(j=0; j<nunlimdims; j++) if (unlimdims[j] == i) is_unlimited = 1;
    if (is_unlimited) {
      sdim->size = SDATIO_UNLIMITED;
      /* We choose the first write to unlimited variables to be a new
       * record, so we set the length to be 1 greater than the current
       * final record.*/
      sdim->start = lengthp;
    }
    else {
      sdim->size = lengthp;
      sdim->start = 0;
    }
    if (strlen(name_tmp)>1){
      sfile->has_long_dim_names = 1;
      /*printf("Dimension names can only be one character long!\n");*/
      /*abort();*/
    }
    sdim->nc_id = i;
    sdim->name = (char *)malloc(sizeof(char)*(strlen(name_tmp)+1));
    strcpy(sdim->name, name_tmp);
    sdatio_append_dimension(sfile, sdim);
  }
  DEBUG_MESS("Finished reading dimensions\n");
  /* Get the number of variables in the file*/
  if ((retval = nc_inq_nvars(sfile->nc_file_id, &nvars))) ERR(retval);
  /* Add each variable to the sfile object*/
  for (i=0; i<nvars; i++){
    if ((retval = nc_inq_varndims(sfile->nc_file_id, i, &ndims))) ERR(retval);
    vardimids = (int*)malloc(sizeof(int)*ndims);
    if ((retval = nc_inq_var(sfile->nc_file_id, i, name_tmp, &vartype,
                            &ndims, vardimids, &dummy))) ERR(retval);
    vartypeint = vartype;
    vartypeint = sdatio_sdatio_variable_type(vartypeint);
    svar = (struct sdatio_variable *) malloc(sizeof(struct sdatio_variable));

    /*Set variable id*/
    svar->nc_id = i;

    /* Set variable name*/
    svar->name = (char *)malloc(sizeof(char)*(strlen(name_tmp)+1));
    strcpy(svar->name, name_tmp);

    /*ndims = strlen(dimension_list);*/
    svar->ndims = ndims;
    DEBUG_MESS("ndims = %d for variable %s\n", ndims, name_tmp);

    /* Set the dimension_ids*/
    svar->dimension_ids = vardimids;

    /*sdatio_get_dimension_ids(sfile, dimension_list, svar);*/
    sdatio_get_dimension_list(sfile, vardimids, svar);
    /*svar->dimension_ids = dimension_ids;*/


    DEBUG_MESS("Setting type for variable %s\n", svar->name);
    switch (vartypeint){
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

    
    svar->type = vartypeint;

    DEBUG_MESS("Allocating manual starts and counts for variable %s; ndims %d\n", svar->name, ndims);
    svar->manual_starts=(int*)malloc(sizeof(int)*ndims);
    svar->manual_counts=(int*)malloc(sizeof(int)*ndims);
    svar->manual_offsets=(int*)malloc(sizeof(int)*ndims);

    DEBUG_MESS("Setting manual starts and counts for variable %s\n", svar->name);
    for (j=0;j<ndims;j++){
      svar->manual_starts[j]=-1;
      svar->manual_counts[j]=-1;
      svar->manual_offsets[j]=-1;
    }

    DEBUG_MESS("Starting sdatio_append_variable\n");

    sdatio_append_variable(sfile, svar);

    DEBUG_MESS("Ending sdatio_append_variable\n");

#ifdef PARALLEL
    if (sfile->is_parallel){
      sdatio_number_of_unlimited_dimensions(sfile, svar->name, &nunlim);
      if (nunlim > 0)
        if ((retval = nc_var_par_access(sfile->nc_file_id, svar->nc_id, NC_COLLECTIVE))) ERR(retval);
    }
#endif
    /*vartypeint = */
  }
  DEBUG_MESS("Finished reading variables\n");
  
  free(unlimdims);
  free(name_tmp);
}




