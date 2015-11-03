#include "stdio.h"
#include "include/simpledataio.h"

int main (int argc, char ** argv){
	struct sdatio_file sdatfile;
	double yvar[2] = {0.1,0.3};
	float floatvar[2] = {0.1,0.3};
	int iy[2] = {1,2};
	double phivar[3][2] = {{0.1,0.3}, {2.0, 4.0}, {-1.0, 3.6}};
	double t;
	double phi_tvar[2];
	double parameter = 0.5;
	int i;
	int j = 4;

	sdatio_debug = 0;

	sdatio_createfile(&sdatfile, "testfile.cdf");

	sdatio_add_dimension(&sdatfile, "x", 3, "The x coordinate", "m");
	sdatio_add_dimension(&sdatfile, "y", 2, "The y coordinate", "m");
	sdatio_add_dimension(&sdatfile, "t", SDATIO_UNLIMITED, "The time coordinate", "s");

	sdatio_print_dimensions(&sdatfile);

	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "phi", "xy", "Some potential", "Vm");
	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "parameter", "", "A scalar parameter", "(none)");
	sdatio_create_variable(&sdatfile, SDATIO_FLOAT, "floatvar", "y", "A single precision variable.", "Vm");
	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "phi_t", "ty", "Some potential as a function of y and time", "Vm");
	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "y", "y", "Values of the y coordinate", "m");
	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "t", "t", "Values of the time coordinate", "m");
	/* This causes a warning */
	sdatio_write_variable(&sdatfile, "parameter", &parameter);

	sdatio_print_variables(&sdatfile);

	sdatio_create_variable(&sdatfile, SDATIO_INT, "iky", "y", "y index values", "(none)");

	sdatio_print_variables(&sdatfile);


	sdatio_write_variable(&sdatfile, "y", &yvar[0]);
	sdatio_write_variable(&sdatfile, "iky", &iy[0]);
	sdatio_write_variable(&sdatfile, "phi", &phivar[0]);
	sdatio_write_variable(&sdatfile, "floatvar", &floatvar[0]);

	for (i=0;i<6;i++){
		t = 0.3 + i;
		phi_tvar[0] = 4 + i/2.0;
		phi_tvar[1] = 6 + i*3.0; 
		sdatio_write_variable(&sdatfile, "t", &t);
		sdatio_write_variable(&sdatfile, "phi_t", &phi_tvar);
		sdatio_increment_start(&sdatfile, "t");
		/*if (i>2) abort();*/
	}

	int idxs[2] = {j,1};
	int idxs2[2] = {2,1};
	double val = 32.9;

	sdatio_write_variable_at_index(&sdatfile, "phi_t", idxs, &val);
	sdatio_write_variable_at_index_fast(&sdatfile, sdatio_find_variable(&sdatfile, "phi_t"), idxs2, &val);
	sdatio_sync(&sdatfile);

	struct sdatio_variable * svar;
	svar = sdatio_find_variable(&sdatfile, "floatvar");
	printf("The name of floatvar is %s\n", svar->name);

	struct sdatio_dimension * sdim;
	sdim = sdatio_find_dimension(&sdatfile, "t");
	printf("The name of t is %s\n", sdim->name);

	printf("This should be 1: %d\n", sdatio_variable_exists(&sdatfile, "t"));
	printf("This should be 1: %d\n", sdatio_variable_exists(&sdatfile, "phi_t"));
	printf("This should be 0: %d\n", sdatio_variable_exists(&sdatfile, "tbbb"));





	sdatio_close(&sdatfile);

	printf("Success!\n");
	return 0;
}
