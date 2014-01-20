#include <cstdlib>
#include "cufft.h"
#include "simpledataio_cuda.h"

int main (int argc, char * argv){
	struct sdatio_file sdatfile;
	double yvar[2] = {0.1,0.3};
	int iy[2] = {1,2};
	double phivar[3][2] = {{0.1,0.3}, {2.0, 4.0}, {-1.0, 3.6}};
	double t;
	double phi_tvar[2];
	int i;

	cuDoubleComplex compvar[3];

	sdatio_debug = 0;

	sdatio_createfile(&sdatfile, "testfile.cdf");

	sdatio_add_dimension(&sdatfile, "r", 2, "Real and imaginary parts", "(none)");
	sdatio_add_dimension(&sdatfile, "x", 3, "The x coordinate", "m");
	sdatio_add_dimension(&sdatfile, "y", 2, "The y coordinate", "m");
	sdatio_add_dimension(&sdatfile, "t", SDATIO_UNLIMITED, "The time coordinate", "s");
	sdatio_print_dimensions(&sdatfile);

	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "comp", "xr", "A complex variable", "(none)");
	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "phi", "xy", "Some potential", "Vm");
	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "phi_t", "ty", "Some potential as a function of y and time", "Vm");
	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "y", "y", "Values of the y coordinate", "m");
	sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "t", "t", "Values of the time coordinate", "m");
	sdatio_create_variable(&sdatfile, SDATIO_INT, "iky", "y", "y index values", "(none)");
	sdatio_print_variables(&sdatfile);


	for (i=0;i<6;i++){
		t = 0.3 + i;
		phi_tvar[0] = 4 + i/2.0;
		phi_tvar[1] = 6 + i*3.0; 
		sdatio_write_variable(&sdatfile, "t", &t);
		sdatio_write_variable(&sdatfile, "phi_t", &phi_tvar);
		sdatio_increment_start(&sdatfile, "t");
		//if (i>2) abort();
	}

	sdatio_write_variable(&sdatfile, "y", &yvar[0]);
	sdatio_write_variable(&sdatfile, "iky", &iy[0]);
	sdatio_write_variable(&sdatfile, "phi", &phivar[0]);

	compvar[0].x = 1.0;
	compvar[0].y = -1.0;
	compvar[1].x = 2.0;
	compvar[1].y = -2.0;
	compvar[2].x = 3.0;
	compvar[2].y = -3.0;

	sdatio_write_variable(&sdatfile, "comp", &compvar[0]);


	sdatio_close(&sdatfile);

	printf("Success!\n");
	return 0;
}
