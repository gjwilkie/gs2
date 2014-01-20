simpledataio
============

A simplified C and Fortran netcdf interface to allow writing to netcdf files with the fewest possible function calls.


Features
--------

1. Simple, easy to use interface
1. Supports C & Fortran
1. Supports parallel I/O using netcdf 4
1. Comprehensive examples (in the test folder)


Examples
--------

An example paints a thousand words: here is the minimal set of C calls to write a 2d parameter:


	#include "include/simpledataio.h"

	int main (int argc, char ** argv){
		double phivar[3][2] = {{0.1,0.3}, {2.0, 4.0}, {-1.0, 3.6}};
		struct sdatio_file sdatfile;

		sdatio_createfile(&sdatfile, "testfile.cdf");
		sdatio_add_dimension(&sdatfile, "x", 3, "The x coordinate", "m");
		sdatio_add_dimension(&sdatfile, "y", 2, "The y coordinate", "m");
		sdatio_create_variable(&sdatfile, SDATIO_DOUBLE, "phi", "xy", "Some potential", "Vm");
		sdatio_write_variable(&sdatfile, "phi", &phivar[0]);
		sdatio_close(&sdatfile);

		return 0;
	}


Here is the minium set of Fortran calls to do the same thing.


	program test
		use simpledataio
		use simpledataio_write
		implicit none
		type (sdatio_file) :: sdatfile
		double precision, dimension(3,2) ::  phivar = reshape((/0.1d0,2.0d0,-1.0d0, 0.3d0,4.0d0, 3.6d0/), (/3,2/))


		call createfile(sdatfile, "test.cdf")
		call add_dimension(sdatfile, "x", 3, "The x coordinate", "m")
		call add_dimension(sdatfile, "y", 2, "The y coordinate", "m")
		call create_variable(sdatfile, SDATIO_DOUBLE, "phi", "xy", "Some potential", "Vm")
		call write_variable(sdatfile, "phi", phivar)
		call closefile(sdatfile)

	end program test


Installing
----------

- Download the source from github

    	git clone https://github.com/edmundhighcock/simpledataio.git

- Configure, make and install...

		./configure
		make 
		make install

Options
-------

Configure options are:


- --enable-parallel ---  build with parallel netcdf (requires netcdf4, hdf5 and mpi)
- --enable-mpi ---  build with mpi version of compilers (automatically switched on by --enable-parallel)

You can also set the C and Fortran compilers used by specifying the FC and CC flags

E.g.

    ./configure --enable-parallel FC=mpifc CC=mpicc
 
If your netCDF is installed in a non-standard location you may need to add

    CFLAGS='-I/path/to/netcdf/include' FCFLAGS='-I/path/to/netcdf/include' LDFLAGS='-L/path/to/netcdf/lib'
