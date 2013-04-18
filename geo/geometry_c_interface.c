/* This file exists to test and exemplify the interface
 * to the geometry module for a C/cuda program.
 * Maintainer: edmundhighcock@sourceforge.net */

#include <stdlib.h>
#include <stdio.h>

extern void geometry_set_inputs_(int * equilibrium_type,
		 											 char * eqfile,
													 int * irho,
		 											 double * rhoc, 
													 int * bishop,
													 int * nperiod,
													 int * ntheta_out);
 
/* These MUST be in the same order as the type declaration in
 * geometry.f90*/
struct advanced_parameters_struct {
	 	int equal_arc;
    double dp_mult;
    double delrho;                                                                
    double rmin;
    double rmax;
    int isym;
    int in_nt;
    int write_lots;
    int itor;
} ;

/* These MUST be in the same order as the type declaration in
 * geometry.f90*/
struct miller_parameters_struct { 
     double rmaj;
     double R_geo;
     double akappa;
     double akappri;
     double tri;
     double tripri;
     double shift;
     double qinp;
     double shat;
     double asym;
     double asympri;
};

struct coefficients_struct {
		 double  grho;   
		 double  bmag;       
		 double  gradpar;    
		 double  cvdrift;    
		 double  cvdrift0;   
		 double  gbdrift;    
		 double  gbdrift0;   
		 double  cdrift;    
		 double  cdrift0;    
		 double  gbdrift_th; 
		 double  cvdrift_th; 
		 double  gds2;       
		 double  gds21;      
		 double  gds22;      
		 double  gds23;      
		 double  gds24;      
		 double  gds24_noq;  
		 double  jacob;      
		 double  Rplot;      
		 double  Zplot;      
		 double  aplot;      
		 double  Rprime;     
		 double  Zprime;     
		 double  aprime;     
		 double  Uk1;        
		 double  Uk2;        
		 double  Bpol;       
};

/*void allocate_coefficients(int ntgrid){*/
/**/
/*}*/

extern void geometry_get_default_advanced_paramters_(
		struct advanced_parameters_struct * advanced_parameters_out);

extern void geometry_set_advanced_paramters_(
		struct advanced_parameters_struct * advanced_parameters_out);

extern void geometry_get_miller_parameters_(
		struct miller_parameters_struct * miller_parameters_out);
		
extern void geometry_set_miller_parameters_(
		struct miller_parameters_struct * miller_parameters_in);

extern void geometry_vary_s_alpha_(
		double * s_hat_input_in, double *beta_prime_input_in);

extern void geometry_get_coefficients_(int * ntheta, struct coefficients_struct * coefficients_out);

void geometry_vary_s_alpha_c(double * s_hat_input_in, double * beta_prime_input_in){
	geometry_vary_s_alpha_(s_hat_input_in, beta_prime_input_in);
}

void geometry_set_inputs_c(int * equilibrium_type,
		 											 char * eqfile,
													 int * irho,
		 											 double * rhoc, 
													 int * bishop,
													 int * nperiod,
													 int * ntheta_out){
	geometry_set_inputs_( equilibrium_type,
		 											  eqfile,
													  irho,
		 											  rhoc, 
													  bishop,
													  nperiod,
													  ntheta_out);
}


