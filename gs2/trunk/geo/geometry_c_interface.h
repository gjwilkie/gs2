/* The .c file exists to get round the worst of the 
 * name mangling issues when using geo.a with CUDA.
 * This header file defines the structs that are
 * used in the function calls to the geometry library. 
 * Maintainer: edmundhighcock@sourceforge.net */

#include <stdlib.h>
#include <stdio.h>
 
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
     double  theta;
     double  theta_eqarc;
     double  theta_prime;
     double  theta_prime_eqarc;
		 double  grho;
		 double  grho_eqarc;   
		 double  bmag;
		 double  bmag_eqarc;       
		 double  gradpar;
		 double  gradpar_eqarc;    
		 double  gradpar_prime;
		 double  gradpar_prime_eqarc;    
		 double  cvdrift;
		 double  cvdrift_eqarc;    
		 double  cvdrift0;
		 double  cvdrift0_eqarc;   
		 double  gbdrift;
		 double  gbdrift_eqarc;    
		 double  gbdrift0;
		 double  gbdrift0_eqarc;   
		 double  cdrift;
		 double  cdrift_eqarc;    
		 double  cdrift0;
		 double  cdrift0_eqarc;    
		 double  gbdrift_th;
		 double  gbdrift_th_eqarc; 
		 double  cvdrift_th;
		 double  cvdrift_th_eqarc; 
		 double  gds2;
		 double  gds2_eqarc;       
		 double  gds21;
		 double  gds21_eqarc;      
		 double  gds22;
		 double  gds22_eqarc;      
		 double  gds23;
		 double  gds23_eqarc;      
		 double  gds24;
		 double  gds24_eqarc;      
		 double  gds24_noq;
		 double  gds24_noq_eqarc;  
		 double  jacob;
		 double  jacob_eqarc;      
		 double  Rplot;
		 double  Rplot_eqarc;      
		 double  Zplot;
		 double  Zplot_eqarc;      
		 double  aplot;
		 double  aplot_eqarc;      
		 double  Rprime;
		 double  Rprime_eqarc;     
		 double  Zprime;
		 double  Zprime_eqarc;     
		 double  aprime;
		 double  aprime_eqarc;     
		 double  Uk1;
		 double  Uk1_eqarc;        
		 double  Uk2;
		 double  Uk2_eqarc;        
		 double  Bpol;
		 double  Bpol_eqarc;       
};

struct constant_coefficients_struct {

    double qsf;
    double rmaj;
    double shat;
    double kxfac;
    double aminor;
    double drhodpsin;
    double bi;
};
