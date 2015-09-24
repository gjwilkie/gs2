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

struct constant_coefficients_struct {

    double qsf;
    double rmaj;
    double shat;
    double kxfac;
    double aminor;
};
