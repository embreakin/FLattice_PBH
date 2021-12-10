#ifndef _LATTICEFIELD_H_
#define _LATTICEFIELD_H_

#include <cmath>
#include "parameters.hpp"



class Field
{
	private:
		double* _faverage;
		double* _fvariance;
        double* _dfaverage;
        double* _dfvariance;
        double* tauvalue;
        double TAU;
        double aTAU;
        int fld;
    
  
    
	public:
    
        Field (): _faverage(new double [num_fields]()), _fvariance(new double [num_fields]()),_dfaverage(new double [num_fields]()), _dfvariance(new double [num_fields]()), tauvalue(new double [num_fields]())  {}
    
        ~Field () {
            delete [] _faverage;
            delete [] _fvariance;
            delete [] _dfaverage;
            delete [] _dfvariance;
            delete [] tauvalue;
        }
    
    void zeromode_initialize();
		double laplacian        ( double* f, int j, int k = 0, int l = 0 ); //You can omit k and l if they are zero
		double gradient_energy  ( double* f ) ;
		double potential_energy ( double** f, double a );
        double f_average  ( double* f, int i );
        double f_variance ( double* f, int i );
        double df_average ( double* df, int i );
        double df_variance ( double* df, int i );
        double V   ( double** f, int i, int idx ) ;
        double dV  ( double** f, int i, int idx );
        double aV  ( double** f, int i, int idx, double a );
        double adV ( double** f, int i, int idx, double a );
        void effective_mass( double* mass_sq, double* field_values);
    
};

void initialize( double**& f, double**& df, Field* field, double**& lattice_var);
void finalize( double**& f, double**& df );


#endif
