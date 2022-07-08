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
        double* f_MPl;
        int fld;
    
  
    
	public:
    
        Field (): _faverage(new double [num_fields]()), _fvariance(new double [num_fields]()),_dfaverage(new double [num_fields]()), _dfvariance(new double [num_fields]()), f_MPl(new double [num_fields]())  {}
    
        ~Field () {
            delete [] _faverage;
            delete [] _fvariance;
            delete [] _dfaverage;
            delete [] _dfvariance;
            delete [] f_MPl;
        }
    
    void zeromode_initialize();
		double laplacian        ( double* f, int j, int k = 0, int l = 0 ); //You can omit k and l if they are zero
        double gradient_energy_eachpoint( double** f ,int i, int idx );
		double gradient_energy  ( double* f ) ;
		double potential_energy ( double** f, double a );
        double f_average  ( double* f, int i );
        double f_variance ( double* f, int i );
        double df_average ( double* df, int i );
        double df_variance ( double* df, int i );
        double V_lattice  ( double** f, int idx, double a = 1 );
        double dV_lattice ( double** f, int i, int idx, double a = 1);
        double mass( int i,  double a );
    
};

double Fk_log_int_calc(int k_int , double**& lattice_var, int num_field);
double dFk_log_int_calc(int k_int , double**& lattice_var, int num_field);
double omega_calc(double distance, double**& lattice_var, int num_field);
void initialize( double**& f, double**& df, Field* field, double radiation_pr, double**& lattice_var);
void finalize( double**& f, double**& df );


#endif
