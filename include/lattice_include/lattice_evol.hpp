#ifndef _LATTICEEVOL_H_
#define _LATTICEEVOL_H_

#include "lattice_field.hpp"
#include "lattice.hpp"

class LeapFrog
{
	private:
        double _a,_da, _dda;
        double sfexponent;
        double sfbase;
        double hubble_init = Hinitial_pr;
    
    //Decay rate for all scalar fields
        double* Gamma_pr = new double [num_fields-1];
    
    //Evolution variables for scalar fields
        double **f_tilde, **df_tilde;

		void evol_fields                 ( double** f, double** df, double h );
        void evol_field_derivs           ( double** f, double** df, Field* field, double h );
        void evol_field_derivs_expansion ( double** f, double** df, Field* field, double h );

		void evol_scale_dderivs        ( Field* field,  double** f,double h );
		void evol_scale_derivs ( double h ){ _da += .5*_dda * h*dt; }
        void evol_scale ( double h ){  _a += _da *h*dt; }
    
	public:
    LeapFrog( Field* field, double** f, double ** df, double rad_pr );
    
    ~LeapFrog() {
        delete [] Gamma_pr;
        delete [] f_tilde[0];
        delete [] df_tilde[0];
        delete [] f_tilde;
        delete [] df_tilde;
    }

        void evolution           ( Field* field, double** f, double** df );
		void evolution_expansion ( Field* field, double** f, double** df, double t );

		double a()  { return _a; }
		double da() { return _da; }
    double hubble() { return _da*rescale_B*pow(_a,-2);}
    double adotdot() {return pw2(rescale_B)*pow(_a,-2)*(_dda-pw2(_da)/_a);}
};


#endif
