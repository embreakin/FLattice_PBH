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
    
    //Decay rate for all scalar fields  (program variables)
        double* Gamma_pr = new double [num_fields-1];
    
    //Scalar fields at the minimum of the potential (program variables)
     double* f_min_MPl = new double [num_fields-1];
    
    //Evolution variables for scalar fields
        double **f_tilde, **df_tilde, **fdf_save;

		void evol_fields                 ( double** f_tilde, double** df_tilde, double h );
        void evol_field_derivs_expansion ( double** f, double** f_tilde, double** df_tilde, Field* field, double t, double h );

        void evol_scale_dderivs( Field* field, double** f, double** f_tilde, double& rho_r, double t, double h);
		void evol_scale_derivs ( double h ){ _da += _dda * h*dt; }
        void evol_scale ( double h ){  _a += _da *h*dt; }
    
    void evol_radiation(Field* field, double** f_tilde, double** df_tilde, double& rho, double t , double h);
    
    void evol_gravpot( double** f, double** df, double h );
    void evol_gravpot_derivs_expansion( double** f, double** df, double** f_tilde, double** df_tilde, Field* field, double t, double h );
    
    void fields_copy( double** f_from, double** f_to);
    void fields_convert( double** f, double** f_tilde, double t, int convert_switch);
    void fields_deriv_convert( double** f, double** df, double** df_tilde, double t, int convert_switch);

    
	public:
    LeapFrog( Field* field, double** f, double** df, double& rad_pr );
    
    ~LeapFrog() {
        delete [] Gamma_pr;
        delete [] f_min_MPl;
        delete [] f_tilde[0];
        delete [] df_tilde[0];
        delete [] fdf_save[0];
        delete [] f_tilde;
        delete [] df_tilde;
        delete [] fdf_save;
    }

		void evolution_expansion ( Field* field, double** f, double** df, double& rad, double t );

		double a()  { return _a; }
        double efolds() { return OSCSTART + log(_a); }
        double da() { return _da; }
        double hubble() { return _da*rescale_B*pow(_a,-2);}
        double adotdot() {return exp(OSCSTART)*pw2(rescale_B)*pow(_a,-2)*(_dda-pw2(_da)/_a);}
};


#endif
