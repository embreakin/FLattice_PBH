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
        double _sfint = 0;
        double _a_save1 = 1;
    
    //Decay rate for all scalar fields  (program variables)
        double* Gamma_pr = new double [num_fields-1];
    
    
    //Evolution variables for scalar fields
        double **f_tilde, **df_tilde, **fdf_save;

        void evol_sfint(double h)
        {
            //Integration by using trapezoidal rule
            _sfint += ( _a + _a_save1 )*h*dt/2;
        }
    
        void evol_fields( double** f_tilde_from, double** f_tilde_to, double** df_tilde, double h );
    
        void evol_field_derivs_expansion(double** df_tilde_from, double** df_tilde_to, double** f, double** f_tilde, Field* field,  double h );
        
    
        void evol_scale_dderivs( Field* field, double** f, double** f_tilde, double& rho_r,  double h);
		void evol_scale_derivs ( double h ){ _da += _dda * h*dt; }
        void evol_scale ( double h ){
            
            _a_save1 = _a;
            _a += _da *h*dt;
            
//            std::cout << "_a = " << _a <<std::endl;
//            std::cout << "_da = " << _da <<std::endl;
//            std::cout << "_dda = " << _dda <<std::endl;
            
        }
    
        void evol_radiation(Field* field, double** f_tilde, double** df_tilde, double& rho,  double h);
    
        void evol_gravpot( double** f, double** df, double h );
        void evol_gravpot_derivs_expansion( double** f, double** df, double** f_tilde, double** df_tilde, Field* field, double h );
    
        void fields_copy( double** f_from, double** f_to);
        void fields_convert( double** f, double** f_tilde, int convert_switch);
        void fields_deriv_convert( double** f, double** df, double** f_tilde, double** df_tilde, int convert_switch);

    
	public:
    LeapFrog( Field* field, double** f, double** df, double& rad_pr );
    
    ~LeapFrog() {
        delete [] Gamma_pr;
        delete [] f_tilde[0];
        delete [] df_tilde[0];
        delete [] fdf_save[0];
        delete [] f_tilde;
        delete [] df_tilde;
        delete [] fdf_save;
    }

		void evolution_expansion ( Field* field, double** f, double** df, double& rad );

		double a()  { return _a; }
        double efolds() { return OSCSTART + log(_a); }
        double da() { return _da; }
        double hubble() { return _da*rescale_B*pow(_a,-2);}
        double adotdot() {return exp(OSCSTART)*pw2(rescale_B)*pow(_a,-2)*(_dda-pw2(_da)/_a);}
};


#endif
