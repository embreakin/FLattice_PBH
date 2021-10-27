#ifndef _EVOULUTION_H_
#define _EVOULUTION_H_

#include "field.hpp"


class LeapFrog
{
	private:
        double _a,_da, _dda;
        const double hubble_init = Hinitial/m;
        double sfexponent;
        double sfbase;

		void evol_fields                 ( double** f, double** df, double h );
        void evol_field_derivs           ( double** f, double** df, Field* field, double h );
        void evol_field_derivs_expansion ( double** f, double** df, Field* field, double h );

		void evol_scale_dderivs        ( Field* field,  double** f,double h );
		void evol_scale_derivs ( double h ){ _da += .5*_dda * h*dt; }
        void evol_scale ( double h ){  _a += _da *h*dt; }
    
	public:
    LeapFrog( Field* field, double** f);

        void evolution           ( Field* field, double** f, double** df );
		void evolution_expansion ( Field* field, double** f, double** df, double t );

		double a()  { return _a; }
		double da() { return _da; }
    double hubble() { return _da*m*pow(_a,-2);}
    double adotdot() {return pw2(m)*pow(_a,-2)*(_dda-pw2(_da)/_a);}
};


#endif
