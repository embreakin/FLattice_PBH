#ifndef _CALCULATION_H_
#define _CALCULATION_H_

//#include <valarray>

#include "field.hpp"
#include "evolution.hpp"


class Calculation
{
    protected:
        double* _average;
        double* _variance;
        double _total_average;
        double _potential_average;
        double  _timederiv_average;
        double _grad_average;
        double _df_average;
        double _f_average;
    public:
        double **value, *total_value;
        Calculation();
       ~Calculation();

        double average  ( int i ) { return _average[i]; }
        double variance ( int i ) { return _variance[i]; }
        double total_average   () { return _total_average*pw2(m)*pow(sqrt(8*M_PI),4); }
        double potential_average () { return _potential_average*pw2(m)*pow(sqrt(8*M_PI),4);}
        double timederiv_average () { return _timederiv_average*pw2(m)*pow(sqrt(8*M_PI),4);}
        double grad_average () { return _grad_average*pw2(m)*pow(sqrt(8*M_PI),4);}
        double df_average () { return _df_average; }
};


class Energy: public Calculation
{
    private:
        
    public:
        Energy( Field* field, LeapFrog* leapfrog, double** f, double** df );
        double gradient_energy( double* f );
        double gradient_energy_eachpoint( double** f ,int i, int idx );
};


#endif
