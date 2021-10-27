#ifndef _CALCULATION_H_
#define _CALCULATION_H_

//#include <valarray>

#include "field.hpp"
#include "evolution.hpp"


class Energy
{
    protected:
        double* _average;
        double* _variance;
        double _total_average;
        double _potential_average;
        double  _timederiv_average;
        double _grad_average;
        double _value_max;
    
    
    public:
    
    double **value;
    
    Energy(): _average(new double [num_fields]()), _variance(new double [num_fields]()), _total_average(),_potential_average(),_timederiv_average(),_grad_average(),_value_max()
    {
        value = new double* [num_fields];
        switch( dim ){
            case 1:
                value[0] = new double [num_fields*N]();//() is for initializing to zero
                for( int i = 0; i < num_fields; ++i ) value[i] = value[0] + i*N;
                break;
            case 2:
                value[0] = new double [num_fields*N*N]();
                for( int i = 0; i < num_fields; ++i ) value[i] = value[0] + i*N*N;
                break;
            case 3:
                value[0] = new double [num_fields*N*N*N]();
                for( int i = 0; i < num_fields; ++i ) value[i] = value[0] + i*N*N*N;
                break;
        }
    }
    
    ~Energy(){
        delete [] value[0];
        delete [] value;
        delete [] _average;
        delete [] _variance;
        
    }

        double average  ( int i ) { return _average[i]; }
        double variance ( int i ) { return sqrt(_variance[i]); }
        double total_average   () { return _total_average*pw2(m)*pow(sqrt(8*M_PI),4); }
        double potential_average () { return _potential_average*pw2(m)*pow(sqrt(8*M_PI),4);}
        double timederiv_average () { return _timederiv_average*pw2(m)*pow(sqrt(8*M_PI),4);}
        double grad_average () { return _grad_average*pw2(m)*pow(sqrt(8*M_PI),4);}
        double energy_max () {return _value_max;}
    
        void energy_calc( Field* field, LeapFrog* leapfrog, double** f, double** df );
        double gradient_energy_eachpoint( double** f ,int i, int idx );

};

#endif
