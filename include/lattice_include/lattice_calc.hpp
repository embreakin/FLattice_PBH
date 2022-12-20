//Doxygen
/**
* @file   lattice_calc.hpp
* @brief    Lattice calclulation header file
* @author   Francis Otani
* @date
* @details
*/

#ifndef _LATTICECALC_H_
#define _LATTICECALC_H_

//#include <valarray>

#include "lattice_field.hpp"
#include "lattice_evol.hpp"


class Energy
{
    protected:
    double _variance;
    double _total_average;
    double _potential_average;
    double _timederiv_average;
    double _grad_average;
    double _rad;
    double _pressure_average;
    double _value_max;
    
    
    public:
    
    double *value;
    
    
    Energy():  _total_average(),_potential_average(),_timederiv_average(),_grad_average(),_rad(),_pressure_average(),_value_max()
    {
        
        switch( dim ){
            case 1:
                value = new double [N]();//() is for initializing to zero
                break;
            case 2:
                value = new double [N*N]();
                break;
            case 3:
                value = new double [N*N*N]();
                break;
        }
    }
    
    ~Energy(){
        delete [] value;
        
    }

        double total_variance() { return sqrt(_variance); }
        double total_average   () { return _total_average*pw2(rescale_B/rescale_A); }
        double potential_average () { return _potential_average*pw2(rescale_B/rescale_A);}
        double timederiv_average () { return _timederiv_average*pw2(rescale_B/rescale_A);}
        double grad_average () { return _grad_average*pw2(rescale_B/rescale_A);}
        double radiation () { return _rad*pw2(rescale_B/rescale_A);}
        double pressure (){ return _pressure_average*pw2(rescale_B/rescale_A); }
        double energy_max () {return _value_max;}
    
        void energy_calc( Field* field, LeapFrog* leapfrog, double** f, double** df, double& rad );
        double gradient_energy_eachpoint( double** f ,int i, int idx );
    
       #pragma omp declare simd
        double kinetic_energy_eachpoint( double** f , double** df, int i, int idx, double a = 1, double da = 0 )
        {   return pow(df[i][idx]*a - f[i][idx]*da, 2.0)/(2*pow(a,2.0)); }

};

#endif
