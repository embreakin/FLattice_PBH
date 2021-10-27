#ifndef _FIELD_H_
#define _FIELD_H_

#include <cmath>
#include "parameter.hpp"



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
        double POTCOEF = pw2(1/m);
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
    
		double laplacian        ( double* f, int j, int k = 0, int l = 0 );
		double gradient_energy  ( double* f ) ;
		double potential_energy ( double** f, double a );
        double f_average  ( double* f, int i );
        double f_variance ( double* f, int i );
        double df_average ( double* df, int i );
        double df_variance ( double* df, int i );

		double V   ( double** f, int i, int idx ) {
            
            TAU = exp((sqrt(2)/sqrt(3))*f[i][idx]*sqrt(8*M_PI))/2;
           
            
            return (POTCOEF*pow(sqrt(8*M_PI),-4)*(exp(-2*aa*TAU)/(pw2(TAU))*aa*pw2(AA)/2
            +exp(-2*aa*TAU)/(TAU)*pw2(aa)*pw2(AA)/6
            +exp(-aa*TAU)/(pw2(TAU))*aa*AA*W0/2
            +1/(pow(TAU,3))*D)
                    );
        }
		double dV  ( double** f, int i, int idx ){
            
            TAU = exp((sqrt(2)/sqrt(3))*f[i][idx]*sqrt(8*M_PI))/2;
          
            
            return (POTCOEF*pow(sqrt(8*M_PI),-3)*(sqrt(2)/sqrt(3))*TAU*(-(7*pw2(aa)*pw2(AA)/(6*pw2(TAU)))*exp(-2*aa*TAU)
                -(aa*pw2(AA)/(pow(TAU,3)))*exp(-2*aa*TAU)
            -(pow(aa,3)*pw2(AA)/(3*TAU))*exp(-2*aa*TAU)
            -(pw2(aa)*AA*W0/(2*pw2(TAU)))*exp(-aa*TAU)
            -(aa*AA*W0/(pow(TAU,3)))*exp(-aa*TAU)
            -3*D/(pow(TAU,4)))
            );
        }
    
    void effective_mass( double mass_sq[], double *field_values){
        for(fld=0;fld<num_fields;fld++)
        {
        tauvalue[fld] = exp((sqrt(2)/sqrt(3))*field_values[fld]*sqrt(8*M_PI))/2;
             std::cout << "tauvalue[0] = " << tauvalue[fld] << std::endl;
       
            mass_sq[fld]=POTCOEF*pow(sqrt(8*M_PI),-2)*pw2(sqrt(2)/sqrt(3))*(
        (tauvalue[fld])*(-(7*pw2(aa)*pw2(AA)/(6*pw2(tauvalue[fld])))*exp(-2*aa*tauvalue[fld])
            -(aa*pw2(AA)/(pow(tauvalue[fld],3)))*exp(-2*aa*tauvalue[fld])
            -(pow(aa,3)*pw2(AA)/(3*tauvalue[fld]))*exp(-2*aa*tauvalue[fld])
            -(pw2(aa)*AA*W0/(2*pw2(tauvalue[fld])))*exp(-aa*tauvalue[fld])
        -(aa*AA*W0/(pow(tauvalue[fld],3)))*exp(-aa*tauvalue[fld])
        -3*D/(pow(tauvalue[fld],4)))  +
       pw2(tauvalue[fld])*(
                           exp(-2*aa*tauvalue[fld])*(
        (2*pow(aa,4)*pw2(AA)/(3*tauvalue[fld]))+
        (8*pow(aa,3)*pw2(AA)/(3*pw2(tauvalue[fld])))+
       (13*pw2(aa)*pw2(AA)/(3*pow(tauvalue[fld],3)))+
        (3*aa*pw2(AA)/(pow(tauvalue[fld],4))))
                           +
        exp(-aa*tauvalue[fld])*(
        (pow(aa,3)*AA*W0/(2*pw2(tauvalue[fld])))+
        (2*pw2(aa)*AA*W0/(pow(tauvalue[fld],3)))+
        (3*aa*AA*W0/(pow(tauvalue[fld],4))))
                           +
         12*D/(pow(tauvalue[fld],5))
        )
                                                                         
        );
        }
    }
    
    double aV  ( double** f, int i, int idx, double a )  {
           // std::cout << "a = " << a << std::endl;
            aTAU = exp((sqrt(2)/sqrt(3))*f[i][idx]*sqrt(8*M_PI)/a)/2;
        
            
            return (POTCOEF*pow(sqrt(8*M_PI),-4)*(exp(-2*aa*aTAU)/(pw2(aTAU))*aa*pw2(AA)/2
                    +exp(-2*aa*aTAU)/(aTAU)*pw2(aa)*pw2(AA)/6
                    +exp(-aa*aTAU)/(pw2(aTAU))*aa*AA*W0/2
                    +1/(pow(aTAU,3))*D)
                    ); }
    
    
		double adV ( double** f, int i, int idx, double a )  {

            aTAU = exp((sqrt(2)/sqrt(3))*f[i][idx]*sqrt(8*M_PI)/a)/2;           

            return (POTCOEF*pow(sqrt(8*M_PI),-3)*(sqrt(2)/sqrt(3))*aTAU*(-(7*pw2(aa)*pw2(AA)/(6*pw2(aTAU)))*exp(-2*aa*aTAU)
                -(aa*pw2(AA)/(pow(aTAU,3)))*exp(-2*aa*aTAU)
                -(pow(aa,3)*pw2(AA)/(3*aTAU))*exp(-2*aa*aTAU)
                -(pw2(aa)*AA*W0/(2*pw2(aTAU)))*exp(-aa*aTAU)
                -(aa*AA*W0/(pow(aTAU,3)))*exp(-aa*aTAU)
                -3*D/(pow(aTAU,4)))
                    ); }
    
};

void initialize( double**& f, double**& df, Field* field);
void finalize( double**& f, double**& df );


#endif
