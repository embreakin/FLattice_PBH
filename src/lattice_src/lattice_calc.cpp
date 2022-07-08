#include <cmath>

#include "lattice_calc.hpp"
#include "parameters.hpp"

void Energy::energy_calc( Field* field, LeapFrog* leapfrog, double** f, double** df, double rad )
{
    double a  = leapfrog->a ();
    double da = leapfrog->da();
    
    _total_average = 0;
    
        
        // These variable need to be declared for OpenMP reduction directive
        double _average_reduc = 0;
        double _potential_average_reduc = 0;
        double _timederiv_average_reduc = 0;
        double _grad_average_reduc = 0;
        double _variance_reduc = 0;
        double _value_max_reduc = 0;
        
        

    #if   dim == 1
    #pragma omp parallel for simd reduction(+:_average_reduc,_potential_average_reduc,_timederiv_average_reduc,_grad_average_reduc) schedule(static) num_threads(num_threads)
    #elif dim >= 2
    #pragma omp parallel for reduction(+:_average_reduc,_potential_average_reduc,_timederiv_average_reduc,_grad_average_reduc) schedule(static) num_threads(num_threads)
    #endif
        for( int j = 0; j < N; ++j ){
            
            #if   dim == 1
                        int idx = j;
            
                        value[idx] = kinetic_energy_eachpoint(f, df, 0, idx, a, da)
                        +kinetic_energy_eachpoint(f, df, 1, idx, a, da)
                        +kinetic_energy_eachpoint(f, df, 2, idx, a, da)
                        +
                        gradient_energy_eachpoint(f, 0, idx)/(pw2(exp(OSCSTART)))
                        +gradient_energy_eachpoint(f, 1, idx)/(pw2(exp(OSCSTART)))
                        +gradient_energy_eachpoint(f, 2, idx)/(pw2(exp(OSCSTART)))
                        + field->V_lattice(f, idx ,a)
                        + rad;
            
                        _average_reduc += value[idx];
            
                        _timederiv_average_reduc +=
                        kinetic_energy_eachpoint(f, df, 0, idx, a, da)
                        +kinetic_energy_eachpoint(f, df, 1, idx, a, da)
                        +kinetic_energy_eachpoint(f, df, 2, idx, a, da);
            
                        _grad_average_reduc +=
                        gradient_energy_eachpoint(f, 0, idx)/(pw2(exp(OSCSTART)))
                        +gradient_energy_eachpoint(f, 1, idx)/(pw2(exp(OSCSTART)))
                        +gradient_energy_eachpoint(f, 2, idx)/(pw2(exp(OSCSTART)));
            
                        _potential_average_reduc += field->V_lattice(f, idx ,a);
            
            #elif   dim == 2
//                #pragma omp simd reduction(+:_average_reduc,_potential_average_reduc,_timederiv_average_reduc,_grad_average_reduc)
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                            
                         
                            value[idx] = kinetic_energy_eachpoint(f, df, 0, idx, a, da)
                            +kinetic_energy_eachpoint(f, df, 1, idx, a, da)
                            +kinetic_energy_eachpoint(f, df, 2, idx, a, da)
                            +
                            gradient_energy_eachpoint(f, 0, idx)/(pw2(exp(OSCSTART)))
                            +gradient_energy_eachpoint(f, 1, idx)/(pw2(exp(OSCSTART)))
                            +gradient_energy_eachpoint(f, 2, idx)/(pw2(exp(OSCSTART)))
                            + field->V_lattice(f, idx ,a)
                            + rad;
                            
                            _average_reduc += value[idx];
                            
                            _timederiv_average_reduc +=
                            kinetic_energy_eachpoint(f, df, 0, idx, a, da)
                            +kinetic_energy_eachpoint(f, df, 1, idx, a, da)
                            +kinetic_energy_eachpoint(f, df, 2, idx, a, da);
                            
                            _grad_average_reduc +=
                            gradient_energy_eachpoint(f, 0, idx)/(pw2(exp(OSCSTART)))
                            +gradient_energy_eachpoint(f, 1, idx)/(pw2(exp(OSCSTART)))
                            +gradient_energy_eachpoint(f, 2, idx)/(pw2(exp(OSCSTART)));
                            
                            _potential_average_reduc += field->V_lattice(f, idx ,a);
                            
                        }
            #elif   dim == 3
            for( int k = 0; k < N; ++k ){
//              #pragma omp simd reduction(+:_average_reduc,_potential_average_reduc,_timederiv_average_reduc,_grad_average_reduc)
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;

                              
                                value[idx] = kinetic_energy_eachpoint(f, df, 0, idx, a, da)
                                +kinetic_energy_eachpoint(f, df, 1, idx, a, da)
                                +kinetic_energy_eachpoint(f, df, 2, idx, a, da)
                                +
                                gradient_energy_eachpoint(f, 0, idx)/(pw2(exp(OSCSTART)))
                                +gradient_energy_eachpoint(f, 1, idx)/(pw2(exp(OSCSTART)))
                                +gradient_energy_eachpoint(f, 2, idx)/(pw2(exp(OSCSTART)))
                                + field->V_lattice(f, idx ,a)
                                + rad;
                                
                                _average_reduc += value[idx];
                                
                                _timederiv_average_reduc +=
                                kinetic_energy_eachpoint(f, df, 0, idx, a, da)
                                +kinetic_energy_eachpoint(f, df, 1, idx, a, da)
                                +kinetic_energy_eachpoint(f, df, 2, idx, a, da);
                                
                                _grad_average_reduc +=
                                gradient_energy_eachpoint(f, 0, idx)/(pw2(exp(OSCSTART)))
                                +gradient_energy_eachpoint(f, 1, idx)/(pw2(exp(OSCSTART)))
                                +gradient_energy_eachpoint(f, 2, idx)/(pw2(exp(OSCSTART)));
                                
                                _potential_average_reduc += field->V_lattice(f, idx ,a);
                            }
                            }

             #endif
            }
        

        
        for( int j = 0; j < dim; ++j ) _potential_average_reduc /= N;
        for( int j = 0; j < dim; ++j ) _timederiv_average_reduc /= N;
        for( int j = 0; j < dim; ++j ) _grad_average_reduc /= N;
        for( int j = 0; j < dim; ++j ) _average_reduc /= N;
    
        
        for( int j = 0; j < N; ++j ){
            switch( dim ){
                case 1:
                    
                    int idx = j;
                    
                    value[idx] /= _average_reduc ;
                    _variance_reduc += pow( value[idx]  - 1 , 2 );
                    
                    if(value[idx] > _value_max_reduc){
                        _value_max_reduc = value[idx];
                    }
                    
                    break;
                case 2:
                    //            #pragma omp simd
                    for( int k = 0; k < N; ++k ){
                        int idx = j*N + k;
                        
                        value[idx] /= _average_reduc ;
                        _variance_reduc += pow( value[idx]  - 1 , 2 );
                        
                        if(value[idx] > _value_max_reduc){
                            _value_max_reduc = value[idx];
                        }
                        
                    }
                    break;
                case 3:
                    for( int k = 0; k < N; ++k ){
                        //                #pragma omp simd
                        for( int l = 0; l < N; ++l ){
                            int idx = (j*N + k)*N + l;
                            value[idx] /= _average_reduc ;
                            _variance_reduc += pow( value[idx]  - 1 , 2 );
                            
                            if(value[idx] > _value_max_reduc){
                                _value_max_reduc = value[idx];
                            }
                        }
                    }
                    break;
                default:
                    break;
            }
        }
        
        for( int j = 0; j < dim; ++j ) _variance_reduc /= N;
        
        //substitute the obtained variables to member variables
        _total_average = _average_reduc*pow(a,-4);
        _potential_average = _potential_average_reduc*pow(a,-4);
        _timederiv_average = _timederiv_average_reduc*pow(a,-4);
        _grad_average = _grad_average_reduc*pow(a,-4);
        _variance = _variance_reduc;
        _value_max = _value_max_reduc;
        
    
}

#pragma omp declare simd
double Energy::gradient_energy_eachpoint( double** f ,int i, int idx )
{
   
#if dim == 1
    int j= idx;
    int jp1 = (j == N-1)?     0: j+1;
    int jp2 = (j >= N-2)? j-N+2: j+2;

    int jm1 = (j == 0)?   N-1: j-1;
    int jm2 = (j <  2)? j+N-2: j-2;
    
    return  pow( ( - f[i][jp2]  + 8*f[i][jp1]  - 8*f[i][jm1] + f[i][jm2] ) / (12*dx), 2 )/2;

#elif dim == 2
    int j = idx / N;
    int jp1 = (j == N-1)?     0: j+1;
    int jp2 = (j >= N-2)? j-N+2: j+2;
    int jm1 = (j == 0)?   N-1: j-1;
    int jm2 = (j <  2)? j+N-2: j-2;
    int k = idx % N;
    int kp1 = (k == N-1)?     0: k+1;
    int kp2 = (k >= N-2)? k-N+2: k+2;
    int km1 = (k ==   0)?   N-1: k-1;
    int km2 = (k <    2)? k+N-2: k-2;
  
//    return (pow( (f[i][jp1*N+k] - f[i][jm1*N+k]) / (2*dx), 2 )
//            + pow( (f[i][j*N+kp1] - f[i][j*N+km1]) / (2*dx), 2 ))/2;
    return (pow( (- f[i][jp2*N+k] + 8*f[i][jp1*N+k] - 8*f[i][jm1*N+k] + f[i][jm2*N+k]) / (12*dx), 2 )
            + pow( (- f[i][j*N+kp2] + 8*f[i][j*N+kp1] - 8*f[i][j*N+km1] + f[i][j*N+km2]) / (12*dx), 2 ))/2;

#elif dim == 3
    int j = idx /(N*N);
    int jp1 = (j == N-1)?     0: j+1;
    int jp2 = (j >= N-2)? j-N+2: j+2;
    int jm1 = (j == 0)?   N-1: j-1;
    int jm2 = (j <  2)? j+N-2: j-2;
    int k =(idx %(N*N))/N;
        int kp1 = (k == N-1)?     0: k+1;
        int kp2 = (k >= N-2)? k-N+2: k+2;
        int km1 = (k ==   0)?   N-1: k-1;
        int km2 = (k <    2)? k+N-2: k-2;
    int l =(idx %(N*N))% N;
            int lp1 = (l == N-1)?     0: l+1;
            int lp2 = (l >= N-2)? l-N+2: l+2;
            int lm1 = (l ==   0)?   N-1: l-1;
            int lm2 = (l <    2)? l+N-2: l-2;
       return  (pow( (- f[i][(jp2*N+k)*N+l] + 8*f[i][(jp1*N+k)*N+l] - 8*f[i][(jm1*N+k)*N+l] + f[i][(jm2*N+k)*N+l]) / (12*dx), 2 )
             + pow( (- f[i][(j*N+kp2)*N+l] + 8*f[i][(j*N+kp1)*N+l] - 8*f[i][(j*N+km1)*N+l] + f[i][(j*N+km2)*N+l]) / (12*dx), 2 )
             + pow( (- f[i][(j*N+k)*N+lp2] + 8*f[i][(j*N+k)*N+lp1] - 8*f[i][(j*N+k)*N+lm1] + f[i][(j*N+k)*N+lm2]) / (12*dx), 2 ))/2;
    
#endif
    
    
}


