#include <cmath>

#include "calculation.hpp"
#include "parameter.hpp"

void Energy::energy_calc( Field* field, LeapFrog* leapfrog, double** f, double** df )
{
    double a  = leapfrog->a ();
    double da = leapfrog->da();
   
    _total_average = 0;
    
    for( int i = 0; i < num_fields; ++i ){

        _average[i] = 0;
        _potential_average = 0;
        _timederiv_average = 0;
        _grad_average = 0;
        _variance[i] = 0;
        _value_max = 0;
        //        #if   dim == 1
        //        #pragma omp parallel for simd schedule(static) num_threads(num_threads)
        //        #elif dim >= 2
        //        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        //        #endif
        //        #pragma omp parallel for schedule(static) num_threads (num_threads)
        
        for( int j = 0; j < N; ++j ){
            if( expansion ) //1,2,3...
            {
                switch( dim ){
                    case 1:
                        int idx = j;
                        value[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6)) + field->aV(f, i, idx ,a)+gradient_energy_eachpoint(f, i, idx)/ pow(a,4);
                        _average[i]+= value[i][idx];
                        _potential_average += field->aV(f, i, idx ,a);
                        _timederiv_average +=  pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6));
                        _grad_average += gradient_energy_eachpoint(f, i, idx)/ pow(a,4);
                        break;
                    case 2:
//                        #pragma omp simd
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                          
                            value[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6)) + field->aV(f, i, idx ,a)+gradient_energy_eachpoint(f, i, idx)/ pow(a,4);
                           
                            _average[i]+= value[i][idx];
                            _potential_average += field->aV(f, i, idx ,a);
                            _timederiv_average +=  pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6));
                            _grad_average += gradient_energy_eachpoint(f, i, idx)/ pow(a,4);
                        
                        }
                        break;
                    case 3:
                        for( int k = 0; k < N; ++k ){
//                            #pragma omp simd
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;
                                value[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6)) + field->aV(f, i, idx ,a)+gradient_energy_eachpoint(f, i, idx)/pow(a,4);
                                _average[i]+= value[i][idx];
                                _potential_average += field->aV(f, i, idx ,a);
                                _timederiv_average +=  pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6));
                                _grad_average += gradient_energy_eachpoint(f, i, idx)/ pow(a,4);
                                
                            }
                        }
                        break;
                }
            }
            else // 0
            {
                switch( dim ){
                    case 1:
                        int idx = j;
                        value[i][idx] = pow(df[i][idx], 2)/2 + field->V(f, i ,idx)+gradient_energy_eachpoint(f, i, idx);
                        _average[i]+= value[i][idx];
                        _potential_average += field->V(f, i ,idx);
                        _timederiv_average +=  pow(df[i][idx], 2)/2;
                        _grad_average += gradient_energy_eachpoint(f, i, idx);
                        break;
                    case 2:
//                        #pragma omp simd
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                            value[i][idx] = pow(df[i][idx], 2)/2 + field->V(f, i ,idx)+gradient_energy_eachpoint(f, i, idx);
                            _average[i] += value[i][idx];
                            _potential_average += field->V(f, i ,idx);
                            _timederiv_average +=  pow(df[i][idx], 2)/2;
                            _grad_average += gradient_energy_eachpoint(f, i, idx);
                        }
                        break;
                    case 3:
                        for( int k = 0; k < N; ++k ){
//                            #pragma omp simd
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;
                                value[i][idx] = pow(df[i][idx], 2)/2 + field->V(f, i ,idx)+gradient_energy_eachpoint(f, i, idx);
                                _average[i]+= value[i][idx];
                                _potential_average += field->V(f, i ,idx);
                                _timederiv_average +=  pow(df[i][idx], 2)/2;
                                _grad_average += gradient_energy_eachpoint(f, i, idx);
                            }
                        }
                        break;
                }
            }
        }// j for loop
    

        for( int j = 0; j < dim; ++j ) _potential_average /= N;
        for( int j = 0; j < dim; ++j ) _timederiv_average /= N;
        for( int j = 0; j < dim; ++j ) _grad_average /= N;
        for( int j = 0; j < dim; ++j ) _average[i] /= N;
        
        _total_average += _average[i];
  
    //devide by average
//    #if   dim == 1
//    #pragma omp parallel for simd schedule(static) num_threads(num_threads)
//    #elif dim >= 2
//    #pragma omp parallel for schedule( static ) num_threads( num_threads )
//    #endif
//    #pragma omp parallel for schedule(static) num_threads (num_threads)
    for( int j = 0; j < N; ++j ){
        switch( dim ){
        case 1:
            
            int idx = j;
            
            value[i][idx] /= _average[i] ;
            _variance[i] += pow( value[i][idx]  - 1 , 2 );
            
            if(value[i][idx] > _value_max){
                _value_max = value[i][idx];
            }
            
            break;
        case 2:
//            #pragma omp simd
            for( int k = 0; k < N; ++k ){
                int idx = j*N + k;
           
                value[i][idx] /= _average[i] ;
                _variance[i] += pow( value[i][idx]  - 1 , 2 );
              
                if(value[i][idx] > _value_max){
                    _value_max = value[i][idx];
                }
                
            }
            break;
        case 3:
            for( int k = 0; k < N; ++k ){
//                #pragma omp simd
                for( int l = 0; l < N; ++l ){
                    int idx = (j*N + k)*N + l;
                    value[i][idx] /= _average[i] ;
                    _variance[i] += pow( value[i][idx]  - 1 , 2 );
                    
                    if(value[i][idx] > _value_max){
                        _value_max = value[i][idx];
                    }
                }
            }
            break;
                   }
        }
        
        for( int j = 0; j < dim; ++j ) _variance[i] /= N;
      
    }
}

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
