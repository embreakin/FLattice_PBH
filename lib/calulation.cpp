#include <cmath>

#include "calculation.hpp"
#include "parameter.hpp"


Calculation::Calculation(): _total_average(),_potential_average(), _average(new double [num_fields]()), _variance(new double [num_fields]())
{
    value = new double* [num_fields];
	switch( dim ){
		case 1:
			value[0] = new double [num_fields*N]();//() is for initializing to zero
          //  total_value = new double [N]();
			for( int i = 0; i < num_fields; ++i ) value[i] = value[0] + i*N;
			break;
		case 2:
			value[0] = new double [num_fields*N*N]();
         //   total_value = new double [N*N]();
			for( int i = 0; i < num_fields; ++i ) value[i] = value[0] + i*N*N;
			break;
		case 3:
			value[0] = new double [num_fields*N*N*N]();
          //  total_value = new double [N*N*N]();
			for( int i = 0; i < num_fields; ++i ) value[i] = value[0] + i*N*N*N;
			break;
	}
}

Calculation::~Calculation()
{
    delete [] value[0];
    delete [] value;
  
}


Energy::Energy( Field* field, LeapFrog* leapfrog, double** f, double** df ): Calculation()
{
    double a  = leapfrog->a ();
    double da = leapfrog->da();
   
   
     //std::cout << "a = " << a << std::endl;
    for( int i = 0; i < num_fields; ++i ){
       // #pragma omp parallel for reduction(+: _potential_average,_average[i] ) schedule( static ) num_threads ( num_threads ) <- causes trouble
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
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                          
                            value[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6)) + field->aV(f, i, idx ,a)+gradient_energy_eachpoint(f, i, idx)/ pow(a,4);
                           
                            _average[i]+= value[i][idx];
                            _potential_average += field->aV(f, i, idx ,a);
                            _timederiv_average +=  pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6));
                            _grad_average += gradient_energy_eachpoint(f, i, idx)/ pow(a,4);
                             _df_average += df[i][idx];
                            _f_average += f[i][idx];
                         //   std::cout << "df[i][idx]*a = " << df[i][idx]*a << std::endl;
                           // std::cout << "f[i][idx] = " << f[i][idx] << std::endl;
                          
                          //  std::cout << "_df_average1 = " << _df_average << std::endl;
                          //  std::cout << "_f_average1 = " << _f_average << std::endl;
                        }
                        break;
                    case 3:
                        for( int k = 0; k < N; ++k ){
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;
                                value[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6)) + field->aV(f, i, idx ,a)+gradient_energy_eachpoint(f, i, idx)/pow(a,4);
                                _average[i]+= value[i][idx];
                                _potential_average += field->aV(f, i, idx ,a);
                                _timederiv_average +=  pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,6));
                                _grad_average += gradient_energy_eachpoint(f, i, idx)/ pow(a,4);
                                _df_average += df[i][idx];
                                _f_average += f[i][idx];
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
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                            value[i][idx] = pow(df[i][idx], 2)/2 + field->V(f, i ,idx)+gradient_energy_eachpoint(f, i, idx);
                            _average[i] += value[i][idx];
                            _potential_average += field->V(f, i ,idx);
                            _timederiv_average +=  pow(df[i][idx], 2)/2;
                            _grad_average += gradient_energy_eachpoint(f, i, idx);
                            _df_average += df[i][idx];
                          //   std::cout << "field->V(f, 0 ,"<< idx << ") =" << field->V(f, i ,idx) << std::endl;
                          //   std::cout << "_potential_average =" << _potential_average << std::endl;
                           // std::cout << "value[i][idx] =" << value[i][idx] << std::endl;
                           // std::cout << "_average[i] =" << _average[i] << std::endl;
                          
                        }
                        break;
                    case 3:
                        for( int k = 0; k < N; ++k ){
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;
                                value[i][idx] = pow(df[i][idx], 2)/2 + field->V(f, i ,idx)+gradient_energy_eachpoint(f, i, idx);
                                _average[i]+= value[i][idx];
                                _potential_average += field->V(f, i ,idx);
                                _timederiv_average +=  pow(df[i][idx], 2)/2;
                                _grad_average += gradient_energy_eachpoint(f, i, idx);
                                _df_average += df[i][idx];
                                _f_average += f[i][idx];
                            }
                        }
                        break;
                }
            }
        }//65536
    }
// std::cout << "a = " << a << std::endl;
   // std::cout << "da = " << da << std::endl;
   // std::cout << "_df_average1 = " <<_df_average << std::endl;
   //std::cout << "_f_average1 = " << _f_average << std::endl;
  //  std::cout << "field->V(f, i ,idx) =" << field->V(f, 0 ,20) << std::endl;
  //  std::cout << "_potential_average =" << _potential_average << std::endl;
   // std::cout << "value[i][idx] =" << value[0][20] << std::endl;
   // std::cout << "_average[i] =" << _average[0] << std::endl;
   for( int i = 0; i < num_fields; ++i ){
       // if( expansion ) _average[i] += gradient_energy(f[i]) / pow(a,4);
      //  else            _average[i] += gradient_energy(f[i]);
        for( int j = 0; j < dim; ++j ) _potential_average /= N;
        for( int j = 0; j < dim; ++j ) _timederiv_average /= N;
        for( int j = 0; j < dim; ++j ) _grad_average /= N;
        for( int j = 0; j < dim; ++j ) _average[i] /= N;
        for( int j = 0; j < dim; ++j ) _df_average /= N;
        for( int j = 0; j < dim; ++j ) _f_average /= N;
     //  std::cout << "_df_average = " <<_df_average << std::endl;
      // std::cout << "_f_average = " << _f_average << std::endl;
   /*    std::cout << "a = " << a << std::endl;
       std::cout << "da = " << da << std::endl;
       std::cout << "_df_average = " <<_df_average << std::endl;
        std::cout << "_f_average = " << _f_average << std::endl;
       _df_average = m*(pow(a,-2)*_df_average-da*pow(a,-3)*_f_average);
       std::cout << "_f_average = " << _f_average << std::endl;*/
        _total_average += _average[i];
   }
   
    //devide by average
    for( int i = 0; i < num_fields; ++i ){
     //  #pragma omp parallel for schedule( static ) num_threads ( num_threads )
        for( int j = 0; j < N; ++j ){
    switch( dim ){
        case 1:
            int idx = j;
         value[i][idx] /= _average[i] ;
            break;
        case 2:
            for( int k = 0; k < N; ++k ){
                int idx = j*N + k;
           //   std::cout << "value[i][" << idx << "]1 =" << value[i][idx] << std::endl;
          //  std::cout << "_average[i] =" << _average[0] << std::endl;
                value[i][idx] /= _average[i] ;
              //  std::cout << "value[i][" << idx << "]2 =" << value[i][idx] << std::endl;
                if( value[i][idx] > 10.){
                     //std::cout << "value[i][" << idx << "] =" << value[i][idx] << std::endl;
                }
            }
            break;
        case 3:
            for( int k = 0; k < N; ++k ){
                for( int l = 0; l < N; ++l ){
                    int idx = (j*N + k)*N + l;
                    value[i][idx] /= _average[i] ;
                }
            }
            break;
                   }
        }
       /* std::cout << "value[i][20] =" << value[i][20] << std::endl;
        std::cout << "value[i][30] =" << value[i][30] << std::endl;
         std::cout << "value[i][100] =" << value[i][100] << std::endl;
       std::cout << "value[i][1200] =" << value[i][1200] << std::endl;
        std::cout << "value[i][12000] =" << value[i][12000] << std::endl;*/
    }
}


double Energy::gradient_energy( double* f )
{
    double gradient_energy = 0;
    
 //  #pragma omp parallel for reduction(+:gradient_energy) schedule( static ) num_threads ( num_threads )
    for( int j = 0; j < N; ++j ){
        int jp1 = (j == N-1)?     0: j+1;
        int jp2 = (j >= N-2)? j-N+2: j+2;
        #ifdef SPHERICAL_SYM
        int jm1 = abs( j-1 );
        int jm2 = abs( j-2 );
        #else
        int jm1 = (j == 0)?   N-1: j-1;
        int jm2 = (j <  2)? j+N-2: j-2;
        #endif
        #if dim == 1
		    #ifdef SPHERICAL_SYM
            int idx = j;
            if( idx == 0 ){
                grad[0] = 0;
            }else if( idx == N-2 ){
                grad[0] = ( f[jp1] - f[jm1] ) / (2*dx);
            }else if( idx == N-1 ){
                grad[0] = ( f[jm2] - 4*f[jm2] + 3*f[idx]  ) / (2*dx);
            }else{
                grad[0] = ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx);
            }
            #else
            gradient_energy +=  pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2 );
		    #endif
	    #elif dim == 2
            for( int k = 0; k < N; ++k ){
                int kp1 = (k == N-1)?     0: k+1;
                int kp2 = (k >= N-2)? k-N+2: k+2;
                int km1 = (k ==   0)?   N-1: k-1;
                int km2 = (k <    2)? k+N-2: k-2;
                gradient_energy += pow( (- f[jp2*N+k] + 8*f[jp1*N+k] - 8*f[jm1*N+k] + f[jm2*N+k]) / (12*dx), 2 );
                gradient_energy += pow( (- f[j*N+kp2] + 8*f[j*N+kp1] - 8*f[j*N+km1] + f[j*N+km2]) / (12*dx), 2 );
            }
        #elif dim == 3
            for( int k = 0; k < N; ++k ){
                int kp1 = (k == N-1)?     0: k+1;
                int kp2 = (k >= N-2)? k-N+2: k+2;
                int km1 = (k ==   0)?   N-1: k-1;
                int km2 = (k <    2)? k+N-2: k-2;
                for( int l = 0; l < N; ++l ){
                    int lp1 = (l == N-1)?     0: l+1;
                    int lp2 = (l >= N-2)? l-N+2: l+2;
                    int lm1 = (l ==   0)?   N-1: l-1;
                    int lm2 = (l <    2)? l+N-2: l-2;                        
                    gradient_energy += pow( (- f[(jp2*N+k)*N+l] + 8*f[(jp1*N+k)*N+l] - 8*f[(jm1*N+k)*N+l] + f[(jm2*N+k)*N+l]) / (12*dx), 2 );
                    gradient_energy += pow( (- f[(j*N+kp2)*N+l] + 8*f[(j*N+kp1)*N+l] - 8*f[(j*N+km1)*N+l] + f[(j*N+km2)*N+l]) / (12*dx), 2 );
                    gradient_energy += pow( (- f[(j*N+k)*N+lp2] + 8*f[(j*N+k)*N+lp1] - 8*f[(j*N+k)*N+lm1] + f[(j*N+k)*N+lm2]) / (12*dx), 2 );
                }
            }
        #endif
    }

	return gradient_energy / 2;
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
  
    return (pow( (f[i][jp1*N+k] - f[i][jm1*N+k]) / (2*dx), 2 )
            + pow( (f[i][j*N+kp1] - f[i][j*N+km1]) / (2*dx), 2 ))/2;
   // return (pow( (- f[i][jp2*N+k] + 8*f[i][jp1*N+k] - 8*f[i][jm1*N+k] + f[i][jm2*N+k]) / (12*dx), 2 )
         //   + pow( (- f[i][j*N+kp2] + 8*f[i][j*N+kp1] - 8*f[i][j*N+km1] + f[i][j*N+km2]) / (12*dx), 2 ))/2;

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
