#include <cmath>
#include <random>
#include "lattice_field.hpp"
#include "utilities.hpp"
#include "equations.hpp"


#pragma omp declare simd
double Field::laplacian( double* f, int j, int k, int l )
{	
    int jp1 = (j == N-1)?     0: j+1;
	int jp2 = (j >= N-2)? j-N+2: j+2;
	#ifdef SPHERICAL_SYM
	if( dim != 1 ){
		Logout( "Error: Dimension must be 1 when you define SPHERICAL_SYM with ABC. \n" );
		exit(1);
	}
	int jm1 = abs( j-1 );
	int jm2 = abs( j-2 );
	#else
	int jm1 = (j == 0)?   N-1: j-1;
	int jm2 = (j <  2)? j+N-2: j-2;
	#endif
	
	#if dim == 1
	    int idx = j;
		#ifdef SPHERICAL_SYM
            if( idx == 0 ){  //1/r*df/dr = 0 when r = 0
                return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*dx*dx);
            }else if( idx == N-2 ){ // Calculated  by 2nd order  central  difference  scheme  at onestep  before  the  boundary
                return (f[jp1] - 2*f[idx] + f[jm1]) / (dx*dx) + 2*gradient(f, 0, i, idx, 0, 0) / (idx*dx);
            }else if( idx == N-1 ){ //Calculated  by 2nd order  backward  difference  scheme  at the  boundary
                return (df[jm2] - 4*df[jm1] + 3*df[idx]) / (2*dx) + df[idx] / (idx*dx);
            }else{
                return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*dx*dx) + 2*gradient(f, 0, i, idx, 0, 0) / (idx*dx);
            }
		#else
		    return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*dx*dx);
		#endif
	#elif dim == 2
		int kp1 = (k == N-1)?     0: k+1;
		int kp2 = (k >= N-2)? k-N+2: k+2;
		int km1 = (k ==   0)?   N-1: k-1;
		int km2 = (k <    2)? k+N-2: k-2;

        int idx = j*N + k;
         return ( (- f[jp2*N+k] + 16*f[jp1*N+k] - 30*f[idx] + 16*f[jm1*N+k] - f[jm2*N+k])
               + (- f[j*N+kp2] + 16*f[j*N+kp1] - 30*f[idx] + 16*f[j*N+km1] - f[j*N+km2]) ) / (12*dx*dx);
//     return ( (f[jp1*N+k] - 2*f[idx] + f[jm1*N+k])
//      + (f[j*N+kp1] - 2*f[idx] + f[j*N+km1]) ) / (dx*dx);
	#elif dim == 3
		int kp1 = (k == N-1)?     0: k+1;
		int kp2 = (k >= N-2)? k-N+2: k+2;
		int km1 = (k ==   0)?   N-1: k-1;
		int km2 = (k <    2)? k+N-2: k-2;
		
		int lp1 = (l == N-1)?     0: l+1;
		int lp2 = (l >= N-2)? l-N+2: l+2;
		int lm1 = (l ==   0)?   N-1: l-1;
		int lm2 = (l <    2)? l+N-2: l-2;

        int idx = (j*N + k)*N + l;
    
        return ( (- f[(jp2*N+k)*N+l] + 16*f[(jp1*N+k)*N+l] - 30*f[idx] + 16*f[(jm1*N+k)*N+l] - f[(jm2*N+k)*N+l])
               + (- f[(j*N+kp2)*N+l] + 16*f[(j*N+kp1)*N+l] - 30*f[idx] + 16*f[(j*N+km1)*N+l] - f[(j*N+km2)*N+l])
               + (- f[(j*N+k)*N+lp2] + 16*f[(j*N+k)*N+lp1] - 30*f[idx] + 16*f[(j*N+k)*N+lm1] - f[(j*N+k)*N+lm2]) ) / (12*dx*dx);
	#endif
}

#pragma omp declare simd
double Field::gradient_energy_eachpoint( double** f ,int i, int idx )
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


double Field::gradient_energy( double* f )
{
    double gradient_energy = 0;
    
  // #pragma omp parallel for reduction (+:gradient_energy) schedule( static ) num_threads ( num_threads )
    for( int j = 0; j < N; ++j ){
        int jp1 = (j == N-1)?     0: j+1;
        int jp2 = (j >= N-2)? j-N+2: j+2;
        #ifdef SPHERICAL_SYM
        int jm1 = abs( j-1 );
        int jm2 = abs( j-2 );
        #else
        int jm1 = (j ==   0)?   N-1: j-1;
        int jm2 = (j <    2)? j+N-2: j-2;
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
//            std::cout << "pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2 ) = " << pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2 ) << "f["<< jp2 << "] = " << f[jp2] << "f["<< jp1 << "] = " << f[jp1] << "f["<< jm1 << "] = " << f[jm1]<< "f["<< jm2 << "] = " << f[jm2] << std::endl;
		    #endif
	    #elif dim == 2
            for( int k = 0; k < N; ++k ){
                int kp1 = (k == N-1)?     0: k+1;
                int kp2 = (k >= N-2)? k-N+2: k+2;
                int km1 = (k ==   0)?   N-1: k-1;
                int km2 = (k <    2)? k+N-2: k-2;
//                gradient_energy +=  pow( (f[jp1*N+k] - f[jm1*N+k]) / (2*dx), 2 );
//                 gradient_energy += pow( ( f[j*N+kp1]  - f[j*N+km1] ) /(2*dx), 2 );
              gradient_energy += pow( ( - f[jp2*N+k] + 8*f[jp1*N+k] - 8*f[jm1*N+k] + f[jm2*N+k] ) / (12*dx), 2 );
                gradient_energy += pow( ( - f[j*N+kp2] + 8*f[j*N+kp1] - 8*f[j*N+km1] + f[j*N+km2] ) / (12*dx), 2 );
            }
        #elif dim == 3
            for( int k = 0; k < N; ++k ){
                int kp1 = (k == N-1)?     0: k+1;
                int kp2 = (k >= N-2)? k-N+2: k+2;
                int km1 = (k ==   0)?   N-1: k-1;
                int km2 = (k <    2)? k+N-2: k-2;
                for( int l = 0; l < N; ++l ){
                    int lp1 = (l == N-1) ?     0: l+1;
                    int lp2 = (l >= N-2) ? l-N+2: l+2;
                    int lm1 = (l ==   0) ?   N-1: l-1;
                    int lm2 = (l <    2) ? l+N-2: l-2;                        
                    gradient_energy += pow( ( - f[(jp2*N+k)*N+l] + 8*f[(jp1*N+k)*N+l] - 8*f[(jm1*N+k)*N+l] + f[(jm2*N+k)*N+l] ) / (12*dx), 2 );
                    gradient_energy += pow( ( - f[(j*N+kp2)*N+l] + 8*f[(j*N+kp1)*N+l] - 8*f[(j*N+km1)*N+l] + f[(j*N+km2)*N+l] ) / (12*dx), 2 );
                    gradient_energy += pow( ( - f[(j*N+k)*N+lp2] + 8*f[(j*N+k)*N+lp1] - 8*f[(j*N+k)*N+lm1] + f[(j*N+k)*N+lm2] ) / (12*dx), 2 );
                }
            }
        #endif
    }

    for( int i = 0; i < dim; ++i ) gradient_energy /= N;
    
	return gradient_energy/2;
    
}

double Field::potential_energy( double** f, double a )
{
	double potential_energy = 0;
    
//#if   dim == 1
//#pragma omp parallel for simd reduction(+:potential_energy) schedule(static) num_threads(num_threads)
//#elif dim >= 2
//#pragma omp parallel for reduction(+:potential_energy) schedule(static) num_threads(num_threads)
//#endif
	for( int j = 0; j < N; ++j ){
        #if dim == 1
            int idx = j;
             potential_energy += V_lattice( f, idx, a );
       // std::cout << "V_lattice( f, idx, a ) = " << V_lattice( f, idx, a ) << std::endl;
	    #elif dim == 2
//        #pragma omp simd reduction(+:potential_energy)
            for( int k = 0; k < N; ++k )
			{
				int idx = j*N + k;
                potential_energy += V_lattice( f, idx, a );
            }
        #elif dim == 3
            for( int k = 0; k < N; ++k ){
//        #pragma omp simd reduction(+:potential_energy)
                for( int l = 0; l < N; ++l ){
					int idx = ( j*N + k)*N + l;
					potential_energy += V_lattice( f, idx, a );
				}
            }
        #endif
	}
    for( int i = 0; i < dim; ++i ) potential_energy /= N;
    
	return potential_energy;
}


double Field::average( double* f, int i )
{
	//Variable needs to be declared for OpenMP reduction directive
    double average = 0;

#if   dim == 1
#pragma omp parallel for simd reduction(+:average) schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for reduction(+:average) schedule(static) num_threads(num_threads)
#endif
	for( int j = 0; j < N; ++j ){
		switch( dim )
		{
			case 1:
				int idx = j;
				average += f[idx];
				break;
			case 2:
            #pragma omp simd reduction(+:average)
				for( int k = 0; k < N; ++k ){
					int idx = j*N + k;
					average += f[idx];
				}
				break;
			case 3:
				for( int k = 0; k < N; ++k ){
            #pragma omp simd reduction(+:average)
					for( int l = 0; l < N; ++l ){
						int idx = (j*N + k)*N + l;
						average += f[idx];
					}
				}
				break;
		}
	}
    for( int j = 0; j < dim; ++j ) average /= N;
    
    //substitute the obtained variable to member variable
    _average[i] = average;
	
    return _average[i];
}


double Field::variance( double* f, int i )
{
    //Variable needs to be declared for OpenMP reduction directive
    double variance = 0;
    
#if   dim == 1
#pragma omp parallel for simd reduction(+:variance) schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for reduction(+:variance) schedule(static) num_threads(num_threads)
#endif

	for( int j = 0; j < N; ++j ){
		switch( dim ){
			case 1:
				int idx = j;
				variance += pow( f[idx] - _average[i], 2 );
				break;
			case 2:
                #pragma omp simd reduction(+:variance)
				for( int k = 0; k < N; ++k ){
					int idx = j*N + k;
                    variance += pow( f[idx] - _average[i], 2 );
                
				}
				break;
			case 3:
				for( int k = 0; k < N; ++k ){
                    #pragma omp simd reduction(+:variance)
					for( int l = 0; l < N; ++l ){
						int idx = (j*N + k)*N + l;
						variance += pow( f[idx] - _average[i], 2 );
					}
				}
				break;
		}
	}
    for( int j = 0; j < dim; ++j ) variance /= N;
    
    //substitute the obtained variable to member variable
    _variance[i] = variance;
   
    return sqrt(_variance[i]);
}

//#pragma omp declare simd
double Field::V_lattice   ( double** f, int idx, double a )  {
    // std::cout << "a = " << a << std::endl;
    for(fld=0;fld< num_fields - 1;fld++)
    {
//        if(fld==1){
//            f_MPl[fld] = FIXPSI - f[fld][idx]/(rescale_A*a);
//        }else{
//            f_MPl[fld] = f[fld][idx]/(rescale_A*a);
//        }
        
        f_MPl[fld] = f[fld][idx]/(rescale_A*a);
    }
//    std::cout << "pow(a,4) = " << pow(a,4) << std::endl;
//    std::cout << "pow(rescale_A/rescale_B,2) = " << pow(rescale_A/rescale_B,2) << std::endl;
//    std::cout << "V(f_MPl[0],f_MPl[1],f_MPl[2]) = " << V(f_MPl[0],f_MPl[1],f_MPl[2]) << std::endl;
    
    return pow(a,4)*pow(rescale_A/rescale_B,2)*V(f_MPl[0],f_MPl[1],f_MPl[2]); }

//#pragma omp declare simd
double Field::dV_lattice ( double** f, int i, int idx, double a )  {
    
    
    for(fld=0;fld< num_fields - 1;fld++)
    {
//        if(fld==1){
//        f_MPl[fld] = FIXPSI - f[fld][idx]/(rescale_A*a);
//        }else{
//        f_MPl[fld] = f[fld][idx]/(rescale_A*a);
//        }
        
        f_MPl[fld] = f[fld][idx]/(rescale_A*a);
    }
    

    switch (i){
        case 0:
          
            
            return pow(a,3)*rescale_A*V_1(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //sigma
            
           
//        case 1: return -pow(a,3)*rescale_A*V_2(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //psi
        case 1: return pow(a,3)*rescale_A*V_2(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //psi
         
        case 2:
            
            
//            std::cout << "pow(a,3)*rescale_A*V_3(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B) = " <<  pow(a,3)*rescale_A*V_3(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B) << std::endl;
            
            return pow(a,3)*rescale_A*V_3(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //phi
           
        default:  Logout( "Parameter 'i' in dV_lattice must be 0 ~ 2. \n" );
                exit(1);
    }
    
}

double Field::ddV_lattice ( double** f, int i, int idx, double a )  {
    
    for( fld=0; fld < num_fields - 1; fld++)
    {
//        if(fld==1){
//            f_MPl[fld] = FIXPSI - f[fld][idx]/(rescale_A*a);
//        }else{
//            f_MPl[fld] = f[fld][idx]/(rescale_A*a);
//        }
        
         f_MPl[fld] = f[fld][idx]/(rescale_A*a);
    }
    
    
    switch (i){
        case 0: return pow(a/rescale_B,2)*V_11(f_MPl[0],f_MPl[1],f_MPl[2]); //sigma
        case 1: return pow(a/rescale_B,2)*V_22(f_MPl[0],f_MPl[1],f_MPl[2]); //psi
        case 2: return pow(a/rescale_B,2)*V_33(f_MPl[0],f_MPl[1],f_MPl[2]); //phi
        default:  Logout( "Parameter 'i' in ddV_lattice must be 0 ~ 2. \n" );
            exit(1);
    }
    
}


void Field::effective_mass(double mass_sq[], double *field_values){

    mass_sq[0] = pow(exp(OSCSTART)/rescale_B,2)*V_11(field_values[0],field_values[1],field_values[2]); //sigma
    mass_sq[1] = pow(exp(OSCSTART)/rescale_B,2)*V_22(field_values[0],field_values[1],field_values[2]); //psi
    mass_sq[2] = pow(exp(OSCSTART)/rescale_B,2)*V_33(field_values[0],field_values[1],field_values[2]); //phi

}

double Field::power_spectrum( double** f, int i, int m)
{
     
    
    
    if (k_lattice_grid_min_MPl < kfrom_MPl_lattice)
    {
        m_start = outrange_num + 1;
    }else
    {
        m_start = 1;
    }
    
    #if  dim==1
    
    
        if(m == m_start){
            
            //subtract zero mode
            for( int j = 0; j < N; ++j ){
                int idx_fluc = j;
               f_fluc[i][idx_fluc] = f[i][idx_fluc] - average(f[i], i);
                
//                std::cout << "f_fluc[" << i << "][" << idx_fluc << "]" << f_fluc[i][idx_fluc] << std::endl;
//                std::cout << "f[" << i << "][" << idx_fluc << "]" << f[i][idx_fluc] << std::endl;
//                std::cout << "_average[" << i << "]" << _average[i] << std::endl;
                
            }
            
            //transform from real space to phase space (Needs to be done only once)
            DFT_r2cD1( f_fluc[i], f_fluc_k[i] );
            
//            for( int m = 0; m < N; ++m ){
//                std::cout << "1: f_fluc_k[" << i << "][" << m << "]" << f_fluc_k[i][m] << std::endl;
//            }
        }
    
         //f[][0] corresponds to Re(f_{i=0})
         //f[][2],f[][3] corresponds to the Re and Im of f_{i=1}
         // ...
         //f[][1] corresponds to Re(f_{i=N/2})
    
//    for( int m = 0; m < N; ++m ){
//        std::cout << "2: f_fluc_k[" << i << "][" << m << "]" << f_fluc_k[i][m] << std::endl;
//    }
        int j = m;
             if(m==0)
             {// zero-frequency (DC)
                 PS[i][m] = m*pow(f_fluc_k[i][0]/N,2);
                 
             }else if(m == N/2){//Nyquist frequency
                 PS[i][m] = m*pow(f_fluc_k[i][1]/N,2);
             }else{
                 
                PS[i][m] = m*( pow(f_fluc_k[i][2*j]/N,2) + pow(f_fluc_k[i][2*j+1]/N,2));
                
             }
    
        PS[i][m] *= 2; //  the conjugate part needs to be taken into account as well, so double the value
     //std::cout << " PS[" << i << "][" << idx << "] = " << PS[i][idx] << std::endl;
        return PS[i][m];
 
    #elif  dim==2
    
    if(m == m_start)
    {
        //Initialize PS to zero for each time step
       std::fill(PS[i].begin(), PS[i].end(), 0);
        
        for( int j = 0; j < N; ++j )
        {
            for( int k = 0; k < N; ++k )
            {
                int idx_fluc = j*N + k;
                f_fluc[i][idx_fluc] = f[i][idx_fluc] - average(f[i], i);
                
                       //         std::cout << "f_fluc[" << i << "][" << idx_fluc << "]" << f_fluc[i][idx_fluc] << std::endl;
                //                std::cout << "f[" << i << "][" << idx_fluc << "]" << f[i][idx_fluc] << std::endl;
             
            }
        }
        
//                        std::cout << "bf f_fluc[" << i << "][20]" << f_fluc[i][20] << std::endl;
//                        std::cout << "bf f[" << i << "][20]" << f[i][20] << std::endl;
//    std::cout << "bf _average[" << i << "]" << _average[i] << std::endl;

                //transform from real space to phase space (Needs to be done only once)
                DFT_r2cD2( f_fluc[i], f_fluc_k[i], f_fluc_k_nyquist[i] );
        
//        std::cout << "af f_fluc[" << i << "][20]" << f_fluc[i][20] << std::endl;
//         std::cout << "af f_fluc_k[" << i << "][20]" << f_fluc_k[i][20] << std::endl;
//        std::cout << "af f_fluc_k_nyquist[" << i << "][20]" << f_fluc_k_nyquist[i][20] << std::endl;
        
//        for( int aaa = 0; aaa < N; ++aaa ){
//                            std::cout << "1: f_fluc_k[" << i << "][" << aaa << "]" << f_fluc_k[i][aaa] << std::endl;
//                        }
    }
    
    
//        for( int bbb = 0; bbb < N; ++bbb ){
//            std::cout << "2: f_fluc_k[" << i << "][" << bbb << "]" << f_fluc_k[i][bbb] << std::endl;
//        }
    
        
      for(int j = 0; j < N; ++j )
      {
          J = (j <= N/2 ? j : j-N);
          
          for(int k = 0; k <= N/2; ++k )
          {
              if ( m-1 < sqrt(pow(J,2)+pow(k,2)) && sqrt(pow(J,2)+pow(k,2)) <= m  )
              {
                  
              if(k==N/2)//Nyquist
              {
                 
                      if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=k, 0<=j <=N/2) doesn't have an equal number of counterparts in the range (0<=k, N/2+1<j<=N-1), we consider the range (1<=k, 1<=j<=N/2) and multiply by 4 and add the Nyquist  (k=N/2, j=0) x 2 and (k=0, j=N/2) x 2.
                      {
                          // This corresponds to (k=N/2, j=0). We multiply the total sum by 4 at the end so we multiply by 0.5 here
                          PS[i][m] += 0.5*m*(pow(f_fluc_k_nyquist[i][2*j],2) + pow(f_fluc_k_nyquist[i][2*j+1],2))/(pow(N,4));
                      }
//                              std::cout << "1: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                      else {//No modes correspond to this branch
                          PS[i][m] += m*(pow(f_fluc_k_nyquist[i][2*j],2) + pow(f_fluc_k_nyquist[i][2*j+1],2))/(pow(N,4));
                      }
                      
                  
                      
                  
                  
              }else//Modes besides Nyquist
              {
                  
//                  std::cout << "m = " << m << std::endl;
//                  std::cout << "j = " << j << std::endl;
//                  std::cout << "J = " << J << std::endl;
//                  std::cout << "k = " << k << std::endl;
//                  std::cout << "sqrt(pow(J,2)+pow(k,2)) = " << sqrt(pow(J,2)+pow(k,2)) << std::endl;
                  
                 
                      if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=k, 0<=j <=N/2) doesn't have an equal number of counterparts in the range (0<=k, N/2+1<j<=N-1), we consider the range (1<=k, 1<=j<=N/2) and multiply by 4 and add the Nyquist  (k=N/2, j=0) x 2 and (k=0, j=N/2) x 2.
                      {
                          if(j==N/2)
                          {// This corresponds to (k=0, j=N/2). We multiply the total sum by 4 at the end so we multiply by 0.5 here
                              PS[i][m] +=
                              0.5*m*(pow(f_fluc_k[i][j*N+2*k],2) + pow(f_fluc_k[i][j*N+2*k+1],2))/(pow(N,4));
                          }else{//No modes correspond to this branch
                              PS[i][m] +=
                              m*(pow(f_fluc_k[i][j*N+2*k],2) + pow(f_fluc_k[i][j*N+2*k+1],2))/(pow(N,4));
                          }
                             
                              
//                               std::cout << "3: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                          
                      }
                      else
                      {
                    PS[i][m] +=  m*(pow(f_fluc_k[i][j*N+2*k],2) + pow(f_fluc_k[i][j*N+2*k+1],2))/(pow(N,4));
                          
//                          std::cout << "4: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                      }
                      
              }//if ( m-1 < sqrt(pow(J,2)+pow(k,2)) && sqrt(pow(J,2)+pow(k,2)) <= m  )
                  
              }
          }//for(int k = 0; k <= N/2; ++k )
          
        
          
      }//for(int j = 0; j < N; ++j )
    
    if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=k, 0<=j <=N/2) doesn't have an equal number of counterparts in the range (0<=k, N/2+1<j<=N-1), we consider the range (1<=k, 1<=j<=N/2) and multiply by 4 and add the Nyquist  (k=N/2, j=0) x 2 and (k=0, j=N/2) x 2.
    {
        PS[i][m] *= 4;
        
//         std::cout << "5: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
    }
    else
    {
         PS[i][m] *= 2; //  the conjugate part needs to be taken into account as well so double the value
        
//         std::cout << "6: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
    }
    
    //std::cout << "7: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
         return PS[i][m];
    
    
    #elif  dim==3
    
    if(m == m_start)
    {
        //Initialize PS to zero for each time step
        std::fill(PS[i].begin(), PS[i].end(), 0);
        
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N; ++k ){
                for( int l = 0; l < N; ++l ){
                    int idx_fluc = (j*N + k)*N + l;
                    f_fluc[i][idx_fluc] = f[i][idx_fluc] - average(f[i], i);

//                             std::cout << "f_fluc[" << i << "][" << idx_fluc << "]" << f_fluc[i][idx_fluc] << std::endl;
//                                    std::cout << "f[" << i << "][" << idx_fluc << "]" << f[i][idx_fluc] << std::endl;
                }
            }
        }
        //std::cout << std::endl;//Not sure why but without this line, we get "Segmentation fault: 11"
        DFT_r2cD3( f_fluc[i], f_fluc_k[i], f_fluc_k_nyquist[i] );
       
    }
    
    
    
   
    for(int j = 0; j < N; ++j )
    {
        J = (j <= N/2 ? j : j-N);
        
        for(int k = 0; k < N; k++)
        {
            K = (k <= N/2 ? k : k-N);
        
        for(int l = 0; l <= N/2; ++l )
        {
            
            if ( m-1 < sqrt(pow(J,2)+pow(K,2)+pow(l,2)) && sqrt(pow(J,2)+pow(K,2)+pow(l,2)) <= m  )
            {
            
            if(l==N/2)//Nyquist
            {
                
                    if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=l, 0<=k<=N/2, 0<=j<=N/2) doesn't have an equal number of counterparts in the range (0<=l && !(0<=k <=N/2, 0<=j <=N/2) ), we consider the range (1<=l, 1<=k<=N/2, 1<=j <=N/2 ) and multiply by 8 and add the Nyquist (l=N/2, k=0,j=0) x 2, (l=0, k=0, j=N/2) x 2 and (l=0, k=N/2, j=0) x 2.
                
                    {
                        // This corresponds to (l=N/2, k=0,j=0). We multiply the total sum by 8 at the end so we multiply by 0.25 here
                        PS[i][m] += 0.25*m*(pow(f_fluc_k_nyquist[i][j][2*k],2) + pow(f_fluc_k_nyquist[i][j][2*k+1],2))/(pow(N,6));
                        
//                                                      std::cout << "1: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                    }
                    else
                    {//No modes correspond to this branch
                        PS[i][m] += m*(pow(f_fluc_k_nyquist[i][j][2*k],2) + pow(f_fluc_k_nyquist[i][j][2*k+1],2))/(pow(N,6));
                        
//                                                   std::cout << "2: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                    }
                    
                
                
            }else//Modes besides Nyquist
            {
                
                //                  std::cout << "m = " << m << std::endl;
                //                  std::cout << "j = " << j << std::endl;
                //                  std::cout << "J = " << J << std::endl;
                //                  std::cout << "k = " << k << std::endl;
                //                  std::cout << "sqrt(pow(J,2)+pow(k,2)) = " << sqrt(pow(J,2)+pow(k,2)) << std::endl;
                
               
                    if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=l, 0<=k<=N/2, 0<=j<=N/2) doesn't have an equal number of counterparts in the range (0<=l && !(0<=k <=N/2, 0<=j <=N/2) ), we consider the range (1<=l, 1<=k<=N/2, 1<=j <=N/2 ) and multiply by 8 and add the Nyquist (l=N/2, k=0,j=0) x 2, (l=0, k=0, j=N/2) x 2 and (l=0, k=N/2, j=0) x 2.
                    {
                        if((j==N/2 && k==0) || (j==0 && k==N/2) )
                        {// This corresponds to (l=0, k=0, j=N/2) and (l=0, k=N/2, j=0). We multiply the total sum by 8 at the end so we multiply by 0.25 here
                            PS[i][m] +=
                            0.25*m*(pow(f_fluc_k[i][(j*N + k)*N + 2*l],2) + pow(f_fluc_k[i][(j*N + k)*N + 2*l+1],2))/(pow(N,6));
                        }else
                        {
                            PS[i][m] += m*(pow(f_fluc_k[i][(j*N + k)*N + 2*l],2) + pow(f_fluc_k[i][(j*N + k)*N + 2*l+1],2))/(pow(N,6));
                            //                           std::cout << "2: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                        }
                        
                        
//                                                       std::cout << "3: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                        
                    }
                    else
                    {
                        PS[i][m] +=
                        m*(pow(f_fluc_k[i][(j*N + k)*N + 2*l],2) + pow(f_fluc_k[i][(j*N + k)*N + 2*l+1],2))/(pow(N,6));
                        
//                                                  std::cout << "4: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                    }
                    
                
            }
                
        } //if ( m-1 < sqrt(pow(J,2)+pow(K,2)+pow(l,2)) && sqrt(pow(J,2)+pow(K,2)+pow(l,2)) <= m  )
            
        }//for(int l = 0; l <= N/2; ++l )
            
        }//for(int k = 0; k <= N/2; ++k )
        
        
        
    }//for(int j = 0; j < N; ++j )
    
    if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=l, 0<=k<=N/2, 0<=j <=N/2) doesn't have an equal number of counterparts in the range (0<=l && !(0<=k <=N/2, 0<=j <=N/2) ), we only consider the range (0<=l, 0<=k<=N/2, 0<=j <=N/2 ) and multiply by 8 at the end.
    {
        PS[i][m] *= 8;
        
//                 std::cout << "5: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
    }
    else
    {
        PS[i][m] *= 2; //  the conjugate part needs to be taken into account as well so double the value
        
//                 std::cout << "6: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
    }
    
//    std::cout << "7: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
    return PS[i][m];
    
    #endif

    

}
