//Doxygen
/**
* @file    lattice_field.cpp
* @brief    Lattice field source file
* @author   Francis Otani
* @date
* @details
*/



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

	int jm1 = (j == 0)?   N-1: j-1;
	int jm2 = (j <  2)? j+N-2: j-2;
	
    switch (dim){
    case 1:
    {
    int idx = j;
    return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*dx*dx);
    }
    case 2:
    {
    int kp1 = (k == N-1)?     0: k+1;
		int kp2 = (k >= N-2)? k-N+2: k+2;
		int km1 = (k ==   0)?   N-1: k-1;
		int km2 = (k <    2)? k+N-2: k-2;

        int idx = j*N + k;
         return ( (- f[jp2*N+k] + 16*f[jp1*N+k] - 30*f[idx] + 16*f[jm1*N+k] - f[jm2*N+k])
               + (- f[j*N+kp2] + 16*f[j*N+kp1] - 30*f[idx] + 16*f[j*N+km1] - f[j*N+km2]) ) / (12*dx*dx);
    }
//     return ( (f[jp1*N+k] - 2*f[idx] + f[jm1*N+k])
//      + (f[j*N+kp1] - 2*f[idx] + f[j*N+km1]) ) / (dx*dx);
    case 3:
    {
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
    }
default:
    std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
    exit(1);
    
    }
    
}

#pragma omp declare simd
double Field::gradient_energy_eachpoint( double** f ,int i, int idx )
{
    
    switch (dim){
    case 1:
    {
    int j= idx;
    int jp1 = (j == N-1)?     0: j+1;
    int jp2 = (j >= N-2)? j-N+2: j+2;
    
    int jm1 = (j == 0)?   N-1: j-1;
    int jm2 = (j <  2)? j+N-2: j-2;
    
    return  pow( ( - f[i][jp2]  + 8*f[i][jp1]  - 8*f[i][jm1] + f[i][jm2] ) / (12*dx), 2.0 )/2;
    }
    case 2:
    {
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
    return (pow( (- f[i][jp2*N+k] + 8*f[i][jp1*N+k] - 8*f[i][jm1*N+k] + f[i][jm2*N+k]) / (12*dx), 2.0 )
            + pow( (- f[i][j*N+kp2] + 8*f[i][j*N+kp1] - 8*f[i][j*N+km1] + f[i][j*N+km2]) / (12*dx), 2.0 ))/2;
    }
case 3:
    {
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
    
    return  (pow( (- f[i][(jp2*N+k)*N+l] + 8*f[i][(jp1*N+k)*N+l] - 8*f[i][(jm1*N+k)*N+l] + f[i][(jm2*N+k)*N+l]) / (12*dx), 2.0 )
             + pow( (- f[i][(j*N+kp2)*N+l] + 8*f[i][(j*N+kp1)*N+l] - 8*f[i][(j*N+km1)*N+l] + f[i][(j*N+km2)*N+l]) / (12*dx), 2.0 )
             + pow( (- f[i][(j*N+k)*N+lp2] + 8*f[i][(j*N+k)*N+lp1] - 8*f[i][(j*N+k)*N+lm1] + f[i][(j*N+k)*N+lm2]) / (12*dx), 2.0 ))/2;
    }
default:
    std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
    exit(1);
    }
}


double Field::gradient_energy( double* f )
{
    double gradient_energy = 0;
    
  // #pragma omp parallel for reduction (+:gradient_energy) schedule( static ) num_threads ( num_threads )
    for( int j = 0; j < N; ++j ){
        int jp1 = (j == N-1)?     0: j+1;
        int jp2 = (j >= N-2)? j-N+2: j+2;
        
        int jm1 = (j ==   0)?   N-1: j-1;
        int jm2 = (j <    2)? j+N-2: j-2;
        
        switch (dim)
        {
        case 1:
         {
             gradient_energy +=  pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2.0 );
 //            std::cout << "pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2 ) = " << pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2 ) << "f["<< jp2 << "] = " << f[jp2] << "f["<< jp1 << "] = " << f[jp1] << "f["<< jm1 << "] = " << f[jm1]<< "f["<< jm2 << "] = " << f[jm2] << std::endl;
        break;
         }
        case 2:
         {
    for( int k = 0; k < N; ++k ){
        int kp1 = (k == N-1)?     0: k+1;
        int kp2 = (k >= N-2)? k-N+2: k+2;
        int km1 = (k ==   0)?   N-1: k-1;
        int km2 = (k <    2)? k+N-2: k-2;
//                gradient_energy +=  pow( (f[jp1*N+k] - f[jm1*N+k]) / (2*dx), 2 );
//                 gradient_energy += pow( ( f[j*N+kp1]  - f[j*N+km1] ) /(2*dx), 2 );
      gradient_energy += pow( ( - f[jp2*N+k] + 8*f[jp1*N+k] - 8*f[jm1*N+k] + f[jm2*N+k] ) / (12*dx), 2.0 );
        gradient_energy += pow( ( - f[j*N+kp2] + 8*f[j*N+kp1] - 8*f[j*N+km1] + f[j*N+km2] ) / (12*dx), 2.0 );
    }
        break;
         }
        case 3:
         {
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
                     gradient_energy += pow( ( - f[(jp2*N+k)*N+l] + 8*f[(jp1*N+k)*N+l] - 8*f[(jm1*N+k)*N+l] + f[(jm2*N+k)*N+l] ) / (12*dx), 2.0 );
                     gradient_energy += pow( ( - f[(j*N+kp2)*N+l] + 8*f[(j*N+kp1)*N+l] - 8*f[(j*N+km1)*N+l] + f[(j*N+km2)*N+l] ) / (12*dx), 2.0 );
                     gradient_energy += pow( ( - f[(j*N+k)*N+lp2] + 8*f[(j*N+k)*N+lp1] - 8*f[(j*N+k)*N+lm1] + f[(j*N+k)*N+lm2] ) / (12*dx), 2.0 );
                 }
             }
        break;
         }
        default:
           std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
                    exit(1);
        }//switch
                
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
        switch (dim)
        {
        case 1:
         {
        int idx = j;
         potential_energy += V_lattice( f, idx, a );
    // std::cout << "V_lattice( f, idx, a ) = " << V_lattice( f, idx, a ) << std::endl;
        break;
         }
        case 2:
         {
//        #pragma omp simd reduction(+:potential_energy)
    for( int k = 0; k < N; ++k )
    {
        int idx = j*N + k;
        potential_energy += V_lattice( f, idx, a );
    }
        break;
         }
        case 3:
         {
    for( int k = 0; k < N; ++k ){
//        #pragma omp simd reduction(+:potential_energy)
        for( int l = 0; l < N; ++l ){
            int idx = ( j*N + k)*N + l;
            potential_energy += V_lattice( f, idx, a );
        }
    }
        break;
         }
        default:
           std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
                    exit(1);
        }//switch
       
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
            {
				int idx = j;
				average += f[idx];
				break;
            }
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
            {
				int idx = j;
				variance += pow( f[idx] - _average[i], 2.0 );
				break;
            }
			case 2:
                #pragma omp simd reduction(+:variance)
				for( int k = 0; k < N; ++k ){
					int idx = j*N + k;
                    variance += pow( f[idx] - _average[i], 2.0 );
                
				}
				break;
			case 3:
				for( int k = 0; k < N; ++k ){
                    #pragma omp simd reduction(+:variance)
					for( int l = 0; l < N; ++l ){
						int idx = (j*N + k)*N + l;
						variance += pow( f[idx] - _average[i], 2.0 );
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
//    std::cout << "pow(a,4.0) = " << pow(a,4.0) << std::endl;
//    std::cout << "pow(rescale_A/rescale_B,2.0) = " << pow(rescale_A/rescale_B,2.0) << std::endl;
//    std::cout << "V(f_MPl[0],f_MPl[1],f_MPl[2]) = " << V(f_MPl[0],f_MPl[1],f_MPl[2]) << std::endl;
    
    return pow(a,4.0)*pow(rescale_A/rescale_B,2.0)*V(f_MPl[0],f_MPl[1],f_MPl[2]); }

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
          
            
            return pow(a,3.0)*rescale_A*V_1(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //sigma
            
           
//        case 1: return -pow(a,3.0)*rescale_A*V_2(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //psi
        case 1: return pow(a,3.0)*rescale_A*V_2(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //psi
         
        case 2:
            
            
//            std::cout << "pow(a,3.0)*rescale_A*V_3(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B) = " <<  pow(a,3.0)*rescale_A*V_3(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B) << std::endl;
            
            return pow(a,3.0)*rescale_A*V_3(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //phi
           
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
        case 0: return pow(a/rescale_B,2.0)*V_11(f_MPl[0],f_MPl[1],f_MPl[2]); //sigma
        case 1: return pow(a/rescale_B,2.0)*V_22(f_MPl[0],f_MPl[1],f_MPl[2]); //psi
        case 2: return pow(a/rescale_B,2.0)*V_33(f_MPl[0],f_MPl[1],f_MPl[2]); //phi
        default:  Logout( "Parameter 'i' in ddV_lattice must be 0 ~ 2. \n" );
            exit(1);
    }
    
}


void Field::effective_mass(double mass_sq[], double *field_values){

    mass_sq[0] = pow(exp(OSCSTART)/rescale_B,2.0)*V_11(field_values[0],field_values[1],field_values[2]); //sigma
    mass_sq[1] = pow(exp(OSCSTART)/rescale_B,2.0)*V_22(field_values[0],field_values[1],field_values[2]); //psi
    mass_sq[2] = pow(exp(OSCSTART)/rescale_B,2.0)*V_33(field_values[0],field_values[1],field_values[2]); //phi

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
    
    
    switch( dim )
    {
        case 1:
        {
        if(m == m_start)
        {
            
            //subtract zero mode
            for( int j = 0; j < N; ++j ){
                int idx_fluc = j;
                
                if(i==4)//Curvature perturbation
                {
                    f_fluc[i-1][idx_fluc] = f[i-1][idx_fluc] - average(f[i-1], i-1);
                    df_fluc[i-1][idx_fluc] = df[i-1][idx_fluc] - average(df[i-1], i-1);

                }
                else{
               f_fluc[i][idx_fluc] = f[i][idx_fluc] - average(f[i], i);
                }
            }
            //transform from real space to phase space (Needs to be done only once)
            if(i==4)//Curvature perturbation
            {
                //Fourier Transform
                DFT_r2cD1( f_fluc[i-1], f_fluc_k[i-1] );
                DFT_r2cD1( df_fluc[i-1], df_fluc_k[i-1] );
            }else{
            DFT_r2cD1( f_fluc[i], f_fluc_k[i] );
            }
        }
        
//            for( int m = 0; m < N; ++m ){
//                std::cout << "1: f_fluc_k[" << i << "][" << m << "]" << f_fluc_k[i][m] << std::endl;
//            }
     
    
         //f[][0] corresponds to Re(f_{i=0})
         //f[][2],f[][3] corresponds to the Re and Im of f_{i=1}
         // ...
         //f[][1] corresponds to Re(f_{i=N/2})
    
//    for( int m = 0; m < N; ++m ){
//        std::cout << "2: f_fluc_k[" << i << "][" << m << "]" << f_fluc_k[i][m] << std::endl;
//    }
        int j = m;
            
            if(i==4)
            {
                if(m==0)
                {// zero-frequency (DC)
                PS[i][m] = m*pow(f_fluc_k[i][0]/N,2.0);
                    
                }else if(m == N/2){//Nyquist frequency
                    PS[i][m] = m*pow(f_fluc_k[i][1]/N,2.0);
                }else{
                    
                   PS[i][m] = m*( pow(f_fluc_k[i][2*j]/N,2.0) + pow(f_fluc_k[i][2*j+1]/N,2.0));
                   
                }
            }
            else
            {
             if(m==0)
             {// zero-frequency (DC)
                 PS[i][m] = m*pow(f_fluc_k[i][0]/N,2.0);
                 
             }else if(m == N/2){//Nyquist frequency
                 PS[i][m] = m*pow(f_fluc_k[i][1]/N,2.0);
             }else{
                 
                PS[i][m] = m*( pow(f_fluc_k[i][2*j]/N,2.0) + pow(f_fluc_k[i][2*j+1]/N,2.0));
                
             }
            }
    
        PS[i][m] *= 2; //  the conjugate part needs to be taken into account as well, so double the value
    

        
//     std::cout << " PS[" << i << "][" << idx << "] = " << PS[i][idx] << std::endl;
        return PS[i][m];
        }
        case 2:
        {
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
                DFT_r2cD2( f_fluc[i], f_fluc_k[i], f_fluc_k_nyquist_2d[i] );
        
//        std::cout << "af f_fluc[" << i << "][20]" << f_fluc[i][20] << std::endl;
//         std::cout << "af f_fluc_k[" << i << "][20]" << f_fluc_k[i][20] << std::endl;
//        std::cout << "af f_fluc_k_nyquist_2d[" << i << "][20]" << f_fluc_k_nyquist_2d[i][20] << std::endl;
        
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
              if ( m-1 < sqrt(pow((double)J,2.0)+pow((double)k,2.0)) && sqrt(pow((double)J,2.0)+pow((double)k,2.0)) <= m  )
              {
                  
              if(k==N/2)//Nyquist
              {
                 
                      if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=k, 0<=j <=N/2) doesn't have an equal number of counterparts in the range (0<=k, N/2+1<j<=N-1), we consider the range (1<=k, 1<=j<=N/2) and multiply by 4 and add the Nyquist  (k=N/2, j=0) x 2 and (k=0, j=N/2) x 2.
                      {
                          // This corresponds to (k=N/2, j=0). We multiply the total sum by 4 at the end so we multiply by 0.5 here
                          PS[i][m] += 0.5*m*(pow(f_fluc_k_nyquist_2d[i][2*j],2.0) + pow(f_fluc_k_nyquist_2d[i][2*j+1],2.0))/(pow((double)N,4.0));
                      }
//                              std::cout << "1: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                      else {//No modes correspond to this branch
                          PS[i][m] += m*(pow(f_fluc_k_nyquist_2d[i][2*j],2.0) + pow(f_fluc_k_nyquist_2d[i][2*j+1],2.0))/(pow((double)N,4.0));
                      }
                      
                  
                      
                  
                  
              }else//Modes besides Nyquist
              {
                  
//                  std::cout << "m = " << m << std::endl;
//                  std::cout << "j = " << j << std::endl;
//                  std::cout << "J = " << J << std::endl;
//                  std::cout << "k = " << k << std::endl;
//                  std::cout << "sqrt(pow(J,2.0)+pow(k,2.0)) = " << sqrt(pow(J,2.0)+pow(k,2.0)) << std::endl;
                  
                 
                      if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=k, 0<=j <=N/2) doesn't have an equal number of counterparts in the range (0<=k, N/2+1<j<=N-1), we consider the range (1<=k, 1<=j<=N/2) and multiply by 4 and add the Nyquist  (k=N/2, j=0) x 2 and (k=0, j=N/2) x 2.
                      {
                          if(j==N/2)
                          {// This corresponds to (k=0, j=N/2). We multiply the total sum by 4 at the end so we multiply by 0.5 here
                              PS[i][m] +=
                              0.5*m*(pow(f_fluc_k[i][j*N+2*k],2.0) + pow(f_fluc_k[i][j*N+2*k+1],2.0))/(pow((double)N,4.0));
                          }else{//No modes correspond to this branch
                              PS[i][m] +=
                              m*(pow(f_fluc_k[i][j*N+2*k],2.0) + pow(f_fluc_k[i][j*N+2*k+1],2.0))/(pow((double)N,4.0));
                          }
                             
                              
//                               std::cout << "3: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                          
                      }
                      else
                      {
                    PS[i][m] +=  m*(pow(f_fluc_k[i][j*N+2*k],2.0) + pow(f_fluc_k[i][j*N+2*k+1],2.0))/(pow((double)N,4.0));
                          
//                          std::cout << "4: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                      }
                      
              }//if ( m-1 < sqrt(pow(J,2.0)+pow(k,2.0)) && sqrt(pow(J,2.0)+pow(k,2.0)) <= m  )
                  
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
    
        }
  case 3:
        {
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
        DFT_r2cD3( f_fluc[i], f_fluc_k[i], f_fluc_k_nyquist_3d[i] );
       
    }
    
    
    
   
    for(int j = 0; j < N; ++j )
    {
        J = (j <= N/2 ? j : j-N);
        
        for(int k = 0; k < N; k++)
        {
            K = (k <= N/2 ? k : k-N);
        
        for(int l = 0; l <= N/2; ++l )
        {
            
            if ( m-1 < sqrt(pow((double)J,2.0)+pow((double)K,2.0)+pow((double)l,2.0)) && sqrt(pow((double)J,2.0)+pow((double)K,2.0)+pow((double)l,2.0)) <= m  )
            {
            
            if(l==N/2)//Nyquist
            {
                
                    if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=l, 0<=k<=N/2, 0<=j<=N/2) doesn't have an equal number of counterparts in the range (0<=l && !(0<=k <=N/2, 0<=j <=N/2) ), we consider the range (1<=l, 1<=k<=N/2, 1<=j <=N/2 ) and multiply by 8 and add the Nyquist (l=N/2, k=0,j=0) x 2, (l=0, k=0, j=N/2) x 2 and (l=0, k=N/2, j=0) x 2.
                
                    {
                        // This corresponds to (l=N/2, k=0,j=0). We multiply the total sum by 8 at the end so we multiply by 0.25 here
                        PS[i][m] += 0.25*m*(pow(f_fluc_k_nyquist_3d[i][j][2*k],2.0) + pow(f_fluc_k_nyquist_3d[i][j][2*k+1],2.0))/(pow((double)N,6.0));
                        
//                                                      std::cout << "1: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                    }
                    else
                    {//No modes correspond to this branch
                        PS[i][m] += m*(pow(f_fluc_k_nyquist_3d[i][j][2*k],2.0) + pow(f_fluc_k_nyquist_3d[i][j][2*k+1],2.0))/(pow((double)N,6.0));
                        
//                                                   std::cout << "2: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                    }
                    
                
                
            }else//Modes besides Nyquist
            {
                
                //                  std::cout << "m = " << m << std::endl;
                //                  std::cout << "j = " << j << std::endl;
                //                  std::cout << "J = " << J << std::endl;
                //                  std::cout << "k = " << k << std::endl;
                //                  std::cout << "sqrt(pow(J,2.0)+pow(k,2.0)) = " << sqrt(pow(J,2.0)+pow(k,2.0)) << std::endl;
                
               
                    if(m==N/2)//Since the outtermost mode (m=N/2) in the range (0<=l, 0<=k<=N/2, 0<=j<=N/2) doesn't have an equal number of counterparts in the range (0<=l && !(0<=k <=N/2, 0<=j <=N/2) ), we consider the range (1<=l, 1<=k<=N/2, 1<=j <=N/2 ) and multiply by 8 and add the Nyquist (l=N/2, k=0,j=0) x 2, (l=0, k=0, j=N/2) x 2 and (l=0, k=N/2, j=0) x 2.
                    {
                        if((j==N/2 && k==0) || (j==0 && k==N/2) )
                        {// This corresponds to (l=0, k=0, j=N/2) and (l=0, k=N/2, j=0). We multiply the total sum by 8 at the end so we multiply by 0.25 here
                            PS[i][m] +=
                            0.25*m*(pow(f_fluc_k[i][(j*N + k)*N + 2*l],2.0) + pow(f_fluc_k[i][(j*N + k)*N + 2*l+1],2.0))/(pow((double)N,6.0));
                        }else
                        {
                            PS[i][m] += m*(pow(f_fluc_k[i][(j*N + k)*N + 2*l],2.0) + pow(f_fluc_k[i][(j*N + k)*N + 2*l+1],2.0))/(pow((double)N,6.0));
                            //                           std::cout << "2: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                        }
                        
                        
//                                                       std::cout << "3: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                        
                    }
                    else
                    {
                        PS[i][m] +=
                        m*(pow(f_fluc_k[i][(j*N + k)*N + 2*l],2.0) + pow(f_fluc_k[i][(j*N + k)*N + 2*l+1],2.0))/(pow((double)N,6.0));
                        
//                                                  std::cout << "4: PS[" << i << "][" << m << "] = " << PS[i][m] << std::endl;
                    }
                    
                
            }
                
        } //if ( m-1 < sqrt(pow(J,2.0)+pow(K,2.0)+pow(l,2.0)) && sqrt(pow(J,2.0)+pow(K,2.0)+pow(l,2.0)) <= m  )
            
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
        }
        default:
            std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
            exit(1);
    }//switch
    

}



void Field::finalize(double** f, double** df, LeapFrog* leapfrog, double radiation_pr, double** lattice_var )
{
    double a = leapfrog->a();
    double da = leapfrog->da();
    
    a_lattice_end = exp(OSCSTART)*a;
    
    
    //Zero modes
    for (int i = 0; i< num_fields - 1; i++){
        
        int j = (num_fields -1) + i;
        
        for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++)
        {
            //
            //    Logout(" lattice_var[%d][%d] = %2.5e \n", lattice_loop, i, lattice_var[lattice_loop][i]);
            //
            lattice_var[lattice_loop][i] = average(f[i], i)/(rescale_A*a) ;//0,1,2
            
            lattice_var[lattice_loop][j] = (rescale_B/rescale_A)*( average(df[i], i)  - (da/a)*average(f[i], i) )/pow(a,2.0);//3,4,5
            
//            Logout("lattice_var[%d][%d] = %2.5e \n",lattice_loop,i, lattice_var[lattice_loop][i] );
//            Logout(" lattice_var[%d][%d] = %2.5e \n",lattice_loop,j,lattice_var[lattice_loop][j] );
            
        }
        
        
        
    }
    
    //Radiation
    for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++)
    {
        
        lattice_var[lattice_loop][6] =  pow((rescale_B/(rescale_A*pow(a,2.0))),2.0)*radiation_pr;//6
        
    }
    
    
    //Perturbations
    
    int m_start, m_end;
    
    if (k_lattice_grid_min_MPl < kfrom_MPl_lattice)
    {
        m_start = outrange_num + 1;
    }else
    {
        m_start = 1;
    }
    
    m_end = N/2;
    
    double Rescale_var;
    double Rescale_var_D3 =sqrt(pow(dx,6.0)/(pow(rescale_B*L,3.0)));
    
    double Re_scalar_pert_pr, Im_scalar_pert_pr;
    double Re_scalar_pert_deriv_pr, Im_scalar_pert_deriv_pr;
    double Re_metric_pert_pr, Im_metric_pert_pr;
    double Re_metric_pert_deriv_pr, Im_metric_pert_deriv_pr;

   
        
    switch (dim)
    {
    case 1:
     {
         for (int i = 0; i < num_fields; i++){
             //subtract zero mode
             for( int j = 0; j < N; ++j ){
                 int idx_fluc = j;
                 f_fluc[i][idx_fluc] = f[i][idx_fluc] - average(f[i], i);
                 df_fluc[i][idx_fluc] = df[i][idx_fluc] - average(df[i], i);
                 
             }
             
             //Fourier Transform
             DFT_r2cD1( f_fluc[i], f_fluc_k[i] );
             DFT_r2cD1( df_fluc[i], df_fluc_k[i] );
             
         }
             
             //f[][0] corresponds to Re(f_{i=0})
             //f[][2],f[][3] corresponds to the Re and Im of f_{i=1}
             // ...
             //f[][1] corresponds to Re(f_{i=N/2})
             
             //    for( int m = 0; m < N; ++m ){
             //        std::cout << "2: f_fluc_k[" << i << "][" << m << "]" << f_fluc_k[i][m] << std::endl;
             //    }
             
             int N_m = 2; //number of fields in the range where m-1 < |m| <= m
         
             for (int m = m_start; m < m_end + 1; ++m  )
             {
                 
             for (int i = 0; i < num_fields; i++)
             {
                 
                 
             int j=m;
             
     //        std::cout << "i = " << i << ", m = " << m  << std::endl;
             
             if(m == N/2)
             {//Nyquist frequency
                 
                 if(i == 3)//Metric Peturbation
                 {
                     for(int num = 0 ; num < 3; num++ )
                     {
                         if(num == 2)
                         {
                             
                             Re_metric_pert_pr = sqrt(N_m*pow(f_fluc_k[i][1],2.0)/(N_m - 1));
                                               
                             Re_metric_pert_deriv_pr = sqrt(N_m*pow(df_fluc_k[i][1],2.0)/(N_m - 1));
                             
                             //27
                             lattice_var[m-1][25+num] = Re_metric_pert_pr/(pow(a,2.0));
                             //51
                             lattice_var[m-1][49+num] = 0;
                             
                             //30
                             lattice_var[m-1][28+num] = rescale_B*(Re_metric_pert_deriv_pr - 2*Re_metric_pert_pr*da/a)/(pow(a,3.0));
                             //54
                             lattice_var[m-1][52+num] = 0;
                             
                         }else{
                             //25,26
                             lattice_var[m-1][25+num] = 0;
                             //49,50
                             lattice_var[m-1][49+num] = 0;
                             
                             //28,29
                             lattice_var[m-1][28+num] = 0;
                             //52,53
                             lattice_var[m-1][52+num] = 0;
                         }
                     }
                     
                 }else//Scalar Peturbation
                 {
                     
                     for(int num = 0 ; num < 3; num++ )
                     {
                         if(num == i)
                         {
                             Re_scalar_pert_pr = sqrt(N_m*pow(f_fluc_k[i][1],2.0)/(N_m - 1));
                             
                             Re_scalar_pert_deriv_pr = sqrt(N_m*pow(df_fluc_k[i][1],2.0)/(N_m - 1));
                             //7,11,15
                             lattice_var[m-1][7 + 3*i + num] = Re_scalar_pert_pr/(rescale_A*a);
                             //31,35,39
                             lattice_var[m-1][31 + 3*i + num] = 0;
                             
                             //16,20,24
                             lattice_var[m-1][16 + 3*i + num] = rescale_B*(Re_scalar_pert_deriv_pr - Re_scalar_pert_pr*da/a)/(rescale_A*pow(a,2.0));
                             //40,44,48
                             lattice_var[m-1][40 + 3*i + num] = 0;
                         }else
                         {
                             //8,9,10,12,13,14
                             lattice_var[m-1][7 + 3*i + num] = 0;
                             //32,33,34,36,37,38
                             lattice_var[m-1][31 + 3*i + num] = 0;
                             
                             //17,18,19,21,22,23
                             lattice_var[m-1][16 + 3*i + num] = 0;
                             //41,42,43,45,46,47
                             lattice_var[m-1][40 + 3*i + num] = 0;
                         }
                     }
                     
                 }
                 
                 
                 
             }else{//Besides Nyquist frequency
                 
                 if(i == 3)//Metric Peturbation
                 {
                     for(int num = 0 ; num < 3; num++ )
                     {
                         if(num == 2)
                         {
                             
                             Re_metric_pert_pr = sqrt(N_m*pow(f_fluc_k[i][2*j],2.0)/(N_m - 1));
                             Im_metric_pert_pr = sqrt(N_m*pow(f_fluc_k[i][2*j+1],2.0)/(N_m - 1));
                                               
                             Re_metric_pert_deriv_pr = sqrt(N_m*pow(df_fluc_k[i][2*j],2.0)/(N_m - 1));
                             Im_metric_pert_deriv_pr = sqrt(N_m*pow(df_fluc_k[i][2*j+1],2.0)/(N_m - 1));
                             
                             //27
                             lattice_var[m-1][25+num] = Re_metric_pert_pr/(pow(a,2.0));
                             //51
                             lattice_var[m-1][49+num] = Im_metric_pert_pr/(pow(a,2.0)) ;
                             
                             //30
                             lattice_var[m-1][28+num] = rescale_B*(Re_metric_pert_deriv_pr - 2*Re_metric_pert_pr*da/a)/(pow(a,3.0));
                             //54
                             lattice_var[m-1][52+num] = rescale_B*(Im_metric_pert_deriv_pr - 2*Im_metric_pert_pr*da/a)/(pow(a,3.0));
                         }else{
                             //25,26
                             lattice_var[m-1][25+num] = 0;
                             //49,50
                             lattice_var[m-1][49+num] = 0;
                             
                             //28,29
                             lattice_var[m-1][28+num] = 0;
                             //52,53
                             lattice_var[m-1][52+num] = 0;
                         }
                     }
                     
                 }else//Scalar Peturbation
                 {
                     
                     for(int num = 0 ; num < 3; num++ )
                     {
                         if(num == i)
                         {
                             
                             Re_scalar_pert_pr = sqrt(N_m*pow(f_fluc_k[i][2*j],2.0)/(N_m - 1));
                             Im_scalar_pert_pr = sqrt(N_m*pow(f_fluc_k[i][2*j+1],2.0)/(N_m - 1));
                             
                             Re_scalar_pert_deriv_pr = sqrt(N_m*pow(df_fluc_k[i][2*j],2.0)/(N_m - 1));
                             Im_scalar_pert_deriv_pr = sqrt(N_m*pow(df_fluc_k[i][2*j+1],2.0)/(N_m - 1));
                             
                             //7,11,15
                             lattice_var[m-1][7 + 3*i + num] = Re_scalar_pert_pr/(rescale_A*a);
                             
     //                        if(m==m_start)
     //                        {
     //                        std::cout << "i = " << i << ": lattice_var[" << m_start-1 << "][" << 7 + 3*i + num << "] = " <<  lattice_var[m-1][7 + 3*i + num] << std::endl;
     //                        }
                             //31,35,39
                             lattice_var[m-1][31 + 3*i + num] = Im_scalar_pert_pr/(rescale_A*a);
                             
                             //16,20,24
                             lattice_var[m-1][16 + 3*i + num] = rescale_B*(Re_scalar_pert_deriv_pr - Re_scalar_pert_pr*da/a)/(rescale_A*pow(a,2.0));
                             //40,44,48
                             lattice_var[m-1][40 + 3*i + num] = rescale_B*(Im_scalar_pert_deriv_pr - Im_scalar_pert_pr*da/a)/(rescale_A*pow(a,2.0));
                             
                         }else
                         {
                             //8,9,10,12,13,14
                             lattice_var[m-1][7 + 3*i + num] = 0;
                             //32,33,34,36,37,38
                             lattice_var[m-1][31 + 3*i + num] = 0;
                             
                             //17,18,19,21,22,23
                             lattice_var[m-1][16 + 3*i + num] = 0;
                             //41,42,43,45,46,47
                             lattice_var[m-1][40 + 3*i + num] = 0;
                         }
                         
                     }
                     
                 }
                 
             }//Besides Nyquist frequency
             
                 
         }//for (int i = 0; i < num_fields; i++)
             
                 
     //           for (int z=0;z<N_pert;z++) Logout("m == m_start end: lattice_var[%d][%d] = %2.5e \n",m-1,z,lattice_var[m-1][z] );
     //
                 
                 
         //Set rescale variables (1D case)
         Rescale_var = Rescale_var_D3*sqrt(2*M_PI)*L/(pow(dx,2.0)*(2*M_PI*m/L));
                 
                 //Rescale program variables back to their original variables
                 for(int q = N_zero; q < N_pert; q++)
                 {
                     lattice_var[m-1][q] *= Rescale_var;
                 }
                 
     //            for (int z=0;z<N_pert;z++) Logout("m == m_start end: lattice_var[%d][%d] = %2.5e \n",m-1,z,lattice_var[m-1][z] );
                 
         }//for (int m = m_start; m < m_end + 1; ++m  )
             
         
         
             
    break;
     }
    case 2:
     {
         for (int i = 0; i < num_fields; i++){

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
         
         DFT_r2cD2( f_fluc[i], f_fluc_k[i], f_fluc_k_nyquist_2d[i] );
             
         }
         
         
    break;
     }
    case 3:
     {
         for (int i = 0; i < num_fields; i++){

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
         
         DFT_r2cD3( f_fluc[i], f_fluc_k[i], f_fluc_k_nyquist_3d[i] );
             
         }
         
    break;
     }
    default:
       std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
                exit(1);
    }//switch
    
        
    delete [] f[0];
    delete [] df[0];
    delete [] f;
    delete [] df;
}

