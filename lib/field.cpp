#include <cmath>
#include <random>
#include "field.hpp"
#include "utilities.hpp"


void initialize( double**& f, double**& df, Field* field)
{   
	f = new double* [num_fields];
	df = new double* [num_fields];

	switch( dim )
	{
		case 1:
			f[0] = new double [num_fields*N];
			df[0] = new double [num_fields*N];
			for( int i = 0; i < num_fields; ++i )
			{
				f[i] = f[0] + i*N;
				df[i] = df[0] + i*N;
			}
			break;
		case 2:
			f[0] = new double [num_fields*N*N];
			df[0] = new double [num_fields*N*N];
			for( int i = 0; i < num_fields; ++i )
			{
				f[i] = f[0] + i*N*N;
				df[i] = df[0] + i*N*N;
			}
			break;
		case 3:
			f[0] = new double [num_fields*N*N*N];
			df[0] = new double [num_fields*N*N*N];
			for( int i = 0; i < num_fields; ++i )
			{
				f[i] = f[0] + i*N*N*N;
				df[i] = df[0] + i*N*N*N;
			}
			break;
		default:
			std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
			exit(1);
	}
 
    double p2;
    double dp2=pw2(2*M_PI/L);
  double* mass_sq = new double [num_fields];
  double* initial_field_values = new double [num_fields];
  double* initial_field_derivs = new double [num_fields];
    
    
    #if  dim==1
      double pz;
    #elif dim==2
    int jconj;
    double py,pz;
  
    
    double* fnyquist = new double [2*N];
    double* fdnyquist = new double [2*N];
    
    #elif dim==3
    int kconj, jconj;
    double px,py,pz;
    
    double** fnyquist = new double* [N];
    double** fdnyquist = new double* [N];
    fnyquist[0] = new double [N*2*N];
    fdnyquist[0] = new double [N*2*N];
    for( int i = 0; i < N; ++i )
    {
        fnyquist[i] = fnyquist[0] + i*2*N;
        fdnyquist[i] = fdnyquist[0] + i*2*N;
    }
     #endif
    
   // Courant condition
    if(dt>dx/sqrt(dim)){
      std::cout <<  "Time step too large. The ratio dt/dx is currently" << dt/dx <<
        "but for stability should never exceed 1/sqrt("<< dim <<") (" << 1/sqrt(dim) << ")"<< std::endl;
        exit(1);
    }
    
    
    for( int i = 0; i < num_fields; ++i ){
        
        initial_field_values[i] = initfield[i];
        initial_field_derivs[i] = initderivs[i];
        field->effective_mass(mass_sq, initial_field_values);
      std::cout << "mass_sq[0] = " << mass_sq[i]*pw2(m)*pow(sqrt(8*M_PI),2) << std::endl;
        
       /* if(expansion)
        {
            aeffective_mass( mass_sq, initial_field_values,a);
        }else{
            effective_mass( mass_sq, initial_field_values);
        }*/
        
  #if  dim==1
      
        //k=N/2
        p2 = dp2*pw2(N/2);
        set_mode(p2, mass_sq[i], &f[i][1], &df[i][1], 1);
        
        //Loop over gridpoints.
        for(int k = 1; k < N/2; k++ ){
            pz = k;
            p2 = dp2*pw2(pz);
             set_mode(p2, mass_sq[i], &f[i][2*k], &df[i][2*k], 0);
        }
        //k=0
        f[i][0] = 0.;
        fd[i][0] = 0.;
#elif dim==2
      
        for(int j = 0; j < N; j++ ){
            py = (j <= N/2 ? j : j-N);
            
            for(int k = 1; k < N/2; k++ ){
                pz = k;
                p2=dp2*(pw2(py)+pw2(pz));
                set_mode(p2, mass_sq[i], &f[i][j*N+2*k], &df[i][j*N+2*k], 0);
        
            }
        
            if(j > N/2)
                
            {
                jconj = N-j;
                
                //k=0
                p2 = dp2*pw2(py);
                set_mode(p2, mass_sq[i], &f[i][j*N], &df[i][j*N], 0);
                f[i][jconj*N] = f[i][j*N];
                f[i][jconj*N+1] = -f[i][j*N+1];
                df[i][jconj*N] = df[i][j*N];
                df[i][jconj*N+1] = -df[i][j*N+1];
                //k=N/2
                p2 = dp2*(pw2(py)+pw2(N/2));
                set_mode(p2, mass_sq[i], &fnyquist[2*j], &fdnyquist[2*j], 0);
                fnyquist[2*jconj] = fnyquist[2*j];
                fnyquist[2*jconj+1] = -fnyquist[2*j+1];
                fdnyquist[2*jconj] = fdnyquist[2*j];
                fdnyquist[2*jconj+1] = -fdnyquist[2*j+1];
                
            }
            else if(j == 0 || j == N/2){
                
                p2 = dp2*pw2(py); //k = 0
                if(p2 > 0.)
                   set_mode(p2, mass_sq[i], &f[i][j*N], &df[i][j*N], 1);
                
                p2 = dp2*(pw2(py)+pw2(N/2));  //k = N/2
                set_mode(p2, mass_sq[i], &fnyquist[2*j], &fdnyquist[2*j], 1);
                
            }
        }  //k = 0, j = 0
        f[i][0] = 0.;
        f[i][1] = 0.;
        df[i][0] = 0.;
        df[i][1] = 0.;
       
       
      /*
        for( int j = 0; j < N; ++j ){
            //#pragma omp parallel for schedule( static ) num_threads( num_threads )
            for( int k = 0; k < N; ++k ){
                int idx = j*N + k;
                std::cout << "f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
                std::cout << "df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
               
            }
        }
        for( int j = 0; j < 2*N; ++j ){
        std::cout << " fnyquist[" << j << "] = " << fnyquist[j] << std::endl;
        std::cout << " fdnyquist[" << j << "] = " << fdnyquist[j] << std::endl;
        }*/
   
#elif dim==3
        for(int j; j < N; j++){
            px = (j <= N/2 ? j : j-N);
            jconj = (j==0 ? 0 : N-j);
            
            for(int k = 0; k < N; k++)
            {
                 py = (k <= N/2 ? k : k-N);
                
                //0<l<N/2
                for(int l = 1; l < N/2; l++){
                    pz = l;
                    p2 = dp2*(pw2(px)+pw2(py)+pw2(pz));
                    set_mode(p2,mass_sq[i], &f[i][(j*N + k)*N + 2*l], &df[i][(j*N + k)*N + 2*l], 0);
                }
                
                if(k > N/2 || (j > N/2 && (k == 0 || k == N/2))){
                    kconj = (k == 0 ? 0 : N-k);
                    //l=0
                     p2 = dp2*(pw2(px)+pw2(py));
                    set_mode(p2,mass_sq[i], &f[i][(j*N + k)*N], &df[i][(j*N + k)*N], 0);
                    f[i][(jconj*N + kconj)*N] = f[i][(j*N + k)*N];
                    f[i][(jconj*N + kconj)*N+1] = -f[i][(j*N + k)*N+1];
                    df[i][(jconj*N + kconj)*N] = df[i][(j*N + k)*N];
                    df[i][(jconj*N + kconj)*N+1] = -df[i][(j*N + k)*N+1];
                    //l=N/2
                    p2 = dp2*(pw2(px)+pw2(py)+pw2(N/2));
                    set_mode(p2, mass_sq[i], &fnyquist[j][2*k], &fdnyquist[j][2*k], 0);
                    fnyquist[jconj][2*kconj] = fnyquist[j][2*k];
                    fnyquist[jconj][2*kconj+1] = -fnyquist[j][2*k+1];
                    fdnyquist[jconj][2*kconj] = fdnyquist[j][2*k];
                    fdnyquist[jconj][2*kconj+1] = -fdnyquist[j][2*k+1];
                    }
                else if((j == 0 || j == N/2) && (k == 0 || k == N/2))
                {
                    p2 = dp2*(pw2(px)+pw2(py)); //l=0
                    if(p2 > 0.)
                      set_mode(p2,mass_sq[i], &f[i][(j*N + k)*N], &df[i][(j*N + k)*N], 1);
                    
                    p2=dp2*(pw2(px)+pw2(py)+pw2(N/2)); //l=N/2
                     set_mode(p2, mass_sq[i], &fnyquist[j][2*k], &fdnyquist[j][2*k], 1);
                }
            }
        }
        f[i][0] = 0.;
        f[i][1] = 0.;
        df[i][0] = 0.;
        df[i][1] = 0.;

#endif
        
      

	

#if  dim==1
    DFT_c2rD1( f[i] ); //transform from phase space to real space
    DFT_c2rD1( df[i] );
        //Add zeromode
       for( int j = 0; j < N; ++j ){
           int idx = j;
           f[i][idx] += initial_field_values[i];
           df[i][idx] += initial_field_derivs[i];
       }
        
#elif dim==2
      
        
         std::cout << "before fourier transform" << std::endl;
        std::cout << "f[0][4]" << f[i][4] << std::endl;
        std::cout << "df[0][4]" << df[i][4] << std::endl;
        std::cout << "f[0][10]" << f[i][10] << std::endl;
        std::cout << "df[0][10]" << df[i][10] << std::endl;
     /*   double   fieldsum, dfieldsum;
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N; ++k ){
                int idx = j*N + k;
                // std::cout << "f[0]["<< idx <<"]" << f[i][idx] << std::endl;
                //  std::cout << "df[0]["<< idx <<"]" << df[i][idx] << std::endl;
                fieldsum +=  f[i][idx];
                dfieldsum += df[i][idx];
                // std::cout << "f[0]["<< idx <<"]" << f[i][idx] << std::endl;
                // std::cout << "df[0]["<< idx <<"]" << df[i][idx] << std::endl;
                
            }
        }
        for( int j = 0; j < dim; ++j ) fieldsum /= N;
        for( int j = 0; j < dim; ++j ) dfieldsum /= N;
        std::cout << "fieldsum = " << fieldsum << std::endl;
        std::cout << "dfieldsum = " << dfieldsum << std::endl;*/
        
    DFT_c2rD2( f[i] , fnyquist ); //transform from phase space to real space
   DFT_c2rD2( df[i], fdnyquist);
        
         std::cout << "after fourier transform" << std::endl;
        std::cout << "f[0][4]" << f[i][4] << std::endl;
        std::cout << "df[0][4]" << df[i][4] << std::endl;
        std::cout << "f[0][10]" << f[i][10] << std::endl;
        std::cout << "df[0][10]" << df[i][10] << std::endl;
         /*
        for( int j = 0; j < N; ++j ){
            //#pragma omp parallel for schedule( static ) num_threads( num_threads )
            for( int k = 0; k < N; ++k ){
                int idx = j*N + k;
                std::cout << "f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
                std::cout << "df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
                
            }
        }*/
    /* double   fieldsum2, dfieldsum2;
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N; ++k ){
                int idx = j*N + k;
                // std::cout << "f[0]["<< idx <<"]" << f[i][idx] << std::endl;
                //  std::cout << "df[0]["<< idx <<"]" << df[i][idx] << std::endl;
                fieldsum2 +=  f[i][idx];
               dfieldsum2 += df[i][idx];
                // std::cout << "f[0]["<< idx <<"]" << f[i][idx] << std::endl;
                // std::cout << "df[0]["<< idx <<"]" << df[i][idx] << std::endl;
                
            }
        }
         for( int j = 0; j < dim; ++j ) fieldsum2 /= N;
        for( int j = 0; j < dim; ++j ) dfieldsum2 /= N;
        std::cout << "fieldsum2 = " << fieldsum2 << std::endl;
         std::cout << "dfieldsum2 = " << dfieldsum2 << std::endl;*/
        //Add zeromode
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N; ++k ){
                int idx = j*N + k;
               // std::cout << "f[0]["<< idx <<"]" << f[i][idx] << std::endl;
              //  std::cout << "df[0]["<< idx <<"]" << df[i][idx] << std::endl;
                f[i][idx] += initial_field_values[i];
                df[i][idx] += initial_field_derivs[i];
               // std::cout << "f[0]["<< idx <<"]" << f[i][idx] << std::endl;
               // std::cout << "df[0]["<< idx <<"]" << df[i][idx] << std::endl;
               
        }
    }
      /*
        std::random_device rnd;
        std::mt19937 mt( rnd() );
        std::uniform_real_distribution<> rand( 0, 1 );
        
        
        
        
        //  #pragma omp parallel for schedule( static ) num_threads ( num_threads )
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N; ++k ){
                int idx = j*N + k;
                f[i][idx] *= ( 1 + pow(10,-4)*rand(mt) );
              //   std::cout << "f[i][" << idx << "]" << f[i][idx]<< std::endl;
               // df[i][idx] = 0;//pow(10,-2)*rand(mt);
            }
        }*/
        
         std::cout << "after adding zeromode" << std::endl;
        std::cout << "f[0][4]" << f[i][4] << std::endl;
        std::cout << "df[0][4]" << df[i][4] << std::endl;
        std::cout << "f[0][10]" << f[i][10] << std::endl;
        std::cout << "df[0][10]" << df[i][10] << std::endl;
    delete [] fnyquist;
    delete [] fdnyquist;
        
#elif dim==3
    DFT_c2rD3( f[i] , fnyquist ); //transform from phase space to real space
    DFT_c2rD3( df[i], fdnyquist);
        //Add zeromode
         for( int j = 0; j < N; ++j ){
             for( int k = 0; k < N; ++k ){
                 for( int l = 0; l < N; ++l ){
                     int idx = (j*N + k)*N + l;
                     f[i][idx] += initial_field_values[i];
                     df[i][idx] += initial_field_derivs[i];
                 }
             }
         }
        
    
    delete [] fnyquist[0];
    delete [] fdnyquist[0];
    delete [] fnyquist;
    delete [] fdnyquist;
#endif
      
        delete[] mass_sq;
        delete[] initial_field_values;
        delete[] initial_field_derivs;
        
    }
    
}

void finalize( double**& f, double**& df )
{
	delete [] f[0];
	delete [] df[0];
    delete [] f;
	delete [] df;
}



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
    
 //   #pragma omp parallel for reduction(+:potential_energy) schedule( static ) num_threads ( num_threads )
	for( int j = 0; j < N; ++j ){
        #if dim == 1
            int idx = j;
             potential_energy += aV( f, 0, idx, a );
	    #elif dim == 2
            for( int k = 0; k < N; ++k )
			{
				int idx = j*N + k;
                potential_energy += aV( f, 0, idx, a);
            }
        #elif dim == 3
            for( int k = 0; k < N; ++k ){
                for( int l = 0; l < N; ++l ){
					int idx = ( j*N + k)*N + l;
					potential_energy += aV( f, 0, idx, a );
				}
            }
        #endif
	}
    for( int i = 0; i < dim; ++i ) potential_energy /= N;
    
	return potential_energy;
}


double Field::f_average( double* f, int i )
{
	_faverage[i] = 0;

	for( int j = 0; j < N; ++j ){
		switch( dim )
		{
			case 1:
				int idx = j;
				_faverage[i] += f[idx];
				break;
			case 2:
				for( int k = 0; k < N; ++k ){
					int idx = j*N + k;
					_faverage[i] += f[idx];
				}
				break;
			case 3:
				for( int k = 0; k < N; ++k ){
					for( int l = 0; l < N; ++l ){
						int idx = (j*N + k)*N + l;
						_faverage[i] += f[idx];
					}
				}
				break;
		}
	}
    for( int j = 0; j < dim; ++j ) _faverage[i] /= N;
	
    return _faverage[i];
}


double Field::f_variance( double* f, int i )
{
	_fvariance[i] = 0;

	for( int j = 0; j < N; ++j ){
		switch( dim ){
			case 1:
				int idx = j;
				_fvariance[i] += pow( f[idx] - _faverage[i], 2 );
				break;
			case 2:
				for( int k = 0; k < N; ++k ){
					int idx = j*N + k;
                    _fvariance[i] += pow( f[idx] - _faverage[i], 2 );
                
				}
				break;
			case 3:
				for( int k = 0; k < N; ++k ){
					for( int l = 0; l < N; ++l ){
						int idx = (j*N + k)*N + l;
						_fvariance[i] += pow( f[idx] - _faverage[i], 2 );
					}
				}
				break;
		}
	}
    for( int j = 0; j < dim; ++j ) _fvariance[i] /= N;
   
	
    return sqrt(_fvariance[i]);
}

double Field::df_average( double* df, int i )
{
    _dfaverage[i] = 0;
    
    for( int j = 0; j < N; ++j ){
        switch( dim )
        {
            case 1:
                int idx = j;
                _dfaverage[i] += df[idx];
                break;
            case 2:
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    _dfaverage[i] += df[idx];
                }
                break;
            case 3:
                for( int k = 0; k < N; ++k ){
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        _dfaverage[i] += df[idx];
                    }
                }
                break;
        }
    }
    for( int j = 0; j < dim; ++j ) _dfaverage[i] /= N;
    
    return _dfaverage[i];
}

double Field::df_variance( double* df, int i )
{
    _dfvariance[i] = 0;
    
    for( int j = 0; j < N; ++j ){
        switch( dim ){
            case 1:
                int idx = j;
                _dfvariance[i] += pow( df[idx] - _dfaverage[i], 2 );
                break;
            case 2:
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    _dfvariance[i] += pow( df[idx] - _dfaverage[i], 2 );
                   
                }
                break;
            case 3:
                for( int k = 0; k < N; ++k ){
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        _dfvariance[i] += pow( df[idx] - _dfaverage[i], 2 );
                    }
                }
                break;
        }
    }
    for( int j = 0; j < dim; ++j ) _dfvariance[i] /= N;
    
    return sqrt(_dfvariance[i]);
    
}
