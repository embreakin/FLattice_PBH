#include <cmath>
#include <random>
#include "lattice_field.hpp"
#include "utilities.hpp"
#include "equations.hpp"

double Fk_log_int_calc(int k_int , double**& lattice_var, int num_field)
{
    double F_k;
    
    switch( num_field )
    {
        case 0: //sigma field
        case 1: //psi field
        case 2:  //phi field
            //Fourier mode before discretization
            F_k =  sqrt(
                        pow(lattice_var[k_int-1][7+3*num_field],2) + pow(lattice_var[k_int-1][31+3*num_field],2)
                        +pow(lattice_var[k_int-1][8+3*num_field],2) + pow(lattice_var[k_int-1][32+3*num_field],2)
                        +pow(lattice_var[k_int-1][9+3*num_field],2) + pow(lattice_var[k_int-1][33+3*num_field],2)
                        );
        
            //            std::cout << "dF_k = " << dF_k << std::endl;
            //            std::cout << "F_k = " << F_k << std::endl;
            //            std::cout << "rescale_B = " << rescale_B << std::endl;
            //            std::cout << "omega = " << omega << std::endl;
            break;
        case 3: //gravitational potential
            //Fourier mode before discretization
            F_k =  sqrt(
                        pow(lattice_var[k_int-1][25],2) + pow(lattice_var[k_int-1][49],2)
                        +pow(lattice_var[k_int-1][26],2) + pow(lattice_var[k_int-1][50],2)
                        +pow(lattice_var[k_int-1][27],2) + pow(lattice_var[k_int-1][51],2)
                        );
            break;
        default:
            Logout( "Parameter 'num_field' must be 0 ~ 3. \n" );
            exit(1);
            
    }
    
    return log(F_k);
}

double dFk_log_int_calc(int k_int , double**& lattice_var, int num_field)
{
    double hubble_lattice_init, dF_k;
    
    //Hubble
    hubble_lattice_init = Fri(lattice_var[k_int-1][0],lattice_var[k_int-1][1],lattice_var[k_int-1][2],lattice_var[k_int-1][3],lattice_var[k_int-1][4],lattice_var[k_int-1][5],lattice_var[k_int-1][6]);
    
    switch( num_field )
    {
        case 0: //sigma field
        case 1: //psi field
        case 2:  //phi field
            //Fourier mode of the derivative before discretization
            dF_k = hubble_lattice_init*sqrt(
                                            pow(lattice_var[k_int-1][16+3*num_field],2) + pow(lattice_var[k_int-1][40+3*num_field],2)
                                            +pow(lattice_var[k_int-1][17+3*num_field],2) + pow(lattice_var[k_int-1][41+3*num_field],2)
                                            +pow(lattice_var[k_int-1][18+3*num_field],2) + pow(lattice_var[k_int-1][42+3*num_field],2)
                                            );
            break;
        case 3: //gravitational potential
            //Fourier mode of the derivative before discretization
            dF_k = hubble_lattice_init*sqrt(
                                            pow(lattice_var[k_int-1][28],2) + pow(lattice_var[k_int-1][52],2)
                                            +pow(lattice_var[k_int-1][29],2) + pow(lattice_var[k_int-1][53],2)
                                            +pow(lattice_var[k_int-1][30],2) + pow(lattice_var[k_int-1][54],2)
                                            );
            break;
        default:
            Logout( "Parameter 'num_field' must be 0 ~ 3. \n" );
            exit(1);
            
    }
    
    return log(dF_k);
}


double omega_calc(double distance, double**& lattice_var, int num_field)
{
    int l;
    double log_Fk, log_Fk2, log_dFk, log_dFk2, omega;

    l = floor(distance);
    if (l < N/2){
        
        log_Fk = Fk_log_int_calc(l, lattice_var, num_field);
        log_Fk2 = Fk_log_int_calc(l+1, lattice_var, num_field);
        log_Fk += (distance - l)*(log_Fk2 - log_Fk);
        
        log_dFk = dFk_log_int_calc(l, lattice_var, num_field);
        log_dFk2 = dFk_log_int_calc(l+1, lattice_var, num_field);
        log_dFk += (distance - l)*(log_dFk2 - log_dFk);
        
        omega = exp(log_dFk)/(exp(log_Fk)*rescale_B);
        
    }else{
        
        log_Fk = Fk_log_int_calc(N/2-1, lattice_var, num_field);
        log_Fk2 = Fk_log_int_calc(N/2, lattice_var, num_field);
        log_Fk += (distance - (N/2-1))*(log_Fk2 - log_Fk);
        
        log_dFk = dFk_log_int_calc(N/2-1, lattice_var, num_field);
        log_dFk2 = dFk_log_int_calc(N/2, lattice_var, num_field);
        log_dFk += (distance - (N/2-1))*(log_dFk2 - log_dFk);
        
        omega = exp(log_dFk)/(exp(log_Fk)*rescale_B);
    }
    
    return omega;
}

void initialize( double**& f, double**& df, Field* field, double radiation_pr, double**& lattice_var)
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
 
   
    
  
//  double* mass_sq = new double [num_fields];
  double* initial_field_values = new double [num_fields];
  double* initial_field_derivs = new double [num_fields];
  double radiation_var;
  double hubble_parameter;
 

    //lattice_var[knum_lattice][j] [j] specifies variables as follows:
    //    j=0-2  : zero modes of inflaton sigma, psi, and phi
    //    j=3-5  : log(a) derivatives of inflaton sigma', psi', phi'
    //    j=6    : energy density of radiation
    //    j=7-15 : mode functions of field perturbation: delta_{sigma,sigma}, delta_{sigma,psi}, delta{sigma,phi}, delta_{psi,sigma}, etc...
    //    j=16-24: log(a) derivatives of mode functions
    //    j=25-27: mode functions of gravitational potential perturbation: delta Phi_{sigma}, delta Phi_{psi}, delta Phi_{phi}
    //    j=28-30: log(a) derivatives of gravitational potential perturbation
    //    j=31-54: complex conjugate of [7]-[30]
 
    //Fields zeromode
    for (int i=0; i< num_fields; i++){

        int j = (num_fields -1) + i;

        initial_field_values[i] = 0.;
        initial_field_derivs[i] = 0.;
//
            if(i < num_fields - 1){ // Zeromode doesn't need to be calculated for gravitational potential and we set it to zero.
                for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){
//
    //    Logout(" lattice_var[%d][%d] = %2.5e \n", lattice_loop, i, lattice_var[lattice_loop][i]);
//
        initial_field_values[i] += lattice_var[lattice_loop][i];


        initial_field_derivs[i] += lattice_var[lattice_loop][j];

                }

        initial_field_values[i] /= (N/2);
        initial_field_derivs[i] /= (N/2);

            }
        
        Logout(" initial_field_values[%d] = %2.5e \n",i,initial_field_values[i]);
        Logout(" initial_field_derivs[%d] = %2.5e \n",i,initial_field_derivs[i]);

        }
    
    
    
    //Calculate Hubble parameter by averaging the hubble parameter for each wave number.
    for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){
    hubble_parameter +=
        Fri(lattice_var[lattice_loop][0],lattice_var[lattice_loop][1],lattice_var[lattice_loop][2],lattice_var[lattice_loop][3],lattice_var[lattice_loop][4],lattice_var[lattice_loop][5],lattice_var[lattice_loop][6]);
    }
    
    hubble_parameter/= (N/2);
    
    Logout("final hubble_parameter = %2.5e \n",hubble_parameter);
    
    
    //Rescaling initial fields and their derivatives to a lattice program variable
    //Only necessary for the scalar fields
    for (int i=0; i< num_fields - 1; i++){
        
        if(i==1){//psi
            initial_field_values[i] = FIXPSI - initial_field_values[i];
        }
        
        initial_field_values[i] *= rescale_A;
        
    initial_field_derivs[i] = (rescale_B/rescale_A)*initial_field_derivs[i] + hubble_parameter*initial_field_values[i];
    
    
        
        Logout("pr initial_field_values[%d] = %2.5e \n",i,initial_field_values[i]);
        Logout("pr initial_field_derivs[%d] = %2.5e \n",i,initial_field_derivs[i]);
    }
    
    Logout("FIXPSI = %2.5e \n",FIXPSI);
    
    //Rescaling hubble parameter to a lattice program variable
    //It doesn't change since a = a_{0}\bar{a}
    Hinitial_pr = hubble_parameter;
    

    //Radiation zeromode
     for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++)
     {

     radiation_var += lattice_var[lattice_loop][6];

     }

    radiation_var /= (N/2);

    Logout("radiation_var = %2.5e \n", radiation_var);
    
    //Rescaling radiation to a lattice program variable
    radiation_pr = pow((rescale_A/rescale_B),2)*radiation_var;
    
    

   
    //Initializing perturbations
    double p2;
    double dp2=pw2(2*M_PI/L);
    double omega;
    double distance;
    
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

     // std::cout << "mass_sq[ " << i <<  "] = " << mass_sq[i] << std::endl;


  #if  dim==1
        
        //k=N/2
        p2 = dp2*pw2(N/2);
        
        //Omega
        omega = omega_calc( N/2 , lattice_var, i);

        //Set mode for lattice simulation
        set_mode(p2, omega, &f[i][1], &df[i][1], 1);

        //Loop over gridpoints.
        for(int k = 1; k < N/2; k++ ){
            pz = k;
            p2 = dp2*pw2(pz);
            //Omega
            omega = omega_calc( k , lattice_var, i);

            //Set mode for lattice simulation
             set_mode(p2, omega, &f[i][2*k], &df[i][2*k], 0);
        }
        
        //k=0
        f[i][0] = 0.;
        df[i][0] = 0.;
        
#elif dim==2
        for(int j = 0; j < N; j++ ){
            py = (j <= N/2 ? j : j-N);

            for(int k = 1; k < N/2; k++ ){
                pz = k;
                p2=dp2*(pw2(py)+pw2(pz));
                distance = sqrt(pw2(py)+pw2(pz));
               omega = omega_calc(distance, lattice_var, i);
                //Set mode for lattice simulation
                set_mode(p2, omega, &f[i][j*N+2*k], &df[i][j*N+2*k], 0);

            }

            if(j > N/2)

            {
                jconj = N-j;

                //k=0
                p2 = dp2*pw2(py);
                omega = omega_calc( jconj , lattice_var, i);
                set_mode(p2, omega, &f[i][j*N], &df[i][j*N], 0);
                f[i][jconj*N] = f[i][j*N];
                f[i][jconj*N+1] = -f[i][j*N+1];
                df[i][jconj*N] = df[i][j*N];
                df[i][jconj*N+1] = -df[i][j*N+1];
                //k=N/2
                distance = sqrt(pw2(py)+pw2(N/2));
                omega = omega_calc(distance, lattice_var, i);
                p2 = dp2*(pw2(py)+pw2(N/2));
                set_mode(p2, omega, &fnyquist[2*j], &fdnyquist[2*j], 0);
                fnyquist[2*jconj] = fnyquist[2*j];
                fnyquist[2*jconj+1] = -fnyquist[2*j+1];
                fdnyquist[2*jconj] = fdnyquist[2*j];
                fdnyquist[2*jconj+1] = -fdnyquist[2*j+1];

            }
            // The 4 "corners" of the lattice are set to real values
            else if(j == 0 || j == N/2){
                
                p2 = dp2*pw2(py); //k = 0
                if(p2 > 0.){
                    omega = omega_calc( N/2 , lattice_var, i);
                    set_mode(p2, omega, &f[i][j*N], &df[i][j*N], 1);
                
                }
                
                p2 = dp2*(pw2(py)+pw2(N/2));  //k = N/2
                distance = sqrt(pw2(py)+pw2(N/2));
                omega = omega_calc(distance, lattice_var, i);
                set_mode(p2, omega, &fnyquist[2*j], &fdnyquist[2*j], 1);
                
            }
        } //k = 0, j = 0
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
                    distance = sqrt(pw2(px)+pw2(py)+pw2(pz));
                    omega = omega_calc(distance, lattice_var, i);

                    set_mode(p2,omega, &f[i][(j*N + k)*N + 2*l], &df[i][(j*N + k)*N + 2*l], 0);
                }

                if(k > N/2 || (j > N/2 && (k == 0 || k == N/2))){
                    kconj = (k == 0 ? 0 : N-k);
                    //l=0
                     p2 = dp2*(pw2(px)+pw2(py));
                    distance = sqrt(pw2(px)+pw2(py));
                    omega = omega_calc(distance, lattice_var, i);
                    set_mode(p2,omega, &f[i][(j*N + k)*N], &df[i][(j*N + k)*N], 0);
                    f[i][(jconj*N + kconj)*N] = f[i][(j*N + k)*N];
                    f[i][(jconj*N + kconj)*N+1] = -f[i][(j*N + k)*N+1];
                    df[i][(jconj*N + kconj)*N] = df[i][(j*N + k)*N];
                    df[i][(jconj*N + kconj)*N+1] = -df[i][(j*N + k)*N+1];
                    //l=N/2
                    p2 = dp2*(pw2(px)+pw2(py)+pw2(N/2));
                    distance = sqrt(pw2(px)+pw2(py)+pw2(N/2));
                    omega = omega_calc(distance, lattice_var, i);
                    set_mode(p2, omega, &fnyquist[j][2*k], &fdnyquist[j][2*k], 0);
                    fnyquist[jconj][2*kconj] = fnyquist[j][2*k];
                    fnyquist[jconj][2*kconj+1] = -fnyquist[j][2*k+1];
                    fdnyquist[jconj][2*kconj] = fdnyquist[j][2*k];
                    fdnyquist[jconj][2*kconj+1] = -fdnyquist[j][2*k+1];
                    }
                else if((j == 0 || j == N/2) && (k == 0 || k == N/2))
                {
                    p2 = dp2*(pw2(px)+pw2(py)); //l=0
                    if(p2 > 0.){
                    distance = sqrt(pw2(px)+pw2(py));
                    omega = omega_calc(distance, lattice_var, i);
                    set_mode(p2, omega, &f[i][(j*N + k)*N], &df[i][(j*N + k)*N], 1);
                    }
                    
                    p2=dp2*(pw2(px)+pw2(py)+pw2(N/2)); //l=N/2
                    distance = sqrt(pw2(px)+pw2(py)+pw2(N/2));
                    omega = omega_calc(distance, lattice_var, i);
                     set_mode(p2, omega, &fnyquist[j][2*k], &fdnyquist[j][2*k], 1);
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

#pragma omp parallel for simd schedule(static) num_threads(num_threads)
       for( int j = 0; j < N; ++j ){
           int idx = j;
           f[i][idx] += initial_field_values[i];
           df[i][idx] += initial_field_derivs[i];
       }

#elif dim==2

         std::cout << "before fourier transform" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
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
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        
//        for( int j = 0; j < N; ++j ){
////            #pragma omp parallel for schedule( static ) num_threads( num_threads )
//            for( int k = 0; k < N; ++k ){
//                int idx = j*N + k;
//                std::cout << "f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
//                std::cout << "df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
//
//            }
//        }
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

#pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
             #pragma omp simd
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
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
    
        
#elif dim==3

        std::cout << "before fourier transform" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;

    DFT_c2rD3( f[i] , fnyquist ); //transform from phase space to real space
    DFT_c2rD3( df[i], fdnyquist);

        std::cout << "after fourier transform" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        //Add zeromode
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
         for( int j = 0; j < N; ++j ){
             for( int k = 0; k < N; ++k ){
                  #pragma omp simd
                 for( int l = 0; l < N; ++l ){
                     int idx = (j*N + k)*N + l;
                     f[i][idx] += initial_field_values[i];
                     df[i][idx] += initial_field_derivs[i];
                 }
             }
         }

        std::cout << "after adding zeromode" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;

    
#endif

    } // for( int i = 0; i < num_fields; ++i ){}

    #if  dim==1
    
    #elif dim==2
    delete [] fnyquist;
    delete [] fdnyquist;
    
    #elif dim==3
    delete [] fnyquist[0];
    delete [] fdnyquist[0];
    delete [] fnyquist;
    delete [] fdnyquist;
    
    #endif
    
    //delete[] mass_sq;
    delete[] initial_field_values;
    delete[] initial_field_derivs;
}

void finalize( double**& f, double**& df )
{
	delete [] f[0];
	delete [] df[0];
    delete [] f;
	delete [] df;
}


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


double Field::f_average( double* f, int i )
{
	//Variable needs to be declared for OpenMP reduction directive
    double faverage = 0;

#if   dim == 1
#pragma omp parallel for simd reduction(+:faverage) schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for reduction(+:faverage) schedule(static) num_threads(num_threads)
#endif
	for( int j = 0; j < N; ++j ){
		switch( dim )
		{
			case 1:
				int idx = j;
				faverage += f[idx];
				break;
			case 2:
            #pragma omp simd reduction(+:faverage)
				for( int k = 0; k < N; ++k ){
					int idx = j*N + k;
					faverage += f[idx];
				}
				break;
			case 3:
				for( int k = 0; k < N; ++k ){
            #pragma omp simd reduction(+:faverage)
					for( int l = 0; l < N; ++l ){
						int idx = (j*N + k)*N + l;
						faverage += f[idx];
					}
				}
				break;
		}
	}
    for( int j = 0; j < dim; ++j ) faverage /= N;
    
    //substitute the obtained variable to member variable
    _faverage[i] = faverage;
	
    return _faverage[i];
}


double Field::f_variance( double* f, int i )
{
    //Variable needs to be declared for OpenMP reduction directive
    double fvariance = 0;
    
#if   dim == 1
#pragma omp parallel for simd reduction(+:fvariance) schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for reduction(+:fvariance) schedule(static) num_threads(num_threads)
#endif

	for( int j = 0; j < N; ++j ){
		switch( dim ){
			case 1:
				int idx = j;
				fvariance += pow( f[idx] - _faverage[i], 2 );
				break;
			case 2:
                #pragma omp simd reduction(+:fvariance)
				for( int k = 0; k < N; ++k ){
					int idx = j*N + k;
                    fvariance += pow( f[idx] - _faverage[i], 2 );
                
				}
				break;
			case 3:
				for( int k = 0; k < N; ++k ){
                    #pragma omp simd reduction(+:fvariance)
					for( int l = 0; l < N; ++l ){
						int idx = (j*N + k)*N + l;
						fvariance += pow( f[idx] - _faverage[i], 2 );
					}
				}
				break;
		}
	}
    for( int j = 0; j < dim; ++j ) fvariance /= N;
    
    //substitute the obtained variable to member variable
    _fvariance[i] = fvariance;
   
    return sqrt(_fvariance[i]);
}

double Field::df_average( double* df, int i )
{
    //Variable needs to be declared for OpenMP reduction directive
    double dfaverage = 0;
    
#if   dim == 1
#pragma omp parallel for simd reduction(+:dfaverage) schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for reduction(+:dfaverage) schedule(static) num_threads(num_threads)
#endif
    for( int j = 0; j < N; ++j ){
        switch( dim )
        {
            case 1:
                int idx = j;
                dfaverage += df[idx];
                break;
            case 2:
#pragma omp simd reduction(+:dfaverage)
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    dfaverage += df[idx];
                }
                break;
            case 3:
                for( int k = 0; k < N; ++k ){
#pragma omp simd reduction(+:dfaverage)
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        dfaverage += df[idx];
                    }
                }
                break;
        }
    }
    for( int j = 0; j < dim; ++j ) dfaverage /= N;
    
    //substitute the obtained variable to member variable
    _dfaverage[i] = dfaverage;
    
    return _dfaverage[i];
}


double Field::df_variance( double* df, int i )
{
    //Variable needs to be declared for OpenMP reduction directive
    double dfvariance = 0;
    
#if   dim == 1
#pragma omp parallel for simd reduction(+:dfvariance) schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for reduction(+:dfvariance) schedule(static) num_threads(num_threads)
#endif
    
    for( int j = 0; j < N; ++j ){
        switch( dim ){
            case 1:
                int idx = j;
                dfvariance += pow( df[idx] - _dfaverage[i], 2 );
                break;
            case 2:
#pragma omp simd reduction(+:dfvariance)
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    dfvariance += pow( df[idx] - _dfaverage[i], 2 );
                    
                }
                break;
            case 3:
                for( int k = 0; k < N; ++k ){
#pragma omp simd reduction(+:dfvariance)
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        dfvariance += pow( df[idx] - _dfaverage[i], 2 );
                    }
                }
                break;
        }
    }
    for( int j = 0; j < dim; ++j ) dfvariance /= N;
    
    //substitute the obtained variable to member variable
    _dfvariance[i] = dfvariance;
    
    return sqrt(_dfvariance[i]);
}

#pragma omp declare simd
double Field::V_lattice   ( double** f, int idx, double a )  {
    // std::cout << "a = " << a << std::endl;
    for(fld=0;fld< num_fields - 1;fld++)
    {
        if(fld==1){
            f_MPl[fld] = FIXPSI - f[fld][idx]/(rescale_A*a);
        }else{
            f_MPl[fld] = f[fld][idx]/(rescale_A*a);
        }
    }
    
    return pow(a,4)*pow(rescale_A/rescale_B,2)*V(f_MPl[0],f_MPl[1],f_MPl[2]); }

#pragma omp declare simd
double Field::dV_lattice ( double** f, int i, int idx, double a )  {
    
    for(fld=0;fld< num_fields - 1;fld++)
    {
        if(fld==1){
        f_MPl[fld] = FIXPSI - f[fld][idx]/(rescale_A*a);
        }else{
        f_MPl[fld] = f[fld][idx]/(rescale_A*a);
        }
    }
    
    
    switch (i){
        case 0: return pow(a,3)*rescale_A*V_1(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //sigma
        case 1: return -pow(a,3)*rescale_A*V_2(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //psi
        case 2: return pow(a,3)*rescale_A*V_3(f_MPl[0],f_MPl[1],f_MPl[2])/(rescale_B*rescale_B); //phi
        default:  Logout( "Parameter 'i' in dV_lattice must be 0 ~ 2. \n" );
                exit(1);
    }
    
}


double Field::mass( int i,  double a  ){
   
    switch (i){
        case 0: return pow(a*exp(OSCSTART),2)*V_11(0,FIXPSI,0)/(rescale_B*rescale_B); //sigma
        case 1: return pow(a*exp(OSCSTART),2)*V_22(0,FIXPSI,0)/(rescale_B*rescale_B); //psi
        case 2: return pow(a*exp(OSCSTART),2)*V_33(0,FIXPSI,0)/(rescale_B*rescale_B); //phi
        default:  Logout( "Parameter 'i' in dV_lattice must be 0 ~ 2. \n" );
                exit(1);
    }
    
}

