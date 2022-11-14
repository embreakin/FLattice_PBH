#include <cmath>
#include "parameters.hpp"
#include "utilities.hpp"
#include "lattice_evol.hpp"

/////////////////////////////////////////////////////
////////////////////CONSTRUCTER//////////////////////
/////////////////////////////////////////////////////

LeapFrog::LeapFrog( Field* field, double** f, double** df, double& rad_pr )
{
    
    //Create arrays for evolution fields and to temporarily save fields
    f_tilde = new double* [num_fields - 1];
    df_tilde = new double* [num_fields - 1];
    fdf_save = new double* [num_fields - 1];
    
    switch( dim )
    {
        case 1:
            f_tilde[0] = new double [(num_fields - 1)*N];
            df_tilde[0] = new double [(num_fields - 1)*N];
            fdf_save[0] = new double [(num_fields - 1)*N];
            for( int i = 0; i < num_fields - 1; ++i )
            {
                f_tilde[i] = f_tilde[0] + i*N;
                df_tilde[i] = df_tilde[0] + i*N;
                fdf_save[i] = fdf_save[0] + i*N;
            }
            break;
        case 2:
            f_tilde[0] = new double [(num_fields - 1)*N*N];
            df_tilde[0] = new double [(num_fields - 1)*N*N];
            fdf_save[0] = new double [(num_fields - 1)*N*N];
            for( int i = 0; i < num_fields - 1; ++i )
            {
                f_tilde[i] = f_tilde[0] + i*N*N;
                df_tilde[i] = df_tilde[0] + i*N*N;
                fdf_save[i] = fdf_save[0] + i*N*N;
            }
            break;
        case 3:
            f_tilde[0] = new double [(num_fields - 1)*N*N*N];
            df_tilde[0] = new double [(num_fields - 1)*N*N*N];
            fdf_save[0] = new double [(num_fields - 1)*N*N*N];
            for( int i = 0; i < num_fields - 1; ++i )
            {
                f_tilde[i] = f_tilde[0] + i*N*N*N;
                df_tilde[i] = df_tilde[0] + i*N*N*N;
                fdf_save[i] = fdf_save[0] + i*N*N*N;
            }
            break;
        default:
            std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
            exit(1);
    }
    
    
    Logout( "rad_pr = %2.5e \n", rad_pr);
    Logout( "df[3][0] = %2.5e \n", df[3][0]);
    
    // Set Decay Rates
    Gamma_pr[0]=Gamma1/rescale_B;
    Gamma_pr[1]=Gamma2/rescale_B;
    Gamma_pr[2]=Gamma3/rescale_B;
    
    
    //Rescale field program variables to field variables necessary for leapfrog evolution
    //We only rescale the scalar fields.
    
    
    fields_convert( f, f_tilde, 0);
    
//    for( int i = 0; i < num_fields-1; ++i ){
//        for( int j = 0; j < N; ++j ){
//            int idx = j;
//            std::cout << " f_tilde[" << i << "][" << idx << "] = " << f_tilde[i][idx] << std::endl;
//            std::cout << " df_tilde[" << i << "][" << idx << "] = " << df_tilde[i][idx] << std::endl;
//        }
//    }
    
    fields_deriv_convert( f, df, f_tilde, df_tilde, 0);
    
//    for( int i = 0; i < num_fields-1; ++i ){
//        for( int j = 0; j < N; ++j ){
//            int idx = j;
//            std::cout << " f_tilde[" << i << "][" << idx << "] = " << f_tilde[i][idx] << std::endl;
//            std::cout << " df_tilde[" << i << "][" << idx << "] = " << df_tilde[i][idx] << std::endl;
//        }
//    }
    
    switch( precision )
    {
        case 2:
        case 4:
            break;
        default:
            Logout( "LeapFrog precision must be 2 or 4. \n" );
            exit(1);
    }
    
    switch( expansion )
    {
        case 0:
        case 1:  // self-consistent evolution
            _a  = 1;
            _da = hubble_init;
            //In the function below, a'' is calculated by converting f_tilde to f.
           evol_scale_dderivs( field, f, f_tilde, rad_pr , 0.); // _dda( Field* field, double** f, double ** f_tilde, double rho_r, double t, double h)
           
            break;
        case 2: // radiation dominant universe
            sfexponent = 1;
            sfbase = 1;
            _a  = 1;
            _da = hubble_init;
            _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
            break;
        case 3: // matter dominant universe
            sfexponent = 2;
            sfbase = 1;
            _a  = 1;
            _da = hubble_init;
            _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
           break;
        default:
            Logout( "Parameter 'expansion' must be 0 ~ 3 when you use class 'LeapFrog'. \n" );
            exit(1);
    }
}


/////////////////////////////////////////////////////
//////////////////PRIVATE FUNCTIONS//////////////////
/////////////////////////////////////////////////////


void LeapFrog::evol_fields( double** f_tilde_from, double** f_tilde_to, double** df_tilde, double h )
{	
    for( int i = 0; i < num_fields-1; ++i ){
    #if   dim == 1
    #pragma omp parallel for simd schedule(static) num_threads(num_threads)
    #elif dim >= 2
    #pragma omp parallel for schedule( static ) num_threads( num_threads )
    #endif
        for( int j = 0; j < N; ++j ){
            #if dim == 1
                int idx = j;
        
                f_tilde_to[i][idx] = f_tilde_from[i][idx] + df_tilde[i][idx]*h*dt;
//            std::cout << " f_tilde[" << i << "][" << idx << "] = " << f_tilde[i][idx] << std::endl;
            #elif dim == 2
            #pragma omp simd
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                     f_tilde_to[i][idx] = f_tilde_from[i][idx] + df_tilde[i][idx]*h*dt;
                }
            #elif dim == 3
                for( int k = 0; k < N; ++k ){
                    #pragma omp simd
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                       f_tilde_to[i][idx] = f_tilde_from[i][idx] + df_tilde[i][idx]*h*dt;
                    }
                }
            #endif
        }
    }
}


void LeapFrog::evol_field_derivs_expansion(double** df_tilde_from, double** df_tilde_to, double** f, double** f_tilde, Field* field,  double h )
{
    
    
    for( int i = 0; i < num_fields - 1; ++i ){
        
        #if   dim == 1
//        #pragma omp parallel for simd schedule(static) num_threads(num_threads)
        #elif dim >= 2
   //     #pragma omp parallel for schedule( static ) num_threads( num_threads )
        #endif
                    for( int j = 0; j < N; ++j ){
        #if dim == 1
                        int idx = j;


                        df_tilde_to[i][idx] = df_tilde_from[i][idx] +
                        (
                         
                         //Laplacian
                         (
                          field->laplacian(f_tilde[i], j)
                          )/pw2(exp(OSCSTART))
                         +//a terms
                         (
                          //- field->mass(i,_a)
                          + _dda/_a + _da*Gamma_pr[i]
                          + pow(Gamma_pr[i]*_a,2)/4
                          ) *f_tilde[i][idx]
                         //dV term
                         - field->dV_lattice(f, i, idx ,_a)*exp(Gamma_pr[i]*_sfint/2)
                         
                         )* h*dt;
                        
        #elif dim == 2
                          //  #pragma omp simd
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;

//                        if(i==2 && idx==0){
//        Logout( "laplacian = %2.5e \n", field->laplacian(f_tilde[i], j, k));
//       Logout( "pw2(exp(OSCSTART)) = %2.5e \n", pw2(exp(OSCSTART)));
//                            Logout( "field->laplacian(f_tilde[i], j)/pw2(exp(OSCSTART)) = %2.5e \n", field->laplacian(f_tilde[i], j, k)
//                                   /pw2(exp(OSCSTART)));
//    Logout( " _dda = %2.5e \n", _dda);
//    Logout( "term2 = %2.5e \n", _dda/_a);
//    Logout( " _da = %2.5e \n", _da);
//    Logout( " _a = %2.5e \n", _a);
//    Logout( "term3 = %2.5e \n",  _da*Gamma_pr[i]);
//    Logout( "term4 = %2.5e \n", pow(Gamma_pr[i],2)/4);
//    Logout( "i = %d, idx = %d \n", i, idx);
//    Logout( "f_tilde[i][idx] = %2.5e \n", f_tilde[i][idx]);
//    Logout( "field->dV_lattice(f, i, idx ,_a) = %2.5e \n", field->dV_lattice(f, i, idx ,_a));
//    Logout( "exp(Gamma_pr[i]*_sfint/2) = %2.5e \n", exp(Gamma_pr[i]*_sfint/2));
//    Logout( "- field->dV_lattice(f, i, idx ,_a)*exp(Gamma_pr[i]*t/2) = %2.5e \n", - field->dV_lattice(f, i, idx ,_a)*exp(Gamma_pr[i]*_sfint/2));
//        Logout( " efolds = %2.5e \n",log(exp(OSCSTART)*_a));
//                            Logout( " addition = %2.5e \n\n",  (
//
//                                                                (
//                                                                 field->laplacian(f_tilde[i], j, k)
//                                                                 )/pw2(exp(OSCSTART))
//                                                                +
//                                                              (
//                                                               //- field->mass(i,_a)
//                                                               + _dda/_a + _da*Gamma_pr[i]
//                                                               + pow(Gamma_pr[i]*_a,2)/4
//                                                               ) *f_tilde[i][idx]
//
//                                                              - field->dV_lattice(f, i, idx ,_a)*exp(Gamma_pr[i]*_sfint/2)
//
//                                                              )* h*dt);
////                            exit(1);
//                        }
//
                        
                            df_tilde_to[i][idx] = df_tilde_from[i][idx] +
                            (
                             
                             //Laplacian
                             (
                              field->laplacian(f_tilde[i], j, k)
                              )/pw2(exp(OSCSTART))
                             +//a terms
                             (
                              //- field->mass(i,_a)
                              + _dda/_a + _da*Gamma_pr[i]
                              + pow(Gamma_pr[i]*_a,2)/4
                              ) *f_tilde[i][idx]
                             //dV term
                             - field->dV_lattice(f, i, idx ,_a)*exp(Gamma_pr[i]*_sfint/2)
                             
                             )* h*dt;
                            
                        }
                        
        #elif dim == 3
                        
                        for( int k = 0; k < N; ++k ){
        //#pragma omp simd
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;
                                
                                df_tilde_to[i][idx] = df_tilde_from[i][idx] +
                                (
                                 //Laplacian
                                 (
                                  field->laplacian(f_tilde[i], j, k, l)
                                  )/pw2(exp(OSCSTART))
                                 +//a terms
                                 (
                                  //- field->mass(i,_a)
                                  + _dda/_a + _da*Gamma_pr[i]
                                  + pow(Gamma_pr[i]*_a,2)/4
                                  ) *f_tilde[i][idx]
                                 //dV term
                                 - field->dV_lattice(f, i, idx ,_a)*exp(Gamma_pr[i]*_sfint/2)
                                 
                                 )* h*dt;
                                
                            }
                        }
        #endif
                    }
        
    }
    
}



void LeapFrog::evol_scale_dderivs( Field* field, double** f, double ** f_tilde, double& rho_r, double h)
{
    //rescale fields from evolution variables to program variables
  
    fields_convert( f, f_tilde, 1);
    
//      for( int i = 0; i < num_fields - 1; ++i ){
//    for( int j = 0; j < N; ++j ){
//        for( int k = 0; k < N; ++k ){
//            int idx = j*N + k;
//            std::cout << " f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
//           // std::cout << " df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
//        }
//    }
//    }
   
//   std::cout << "pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2 ) = " << pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2 ) << "f["<< jp2 << "] = " << f[jp2] << "f["<< jp1 << "] = " << f[jp1] << "f["<< jm1 << "] = " << f[jm1]<< "f["<< jm2 << "] = " << f[jm2] << std::endl;
//
   if(h==0){
       
       _dda = -pw2(_da)/_a
       +(1/(_a*rescale_A*rescale_A))*
       (
        2.*(field->gradient_energy(f[0])
            +field->gradient_energy(f[1])
            +field->gradient_energy(f[2])
            )/(3.*pw2(exp(OSCSTART)))
        + field->potential_energy(f, _a)
        + rho_r/3
        );
       
    }else{
    
    double C = 0;
    C += 2*(field->gradient_energy(f[0])
            +field->gradient_energy(f[1])
            +field->gradient_energy(f[2])
            )/(3.*pw2(exp(OSCSTART)))
        + field->potential_energy( f, _a )
        +rho_r/3;
    
//       std::cout << "field->gradient_energy(f[0]) = " << field->gradient_energy(f[0]) << ",field->gradient_energy(f[1]) = " << field->gradient_energy(f[1]) << ", field->gradient_energy(f[2]) = " << field->gradient_energy(f[2]) << std::endl;
//        std::cout << " field->potential_energy( f, _a ) = " << field->potential_energy( f, _a ) << " rho_r = \n" << rho_r << std::endl;
//std::cout << " pow(h*dt,2.)*C/pow(rescale_A*_a,2) = " << pow(h*dt,2.)*C/pow(rescale_A*_a,2) << " 2*h*dt*_da/_a = \n" << 2*h*dt*_da/_a << std::endl;
//        std::cout << " sqrt( - ) = \n" << sqrt( 1 + 2*h*dt*_da/_a+ pow(h*dt,2.)*C/pow(rescale_A*_a,2)  ) << " _a/(h*dt)*(1-sqrt(-)) = " << _a/(h*dt) *
//        ( 1 -
//         sqrt( 1 + 2*h*dt*_da/_a
//              + pow(h*dt,2.)*C/pow(rescale_A*_a,2)  )
//         ) << std::endl;
//
//        std::cout << " _da + _a/(h*dt)*(1-sqrt(-)) = " <<
//        _da +
//        _a/(h*dt) *
//        ( 1 -
//         sqrt( 1 + 2*h*dt*_da/_a
//              + pow(h*dt,2.)*C/pow(rescale_A*_a,2)  )
//         ) << std::endl;
//
//        std::cout << " - 2./(h*dt) = " <<
//        - 2./(h*dt) << std::endl;
//
    _dda = - 2./(h*dt) * (
                          _da +
                          _a/(h*dt) *
                          ( 1 -
                           sqrt( 1 + 2*h*dt*_da/_a
                                + pow(h*dt,2.)*C/pow(rescale_A*_a,2)  )
                           )
                          );
       //std::cout << "1-1: _a = "<< _a << ", _da = " << _da << ", _dda = " << _dda << ", C = " << C << std::endl;
        _da += .5*_dda*h*dt;
      // std::cout << "1-2: _a = "<< _a << ", _da = " << _da << ", _dda = " << _dda << ", C = " << C << std::endl;
       
    }
    
    
}

void LeapFrog::evol_radiation(Field* field, double** f_tilde, double** df_tilde, double& rho, double h)
{
    double drho = 0;
   
        for( int i = 0; i < num_fields - 1; ++i ){
            
         drho += _a*Gamma_pr[i]*exp(-Gamma_pr[i]*_sfint)*pow(field->average(df_tilde[i], i)-(_a*Gamma_pr[i]/2 + _da/_a)*field->average(f_tilde[i], i),2);
            
//        std::cout << "f_tilde[" << i << "][0] = " << f_tilde[i][0] << std::endl;
//         std::cout << "df_tilde[" << i << "][0] = " << df_tilde[i][0] << std::endl;
//             std::cout << "_a*Gamma_pr[" << i << "]*exp(-Gamma_pr[" << i << "]*_sfint) = " << _a*Gamma_pr[i]*exp(-Gamma_pr[i]*_sfint) << std::endl;
//             std::cout << "_a*(Gamma_pr[i]/2 + _da/_a)*field->average(f_tilde[i], i) = " << (_a*Gamma_pr[i]/2 + _da/_a)*field->average(f_tilde[i], i) << std::endl;
//            std::cout << "pow(field->average(df_tilde[i], i)-(_a*Gamma_pr[i]/2 + _da/_a)*field->average(f_tilde[i], i),2) = " << pow(field->average(df_tilde[i], i)-(_a*Gamma_pr[i]/2 + _da/_a)*field->average(f_tilde[i], i),2) << std::endl;
//            std::cout <<  std::endl;
        }
    
//   std::cout << "drho*h*dt = " << drho*h*dt << std::endl;
//    std::cout << "rho = " << rho << std::endl;
    
    rho = rho + drho*h*dt;
    
//    std::cout << "rho = " << rho << std::endl;
    
//    exit(1);
}

void LeapFrog::evol_gravpot( double** f, double** df, double h )
{
    
    #if   dim == 1
    #pragma omp parallel for simd schedule(static) num_threads(num_threads)
    #elif dim >= 2
    #pragma omp parallel for schedule( static ) num_threads( num_threads )
    #endif
            for( int j = 0; j < N; ++j ){
    #if dim == 1
                int idx = j;
                f[3][idx] += df[3][idx] * h*dt;
//                 if(idx==3){
//                Logout( "f[3][idx] = %2.5e \n", f[3][idx]);
//                Logout( "df[3][idx] = %2.5e \n", df[3][idx] );
//                 }
    #elif dim == 2
    #pragma omp simd
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f[3][idx] += df[3][idx] * h*dt;
                }
    #elif dim == 3
                for( int k = 0; k < N; ++k ){
    #pragma omp simd
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f[3][idx] += df[3][idx] * h*dt;
                    }
                }
    #endif
            }
    
}


void LeapFrog::evol_gravpot_derivs_expansion( double** f, double** df, double** f_tilde, double** df_tilde, Field* field, double h )
{
    //rescale fields from evolution variables to program variables
    fields_convert( f, f_tilde, 1);
    fields_deriv_convert( f, df, f_tilde, df_tilde, 1);
    
    #if   dim == 1
//    #pragma omp parallel for simd schedule(static) num_threads(num_threads)
    #elif dim >= 2
    #pragma omp parallel for schedule( static ) num_threads( num_threads )
    #endif
            for( int j = 0; j < N; ++j ){
    #if dim == 1
                int idx = j;
                
//                                        if(idx==3){
//                    Logout( "h*dt = %2.5e \n", h*dt);
//                    Logout( "_dda*f[3][idx]/_a = %2.5e \n", _dda*f[3][idx]/_a);
//                    Logout( " -(_a + 2*f[3][idx])*(_dda/_a -pow(_da,2)/pow(_a,2)) = %2.5e \n", -(_a + 2*f[3][idx])*(_dda/_a -pow(_da,2)/pow(_a,2)));
//
//                    Logout( "- field->laplacian(f[3], j)/(3*pw2(exp(OSCSTART))) = %2.5e \n", - field->laplacian(f[3], j)/(3*pw2(exp(OSCSTART))));
//
//                                            Logout( " term4 = %2.5e \n",-(rho_rad -field->V_lattice(f, idx ,_a))/(3*rescale_A*rescale_A*_a));
//                                            Logout( "term5 = %2.5e \n",   - ( _a + 2*f[3][idx] )*(
//                                                                                                  pow(df[0][idx]/_a - _da*f[0][idx]/pow(_a,2),2)
//                                                                                                  +
//                                                                                                  pow(df[1][idx]/_a - _da*f[1][idx]/pow(_a,2),2)
//                                                                                                  +
//                                                                                                  pow(df[2][idx]/_a - _da*f[2][idx]/pow(_a,2),2)
//                                                                                                  )
//                                                   /(3*rescale_A*rescale_A));
//                                            Logout( "term6 = %2.5e \n",  - ( _a - 2*f[3][idx] )*2*(field->gradient_energy_eachpoint(f , 0, idx)
//                                                                                                   +
//                                                                                                   field->gradient_energy_eachpoint(f , 1, idx)
//                                                                                                   +
//                                                                                                   field->gradient_energy_eachpoint(f , 2, idx)
//
//                                                                                                   )/(3*rescale_A*rescale_A*pw2(exp(OSCSTART)*_a)));
//
//                                            Logout( " addition = %2.5e \n\n",   (
//                                                                                 _dda*f[3][idx]/_a
//                                                                                 -(_a + 2*f[3][idx])*(_dda/_a -pow(_da,2)/pow(_a,2))
//                                                                                 - field->laplacian(f[3], j)/(3*pw2(exp(OSCSTART)))
//                                                                                 -
//                                                                                 (rho_rad -field->V_lattice(f, idx ,_a))/(3*rescale_A*rescale_A*_a)
//                                                                                 - ( _a + 2*f[3][idx] )*(
//                                                                                                         pow(df[0][idx]/_a - _da*f[0][idx]/pow(_a,2),2)
//                                                                                                         +
//                                                                                                         pow(df[1][idx]/_a - _da*f[1][idx]/pow(_a,2),2)
//                                                                                                         +
//                                                                                                         pow(df[2][idx]/_a - _da*f[2][idx]/pow(_a,2),2)
//                                                                                                         )
//                                                                                 /(3*rescale_A*rescale_A)
//                                                                                 - ( _a - 2*f[3][idx] )*2*(field->gradient_energy_eachpoint(f , 0, idx)
//                                                                                                           +
//                                                                                                           field->gradient_energy_eachpoint(f , 1, idx)
//                                                                                                           +
//                                                                                                           field->gradient_energy_eachpoint(f , 2, idx)
//
//                                                                                                           )/(3*rescale_A*rescale_A*pw2(exp(OSCSTART)*_a))
//                                                                                 ) * h*dt);
//
//                                        }
//
                
                
                df[3][idx] += (
                              2*(pow(_da,2)/pow(_a,2)+_dda/_a)*f[3][idx]
                               -_a*_dda
                               + field->laplacian(f[3], j)/(3*pw2(exp(OSCSTART)))
                               +
                               2*field->V_lattice(f, idx ,_a)/(3*rescale_A*rescale_A)
                               - ( pow(_a,2) + 2*f[3][idx] )*(
                                            pow(df[0][idx]/_a - _da*f[0][idx]/pow(_a,2),2)
                                            +
                                            pow(df[1][idx]/_a - _da*f[1][idx]/pow(_a,2),2)
                                            +
                                            pow(df[2][idx]/_a - _da*f[2][idx]/pow(_a,2),2)
                                            )
                                            /(6*rescale_A*rescale_A)
                               - ( pow(_a,2) - 2*f[3][idx] )*2*(field->gradient_energy_eachpoint(f , 0, idx)
                                                     +
                                                     field->gradient_energy_eachpoint(f , 1, idx)
                                                     +
                                                     field->gradient_energy_eachpoint(f , 2, idx)
                                                     
                                                     )/(6*rescale_A*rescale_A*pw2(exp(OSCSTART)*_a))
                               ) * h*dt;
                
    #elif dim == 2
    //#pragma omp simd
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    df[3][idx] += (
                                   2*(pow(_da,2)/pow(_a,2)+_dda/_a)*f[3][idx]
                                   -_a*_dda
                                   + field->laplacian(f[3], j, k)/(3*pw2(exp(OSCSTART)))
                                   +
                                   2*field->V_lattice(f, idx ,_a)/(3*rescale_A*rescale_A)
                                   - ( pow(_a,2) + 2*f[3][idx] )*(
                                                                  pow(df[0][idx]/_a - _da*f[0][idx]/pow(_a,2),2)
                                                                  +
                                                                  pow(df[1][idx]/_a - _da*f[1][idx]/pow(_a,2),2)
                                                                  +
                                                                  pow(df[2][idx]/_a - _da*f[2][idx]/pow(_a,2),2)
                                                                  )
                                   /(6*rescale_A*rescale_A)
                                   - ( pow(_a,2) - 2*f[3][idx] )*2*(field->gradient_energy_eachpoint(f , 0, idx)
                                                                    +
                                                                    field->gradient_energy_eachpoint(f , 1, idx)
                                                                    +
                                                                    field->gradient_energy_eachpoint(f , 2, idx)
                                                                    
                                                                    )/(6*rescale_A*rescale_A*pw2(exp(OSCSTART)*_a))
                                   ) * h*dt;
                }
                
    #elif dim == 3
                
                for( int k = 0; k < N; ++k ){
   // #pragma omp simd
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        df[3][idx] += (
                                       2*(pow(_da,2)/pow(_a,2)+_dda/_a)*f[3][idx]
                                       -_a*_dda
                                       + field->laplacian(f[3], j, k, l)/(3*pw2(exp(OSCSTART)))
                                       +
                                       2*field->V_lattice(f, idx ,_a)/(3*rescale_A*rescale_A)
                                       - ( pow(_a,2) + 2*f[3][idx] )*(
                                                                      pow(df[0][idx]/_a - _da*f[0][idx]/pow(_a,2),2)
                                                                      +
                                                                      pow(df[1][idx]/_a - _da*f[1][idx]/pow(_a,2),2)
                                                                      +
                                                                      pow(df[2][idx]/_a - _da*f[2][idx]/pow(_a,2),2)
                                                                      )
                                       /(6*rescale_A*rescale_A)
                                       - ( pow(_a,2) - 2*f[3][idx] )*2*(field->gradient_energy_eachpoint(f , 0, idx)
                                                                        +
                                                                        field->gradient_energy_eachpoint(f , 1, idx)
                                                                        +
                                                                        field->gradient_energy_eachpoint(f , 2, idx)
                                                                        
                                                                        )/(6*rescale_A*rescale_A*pw2(exp(OSCSTART)*_a))
                                       ) * h*dt;
                    }
                        
                    }
                
    #endif
            }
            
    
}

void LeapFrog::fields_copy( double** f_from, double** f_to){
    
    
                    for( int i = 0; i < num_fields - 1; ++i ){
        #if   dim == 1
        #pragma omp parallel for simd schedule(static) num_threads(num_threads)
        #elif dim >= 2
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        #endif
                        for( int j = 0; j < N; ++j ){
        #if dim == 1
                            int idx = j;
                            f_to[i][idx] = f_from[i][idx];
        #elif dim == 2
        #pragma omp simd
                            for( int k = 0; k < N; ++k ){
                                int idx = j*N + k;
                                f_to[i][idx] = f_from[i][idx];
                            }
        #elif dim == 3
                            for( int k = 0; k < N; ++k ){
        #pragma omp simd
                                for( int l = 0; l < N; ++l ){
                                    int idx = (j*N + k)*N + l;
                                    f_to[i][idx] = f_from[i][idx];
                                }
                            }
        #endif
                        }
                    }
    
    
    

}

void LeapFrog::fields_convert( double** f, double** f_tilde, int convert_switch){

    switch (convert_switch) {
            //convert from program variables to evolution variables
        case 0:

            for( int i = 0; i < num_fields - 1; ++i ){
#if   dim == 1
#pragma omp parallel for simd schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for schedule( static ) num_threads( num_threads )
#endif
                for( int j = 0; j < N; ++j ){
#if dim == 1
                    int idx = j;
                    f_tilde[i][idx] = exp(Gamma_pr[i]*_sfint/2)*f[i][idx];

#elif dim == 2
#pragma omp simd
                    for( int k = 0; k < N; ++k ){
                        int idx = j*N + k;
                        
                        f_tilde[i][idx] = exp(Gamma_pr[i]*_sfint/2)*f[i][idx];
                      
                    }
#elif dim == 3
                    for( int k = 0; k < N; ++k ){
#pragma omp simd
                        for( int l = 0; l < N; ++l ){
                            int idx = (j*N + k)*N + l;
                            f_tilde[i][idx] = exp(Gamma_pr[i]*_sfint/2)*f[i][idx];

                        }
                    }
#endif
                }
            }

            break;

        case 1:
            //rescale fields from evolution variables to program variables
            for( int i = 0; i < num_fields - 1; ++i ){
#if   dim == 1
#pragma omp parallel for simd schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for schedule( static ) num_threads( num_threads )
#endif
                for( int j = 0; j < N; ++j ){
#if dim == 1
                    int idx = j;
                    f[i][idx] = f_tilde[i][idx]*exp(-Gamma_pr[i]*_sfint/2);
                
#elif dim == 2
#pragma omp simd
                    for( int k = 0; k < N; ++k ){
                        int idx = j*N + k;
                        f[i][idx] = f_tilde[i][idx]*exp(-Gamma_pr[i]*_sfint/2);
                       
                    }
#elif dim == 3
                    for( int k = 0; k < N; ++k ){
#pragma omp simd
                        for( int l = 0; l < N; ++l ){
                            int idx = (j*N + k)*N + l;
                            f[i][idx] = f_tilde[i][idx]*exp(-Gamma_pr[i]*_sfint/2);
                            
                        }
                    }
#endif
                }
            }
            break;

        default:
            break;
    }

}

void LeapFrog::fields_deriv_convert( double** f, double** df, double** f_tilde, double** df_tilde, int convert_switch){
    
    switch (convert_switch) {
            
            //convert from program variables to evolution variables
        case 0:
            for( int i = 0; i < num_fields - 1; ++i ){
#if   dim == 1
#pragma omp parallel for simd schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for schedule( static ) num_threads( num_threads )
#endif
                for( int j = 0; j < N; ++j ){
#if dim == 1
                    int idx = j;
//                    std::cout << "1: df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
//                    std::cout << "1: f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
                    df_tilde[i][idx] = exp(Gamma_pr[i]*_sfint/2)*(df[i][idx] +f[i][idx]*Gamma_pr[i]*_a/2);
//                    std::cout << "2: df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
//                    std::cout << "2: f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
#elif dim == 2
#pragma omp simd
                    for( int k = 0; k < N; ++k ){
                        int idx = j*N + k;
                        df_tilde[i][idx] = exp(Gamma_pr[i]*_sfint/2)*(df[i][idx] +f[i][idx]*Gamma_pr[i]*_a/2);
                    }
#elif dim == 3
                    for( int k = 0; k < N; ++k ){
#pragma omp simd
                        for( int l = 0; l < N; ++l ){
                            int idx = (j*N + k)*N + l;
                             df_tilde[i][idx] = exp(Gamma_pr[i]*_sfint/2)*(df[i][idx] +f[i][idx]*Gamma_pr[i]*_a/2);
                        }
                    }
#endif
                }
            }
            
            break;
            
        case 1:
            //rescale fields from evolution variables to program variables
            for( int i = 0; i < num_fields - 1; ++i ){
#if   dim == 1
#pragma omp parallel for simd schedule(static) num_threads(num_threads)
#elif dim >= 2
#pragma omp parallel for schedule( static ) num_threads( num_threads )
#endif
                for( int j = 0; j < N; ++j ){
#if dim == 1
                    int idx = j;
                    df[i][idx] = (df_tilde[i][idx] - Gamma_pr[i]*f_tilde[i][idx]*_a/2)*exp(-Gamma_pr[i]*_sfint/2);
#elif dim == 2
#pragma omp simd
                    for( int k = 0; k < N; ++k ){
                        int idx = j*N + k;
                        df[i][idx] = (df_tilde[i][idx] - Gamma_pr[i]*f_tilde[i][idx]*_a/2)*exp(-Gamma_pr[i]*_sfint/2);
                    }
#elif dim == 3
                    for( int k = 0; k < N; ++k ){
#pragma omp simd
                        for( int l = 0; l < N; ++l ){
                            int idx = (j*N + k)*N + l;
                            df[i][idx] = (df_tilde[i][idx] - Gamma_pr[i]*f_tilde[i][idx]*_a/2)*exp(-Gamma_pr[i]*_sfint/2);
                        }
                    }
#endif
                }
            }
            break;
            
        default:
            break;
    }
    
}
/////////////////////////////////////////////////////
//////////////////PUBLIC FUNCTIONS///////////////////
/////////////////////////////////////////////////////

void LeapFrog::evolution_expansion( Field* field, double** f, double** df, double& rad)
{
    
    switch ( expansion )
    {   
       
        case 1: // self-consistent evolution
            switch( precision )
            {
                case 2:
//                    std::cout << "before evolution expansion" << std::endl;
//                    std::cout << "f[0][4] = " << f[0][4] << std::endl;
//                    std::cout << "df[0][4] = " <<  df[0][4] << std::endl;
//                    std::cout << "f[0][10] = " <<  f[0][10] << std::endl;
//                    std::cout << "df[0][10] = " <<  df[0][10] << std::endl;
//                    std::cout << "f_tilde[0][4] = " << f_tilde[0][4] << std::endl;
//                    std::cout << "df_tilde[0][4] = " <<  df_tilde[0][4] << std::endl;
//                    std::cout << "f_tilde[0][10] = " <<  f_tilde[0][10] << std::endl;
//                    std::cout << "df_tilde[0][10] = " <<  df_tilde[0][10] << std::endl;
                    
                    evol_radiation(field, f_tilde, df_tilde, rad, 0.5);
                   
                    evol_fields( f_tilde, f_tilde, df_tilde, 0.5 );
                    
                    evol_scale(0.5);
                    
                    evol_gravpot( f, df,  0.5 ); //gravitational potential is stored in f[3][idx], df[3][idx]
                    
                    //t += 0.5*dt;
                    evol_sfint(0.5);
                    
                    for( int i = 0; i < st_output_step; ++i ){
//                        std::cout << "step 1" << std::endl;
                        evol_scale_dderivs( field, f, f_tilde, rad, 1);
//                        std::cout << "step 2" << std::endl;
//                        std::cout << "1 f_tilde[0][0] = " << f_tilde[0][0] << std::endl;
//                        std::cout << "1 df_tilde[0][0] = " << df_tilde[0][0] << std::endl;
//                        std::cout << "1 fdf_save[0][0] = " << fdf_save[0][0] << std::endl;
                        fields_copy( df_tilde, fdf_save);//Temporarily save df_tilde data to fdf_save
//                        std::cout << "2 f_tilde[0][0] = " << f_tilde[0][0] << std::endl;
//                        std::cout << "2 df_tilde[0][0] = " << df_tilde[0][0] << std::endl;
//                        std::cout << "2 fdf_save[0][0] = " << fdf_save[0][0] << std::endl;
                        evol_field_derivs_expansion( df_tilde, df_tilde, f, f_tilde, field,  0.5 );
//                         std::cout << "step 4" << std::endl;
//                        std::cout << "3 f_tilde[0][0] = " << f_tilde[0][0] << std::endl;
//                        std::cout << "3 df_tilde[0][0] = " << df_tilde[0][0] << std::endl;
//                        std::cout << "3 fdf_save[0][0] = " << fdf_save[0][0] << std::endl;
                        evol_gravpot_derivs_expansion( f, df, f_tilde, df_tilde, field,  1);
//                        std::cout << "step 5" << std::endl;
//                        std::cout << "4 f_tilde[0][0] = " << f_tilde[0][0] << std::endl;
//                        std::cout << "4 df_tilde[0][0] = " << df_tilde[0][0] << std::endl;
//                        std::cout << "4 fdf_save[0][0] = " << fdf_save[0][0] << std::endl;
                        evol_field_derivs_expansion( fdf_save, df_tilde, f, f_tilde, field, 1 );//Use data in fdf_save for leapfrog
//                        std::cout << "step 6" << std::endl;
//                        std::cout << "5 f_tilde[0][0] = " << f_tilde[0][0] << std::endl;
//                        std::cout << "5 df_tilde[0][0] = " << df_tilde[0][0] << std::endl;
//                        std::cout << "5 fdf_save[0][0] = " << fdf_save[0][0] << std::endl;
                        evol_gravpot( f, df, 1 ); //gravitational potential is stored in f[3][idx], df[3][idx]
//                        std::cout << "6 f_tilde[0][0] = " << f_tilde[0][0] << std::endl;
//                        std::cout << "6 df_tilde[0][0] = " << df_tilde[0][0] << std::endl;
//                        std::cout << "6 fdf_save[0][0] = " << fdf_save[0][0] << std::endl;
                        fields_copy( f_tilde, fdf_save);//Temporarily save f_tilde data to fdf_save
//                        std::cout << "7 f_tilde[0][0] = " << f_tilde[0][0] << std::endl;
//                        std::cout << "7 df_tilde[0][0] = " << df_tilde[0][0] << std::endl;
//                        std::cout << "7 fdf_save[0][0] = " << fdf_save[0][0] << std::endl;
                        evol_scale_derivs(0.5);
                        
                        evol_fields( f_tilde, f_tilde, df_tilde, 0.5 );
                        
                        evol_scale(0.5);
                        
                        //t += 0.5*dt;
                        evol_sfint(0.5);
                        
                        evol_radiation(field, f_tilde, df_tilde, rad, 1);
                        evol_fields(fdf_save, f_tilde, df_tilde, 1 );//Use data in fdf_save for leapfrog
                        evol_scale(0.5);
                        
                        //t += 0.5*dt;
                        evol_sfint(0.5);
                        
                        if( i == st_output_step - 1 ){
                            
                            evol_gravpot( f, df, -0.5 );//gravitational potential is stored in f[3][idx], df[3][idx]
                            evol_fields( f_tilde, f_tilde, df_tilde, -0.5 );
                            evol_scale(-0.5);
                            
                            //t -= 0.5*dt;
                            evol_sfint(-0.5);
                            
                            evol_radiation(field, f_tilde, df_tilde, rad, -0.5);
                            
                                                    }
                        
                    }
                    
                    evol_scale_dderivs( field, f, f_tilde, rad, 0); //scalar fields in f are also updated
                    
                    //set derivative of scalar fields back to program variables before leaving the evolution_expansion function
                    fields_deriv_convert( f, df, f_tilde, df_tilde, 1);
                    
                    
//                    std::cout << "after evolution expansion" << std::endl;
//                    std::cout << "f[0][4] = " << f[0][4] << std::endl;
//                    std::cout << "df[0][4] = " <<  df[0][4] << std::endl;
//                    std::cout << "f[0][10] = " <<  f[0][10] << std::endl;
//                    std::cout << "df[0][10] = " <<  df[0][10] << std::endl;
//                    std::cout << "f_tilde[0][4] = " << f_tilde[0][4] << std::endl;
//                    std::cout << "df_tilde[0][4] = " <<  df_tilde[0][4] << std::endl;
//                    std::cout << "f_tilde[0][10] = " <<  f_tilde[0][10] << std::endl;
//                    std::cout << "df_tilde[0][10] = " <<  df_tilde[0][10] << std::endl;
//                   //a(0) adot(0) addot(0)
//                   evol_fields( f, df, 0.5 );  //a(0) adot(0) addot(0)
//                     evol_scale(0.5); //a(0.5dt) adot(0) addot(0)
//                    for( int i = 0; i < st_output_step; ++i ){ //st_output_step=3
//                 //       evol_scale_dderivs( field, f, 1. ); //a(0.5dt) adot(0.5dt) addot(0.5dt) //a(1.5dt) adot(1.5dt) addot(1.5dt)//a(2.5dt) adot(2.5dt) addot(2.5dt)
//                       evol_field_derivs_expansion( f, df, field, 1.0 ); //a(0.5dt) adot(0.5dt) addot(0.5dt) //a(1.5dt) adot(1.5dt) addot(1.5dt)//a(2.5dt) adot(2.5dt) addot(2.5dt)
//                         evol_scale_derivs( 1.0 ); //a(0.5dt) adot(dt) addot(0.5dt) //a(1.5dt) adot(2dt) addot(1.5dt)//a(2.5dt) adot(3dt) addot(2.5dt)
//                        evol_fields( f, df, 1.0 ); //a(0.5dt) adot(dt) addot(0.5dt)//a(1.5dt) adot(2dt) addot(1.5dt)//a(2.5dt) adot(3dt) addot(2.5dt)
//                        evol_scale( 1. ); //a(1.5dt) adot(dt) addot(0.5dt)//a(2.5dt) adot(2dt) addot(1.5dt)//a(3.5dt) adot(3dt) addot(2.5dt)
//                      if( i == st_output_step - 1 )
//                        {
//                              evol_fields( f, df, -0.5 );  //a(3.5dt) adot(3dt) addot(2.5dt)
//                            evol_scale( -0.5 );  //a(3dt) adot(3dt) addot(2.5dt)*/
//                        }
//                    }
//                //     evol_scale_dderivs( field, f, 0.);//a(3dt) adot(3dt) addot(3dt)
                    break;
              
                case 4:
                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
                   /* for( int i = 0; i < st_output_step; ++i )
                    {
                        evol_fields( f, df, C[0] );
                        evol_scale( field, f, C[0] );
                        for( int p = 0; p < 3; ++p ){
                            evol_field_derivs_expansion( f, df, field, D[p] );
                            evol_scale_derivs( D[p] );
                            evol_fields( f, df, C[p+1] );
                            evol_scale( field, f, C[p+1] );
                        }
                    }*/
                    break;
            }
            break;

//        case 2: // radiation dominant universe
//
//            switch( precision )
//            {
//
//                case 2:
//                    evol_fields( f, df, 0.5 );
//                    t +=  0.5*dt;
//                    sfbase = t*hubble_init/sfexponent +1;
//                    _a = pow( sfbase, sfexponent );
//
//                    _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
//                    for( int i = 0; i < st_output_step; ++i ){
//                        evol_field_derivs_expansion( f, df, field, 1.0 );
//                        if( i == st_output_step - 1 )
//                        {
//                            evol_fields( f, df, 0.5 );
//                            t +=  0.5*dt;
//                        }
//                        else
//                        {
//                            evol_fields( f, df, 1.0 );
//                            t +=  dt;
//                        }
//                        sfbase = t*hubble_init/sfexponent +1;
//                        _a = pow( sfbase, sfexponent );
//                        _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
//                    }
//                    _da = hubble_init/sfbase*_a;
//                    break;
//
//                case 4:
//                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
//                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
//                    for( int i = 0; i < st_output_step; ++i ){
//                        evol_fields( f, df, C[0] );
//                        for( int p = 0; p < 3; ++p ) {
//                            _a +=  C[p]*dt;
//                            evol_field_derivs_expansion( f, df, field, D[p] );
//                            evol_fields( f, df, C[p+1] );
//                        }
//                        _a +=  C[3]*dt;
//                    }
//                    break;
//            }
//            break;
//
//        case 3: // matter dominant universe
//
//           // std::cout << _a << sfbase << sfexponent << std::endl;
//            switch( precision )
//            {
//
//                case 2:
//                  // std::cout << "_a = " << _a << "_da = " << _da << std::endl;
//                    evol_fields( f, df, 0.5 );
//                    t +=  0.5*dt;
//                     sfbase = t*hubble_init/sfexponent +1;
//                    _a = pow( sfbase, sfexponent );
//                  //     std::cout << hubble_init << std::endl;
//                //    std::cout << _a << sfbase << sfexponent << std::endl;
//                    _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
//
//                    for( int i = 0; i < st_output_step; ++i ){
//                        evol_field_derivs_expansion( f, df, field, 1.0 );
//                        if( i == st_output_step - 1 )
//                        {
//                            evol_fields( f, df, 0.5 );
//                            t +=  0.5*dt;
//                        }
//                        else
//                        {
//                            evol_fields( f, df, 1.0 );
//                            t +=  dt;
//                        }
//                         sfbase = t*hubble_init/sfexponent +1;
//                        _a = pow( sfbase, sfexponent );
//                        _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
//                    }
//                    _da = hubble_init/sfbase*_a;
//                    break;
//
//                case 4:
//                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
//                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
//                    for( int i = 0; i < st_output_step; ++i ){
//                        evol_fields( f, df, C[0] );
//                        for( int p = 0; p < 3; ++p ) {
//                            t +=  C[p]*dt;
//                            _a = pow( t, 2 );
//                            evol_field_derivs_expansion( f, df, field, D[p] );
//                            evol_fields( f, df, C[p+1] );
//                        }
//                        t +=  C[3]*dt;
//                        _a = pow( t, 2 );
//                    }
//                    _da = 2*t;
//                    break;
//            }
//            break;
//
        default:
            Logout( "Error: Parameter 'expansion' must be 0 ~ 3. \n" );
            exit(1);
    }
}
