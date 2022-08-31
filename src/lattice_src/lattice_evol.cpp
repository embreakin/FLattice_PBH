#include <cmath>
#include "parameters.hpp"
#include "utilities.hpp"
#include "lattice_evol.hpp"

/////////////////////////////////////////////////////
////////////////////CONSTRUCTER//////////////////////
/////////////////////////////////////////////////////

LeapFrog::LeapFrog( Field* field, double** f, double ** df, double rad_pr )
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
    
    // Set Decay Rates
        Gamma_pr[0]=Gamma1/rescale_B;
        Gamma_pr[1]=Gamma2/rescale_B;
        Gamma_pr[2]=Gamma3/rescale_B;
    
    // Set scalar fields at the minimum of the potential
    f_min_MPl[0]=0;
    f_min_MPl[1]=FIXPSI;
    f_min_MPl[2]=0;
    
    //Rescale field program variables to field variables necessary for leapfrog evolution
    //We only rescale the scalar fields.
    for( int i = 0; i < num_fields - 1; ++i ){
    #if   dim == 1
    #pragma omp parallel for simd schedule(static) num_threads(num_threads)
    #elif dim >= 2
    #pragma omp parallel for schedule( static ) num_threads( num_threads )
    #endif
        for( int j = 0; j < N; ++j ){
            #if dim == 1
                int idx = j;
            f_tilde[i][idx] = exp(Gamma_pr[i]*t0/2)*f[i][idx];
            df_tilde[i][idx] = exp(Gamma_pr[i]*t0/2)*(df[i][idx] +f[i][idx]*Gamma_pr[i]/2);
            #elif dim == 2
            #pragma omp simd
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f_tilde[i][idx] = exp(Gamma_pr[i]*t0/2)*f[i][idx];
                    df_tilde[i][idx] = exp(Gamma_pr[i]*t0/2)*(df[i][idx] +f[i][idx]*Gamma_pr[i]/2);
                }
            #elif dim == 3
                for( int k = 0; k < N; ++k ){
                    #pragma omp simd
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f_tilde[i][idx] = exp(Gamma_pr[i]*t0/2)*f[i][idx];
                        df_tilde[i][idx] = exp(Gamma_pr[i]*t0/2)*(df[i][idx] +f[i][idx]*Gamma_pr[i]/2);
                    }
                }
            #endif
        }
    }
    
    
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
            _da = hubble_init/rescale_B;
            //In the function below, a'' is calculated by converting f_tilde to f.
           evol_scale_dderivs( field, f, f_tilde, rad_pr ,t0, 0.); // _dda( Field* field, double** f, double ** f_tilde, double rho_r, double t, double h)
            break;
        case 2: // radiation dominant universe
            sfexponent = 1;
            sfbase = 1;
            _a  = 1;
            _da = hubble_init/rescale_B;
            _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
            break;
        case 3: // matter dominant universe
            sfexponent = 2;
            sfbase = 1;
            _a  = 1;
            _da = hubble_init/rescale_B;
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

void LeapFrog::evol_fields( double** f, double** df, double h )
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
                f[i][idx] += df[i][idx] * h*dt;
            #elif dim == 2
            #pragma omp simd
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f[i][idx] += df[i][idx] * h*dt;
                }
            #elif dim == 3
                for( int k = 0; k < N; ++k ){
                    #pragma omp simd
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f[i][idx] += df[i][idx] * h*dt;
                    }
                }
            #endif
        }
    }
}


void LeapFrog::evol_field_derivs_expansion( double** f, double** df, Field* field, double h )
{
    
    for( int i = 0; i < num_fields - 1; ++i ){
        
        #if   dim == 1
        #pragma omp parallel for simd schedule(static) num_threads(num_threads)
        #elif dim >= 2
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        #endif
                    for( int j = 0; j < N; ++j ){
        #if dim == 1
                        int idx = j;
                        df[i][idx] += (
                                       (
                                        field->laplacian(f[i], j) - field->mass(i,_a)
                                        )/pw2(exp(OSCSTART))
                                       + _dda/_a + (_da/_a)*Gamma_pr[i]
                                       + pow(Gamma_pr[i],2)
                                       ) *f[i][idx]* h*dt;
                        
        #elif dim == 2
                        //    #pragma omp simd
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                            df[i][idx] += (
                                        
                                           (
                                            field->laplacian(f[i], j, k) - field->mass(i,_a)
                                            )/pw2(exp(OSCSTART))
                                           + _dda/_a + (_da/_a)*Gamma_pr[i]
                                           + pow(Gamma_pr[i],2)
                                           ) *f[i][idx]* h*dt;
                            
                        }
                        
        #elif dim == 3
                        
                        for( int k = 0; k < N; ++k ){
        #pragma omp simd
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;
                                df[i][idx] +=(
                                    
                                              (
                                               field->laplacian(f[i], j, k, l) - field->mass(i,_a)
                                               )/pw2(exp(OSCSTART))
                                              + _dda/_a + (_da/_a)*Gamma_pr[i]
                                              + pow(Gamma_pr[i],2)
                                              ) *f[i][idx]* h*dt;
                                
                            }
                        }
        #endif
                    }
        
    }
    
}


void LeapFrog::evol_scale_dderivs( Field* field, double** f, double ** f_tilde, double rho_r, double t, double h)
{
    
    
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
            f[i][idx] = f_tilde[i][idx]*exp(-Gamma_pr[i]*t/2);
            #elif dim == 2
            #pragma omp simd
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f[i][idx] = f_tilde[i][idx]*exp(-Gamma_pr[i]*t/2);
                }
            #elif dim == 3
                for( int k = 0; k < N; ++k ){
                    #pragma omp simd
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f[i][idx] = f_tilde[i][idx]*exp(-Gamma_pr[i]*t/2);
                    }
                }
            #endif
        }
    }
    
   if(h==0){
        _dda = -pw2(_da)/_a
       +(1/(_a*rescale_A*rescale_A))*
       (
        2.*(field->gradient_energy(f[0])
            +field->gradient_energy(f[1])
            +field->gradient_energy(f[2])
            )/(3.*pw2(exp(OSCSTART)))
        + field->potential_energy(f, _a)
        + rho_r
        );
    }else{
    
    double C = 0;
    C += 2*(field->gradient_energy(f[0])
            +field->gradient_energy(f[1])
            +field->gradient_energy(f[2])
            )/(3.*pw2(exp(OSCSTART)))
        + field->potential_energy( f, _a )
        +rho_r;
    

    _dda = - 2./(h*dt) * (
                          _da +
                          _a/(h*dt) *
                          ( 1 -
                           sqrt( 1 + 2*h*dt*_da/_a
                                + pow(h*dt,2.)*C/pow(rescale_A*_a,2)  )
                           )
                          );
       
        _da += .5*_dda*h*dt;
        
    }
}

void LeapFrog::evol_radiation(Field* field, double** f, double** df, double rho, double t , double h)
{
    double drho = 0;
    
        for( int i = 0; i < num_fields - 1; ++i ){
            
         drho += Gamma_pr[i]*exp(-Gamma_pr[i]*t)*pow(field->df_average(df[i], i)-(Gamma_pr[i]/2 + _da/_a)*field->f_average(f[i], i),2)* h*dt;
            
        }
    
    rho = rho + drho*h*dt;
    
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


void LeapFrog::evol_gravpot_derivs_expansion( double** f, double** df, double** f_tilde, double** df_tilde, Field* field, double t, double h )
{
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
                f[i][idx] = f_tilde[i][idx]*exp(-Gamma_pr[i]*t/2);
                df[i][idx] = (df_tilde[i][idx] - Gamma_pr[i]*f_tilde[i][idx]/2)*exp(-Gamma_pr[i]*t/2);
    #elif dim == 2
    #pragma omp simd
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f[i][idx] = f_tilde[i][idx]*exp(-Gamma_pr[i]*t/2);
                    df[i][idx] = (df_tilde[i][idx] - Gamma_pr[i]*f_tilde[i][idx]/2)*exp(-Gamma_pr[i]*t/2);
                }
    #elif dim == 3
                for( int k = 0; k < N; ++k ){
    #pragma omp simd
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f[i][idx] = f_tilde[i][idx]*exp(-Gamma_pr[i]*t/2);
                        df[i][idx] = (df_tilde[i][idx] - Gamma_pr[i]*f_tilde[i][idx]/2)*exp(-Gamma_pr[i]*t/2);
                    }
                }
    #endif
            }
        }
    
    #if   dim == 1
    #pragma omp parallel for simd schedule(static) num_threads(num_threads)
    #elif dim >= 2
    #pragma omp parallel for schedule( static ) num_threads( num_threads )
    #endif
            for( int j = 0; j < N; ++j ){
    #if dim == 1
                int idx = j;
                df[3][idx] += (
                               - field->laplacian(f[3], j)/(3*pw2(exp(OSCSTART)))
                               + f[3][idx]*(
                                            pow(df[0][idx]/_a - _da*f[0][idx]/pow(_a,2),2)
                                            +
                                            pow(df[1][idx]/_a - _da*f[1][idx]/pow(_a,2),2)
                                            +
                                            pow(df[2][idx]/_a - _da*f[2][idx]/pow(_a,2),2)
                                            )
                                            /(6*rescale_A*rescale_A)
                               + (1+f[3][idx]/_a)*2*(field->gradient_energy_eachpoint(f , 0, idx)
                                                     +
                                                     field->gradient_energy_eachpoint(f , 1, idx)
                                                     +
                                                     field->gradient_energy_eachpoint(f , 2, idx)
                                                     
                                                     )/(6*rescale_A*rescale_A*pw2(exp(OSCSTART))*_a)
                               +
                               field->V_lattice(f, idx ,_a)/(6*rescale_A*rescale_A*_a)
                               -
                               (_dda -pow(_da,2)/_a)/2
                               ) * h*dt;
                
    #elif dim == 2
    //#pragma omp simd
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    df[3][idx] += (
                                   - field->laplacian(f[3], j, k)/(3*pw2(exp(OSCSTART)))
                                   + f[3][idx]*(
                                                pow(df[0][idx]/_a - _da*f[0][idx]/pow(_a,2),2)
                                                +
                                                pow(df[1][idx]/_a - _da*f[1][idx]/pow(_a,2),2)
                                                +
                                                pow(df[2][idx]/_a - _da*f[2][idx]/pow(_a,2),2)
                                                )
                                   /(6*rescale_A*rescale_A)
                                   + (1+f[3][idx]/_a)*2*(field->gradient_energy_eachpoint(f , 0, idx)
                                                         +
                                                         field->gradient_energy_eachpoint(f , 1, idx)
                                                         +
                                                         field->gradient_energy_eachpoint(f , 2, idx)
                                                         
                                                         )/(6*rescale_A*rescale_A*pw2(exp(OSCSTART))*_a)
                                   +
                                   field->V_lattice(f, idx ,_a)/(6*rescale_A*rescale_A*_a)
                                   -
                                   (_dda -pow(_da,2)/_a)/2
                                   ) * h*dt;
                    
                }
                
    #elif dim == 3
                
                for( int k = 0; k < N; ++k ){
    #pragma omp simd
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        df[3][idx] += (
                                       - field->laplacian(f[3], j, k, l)/(3*pw2(exp(OSCSTART)))
                                       + f[3][idx]*(
                                                    pow(df[0][idx]/_a - _da*f[0][idx]/pow(_a,2),2)
                                                    +
                                                    pow(df[1][idx]/_a - _da*f[1][idx]/pow(_a,2),2)
                                                    +
                                                    pow(df[2][idx]/_a - _da*f[2][idx]/pow(_a,2),2)
                                                    )
                                       /(6*rescale_A*rescale_A)
                                       + (1+f[3][idx]/_a)*2*(field->gradient_energy_eachpoint(f , 0, idx)
                                                             +
                                                             field->gradient_energy_eachpoint(f , 1, idx)
                                                             +
                                                             field->gradient_energy_eachpoint(f , 2, idx)
                                                             
                                                             )/(6*rescale_A*rescale_A*pw2(exp(OSCSTART))*_a)
                                       +
                                       field->V_lattice(f, idx ,_a)/(6*rescale_A*rescale_A*_a)
                                       -
                                       (_dda -pow(_da,2)/_a)/2
                                       ) * h*dt;
                        
                    }
                        
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


/////////////////////////////////////////////////////
//////////////////PUBLIC FUNCTIONS///////////////////
/////////////////////////////////////////////////////

void LeapFrog::evolution_expansion( Field* field, double** f, double** df, double rad, double t )
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
                    evol_radiation(field, f_tilde, df_tilde, rad, t , 0.5);
                    evol_fields( f_tilde, df_tilde, 0.5 );
                    evol_scale(0.5);
                    evol_gravpot( f, df,  0.5 ); //gravitational potential is stored in f[3][idx], df[3][idx]
                    
                    t += 0.5*dt;
                    
                    for( int i = 0; i < st_output_step; ++i ){
                        
                        evol_scale_dderivs( field, f, f_tilde, rad, t, 1);
                        fields_copy( df_tilde, fdf_save);//Temporarily save df_tilde data to fdf_save
                        evol_field_derivs_expansion( f_tilde, df_tilde, field,  0.5 );
                        evol_gravpot_derivs_expansion( f, df, f_tilde, df_tilde, field, t, 1);
                        evol_field_derivs_expansion( f_tilde, fdf_save, field, 1 );//Use data in fdf_save for leapfrog
                        evol_gravpot( f, df, 1 ); //gravitational potential is stored in f[3][idx], df[3][idx]
                        fields_copy( f_tilde, fdf_save);//Temporarily save f_tilde data to fdf_save
                        evol_scale_derivs(0.5);
                        evol_fields( f_tilde, df_tilde, 0.5 );
                        evol_scale(0.5);
                        t += 0.5*dt;
                        evol_radiation(field, f_tilde, df_tilde, rad, t , 1);
                        evol_fields( fdf_save, df_tilde, 1 );//Use data in fdf_save for leapfrog
                        evol_scale(0.5);
                       t += 0.5*dt;
                        
                        if( i == st_output_step - 1 ){
                            t -= 0.5*dt;
                            evol_gravpot( f, df, -0.5 );//gravitational potential is stored in f[3][idx], df[3][idx]
                            evol_fields( f_tilde, df_tilde, -0.5 );
                            evol_scale( -0.5 );
                            evol_radiation(field, f_tilde, df_tilde, rad, t , -0.5);
                            
                                                    }
                        
                    }
                    
                    evol_scale_dderivs( field, f, f_tilde, rad, t, 0); //scalar fields in f are also updated
                    fields_copy( df_tilde, df); //scalar fields in df are updated
                    //f,df are also updated
                    
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
