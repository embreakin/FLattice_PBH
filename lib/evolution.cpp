#include <cmath>
#include "parameter.hpp"
#include "utilities.hpp"
#include "evolution.hpp"

//double LeapFrog::_a  = 1;
//double LeapFrog::_da = 3.99076*pow(10,-7)/sqrt(6.80342*pow(10,-10));

LeapFrog::LeapFrog(): _dda()
{
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
        case 1:
           // _a  = 1;
           // _da = hubble_init;
            break;
        case 2:
           // _a  = 1;
          //  _da = hubble_init;
           // _da  = 1;
           // _dda = 0;
            break;
        case 3:
        _a  = 1;
       _da = hubble_init;//3.99076*pow(10,-7)/sqrt(6.80342*pow(10,-10));//hubble_init;
           // _dda = 2;
           // std::cout << "constructor _a = " << _a << "_da = " << _da << std::endl;
            break;
        default:
            Logout( "Parameter 'expansion' must be 0 ~ 3 when you use class 'LeapFrog'. \n" );
            exit(1);
    }
}

void LeapFrog::evol_fields( double** f, double** df, double h )
{	
    for( int i = 0; i < num_fields; ++i ){
  // #pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            #if dim == 1
                int idx = j;
                f[i][idx] += df[i][idx] * h*dt;
            #elif dim == 2
  
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f[i][idx] += df[i][idx] * h*dt;
                }
            #elif dim == 3
  
                for( int k = 0; k < N; ++k ){
  // #pragma omp parallel for schedule( static ) num_threads( num_threads )
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f[i][idx] += df[i][idx] * h*dt;
                    }
                }
            #endif
        }
    }
}

void LeapFrog::evol_field_derivs( double** f, double** df, Field* field, double h )
{
    for( int i = 0; i < num_fields; ++i ){
        #ifdef SPHERICAL_SYM
            df[i][N-1] -= ( field->laplacian(f[i], df[i], N-1, 0, 0) + f[i][N-1]/2) * h*dt;
        //    #pragma omp parallel for schedule( static ) num_threads(num_threads)
            for( int j = 0; j < N-1; ++j ){
		#else
      //    #pragma omp parallel for schedule( static ) num_threads( num_threads )
            for( int j = 0; j < N; ++j ){
        #endif
            #if dim == 1
                int idx = j;
                df[i][idx] += ( field->laplacian(f[i], j) - field->dV(f, i, idx) ) * h*dt;
            #elif dim == 2
   
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    df[i][idx] += ( field->laplacian(f[i], j, k) - field->dV(f, i, idx) ) * h*dt;
                }
            #elif dim == 3
       
                for( int k = 0; k < N; ++k ){
   // #pragma omp parallel for schedule( static ) num_threads( num_threads )
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        df[i][idx] += ( field->laplacian(f[i], j, k, l) - field->dV(f, i, idx) ) * h*dt;
                    }
                }
            #endif
        }
    }
}

void LeapFrog::evol_field_derivs_expansion( double** f, double** df, Field* field, double h )
{
    for( int i = 0; i < num_fields; ++i ){
  
        for( int j = 0; j < N; ++j ){
            #if dim == 1
                int idx = j;
                df[i][idx] += ( field->laplacian(f[i], j, 0, 0) + _dda*f[i][idx]/_a - pow(_a,3)*field->adV(f, i, idx, _a) ) * h*dt;
            #elif dim == 2
  //#pragma omp parallel for schedule( static ) num_threads( num_threads )
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    //df[i][idx] += ( field->laplacian(f[i], j, k, 0) + _dda*f[i][idx]/_a - pow(_a,3)* (field->*(field->adV[i]))(f, _a, i, idx) ) * h*dt;
                  
                    df[i][idx] += ( field->laplacian(f[i], j, k, 0) + _dda*f[i][idx]/_a - pow(_a,3)*field->adV(f, i, idx, _a) ) * h*dt;
                }
           
            #elif dim == 3
      
                for( int k = 0; k < N; ++k ){
    // #pragma omp parallel for schedule( static ) num_threads( num_threads )
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        df[i][idx] += ( field->laplacian(f[i], j, k, l) + _dda*f[i][idx]/_a - pow(_a,3)*field->adV(f, i, idx, _a) ) * h*dt;
                    }
                }
            #endif
        }
      /*  std::cout << "a = " << _a << std::endl;
       std::cout << "field->laplacian(f[i], 0, 10, 0) = " << field->laplacian(f[i], 0, 10, 0) << std::endl;
        std::cout << "_dda*f[i][10]/_a = " << _dda*f[i][10]/_a << std::endl;
        std::cout << "pow(_a,3)*field->adV(f, i, 10, _a) = " << pow(_a,3)*field->adV(f, i, 10, _a) << std::endl;
        std::cout << "field->laplacian(f[i], 0, 100, 0) = " << field->laplacian(f[i], 0, 100, 0) << std::endl;
        std::cout << "_dda*f[i][100]/_a = " << _dda*f[i][100]/_a << std::endl;
        std::cout << "pow(_a,3)*field->adV(f, i, 100, _a) = " << pow(_a,3)*field->adV(f, i, 100, _a) << std::endl;*/
    }
}


void LeapFrog::evol_scale_dderivs( Field* field, double** f, double h )
{
    double sfev1, sfev2, sfev3;
     sfev1 = 1.;
    sfev2 = 2.;
    sfev3 = 4.;
   if(h==0){
        _dda = -sfev1*pw2(_da)/_a + 8.*M_PI/pow(_a,sfev2-1)*(2.*field->gradient_energy(f[0])/3. + pow(_a,sfev3)*field->potential_energy(f, _a));
    }else{
    
    double C = 0;
    C += 2*field->gradient_energy(f[0])/3 + pow(_a,sfev3)*field->potential_energy( f, _a );
    

    _dda = - 2./(h*dt) * ( _da + _a/(sfev1*h*dt) * ( 1 - sqrt( 1 + 2*sfev1*h*dt*_da/_a + sfev1*8*M_PI*pow(h*dt,2.)*C/pow(_a,sfev2)  ) ) );
       
        _da += .5*_dda*h*dt;
        
    }
}


void LeapFrog::evolution( Field* field, double** f, double** df )
{
    if( expansion != 0 )
    {
        Logout( "Error: Parameter expansion must be 0 when you use member-function 'evolution'. \n" );
        Logout( "Use expansion = 1 ~ 3 and member-function 'evolution_expansion' when you simulate the expanding universe. \n" );
        exit(1);
    }
    
    switch( precision )
    {
        case 2:
            evol_fields( f, df, 0.5 );
            for( int i = 0; i < st_output_step; ++i )
            {
                evol_field_derivs( f, df, field, 1.0 );
                if( i == st_output_step - 1 ) evol_fields( f, df, 0.5 );
                else evol_fields( f, df, 1.0 );
            }
            break;

        case 4:
            const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
            const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
            for( int i = 0; i < st_output_step; ++i )
            {
                evol_fields( f, df, C[0] );
                for( int p = 0; p < 3; ++p )
                {
                    evol_field_derivs( f, df, field, D[p] );
                    evol_fields( f, df, C[p+1] );
                }
            }
            break;
    }
}

void LeapFrog::evolution_expansion( Field* field, double** f, double** df, double t )
{
    
    double sfexponent, sfbase;
    
    switch ( expansion )
    {   
        /* // The code is not well optimized when include member-function 'evol_field_derivs'
        case 0: // no expansion
            switch( precision )
            {
                case 2:
                    evol_fields( f, df, 0.5 );
                    for( int i = 0; i < output_step; ++i )
                    {
                        evol_field_derivs( f, df, field, 1.0 );
                        if( i == output_step - 1 ) evol_fields( f, df, 0.5 );
                        else evol_fields( f, df, 1.0 );
                    }
                    break;

                case 4:
                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
                    for( int i = 0; i < output_step; ++i )
                    {
                        evol_fields( f, df, C[0] );
                        for( int p = 0; p < 3; ++p ){
                            evol_field_derivs( f, df, field, D[p] );
                            evol_fields( f, df, C[p+1] );
                        }
                    }
                    break;
            }
            break;
        */
        case 1: // self-consistent evolution
           // _da = sqrt( ( field->gradient_energy(f[0]) + 2*pow(_a,4)*field->potential_energy(f, _a) )/6 );
          //  for( int i = 0; i < dim; ++i ) _da /= sqrt(N);
            switch( precision )
            {
                case 2:
                    evol_scale_dderivs( field, f, 0.);//a(0) adot(0) addot(0)
                   evol_fields( f, df, 0.5 );  //a(0) adot(0) addot(0)
                     evol_scale(0.5); //a(0.5dt) adot(0) addot(0)
                    for( int i = 0; i < st_output_step; ++i ){ //st_output_step=3
                        evol_scale_dderivs( field, f, 1. ); //a(0.5dt) adot(0.5dt) addot(0.5dt) //a(1.5dt) adot(1.5dt) addot(1.5dt)//a(2.5dt) adot(2.5dt) addot(2.5dt)
                       evol_field_derivs_expansion( f, df, field, 1.0 ); //a(0.5dt) adot(0.5dt) addot(0.5dt) //a(1.5dt) adot(1.5dt) addot(1.5dt)//a(2.5dt) adot(2.5dt) addot(2.5dt)
                         evol_scale_derivs( 1.0 ); //a(0.5dt) adot(dt) addot(0.5dt) //a(1.5dt) adot(2dt) addot(1.5dt)//a(2.5dt) adot(3dt) addot(2.5dt)
                        evol_fields( f, df, 1.0 ); //a(0.5dt) adot(dt) addot(0.5dt)//a(1.5dt) adot(2dt) addot(1.5dt)//a(2.5dt) adot(3dt) addot(2.5dt)
                        evol_scale( 1. ); //a(1.5dt) adot(dt) addot(0.5dt)//a(2.5dt) adot(2dt) addot(1.5dt)//a(3.5dt) adot(3dt) addot(2.5dt)
                      if( i == st_output_step - 1 )
                        {
                              evol_fields( f, df, -0.5 );  //a(3.5dt) adot(3dt) addot(2.5dt)
                            evol_scale( -0.5 );  //a(3dt) adot(3dt) addot(2.5dt)*/
                        }
                    }
                     evol_scale_dderivs( field, f, 0.);//a(3dt) adot(3dt) addot(3dt)
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

        case 2: // radiation dominant universe
            sfexponent = 1;
            switch( precision )
            {
                    
                case 2:
                    evol_fields( f, df, 0.5 );
                    t +=  0.5*dt;
                    sfbase = t*hubble_init/sfexponent +1;
                    _a = pow( sfbase, sfexponent );
                   
                    _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
                    for( int i = 0; i < st_output_step; ++i ){
                        evol_field_derivs_expansion( f, df, field, 1.0 );
                        if( i == st_output_step - 1 )
                        {
                            evol_fields( f, df, 0.5 );
                            t +=  0.5*dt;
                        }
                        else
                        {
                            evol_fields( f, df, 1.0 );
                            t +=  dt;
                        }
                        sfbase = t*hubble_init/sfexponent +1;
                        _a = pow( sfbase, sfexponent );
                        _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
                    }
                    _da = hubble_init/sfbase*_a;
                    break;

                case 4:
                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
                    for( int i = 0; i < st_output_step; ++i ){
                        evol_fields( f, df, C[0] );
                        for( int p = 0; p < 3; ++p ) {
                            _a +=  C[p]*dt;
                            evol_field_derivs_expansion( f, df, field, D[p] );
                            evol_fields( f, df, C[p+1] );
                        }
                        _a +=  C[3]*dt;
                    }
                    break;
            }
            break;

        case 3: // matter dominant universe
            sfexponent = 2;
           // std::cout << _a << sfbase << sfexponent << std::endl;
           
           // std::cout << _a << sfbase << sfexponent << std::endl;
            switch( precision )
            {
                  
                    
                case 2:
                  // std::cout << "_a = " << _a << "_da = " << _da << std::endl;
                    evol_fields( f, df, 0.5 );
                    t +=  0.5*dt;
                     sfbase = t*hubble_init/sfexponent +1;
                    _a = pow( sfbase, sfexponent );
                  //     std::cout << hubble_init << std::endl;
                //    std::cout << _a << sfbase << sfexponent << std::endl;
                    _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
                    for( int i = 0; i < st_output_step; ++i ){
                        evol_field_derivs_expansion( f, df, field, 1.0 );
                        if( i == st_output_step - 1 )
                        {
                            evol_fields( f, df, 0.5 );
                            t +=  0.5*dt;
                        }
                        else
                        {
                            evol_fields( f, df, 1.0 );
                            t +=  dt;
                        }
                         sfbase = t*hubble_init/sfexponent +1;
                        _a = pow( sfbase, sfexponent );
                        _dda = (sfexponent -1)/sfexponent*pw2(hubble_init)/pw2(sfbase)*_a;
                    }
                    _da = hubble_init/sfbase*_a;
                    break;

                case 4:
                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
                    for( int i = 0; i < st_output_step; ++i ){
                        evol_fields( f, df, C[0] );
                        for( int p = 0; p < 3; ++p ) {
                            t +=  C[p]*dt;
                            _a = pow( t, 2 );
                            evol_field_derivs_expansion( f, df, field, D[p] );
                            evol_fields( f, df, C[p+1] );
                        }
                        t +=  C[3]*dt;
                        _a = pow( t, 2 );
                    }
                    _da = 2*t;
                    break;
            }
            break;
        
        default:
            Logout( "Error: Parameter 'expansion' must be 0 ~ 3. \n" );
            exit(1);
    }
}
