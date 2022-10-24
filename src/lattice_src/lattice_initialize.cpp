
#include "lattice_initialize.hpp"
#include "utilities.hpp"
#include <random>

double Fk_log_int_calc(int k_int , double** lattice_var, int num_field)
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
            // std::cout << "pow(lattice_var[k_int-1][7+3*num_field],2) = " << pow(lattice_var[k_int-1][7+3*num_field],2) << std::endl;
            //             std::cout << "pow(lattice_var[k_int-1][31+3*num_field],2) = " << pow(lattice_var[k_int-1][31+3*num_field],2) << std::endl;
            //             std::cout << "pow(lattice_var[k_int-1][8+3*num_field],2) = " << pow(lattice_var[k_int-1][8+3*num_field],2) << std::endl;
            //             std::cout << "pow(lattice_var[k_int-1][32+3*num_field],2) = " << pow(lattice_var[k_int-1][32+3*num_field],2) << std::endl;
            //             std::cout << "pow(lattice_var[k_int-1][9+3*num_field],2) = " << pow(lattice_var[k_int-1][9+3*num_field],2) << std::endl;
            //             std::cout << "pow(lattice_var[k_int-1][33+3*num_field],2) = " << pow(lattice_var[k_int-1][33+3*num_field],2) << std::endl;
            //            std::cout << "dF_k = " << dF_k << std::endl;
            //           std::cout << "F_k = " << F_k << std::endl;
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
    
    return log10(F_k);
}

double dFk_log_int_calc(int k_int , double** lattice_var, int num_field)
{
    double dF_k;
    
    
    switch( num_field )
    {
        case 0: //sigma field
        case 1: //psi field
        case 2:  //phi field
            //Fourier mode of the derivative before discretization
            dF_k = sqrt(pow(lattice_var[k_int-1][16+3*num_field],2) + pow(lattice_var[k_int-1][40+3*num_field],2)
                        +pow(lattice_var[k_int-1][17+3*num_field],2) + pow(lattice_var[k_int-1][41+3*num_field],2)
                        +pow(lattice_var[k_int-1][18+3*num_field],2) + pow(lattice_var[k_int-1][42+3*num_field],2)
                        );
            
            //std::cout << "dF_k = " << dF_k << std::endl;
            
            break;
        case 3: //gravitational potential
            //Fourier mode of the derivative before discretization
            dF_k = sqrt(
                        pow(lattice_var[k_int-1][28],2) + pow(lattice_var[k_int-1][52],2)
                        +pow(lattice_var[k_int-1][29],2) + pow(lattice_var[k_int-1][53],2)
                        +pow(lattice_var[k_int-1][30],2) + pow(lattice_var[k_int-1][54],2)
                        );
            break;
        default:
            Logout( "Parameter 'num_field' must be 0 ~ 3. \n" );
            exit(1);
            
    }
    
    return log10(dF_k);
}


void fdf_calc(double distance, double** lattice_var, double *field, double *deriv, int num_field)
{
    int l;
    double log10_Fk, log10_Fk2, log10_dFk, log10_dFk2, omega;
    l = floor(distance);
    
    if (k_lattice_grid_min_MPl < kfrom_MPl_lattice)
    {
        if(l < outrange_num + 1){
            
            log10_Fk = Fk_log_int_calc(outrange_num+1, lattice_var, num_field);
            
            log10_Fk2 = Fk_log_int_calc(outrange_num+2, lattice_var, num_field);
            
            log10_Fk += (log10(distance) - log10(outrange_num+1))*(log10_Fk2 - log10_Fk)/(log10(outrange_num+2) - log10(outrange_num+1));
            
            log10_dFk = dFk_log_int_calc(outrange_num+1, lattice_var, num_field);
            
            log10_dFk2 = dFk_log_int_calc(outrange_num+2, lattice_var, num_field);
            
            log10_dFk += (log10(distance) - log10(outrange_num+1))*(log10_dFk2 - log10_dFk)/(log10(outrange_num+2) - log10(outrange_num+1));
            
            
            
            field[0] = pow(10,log10_Fk);
            deriv[0] = pow(10,log10_dFk);
            
        }else if (l < N/2)
        {
            log10_Fk = Fk_log_int_calc(l, lattice_var, num_field);
            log10_Fk2 = Fk_log_int_calc(l+1, lattice_var, num_field);
            log10_Fk += (log10(distance) - log10(l))*(log10_Fk2 - log10_Fk)/(log10(l+1) - log10(l));
            
            log10_dFk = dFk_log_int_calc(l, lattice_var, num_field);
            log10_dFk2 = dFk_log_int_calc(l+1, lattice_var, num_field);
            
            log10_dFk += (log10(distance) - log10(l))*(log10_dFk2 - log10_dFk)/(log10(l+1) - log10(l));
            field[0] = pow(10,log10_Fk);
            deriv[0] = pow(10,log10_dFk);
        }
        else
        {
            
            log10_Fk = Fk_log_int_calc(N/2-1, lattice_var, num_field);
            
            log10_Fk2 = Fk_log_int_calc(N/2, lattice_var, num_field);
            
            log10_Fk += (log10(distance) - log10(N/2-1))*(log10_Fk2 - log10_Fk)/(log10(N/2) - log10(N/2-1));
            
            log10_dFk = dFk_log_int_calc(N/2-1, lattice_var, num_field);
            
            log10_dFk2 = dFk_log_int_calc(N/2, lattice_var, num_field);
            
            log10_dFk += (log10(distance) - log10(N/2-1))*(log10_dFk2 - log10_dFk)/(log10(N/2) - log10(N/2-1));
            
            
            
            field[0] = pow(10,log10_Fk);
            deriv[0] = pow(10,log10_dFk);
        }
        
    }
    else
    {
        if (l < N/2)
        {
            
            
            
            log10_Fk = Fk_log_int_calc(l, lattice_var, num_field);
            log10_Fk2 = Fk_log_int_calc(l+1, lattice_var, num_field);
            log10_Fk += (log10(distance) - log10(l))*(log10_Fk2 - log10_Fk)/(log10(l+1) - log10(l));
            
            log10_dFk = dFk_log_int_calc(l, lattice_var, num_field);
            log10_dFk2 = dFk_log_int_calc(l+1, lattice_var, num_field);
            
            log10_dFk += (log10(distance) - log10(l))*(log10_dFk2 - log10_dFk)/(log10(l+1) - log10(l));
            field[0] = pow(10,log10_Fk);
            deriv[0] = pow(10,log10_dFk);
            
            //        std::cout << "field[0] = " << field[0] << std::endl;
            //        std::cout << "deriv[0] = " << deriv[0] << std::endl;
            //
            
            //omega = abs(pow(10,log10_dFk)/(pow(10, log10_Fk)/rescale_B));
            
        }
        else
        {
            
            log10_Fk = Fk_log_int_calc(N/2-1, lattice_var, num_field);
            
            log10_Fk2 = Fk_log_int_calc(N/2, lattice_var, num_field);
            
            log10_Fk += (log10(distance) - log10(N/2-1))*(log10_Fk2 - log10_Fk)/(log10(N/2) - log10(N/2-1));
            
            log10_dFk = dFk_log_int_calc(N/2-1, lattice_var, num_field);
            
            log10_dFk2 = dFk_log_int_calc(N/2, lattice_var, num_field);
            
            log10_dFk += (log10(distance) - log10(N/2-1))*(log10_dFk2 - log10_dFk)/(log10(N/2) - log10(N/2-1));
            
            
            
            field[0] = pow(10,log10_Fk);
            deriv[0] = pow(10,log10_dFk);
            
            
            //omega =  abs(pow(10,log10_dFk)/(pow(10, log10_Fk)/rescale_B));
        }
    }
    
}

double rand_uniform(void)
{
    std::random_device rnd;
    std::mt19937 mt( rnd() );
    std::uniform_real_distribution<> rand( 0, 1 );
    //   std::cout << rand(mt) << "\n";
    return (rand(mt));
}

void set_mode(double p2, double m2, double *field, double *deriv, int i, int real)
{
    double phase, amplitude, norm, omega;
    double re_f_left, im_f_left, re_f_right, im_f_right;
    
    // std::cout << "field i = " << i << std::endl;
    
    if (i==3){//Gravitational perturbation has to be calculated by case 1
        fluc_calc_switch = 1;
    }
    
    switch (fluc_calc_switch) {
        case 0: //Latice Easy case (when amplitudes of fluctuations are not predetermined)
            //   std::cout << "p2 = " << p2 << std::endl;
            //  std::cout << "m2 = " << m2 << std::endl;
            static int tachyonic=0; // Avoid printing the same error repeatedly
            
            if(p2+m2>0.) // Check to avoid floating point errors
                omega=sqrt(p2+m2); // Omega = Sqrt(p^2 + m^2)
            else
            {
                if(tachyonic==0)
                    printf("Warning: Tachyonic mode(s) may be initialized inaccurately\n");
                omega=sqrt(p2); // If p^2 + m^2 < 0 use m^2=0
                tachyonic=1;
            }
            
            if(omega>0.) {// Avoid dividing by zero
#if  dim==1
                norm = rescale_A*rescale_B*pow(L/(dx*dx),.5)/(exp(OSCSTART)*sqrt(4*M_PI*omega));
#elif  dim==2
                norm =  rescale_A*rescale_B*pow(L/(dx*dx),1)/(exp(OSCSTART)*sqrt(2*M_PI*omega));
#elif  dim==3
                norm =  rescale_A*rescale_B*pow(L/(dx*dx),1.5)/(exp(OSCSTART)*sqrt(2*omega));
#endif
            }else{
                norm = 0.0;
            }
            
            // std::cout << "norm = " << norm << std::endl;
            // std::cout << "|F_k| = " << 1/(exp(OSCSTART)*sqrt(2*omega)) << std::endl;
            //Amplitude = RMS amplitude x Rayleigh distributed random number
            // The same amplitude is used for left and right moving waves to generate standing waves. The extra 1/sqrt(2) normalizes the initial occupation number correctly.
            
            amplitude = norm*sqrt(log(1./rand_uniform()))*pow(p2,.75-(double)dim/4.)/sqrt(2);
            phase = 2*M_PI*rand_uniform();
            //   std::cout << "norm = " << norm << std::endl;
            // std::cout << "omega = " << omega << std::endl;
            //    std::cout << "norm/sqrt(2*omega) = " << norm/sqrt(2*omega) << std::endl;
            //    std::cout << "norm/sqrt(2*omega)*sqrt(log(1./rand_uniform())) = " << norm/sqrt(2*omega)*sqrt(log(1./rand_uniform())) << std::endl;
            //  std::cout << "amplitude = " << amplitude << std::endl;
            
            break;
            
        case 1: //when amplitudes of fluctuations are predetermined
            
            omega = abs(deriv[0]/(field[0]*rescale_B))*exp(OSCSTART);
            
#if  dim==1
            norm = (dx*dx/(L*sqrt(2*M_PI)))*field[0]*rescale_A*sqrt(pow(rescale_B*L/(dx*dx),3));
#elif  dim==2
            norm = (dx/sqrt(L*M_PI))*field[0]*rescale_A*sqrt(pow(rescale_B*L/(dx*dx),3));
#elif  dim==3
            
            norm =  field[0]*rescale_A*sqrt(pow(rescale_B*L/(dx*dx),3));
#endif
            
            // std::cout << "dx*dx/(L*sqrt(2*M_PI)) = " << dx*dx/(L*sqrt(2*M_PI)) << std::endl;
            //  std::cout << "field[0]*rescale_A*sqrt(pow(rescale_B*L/(dx*dx),3)) = " << field[0]*rescale_A*sqrt(pow(rescale_B*L/(dx*dx),3)) << std::endl;
            //   std::cout << "norm = " << norm << std::endl;
            //Amplitude = RMS amplitude x Rayleigh distributed random number
            // The same amplitude is used for left and right moving waves to generate standing waves. The extra 1/sqrt(2) normalizes the initial occupation number correctly.
            
            amplitude = norm*sqrt(log(1./rand_uniform()))*pow(p2,.75-(double)dim/4.)/sqrt(2);
            phase = 2*M_PI*rand_uniform();
            //   std::cout << "norm = " << norm << std::endl;
            // std::cout << "omega = " << omega << std::endl;
            //    std::cout << "norm/sqrt(2*omega) = " << norm/sqrt(2*omega) << std::endl;
            //    std::cout << "norm/sqrt(2*omega)*sqrt(log(1./rand_uniform())) = " << norm/sqrt(2*omega)*sqrt(log(1./rand_uniform())) << std::endl;
            //   std::cout << "amplitude = " << amplitude << std::endl;
            
            break;
            
        default:
            Logout( "Parameter 'fluc_calc_switch' must be 0 ~ 1. \n" );
            break;
    }
    
    
    
    //Left moving component
    re_f_left = amplitude * cos( phase );
    im_f_left = amplitude * sin( phase );
    //Right moving component
    phase = 2*M_PI*rand_uniform();
    //  std::cout << "phase2 " << phase/(2*M_PI) << std::endl;
    re_f_right = amplitude * cos( phase );
    im_f_right = amplitude * sin( phase );
    
    field[0] = re_f_left + re_f_right;
    field[1] = im_f_left + im_f_right;
    
    deriv[0] = (omega/exp(OSCSTART))*(im_f_left - im_f_right);
    deriv[1] = -(omega/exp(OSCSTART))*(re_f_left - re_f_right);
    
    if(real==1)
    {
        field[1]=0;
        deriv[1]=0;
    }
    
    //            std::cout << "field[0] = " << field[0] << std::endl;
    //            std::cout << "deriv[0] = " << deriv[0] << std::endl;
    
    return;
}




void initialize_perturb(double** f, double** df, double** lattice_var, double mass_sq[])
{
    
    
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
        
        // std::cout << "k = " << N/2 << ", f[i][1] = "<< f[i][1] << std::endl;
        
        fdf_calc( N/2 , lattice_var, &f[i][1], &df[i][1], i);
        
        //Set mode for lattice simulation
        set_mode(p2, mass_sq[i], &f[i][1], &df[i][1], i, 1);
        
        //Loop over gridpoints.
        for(int k = 1; k < N/2; k++ ){
            pz = k;
            p2 = dp2*pw2(pz);
            // std::cout << "k = " << k << ", f[i][" << 2*k << "] = " << f[i][2*k] << ", f[i][" << 2*k+1 << "] = " << f[i][2*k+1] << std::endl;
            // std::cout << "k = " << k << ", df[i][" << 2*k << "] = " << df[i][2*k] << ", df[i][" << 2*k+1 << "] = " << df[i][2*k+1] << std::endl;
            
            fdf_calc( k , lattice_var, &f[i][2*k], &df[i][2*k], i);
            
            //Set mode for lattice simulation
            set_mode(p2, mass_sq[i], &f[i][2*k], &df[i][2*k], i, 0);
        }
        // std::cout << "k = " << 0 << std::endl;
        //k=0
        f[i][0] = 0.;
        df[i][0] = 0.;
        
        std::cout << "before fourier transform" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        
        //transform from phase space to real space
        DFT_c2rD1( f[i] );
        DFT_c2rD1( df[i] );
        
        std::cout << "after fourier transform" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        
#elif dim==2
        for(int j = 0; j < N; j++ ){
            py = (j <= N/2 ? j : j-N);
            
            for(int k = 1; k < N/2; k++ ){
                pz = k;
                p2=dp2*(pw2(py)+pw2(pz));
                distance = sqrt(pw2(py)+pw2(pz));
                
                fdf_calc( distance , lattice_var, &f[i][j*N+2*k], &df[i][j*N+2*k], i);
                
                //Set mode for lattice simulation
                set_mode(p2, mass_sq[i], &f[i][j*N+2*k], &df[i][j*N+2*k], i, 0);
                
            }
            
            if(j > N/2)
                
            {
                jconj = N-j;
                
                //k=0
                p2 = dp2*pw2(py);
                fdf_calc( jconj , lattice_var, &f[i][j*N], &df[i][j*N], i);
                set_mode(p2, mass_sq[i], &f[i][j*N], &df[i][j*N], i, 0);
                f[i][jconj*N] = f[i][j*N];
                f[i][jconj*N+1] = -f[i][j*N+1];
                df[i][jconj*N] = df[i][j*N];
                df[i][jconj*N+1] = -df[i][j*N+1];
                //k=N/2
                distance = sqrt(pw2(py)+pw2(N/2));
                fdf_calc( distance , lattice_var, &fnyquist[2*j], &fdnyquist[2*j], i);
                p2 = dp2*(pw2(py)+pw2(N/2));
                set_mode(p2, mass_sq[i], &fnyquist[2*j], &fdnyquist[2*j], i, 0);
                fnyquist[2*jconj] = fnyquist[2*j];
                fnyquist[2*jconj+1] = -fnyquist[2*j+1];
                fdnyquist[2*jconj] = fdnyquist[2*j];
                fdnyquist[2*jconj+1] = -fdnyquist[2*j+1];
                
            }
            // The  "corners" of the lattice are set to real values
            else if(j == 0 || j == N/2){
                
                p2 = dp2*pw2(py); //k = 0
                if(p2 > 0.){
                    fdf_calc( N/2 , lattice_var, &f[i][j*N], &df[i][j*N], i);
                    set_mode(p2, mass_sq[i], &f[i][j*N], &df[i][j*N], i, 1);
                    
                }
                
                p2 = dp2*(pw2(py)+pw2(N/2));  //k = N/2
                distance = sqrt(pw2(py)+pw2(N/2));
                fdf_calc(distance, lattice_var, &fnyquist[2*j], &fdnyquist[2*j], i);
                set_mode(p2, mass_sq[i], &fnyquist[2*j], &fdnyquist[2*j], i, 1);
                
            }
        } //k = 0, j = 0
        f[i][0] = 0.;
        f[i][1] = 0.;
        df[i][0] = 0.;
        df[i][1] = 0.;
        
        std::cout << "before fourier transform" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        
        //transform from phase space to real space
        DFT_c2rD2( f[i] , fnyquist );
        DFT_c2rD2( df[i], fdnyquist);
        
        std::cout << "after fourier transform" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        
        
        
#elif dim==3
        for(int j = 0; j < N; j++){
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
                    fdf_calc(distance, lattice_var, &f[i][(j*N + k)*N + 2*l], &df[i][(j*N + k)*N + 2*l], i);
                    
                    set_mode(p2,mass_sq[i], &f[i][(j*N + k)*N + 2*l], &df[i][(j*N + k)*N + 2*l], i, 0);
                }
                
                if(k > N/2 || (j > N/2 && (k == 0 || k == N/2))){
                    kconj = (k == 0 ? 0 : N-k);
                    //l=0
                    p2 = dp2*(pw2(px)+pw2(py));
                    distance = sqrt(pw2(px)+pw2(py));
                    fdf_calc(distance, lattice_var, &f[i][(j*N + k)*N], &df[i][(j*N + k)*N], i);
                    set_mode(p2,mass_sq[i], &f[i][(j*N + k)*N], &df[i][(j*N + k)*N], i, 0);
                    f[i][(jconj*N + kconj)*N] = f[i][(j*N + k)*N];
                    f[i][(jconj*N + kconj)*N+1] = -f[i][(j*N + k)*N+1];
                    df[i][(jconj*N + kconj)*N] = df[i][(j*N + k)*N];
                    df[i][(jconj*N + kconj)*N+1] = -df[i][(j*N + k)*N+1];
                    //l=N/2
                    p2 = dp2*(pw2(px)+pw2(py)+pw2(N/2));
                    distance = sqrt(pw2(px)+pw2(py)+pw2(N/2));
                    fdf_calc(distance, lattice_var, &fnyquist[j][2*k], &fdnyquist[j][2*k], i);
                    set_mode(p2, mass_sq[i], &fnyquist[j][2*k], &fdnyquist[j][2*k], i, 0);
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
                        fdf_calc(distance, lattice_var, &f[i][(j*N + k)*N], &df[i][(j*N + k)*N], i);
                        set_mode(p2, mass_sq[i], &f[i][(j*N + k)*N], &df[i][(j*N + k)*N], i, 1);
                    }
                    
                    p2=dp2*(pw2(px)+pw2(py)+pw2(N/2)); //l=N/2
                    distance = sqrt(pw2(px)+pw2(py)+pw2(N/2));
                    fdf_calc(distance, lattice_var, &fnyquist[j][2*k], &fdnyquist[j][2*k], i);
                    set_mode(p2, mass_sq[i], &fnyquist[j][2*k], &fdnyquist[j][2*k], i, 1);
                }
            }
        }
        f[i][0] = 0.;
        f[i][1] = 0.;
        df[i][0] = 0.;
        df[i][1] = 0.;
        
        std::cout << "before fourier transform" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        
        //transform from phase space to real space
        DFT_c2rD3( f[i] , fnyquist );
        DFT_c2rD3( df[i], fdnyquist);
        
        std::cout << "after fourier transform" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        
        
        
#endif
        
    } //for( int i = 0; i < num_fields; ++i ){
    
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
    
}

void initialize( double**& f, double**& df, Field* field, double &radiation_pr, double**& lattice_var)
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
                
                
                //Initialize all elements to 0 first. If we don't do this, a very small (or in some cases very large) value will be assigned instead and it can cause errors. Also, for some reason std::fill() doesn't work for some cases, so I've avoided using it. //
                for( int j = 0; j < N; ++j ){
                    int idx = j;
                    
                    f[i][idx] = 0;
                    // std::cout << " f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
                    
                    df[i][idx]  = 0;
                    // std::cout << " df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
                    
                }
                
            }
            
            break;
        case 2:
            f[0] = new double [num_fields*N*N];
            df[0] = new double [num_fields*N*N];
            for( int i = 0; i < num_fields; ++i )
            {
                f[i] = f[0] + i*N*N;
                df[i] = df[0] + i*N*N;
                
                for( int j = 0; j < N; ++j ){
                    for( int k = 0; k < N; ++k ){
                        int idx = j*N + k;
                        f[i][idx] = 0;
                        df[i][idx] = 0;
                    }
                }
                
            }
            break;
        case 3:
            f[0] = new double [num_fields*N*N*N];
            df[0] = new double [num_fields*N*N*N];
            for( int i = 0; i < num_fields; ++i )
            {
                f[i] = f[0] + i*N*N*N;
                df[i] = df[0] + i*N*N*N;
                
                for( int j = 0; j < N; ++j ){
                    for( int k = 0; k < N; ++k ){
                        for( int l = 0; l < N; ++l ){
                            int idx = (j*N + k)*N + l;
                            f[i][idx] = 0;
                            df[i][idx] = 0;
                        }
                    }
                }
            }
            break;
        default:
            std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
            exit(1);
    }
    
    
    double* initial_field_values = new double [num_fields];
    double* initial_field_derivs = new double [num_fields];
    double* mass_sq = new double [num_fields]; // Effective mass squared of fields
    double radiation_var = 0;
    double hubble_parameter = 0;
    
    
    //lattice_var[knum_lattice][j] [j] specifies variables as follows:
    //    j=0-2  : zero modes of inflaton sigma, psi, and phi
    //    j=3-5  : log(a) derivatives of inflaton sigma', psi', phi'
    //    j=6    : energy density of radiation
    //    j=7-15 : mode functions of field perturbation: delta_{sigma,sigma}, delta_{sigma,psi}, delta{sigma,phi}, delta_{psi,sigma}, etc...
    //    j=16-24: log(a) derivatives of mode functions
    //    j=25-27: mode functions of gravitational potential perturbation: delta Phi_{sigma}, delta Phi_{psi}, delta Phi_{phi}
    //    j=28-30: log(a) derivatives of gravitational potential perturbation
    //    j=31-54: complex conjugate of [7]-[30]
    
//        for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){
//
//
//            for (int i=0;i<N_pert;i++) Logout("lattice_var[%d][%d] = %2.5e \n",lattice_loop, i , lattice_var[lattice_loop][i] );
//
//        }
    
    
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
    
    //Calculate effective mass
    field->effective_mass(mass_sq, initial_field_values);
    
    //Calculate Hubble parameter by averaging the hubble parameter for each wave number.
    for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){
        hubble_parameter +=
        Fri(lattice_var[lattice_loop][0],lattice_var[lattice_loop][1],lattice_var[lattice_loop][2],lattice_var[lattice_loop][3],lattice_var[lattice_loop][4],lattice_var[lattice_loop][5],lattice_var[lattice_loop][6]);
    }
    
    hubble_parameter/= (N/2);
    
    Logout("hubble_parameter = %2.5e \n",hubble_parameter);
    
    
    //Rescaling initial fields and their derivatives to  lattice program variables
    //Only necessary for the scalar fields
    for (int i=0; i< num_fields - 1; i++){
        
        initial_field_values[i] *= rescale_A;
        
        initial_field_derivs[i] = (rescale_A/rescale_B)*initial_field_derivs[i] + hubble_parameter*initial_field_values[i]/rescale_B;
        
        
        Logout("pr initial_field_values[%d] = %2.5e \n",i,initial_field_values[i]);
        Logout("pr initial_field_derivs[%d] = %2.5e \n",i,initial_field_derivs[i]);
    }
    
    Logout("FIXPSI = %2.5e \n",FIXPSI);
    
    //    double Kinetic1;
    //
    //    //Calculate Kinetic Energy
    //    Kinetic1 = (initial_field_derivs[0]*initial_field_derivs[0]+initial_field_derivs[1]*initial_field_derivs[1]+initial_field_derivs[2]*initial_field_derivs[2])/2;
    //
    //    Logout(" Kinetic1 = %2.5e \n",Kinetic1);
    
    //    double Kinetic2;
    //
    //    //Calculate Kinetic Energy
    //
    //    Kinetic2 = pw2(rescale_B/rescale_A)*(
    //pw2(initial_field_derivs[0]-hubble_parameter*initial_field_values[0]/rescale_B)
    //                +pw2(initial_field_derivs[1]-hubble_parameter*initial_field_values[1]/rescale_B)
    //                +pw2(initial_field_derivs[2]-hubble_parameter*initial_field_values[2]/rescale_B))/2;
    //
    //    Logout(" Kinetic2 = %2.5e \n",Kinetic2);
    
    //Rescaling hubble parameter to a lattice program variable
    Hinitial_pr = hubble_parameter/rescale_B;
    Logout("Hinitial_pr = %2.5e \n",Hinitial_pr);
    
    //Radiation zeromode
    for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++)
    {
        
        radiation_var += lattice_var[lattice_loop][6];
        
    }
    
    radiation_var /= (N/2);
    
    Logout("radiation_var = %2.5e \n", radiation_var);
    
    //Rescaling radiation to a lattice program variable
    radiation_pr = pow((rescale_A/rescale_B),2)*radiation_var;
    
    Logout("radiation_pr = %2.5e \n", radiation_pr);
    
    //------------------------------
    //  Initializing perturbations
    //------------------------------
    if (initialize_perturb_switch){
    initialize_perturb(f, df, lattice_var, mass_sq);
    }
    
    //------------------------------------------------
    //       Adding zeromodes for scalar fields
    //------------------------------------------------
    for( int i = 0; i < num_fields-1; ++i ){
        
#if  dim==1
        
        //#pragma omp parallel for simd schedule(static) num_threads(num_threads)
        
        for( int j = 0; j < N; ++j ){
            int idx = j;
            //std::cout << " f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
            f[i][idx] += initial_field_values[i];
           // std::cout << " f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
            
            //  std::cout << " df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
            
           // std::cout << " df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
            df[i][idx] += initial_field_derivs[i];
           // std::cout << " df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
        }
        
        std::cout << "after adding zeromode" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        
#elif dim==2
        
#pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            // #pragma omp simd
            for( int k = 0; k < N; ++k ){
                int idx = j*N + k;
              // std::cout << "f[" << i << "]["<< idx <<"]" << f[i][idx] << std::endl;
//                std::cout << "df[0]["<< idx <<"]" << df[i][idx] << std::endl;
                f[i][idx] += initial_field_values[i];
                df[i][idx] += initial_field_derivs[i];
//                std::cout << "f[0]["<< idx <<"]" << f[i][idx] << std::endl;
//                std::cout << "df[0]["<< idx <<"]" << df[i][idx] << std::endl;
                
            }
        }
        
//                 for( int j = 0; j < N; ++j ){
//                     for( int k = 0; k < N; ++k ){
//                         int idx = j*N + k;
//                         std::cout << " f[" << i << "][" << idx << "] = " << f[i][idx] << std::endl;
//                         std::cout << " df[" << i << "][" << idx << "] = " << df[i][idx] << std::endl;
//                     }
//                 }
        
        std::cout << "after adding zeromode" << std::endl;
        std::cout << "f[" << i << "][4] = " << f[i][4] << std::endl;
        std::cout << "df[" << i << "][4] = " <<  df[i][4] << std::endl;
        std::cout << "f[" << i << "][10] = " <<  f[i][10] << std::endl;
        std::cout << "df[" << i << "][10] = " <<  df[i][10] << std::endl;
        
        
#elif dim==3
        
        //Add zeromode
#pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N; ++k ){
                // #pragma omp simd
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
    }
    
   
    delete[] initial_field_values;
    delete[] initial_field_derivs;
    delete[] mass_sq;
    
    //    exit(1);
}

void finalize( double** f, double** df )
{
    delete [] f[0];
    delete [] df[0];
    delete [] f;
    delete [] df;
}


