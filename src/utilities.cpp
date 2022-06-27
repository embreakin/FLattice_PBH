#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <sstream>
#include <fftw3.h>
#include "parameters.hpp"
#include "utilities.hpp"
#include <boost/filesystem.hpp>
#include "equations.hpp"


//--------------------------------
// Directory & File Management
//--------------------------------

namespace fs = boost::filesystem;

void dir_manage(const std::string exist_dir, const std::string new_dir )
{
    
    //Path of the existing data directory
    const fs::path existpath("../" + exist_dir );
    
    //Tries to remove all files in there. If it fails, it throws an error.
    try {
        fs::remove_all(existpath);
    }
    catch (fs::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }
    
    //Path of the newly set directory
    const fs::path newpath("../" + new_dir );
    
    //Tries to create a directory. If it fails, it throws an error.
    boost::system::error_code error;
    const bool result = fs::create_directory(newpath, error);
    if (!result || error) {
        std::cout << "failed to create directory" << std::endl;
    }
    
}

void file_manage(const std::string exist_file)
{
    //Path of the existing status file
    const fs::path statuspath("../" + exist_file);
    
    //Tries to remove the file. If it fails, it throws an error.
    try {
        fs::remove(statuspath);
    }
    catch (fs::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }
    
}

//--------------------------------
// Double Inflation Subroutines
//--------------------------------

//subroutine for zeromode output
void zeromode_output(const std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount){
    static int output_timecount = 0; ++output_timecount;
    std::stringstream ss;
    std::ofstream zeromode_output;
    ss << "../" << file;
    
    if( output_timecount == 1 )
    {
        zeromode_output.open(ss.str().c_str(),std::ios::out);
    }else{
        zeromode_output.open(ss.str().c_str(),std::ios::app);
    }
    
    DP H,la,rho,rhop,w,a;
    Vec_DP tr(N_zero);
    int i,j;
    
    for (j=0;j<timecount;j++) {
    for (i=0;i<N_zero;i++) tr[i]=yp[i][j];
    la=xx[j];
    rho=rho_tot(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    rhop=rhoandp(tr[3],tr[4],tr[5],tr[6]);
    H=Fri(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    w=log10(H);
    a=exp(la);
    //output data
    zeromode_output << std::setw(6) << (la-Ini)*msigma/H << " " //t
    << std::setw(10) << la << " " //log(a)
    << std::setw(10) << tr[0] << " "  //buffer for yp_p[i][0] sigma
    << std::setw(20) << std::setprecision(20) << tr[1] << " " //buffer for yp_p[i][1] psi
    << std::setw(10) << tr[2] << " " //buffer for yp_p[i][1] phi
    << std::setw(10) << w << " "     //log10(H)
    << std::setw(10) << rho << " "
    << std::setw(10) << V_11(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << V_12(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << V_13(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << V_22(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << V_23(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << V_33(tr[0],tr[1],tr[2]) << " ";
        for(int loop = 0; loop < knum_zero.size(); loop++ ){
            if(loop ==  knum_zero.size()-1){
    zeromode_output << std::setw(10) << log10(UC::knum_to_kMPl(knum_zero[loop])/(a*H)) << "\n";
            }else{
                zeromode_output << std::setw(10) << log10(UC::knum_to_kMPl(knum_zero[loop])/(a*H)) << " ";
            };
        };
        
    };
}

//subroutine for k-analyze output
void kanalyze_output(const std::string dir, std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount, int knum, const DP k_comoving ){
    //output results
    static int output_timecount = 0;
    static int knum_static = knum;
    if(knum_static == knum){
        ++output_timecount;
    }else{
        knum_static = knum;
        output_timecount = 1;
    };
    
//    std::cout << knum_static << "\n";
//    std::cout << knum << "\n";
//    std::cout << output_timecount << "\n";
    
    std::stringstream ss;
    std::ofstream k_output;
    
    ss << "../" << dir << "/" << file << "_" << std::setw(4) << std::setfill('0') << knum <<".txt";
    
    if( output_timecount == 1 )
    {
        k_output.open(ss.str().c_str(),std::ios::out);
    }else{
        k_output.open(ss.str().c_str(),std::ios::app);
    }
    
    DP H,la,rho,rhop,w,a,Pzeta,Pzeta_raw,PPot,P_sigma,P_psi,P_phi;

    Vec_DP zeta(6);
    Vec_DP tr(N_pert);
    int i,j;
    for (j=0;j<timecount;j++) {
        for (i=0;i<N_pert;i++) tr[i]=yp[i][j];
        la=xx[j];
//        if(j==0){std::cout << la << "\n"; };
        rho=rho_tot(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
        rhop=rhoandp(tr[3],tr[4],tr[5],tr[6]);
        H=Fri(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
        w=log10(H);
        a=exp(la);
        for (i=0;i<3;i++) zeta[i]=2*rho*(tr[i+25] + tr[i+28]/H)/(rhop) + (1 + 2*k_comoving*k_comoving*rho/(9*a*a*H*H*rhop))*3*tr[i+25];
        for (i=0;i<3;i++) zeta[i+3]=2*rho*(tr[i+49] + tr[i+52]/H)/(rhop) + (1 + 2*k_comoving*k_comoving*rho/(9*a*a*H*H*rhop))*3*tr[i+49];
        
        Pzeta=0;
        for (i=0;i<6;i++) Pzeta = Pzeta + zeta[i]*zeta[i];
        Pzeta_raw = Pzeta;
        Pzeta = Pzeta/(2*M_PI*M_PI*9);
        Pzeta = log10(Pzeta) + 3*log10(k_comoving);
        
        PPot = 0;
        for (i=0;i<3;i++) PPot = PPot + tr[i+25]*tr[i+25];
        for (i=0;i<3;i++) PPot = PPot + tr[i+49]*tr[i+49];
        PPot = PPot/(2*M_PI*M_PI);
        PPot = log10(PPot) + 3*log10(k_comoving);
        
        P_sigma = 0;
        for (i=0;i<3;i++) P_sigma = P_sigma + tr[i+7]*tr[i+7];
        for (i=0;i<3;i++) P_sigma = P_sigma + tr[i+31]*tr[i+31];
        P_sigma = P_sigma/(2*M_PI*M_PI);
        P_sigma = log10(P_sigma) + 3*log10(k_comoving);
        
        P_psi = 0;
        for (i=0;i<3;i++) P_psi = P_psi + tr[i+10]*tr[i+10];
        for (i=0;i<3;i++) P_psi = P_psi + tr[i+34]*tr[i+34];
        P_psi = P_psi/(2*M_PI*M_PI);
        P_psi = log10(P_psi) + 3*log10(k_comoving);
        
        P_phi = 0;
        for (i=0;i<3;i++) P_phi = P_phi + tr[i+13]*tr[i+13];
        for (i=0;i<3;i++) P_phi = P_phi + tr[i+37]*tr[i+37];
        P_phi = P_phi/(2*M_PI*M_PI);
        P_phi = log10(P_phi) + 3*log10(k_comoving);
        
//        if(j==0){std::cout << la << "\n"; };
        k_output << std::setw(6) << la << " "               //log(a)
        << std::setw(10) << w << " "                        //log(H)
        << std::setw(10) << Pzeta << " "                    //log(Zeta)
        << std::setw(10) << PPot << " "                    //log(Gravitational Potential)
        << std::setw(10) << log10(rhop/rho) << " "        //log(p/rho)
        << std::setw(10) << log10(H/a) << " "            //log(H/a)
        << std::setw(10) << log10(k_comoving/(a*H)) << " "        //log(k/(a*H))
        << std::setw(10) << P_sigma << " "               //log(P_sigma)
        << std::setw(10) << P_psi << " "                 //log(P_psi)
        << std::setw(10) << P_phi << " "                 //log(P_phi)
        << std::setw(10) << zeta[0]*zeta[0]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[1]*zeta[1]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[2]*zeta[2]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[3]*zeta[3]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[4]*zeta[4]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[5]*zeta[5]/Pzeta_raw << " "                 //
        << "\n";
    };

}

//subroutine for spectrum output before oscillation period
void spectrum_bfosc_output(const std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount, int knum, const DP k_comoving){
    static int output_timecount = 0; ++output_timecount;
    std::stringstream ss;
    std::ofstream sp_output;
    ss << "../" << file;
    
    if( output_timecount == 1 )
    {
        sp_output.open(ss.str().c_str(),std::ios::out);
    }else{
        sp_output.open(ss.str().c_str(),std::ios::app);
    }
    
    DP H,a,la,rho,rhop,Pzeta,Pzeta_raw,PPot,P_sigma,P_psi,P_phi,PSTR;

    Vec_DP zeta(6);
    Vec_DP tr(N_pert);

    int i;
    
    la=xx[timecount-1];
    a=exp(la);

    for (i=0;i<N_pert;i++) tr[i] = yp[i][timecount-1];
    rho=rho_tot(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    rhop=rhoandp(tr[3],tr[4],tr[5],tr[6]);
    H=Fri(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    
    for (i=0;i<3;i++) zeta[i]=2*rho*(tr[i+25] + tr[i+28]/H)/(rhop) + (1 + 2*k_comoving*k_comoving*rho/(9*a*a*H*H*rhop))*3*tr[i+25];
    for (i=0;i<3;i++) zeta[i+3]=2*rho*(tr[i+49] + tr[i+52]/H)/(rhop) + (1 + 2*k_comoving*k_comoving*rho/(9*a*a*H*H*rhop))*3*tr[i+49];
    Pzeta=0;
    for (i=0;i<6;i++) Pzeta = Pzeta + zeta[i]*zeta[i];
    Pzeta = Pzeta/(2*M_PI*M_PI*9);
    Pzeta = log10(Pzeta) + 3*log10(k_comoving);
    PPot = 0;
    for (i=0;i<3;i++) PPot = PPot + tr[i+25]*tr[i+25];
    for (i=0;i<3;i++) PPot = PPot + tr[i+49]*tr[i+49];
    PSTR = PPot;
    PPot = PPot/(2*M_PI*M_PI);
    PPot = log10(PPot) + 3*log10(k_comoving);
    //output log(Gravitational Potential) and log(Zeta), along with k and knum.
    sp_output << std::setprecision (10) << std::setw(10) << k_comoving << " "
    << std::setw(5) << knum << " "
    << std::setw(10) << UC::knum_to_kMpc(knum) << " "
    << std::setw(10) << PPot << " "
    << std::setw(10) << Pzeta <<  " "
    << std::setw(10) << PSTR*sqrt(k_comoving*k_comoving*k_comoving) << " ";
    for (i=7;i<N_pert;i++)
    {
        if(i == N_pert -1){
            sp_output << std::setw(10) << tr[i] << "\n";
            }else{
        sp_output<< std::setw(10) << tr[i] << " " ;
            }        
    }
}

//subroutine for spectrum output
void spectrum_output(const std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount, int knum, const DP k_comoving){
    static int output_timecount = 0; ++output_timecount;
    std::stringstream ss;
    std::ofstream sp_output;
    ss << "../" << file;
    
    if( output_timecount == 1 )
    {
        sp_output.open(ss.str().c_str(),std::ios::out);
    }else{
        sp_output.open(ss.str().c_str(),std::ios::app);
    }
    
    DP H,a,la,rho,rhop,Pzeta,Pzeta_raw,PPot,P_sigma,P_psi,P_phi,PSTR;
    
    Vec_DP zeta(6);
    Vec_DP tr(N_pert);
    int i;
    
    la=xx[timecount-1];
    a=exp(la);
    
    for (i=0;i<N_pert;i++) tr[i] = yp[i][timecount-1];

    rho=rho_tot(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    rhop=rhoandp(tr[3],tr[4],tr[5],tr[6]);
    H=Fri(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    
    for (i=0;i<3;i++) zeta[i]=2*rho*(tr[i+25] + tr[i+28]/H)/(rhop) + (1 + 2*k_comoving*k_comoving*rho/(9*a*a*H*H*rhop))*3*tr[i+25];
    for (i=0;i<3;i++) zeta[i+3]=2*rho*(tr[i+49] + tr[i+52]/H)/(rhop) + (1 + 2*k_comoving*k_comoving*rho/(9*a*a*H*H*rhop))*3*tr[i+49];
    Pzeta=0;
    for (i=0;i<6;i++) Pzeta = Pzeta + zeta[i]*zeta[i];
    Pzeta = Pzeta/(2*M_PI*M_PI*9);
    Pzeta = log10(Pzeta) + 3*log10(k_comoving);
    PPot = 0;
    for (i=0;i<3;i++) PPot = PPot + tr[i+25]*tr[i+25];
    for (i=0;i<3;i++) PPot = PPot + tr[i+49]*tr[i+49];
    PSTR = PPot;
    PPot = PPot/(2*M_PI*M_PI);
    PPot = log10(PPot) + 3*log10(k_comoving);
    //output log(Gravitational Potential) and log(Zeta), along with k and knum.
    sp_output << std::setprecision (10) << std::setw(10) << k_comoving << " " << std::setw(5) << knum << " " << std::setw(10) << UC::knum_to_kMpc(knum) << " " << std::setw(10) << PPot << " " << std::setw(10) << Pzeta <<  " " << std::setw(10) << PSTR*sqrt(k_comoving*k_comoving*k_comoving) << "\n" << std::flush;
}
//--------------------------------
// Lattice Simulation Subroutines
//--------------------------------

double rand_uniform(void)
{
    std::random_device rnd;
    std::mt19937 mt( rnd() );
    std::uniform_real_distribution<> rand( 0, 1 );
 //   std::cout << rand(mt) << "\n";
    return (rand(mt));
}

void set_mode(double p2, double omega, double *field, double *deriv, int real)
{

    double phase, amplitude, rms_amplitude;
    double re_f_left, im_f_left, re_f_right, im_f_right;
#if  dim==1
    static double norm = rescale_A*rescale_B*pow(L/(dx*dx),.5)/(exp(OSCSTART)*sqrt(4*M_PI));
#elif  dim==2
    static double norm =  rescale_A*rescale_B*pow(L/(dx*dx),1)/(exp(OSCSTART)*sqrt(2*M_PI));
#elif  dim==3
    static double norm =  rescale_A*rescale_B*pow(L/(dx*dx),1.5)/(exp(OSCSTART)*sqrt(2));
#endif

        //Amplitude = RMS amplitude x Rayleigh distributed random number
        // The same amplitude is used for left and right moving waves to generate standing waves. The extra 1/sqrt(2) normalizes the initial occupation number correctly.
        
    amplitude = norm/sqrt(2*omega)*sqrt(log(1./rand_uniform()))*pow(p2,.75-(double)dim/4.);
    phase = 2*M_PI*rand_uniform();
//   std::cout << "norm = " << norm << std::endl;
//    std::cout << "omega = " << omega << std::endl;
//    std::cout << "norm/sqrt(2*omega) = " << norm/sqrt(2*omega) << std::endl;
//    std::cout << "norm/sqrt(2*omega)*sqrt(log(1./rand_uniform())) = " << norm/sqrt(2*omega)*sqrt(log(1./rand_uniform())) << std::endl;
//    std::cout << "amplitude = " << amplitude << std::endl;

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

    deriv[0] = omega*(im_f_left - im_f_right);
    deriv[1] = -omega*(re_f_left - re_f_right);

    if(real==1)
    {
        field[1]=0;
        deriv[1]=0;
    }
    return;
   }


void DFT_c2rD1( double* f)
{
    fftw_plan p;
    double* out;
    fftw_complex *in;
    size_t in_size;
    
    in_size = sizeof(fftw_complex) * (N/2+1);
    in  = (fftw_complex*)fftw_malloc( in_size );
    out = new double [N]();
    p = fftw_plan_dft_c2r_1d( N, in, out, FFTW_ESTIMATE );
    
        for( int j = 0; j < N/2+1; ++j ){
            int idx = j;
            if(idx==0)
            {
                in[idx][0] = f[idx]  ;
                in[idx][1] = 0  ;
                
            }else if(idx==N/2){
                in[idx][0] = f[1]  ;
                in[idx][1] = 0  ;
            }else{
                in[idx][0] = f[2*idx]  ;
                in[idx][1] = f[2*idx+1]  ;
            }
            
        }
    
    fftw_execute(p);
    
    // Set output data
    for( int j = 0; j < N; ++j ){
        int idx = j;
        f[idx] = out[idx]/N;
    }
    
    if( p ) fftw_destroy_plan(p);
    if( in ) fftw_free(in);
    delete[] out;
}
/*
void DFT_c2rD2d( double* df,double* fdnyquist )
{
    fftw_plan p;
    double* out;
    fftw_complex *in;
    size_t in_size;
    
    in_size = sizeof(fftw_complex) * N * (N/2+1);
    in  = (fftw_complex*)fftw_malloc( in_size );
    out = new double [N*N]();
    p = fftw_plan_dft_c2r_2d( N, N, in, out, FFTW_ESTIMATE );
    
    
    
    for( int j = 0; j < N; ++j ){
        //#pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int k = 0; k < N/2+1; ++k ){
            int idx = (N/2+1)*j + k;
            // int idx = j*N + k;
            if(k == N/2){
                in[idx][0] = (double)fdnyquist[2*j] ;
                in[idx][1] = (double)fdnyquist[2*j+1] ;
               // std::cout << "in[" << idx << "][0] = " << in[idx][0] << std::endl;
               // std::cout << "in[" << idx << "][1] = " << in[idx][1] << std::endl;
            }else {
                in[idx][0] =  (double)df[j*N + 2*k] ;
                in[idx][1] = (double)df[j*N + 2*k+1] ;
               // std::cout << "in[" << idx << "][0] = " << in[idx][0] << std::endl;
               // std::cout << "in[" << idx << "][1] = " << in[idx][1] << std::endl;
            }
        }
    }
    fftw_execute(p);
    
    // Set output data
    
    for( int j = 0; j < N; ++j ){
        //#pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int k = 0; k < N; ++k ){
            int idx = j*N + k;
            std::cout << "out[" << idx << "] = " << out[idx] << std::endl;
            df[idx] = out[idx]/(N*N);
        }
    }
    
    if( p ) fftw_destroy_plan(p);
    if( in ) fftw_free(in);
    delete[] out;
}
*/
void DFT_c2rD2( double* f,double* fnyquist )
{
    fftw_plan p;
    double* out;
    fftw_complex *in;
    size_t in_size;
    
    in_size = sizeof(fftw_complex) * N * (N/2+1);
    in  = (fftw_complex*)fftw_malloc( in_size );
    out = new double [N*N]();
    p = fftw_plan_dft_c2r_2d( N, N, in, out, FFTW_ESTIMATE );
    
    
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N/2+1; ++k ){
                int idx = (N/2+1)*j + k;
              //  int idx = j*N + k;
                if(k == N/2){
                    in[idx][0] =  fnyquist[2*j] ;
                    in[idx][1] =  fnyquist[2*j+1] ;
                }else {
                    in[idx][0] =  f[j*N + 2*k]; 
                    in[idx][1] =  f[j*N + 2*k+1] ;
                }
            }
        }
        fftw_execute(p);
        
        // Set output data

     for( int j = 0; j < N; ++j ){
        for( int k = 0; k < N; ++k ){
            int idx = j*N + k;
            f[idx] = out[idx]/(N*N);
        }
     }
    
    if( p ) fftw_destroy_plan(p);
    if( in ) fftw_free(in);
    delete[] out;
}




void DFT_c2rD3( double* f,double** fnyquist )
{
    fftw_plan p;
    double* out;
    fftw_complex *in;
    size_t in_size;
    
    in_size = sizeof(fftw_complex) * N * N * (N/2+1);
    in  = (fftw_complex*)fftw_malloc( in_size );
    out = new double [N*N*N]();
    p = fftw_plan_dft_c2r_3d( N, N, N, in, out, FFTW_ESTIMATE );
    
   
//#pragma omp parallel for simd collapse(3) schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N; ++k ){
                for( int l = 0; l < N/2+1; ++l ){
                    int idx = (j*N + k)*(N/2+1) + l;
                    if(l == N/2){
                        in[idx][0] = fnyquist[j][2*k] ;
                        in[idx][1] = fnyquist[j][2*k+1] ;
                    }
                    else{
                        in[idx][0] =  f[(j*N + k)*N + 2*l] ;
                        in[idx][1] =  f[(j*N + k)*N + 2*l+1] ;
                    }
                }
            }
        }
        fftw_execute(p);
        
        // Set output data


         for( int j = 0; j < N; ++j ){
        for( int k = 0; k < N; ++k ){
            for( int l = 0; l < N; ++l ){
                int idx = (j*N + k)*N + l;
                f[idx] = out[idx]/(N*N*N);
            }
        }
    }
    
    if( p ) fftw_destroy_plan(p);
    if( in ) fftw_free(in);
    delete[] out;
}

/*
void DFT_c2r( double** f,double* fnyquist )
{
	  fftw_plan p;
		double* out;
    fftw_complex *in;
		size_t in_size;
	  
		switch( dim )
		{
			case 1:
				in_size = sizeof(fftw_complex) * (N/2+1);
				in  = (fftw_complex*)fftw_malloc( in_size );
				out = new double [N]();
				p = fftw_plan_dft_c2r_1d( N, in, out, FFTW_ESTIMATE );
				break;
			case 2:
				in_size = sizeof(fftw_complex) * N * (N/2+1);
				in  = (fftw_complex*)fftw_malloc( in_size );
				out = new double [N*N]();
				p = fftw_plan_dft_c2r_2d( N, N, in, out, FFTW_ESTIMATE );
				break;
			case 3:
				in_size = sizeof(fftw_complex) * N * N * (N/2+1);
				in  = (fftw_complex*)fftw_malloc( in_size );
				out = new double [N*N*N]();
				p = fftw_plan_dft_c2r_3d( N, N, N, in, out, FFTW_ESTIMATE );
				break;
		}

		for( int i = 0; i < num_fields; ++i )
		{
				
				// Create input data
				switch( dim ){
					case 1:
						#pragma omp parallel for schedule( static ) num_threads( num_threads )
						for( int j = 0; j < N/2+1; ++j ){
							int idx = j;
							if(idx==0)
                            {
                                in[idx][0] = f[i][idx]  ;
                                in[idx][1] = 0  ;
                                
                            }else if(idx==N/2){
                                in[idx][0] = f[i][1]  ;
                                in[idx][1] = 0  ;
                            }else{
                                in[idx][0] = f[i][2*idx]  ;
                                in[idx][1] = f[i][2*idx+1]  ;
                            }
							
						}
						break;
					case 2:
						#pragma omp parallel for schedule( static ) num_threads( num_threads )
						for( int j = 0; j < N; ++j ){
							for( int k = 0; k < N/2+1; k=k+2 ){
									int idx = j*N + k;
                                if(k == N/2){
                                    in[idx][0] = fnyquist[2*j] ;
                                    in[idx][1] = fnyquist[2*j+1] ;
                                }else {
                                    in[idx][0] =  f[i][idx] ;
                                    in[idx][1] =  f[i][idx+1] ;
                                }
							}
						}
						break;
					case 3:
						#pragma omp parallel for schedule( static ) num_threads( num_threads )
						for( int j = 0; j < N; ++j ){
							for( int k = 0; k < N; ++k ){
									for( int l = 0; l < N/2+1; l=l+2 ){
											int idx = (j*N + k)*N + l;
                                        if(l == N/2){
                                        in[idx][0] = fnyquist[j][2*k] ;
                                        in[idx][1] = fnyquist[j][2*k+1] ;
                                        }
                                        else{
                                            in[idx][0] =  f[i][idx] ;
                                            in[idx][1] =  f[i][idx+1] ;
                                        }
                                        }
							}
						}
						break;
				}
		
				fftw_execute(p);
	
				// Set output data
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            switch( dim ){
							case 1:
                int idx = j;
                f[i][idx] = out[idx]/N;
								break;
            	case 2:
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f[i][idx] = out[idx]/(N*N);
                }
								break;
            	case 3:
                for( int k = 0; k < N; ++k ){
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f[i][idx] = out[idx]/(N*N*N);
                    }
                }
								break;
						}
        }
    }

		if( p ) fftw_destroy_plan(p);
		if( in ) fftw_free(in);
		delete[] out;
}
*/

void write_VTK_f( const std::string dir_f, double* f, std::string str, int loop )
{
   // double a = leapfrog->a();
	unsigned int size;
	std::stringstream ss;
	std::ofstream fout;
    
    
	
	if( dim == 1 ){
		ss << "../" << dir_f << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".txt";
    	fout.open( ss.str().c_str() );	
 
    	for( int j = 0; j < N; j++ ){
			int idx = j;
			fout << idx*dx << " " << f[idx] << std::endl;
		}
    }else{
    	ss << "../" << dir_f << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".vti";
    	fout.open( ss.str().c_str() );
    
  	  	fout << "<?xml version=\"1.0\"?>" << std::endl;
  		fout << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl <<std::endl;
    	switch( dim ){
    		case 2:
				size = sizeof(double) * pow(N, 2);//8byte*pow(N,2)
				fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
				fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\">" << std::endl;
    			break;
			case 3:
				size = sizeof(double) * pow(N, 3);
				fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 <<" 0 " << N-1 << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
				fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << N-1 << "\">" << std::endl;
				break;
    	}
    	fout << "<PointData Scalars=\"field\">" << std::endl;
    	fout << "<DataArray type=\"Float64\" Name=\"field\" format=\"appended\" offset=\"0\" />" << std::endl;
    	fout << "</PointData>" << std::endl;
	    fout << "<CellData>" << std::endl;
	    fout << "</CellData>" << std::endl;
	    fout << "</Piece>" << std::endl << std::endl;
	    
	    fout << "</ImageData>" << std::endl << std::endl;
	    fout << "<AppendedData encoding=\"raw\">" << std::endl;
	    fout << "_" ;
	    fout.close();
	    
	    fout.open( ss.str().c_str(), std::ios::binary | std::ios::app);
	    fout.write( (char*) &size, sizeof(unsigned int) );//4byte
       
           fout.write( (char*) f, size );
	    fout.close();	
	    fout.open( ss.str().c_str(), std::ios::app);
	    fout << std::endl << "</AppendedData>" << std::endl;
		fout << "</VTKFile>" ;
	}
	fout.close();
}

void write_VTK_ed( const std::string dir_ed, double* f, std::string str, int loop )
{
  //  double a = leapfrog->a();
   // std::cout << "aaa = " <<  a << std::endl;
    unsigned int size;
    std::stringstream ss;
    std::ofstream fout;
    
    
    
    if( dim == 1 ){
        ss << "../" << dir_ed << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".txt";
        fout.open( ss.str().c_str() );
        
        for( int j = 0; j < N; j++ ){
            int idx = j;
            fout << idx*dx << " " << f[idx] << std::endl;
        }
    }else{
        ss << "../" << dir_ed << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".vti";
        fout.open( ss.str().c_str() );
        
        fout << "<?xml version=\"1.0\"?>" << std::endl;
        fout << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl <<std::endl;
        switch( dim ){
            case 2:
                size = sizeof(double) * pow(N, 2);//8byte*pow(N,2)
                fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
                fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\">" << std::endl;
                break;
            case 3:
                size = sizeof(double) * pow(N, 3);
                fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 <<" 0 " << N-1 << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
                fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << N-1 << "\">" << std::endl;
                break;
        }
        fout << "<PointData Scalars=\"energy\">" << std::endl;
        fout << "<DataArray type=\"Float64\" Name=\"energy density\" format=\"appended\" offset=\"0\" />" << std::endl;
        fout << "</PointData>" << std::endl;
        fout << "<CellData>" << std::endl;
        fout << "</CellData>" << std::endl;
        fout << "</Piece>" << std::endl << std::endl;
        
        fout << "</ImageData>" << std::endl << std::endl;
        fout << "<AppendedData encoding=\"raw\">" << std::endl;
        fout << "_" ;
        fout.close();
        
        fout.open( ss.str().c_str(), std::ios::binary | std::ios::app);
        fout.write( (char*) &size, sizeof(unsigned int) );//4byte
        
            fout.write( (char*) f, size );
        
        fout.close();
        fout.open( ss.str().c_str(), std::ios::app);
        fout << std::endl << "</AppendedData>" << std::endl;
        fout << "</VTKFile>" ;
    }
    fout.close();
}


void write_status( const std::string status_file, Field* field, LeapFrog* leapfrog, Energy* energy, double** f, double t )
{
	double a = leapfrog->a();
	std::ofstream ofs;
	
	if( t == t0 )
	{
		ofs.open( "../" + status_file, std::ios::trunc );

		ofs << std::setw(3) << std::right << "  t ";
		if( expansion ) ofs << "  a ";
		for( int i = 0; i < num_fields; ++i ) ofs << "field_ave["  << i << "] ";
		for( int i = 0; i < num_fields; ++i ) ofs << "field_var["  << i << "] ";
        for( int i = 0; i < num_fields; ++i ) ofs << "field_deriv_ave["  << i << "] ";
        for( int i = 0; i < num_fields; ++i ) ofs << "field_deriv_var["  << i << "] ";
		for( int i = 0; i < num_fields; ++i ) ofs << "energy_ave[" << i << "] ";
        for( int i = 0; i < num_fields; ++i ) ofs << "energy_var[" << i << "] ";
        ofs << "total_energy_ave ";
         ofs << "time_deriv_ave ";
         ofs << "gradient_ave ";
        ofs << "potential_ave ";
        ofs << "hubble ";
        ofs << "adotdot ";
        ofs << "energy_max" << std::endl;
	}
	else ofs.open( "../" + status_file, std::ios::app );
	
	ofs << std::setw(3) << std::right << t << " ";
	if( expansion )
	{
		ofs << std::setw(3) << std::right << a << " ";
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->f_average(f[i], i)/a << " "; //Reduced Plank units
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->f_variance(f[i], i)/a << " ";//Reduced Plank units
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->df_average(f[i], i) << " ";//Programming variable
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->df_variance(f[i], i) << " ";//Programming variable
	}
	else
	{
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->f_average(f[i], i) << " ";//Reduced Plank units
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->f_variance(f[i], i) << " ";//Reduced Plank units
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->df_average(f[i], i) << " ";//Programming variable
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->df_variance(f[i], i) << " ";//Programming variable
	}
	for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << energy->average(i) << " ";
    for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << energy->variance(i) << " ";
	ofs << std::showpos << std::scientific << std::setprecision(4) << energy->total_average() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->timederiv_average () << " ";
     ofs << std::showpos << std::scientific << std::setprecision(4) << energy->grad_average ()  << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->potential_average () << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << leapfrog->hubble() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << leapfrog->adotdot() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->energy_max() << std::endl;
    
    
}



