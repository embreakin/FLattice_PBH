//Doxygen
/**
* @file    utilities.cpp
* @brief      Utilities source file
* @author   Francis Otani
* @date
* @details
*/


#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
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
    static int par_set = 0; par_set++;
    
    if (par_set==1){
    const fs::path p = fs::current_path();
    std::cout << "Current path is " << p << std::endl<< std::endl;
    }
    
    //Path of the existing parameter set directory
    const fs::path existpath_par_set("../" + par_set_name_rm );
    
    //Check if the existing parameter set directory exists
    
    boost::system::error_code error;
        const bool result = fs::exists(existpath_par_set, error);
        if (!result || error) {
            
            if(par_set==1){
            std::cout << "Parameter set directory that you want to remove doesn't exist." << std::endl;
            std::cout << "Skip removing process..." << std::endl << std::endl;
            }
        }
        else
        {
            if(par_set==1){
            std::cout << "Parameter set directory that you want to remove exists." << std::endl << std::endl;
            }
                
                if (exist_par_set_rmall_switch && par_set==1)//This process only needs to be done once
                {
                   
                    
                    //Tries to remove all files in there. If it fails, it throws an error.
                    try {
                        fs::remove_all(existpath_par_set);
                        std::cout << "Removing existing parameter set directory " << par_set_name_rm << " ..." << std::endl << std::endl;
                        std::this_thread::sleep_for(std::chrono::seconds(1));
                        std::cout << "Existing parameter set directory " << par_set_name_rm << " removed successfully" << std::endl<< std::endl;
                    }
                    catch (fs::filesystem_error& ex) {
                        std::cout << ex.what() << std::endl;
                        std::cout << "Failed to remove existing parameter set directory " << par_set_name_rm << std::endl<< std::endl;
                        throw;
                    }
                }
                else
                {
                    if ( (k_switch_rm && exist_dir == exist_dirname_k) ||
                    (k_lattice_switch_rm && exist_dir == exist_dirname_k_lattice) )
                    {
                        
                        //Path of the existing data directory
                        const fs::path existpath("../" + par_set_name_rm + "/" + exist_dir );
                        
                        //Tries to remove all files in there. If it fails, it throws an error.
                        try {
                            fs::remove_all(existpath);
                            std::cout << "Removing existing data directory " << exist_dir << " in existing parameter set directory " << par_set_name_rm << " ..." << std::endl<< std::endl;
                            std::this_thread::sleep_for(std::chrono::seconds(1));
                            std::cout << "Existing data directory " << exist_dir << " in existing parameter set directory " << par_set_name_rm << " removed successfully" << std::endl<< std::endl;
                        }
                        catch (fs::filesystem_error& ex) {
                            std::cout << ex.what() << std::endl;
                            std::cout << "Failed to remove existing data directory " << exist_dir << " in existing parameter set directory " << par_set_name_rm << std::endl<< std::endl;
                            throw;
                        }
                    }
                }
        }
    
    
    
        
        
        
    //Creating Parameter Set Directory
    if(par_set==1)//This process only needs to be done once
        {
    //Path of the new parameter set directory
    const fs::path newpath_par_set("../" + par_set_name );
    
    //Tries to create the directory. If it fails, it throws an error.
    boost::system::error_code error_par_set;
    const bool result_par_set = fs::create_directory(newpath_par_set, error_par_set);
            
    int result_par_set_time = 0;
    std::cout << "Creating parameter set directory " << par_set_name << " ..." << std::endl<< std::endl;
           
        //Give it some time to create the directory (max 30s)
        while(!result_par_set){
        std::this_thread::sleep_for(std::chrono::seconds(1));
            result_par_set_time++;
            
            if(result_par_set_time >= 30)
            {
                break;
            }
        }
     
    if (!result_par_set || error_par_set) {
        std::cout << "Failed to create parameter set directory " << par_set_name << std::endl<< std::endl;
            }else
            {
                std::cout << "Parameter set directory " << par_set_name << " created successfully " << std::endl<< std::endl;
                std::cout << "Elapsed Time: " << result_par_set_time  << std::endl<< std::endl;
            }
        }
    
    std::cout << "Directory::: " << new_dir << std::endl << std::endl << std::endl;
    //Creating Directories in Parameter Set Directory
    // perturbation_switch must be turned on
    if(perturbation_switch)
    {
        
    if((!latticerange_switch && new_dir == new_dirname_k)//non-lattice simulation case
        ||(latticerange_switch && ((kfrom_knum < kfrom_knum_lattice) || (kstart_knum < kto_knum))
           ) // lattice simulation case
       //if ((kfrom_knum < kfrom_knum_lattice) || (kstart_knum < kto_knum)) holds then output all four directories
        ||(latticerange_switch && !((kfrom_knum < kfrom_knum_lattice) || (kstart_knum < kto_knum)) && !(new_dir == new_dirname_k)
       )//lattice simulation case
       // if ((kfrom_knum < kfrom_knum_lattice) && (kstart_knum < kto_knum)) doesn't hold then output all directories expect for new_dirname_k
       )
        {
            
            std::cout << "Directory:::::: " << new_dir << std::endl<< std::endl<< std::endl;
       
            //Path of the newly set directory in the new parameter set directory
            const fs::path newpath("../" + par_set_name + "/"  + new_dir );
            
            //Tries to create the directory. If it fails, it throws an error.
            boost::system::error_code error;
            const bool result = fs::create_directory(newpath, error);

            int result_time = 0;
            std::cout << "Creating directory " << new_dir << " in the newly created parameter set directory " << par_set_name << " ..." << std::endl<< std::endl;
            
            //Give it some time to create the directory (max 30s)
            while(!result){
            std::this_thread::sleep_for(std::chrono::seconds(1));
                result_time++;
                if(result_time >= 30)
                {
                    break;
                }
            }
                        
            if (!result || error) {
                std::cout << "Failed to create directory " << new_dir << " in the newly created parameter set directory " << par_set_name << std::endl<< std::endl;
            }else
            {
                std::cout << "Directory " << new_dir << " in parameter set directory " << par_set_name << " created successfully " << std::endl<< std::endl;
                std::cout << "Elapsed Time: " << result_time  << std::endl<< std::endl;
            }
        }
       }
    
}

void file_manage(const std::string exist_file)
{
    //Path of the existing status file
    const fs::path statuspath("../" + par_set_name_rm + "/" + exist_file);
    
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
// Elapsed Time Calculation
//--------------------------------
void time_calc(std::chrono::system_clock::time_point time_start, std::chrono::system_clock::time_point time_end, std::string time_name)
{
    int time_hours = std::chrono::duration_cast<std::chrono::hours>(time_end - time_start).count();
    int time_minutes = std::chrono::duration_cast<std::chrono::minutes>(time_end - time_start).count() - time_hours*60;
    int time_seconds = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start).count() - time_hours*60*60 - time_minutes*60;
    int time_milliseconds =
    std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count() - time_hours*60*60*1000 - time_minutes*60*1000 - time_seconds*1000;
    
    int time_days = time_hours / 24;
    
    time_hours = time_hours % 24;
    
    Logout("-----------------------------------------------------\n\n");
    Logout( "%s: %d d %d h %d m %d s %d ms\n\n",time_name.c_str(),time_days,time_hours,time_minutes,time_seconds,time_milliseconds);
    Logout("-----------------------------------------------------\n");
}

//--------------------------------
// Double Inflation Subroutines
//--------------------------------

//subroutine for zeromode output
void zeromode_output(const std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount){
    static int output_timecount = 0; ++output_timecount;
    std::stringstream ss;
    std::ofstream zeromode_output;
    ss << "../" << par_set_name << "/" << file;
    
    if( output_timecount == 1 )
    {
        zeromode_output.open(ss.str().c_str(),std::ios::out);
    }else{
        zeromode_output.open(ss.str().c_str(),std::ios::app);
    }
    
    DP H,la,rho,w,a,rho_rad,Kinetic, dda, pressure, epsilon;
    Vec_DP tr(N_zero);
    int i,j;
    double k_comoving = UC::knum_to_kMPl(k_target); //knum 200
    
    
    for (j=0;j<timecount;j++) {
    for (i=0;i<N_zero;i++) tr[i]=yp[i][j];
    la=xx[j];
    rho=rho_tot(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    H=Fri(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    w=log10(H);
    a=exp(la);
    rho_rad = tr[6];
    Kinetic = rho - V(tr[0],tr[1],tr[2]) - rho_rad;
    pressure = p_tot(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    dda = -a*(rho+3*pressure)/6;
    epsilon = 1 - dda/(pw2(H)*a);
    double k_horizoncrossing = a*H;
    double k_horizoncrossing2 = 0.7*a*H;
        
        if (epsilon > 1){
            static int epsilon_count=0; ++epsilon_count;
            
            if(epsilon_count == 1)
            {
                OSCSTART = la; //Set oscillation start to the ln(a) when epsilon=1
                Logout("-----------------------------------------------------\n\n");
                std::cout << "OSCSTART (Zeromode) = " << OSCSTART << std::endl ;
                std::cout << "Variable values at OSCSTART (Zeromode class)" << std::endl ;
                std::cout << "sigma = " << tr[0]<< ", psi = " << tr[1] << ", phi = " << tr[2]  << std::endl ;
                std::cout << "sigma_dot = " << tr[3]<< ", psi_dot = " << tr[4] << ", phi_dot = " << tr[5]  << std::endl ;
                std::cout << "rho_rad = " << tr[6] << std::endl ;
                std::cout << "dda = " << dda << std::endl ;
                std::cout << "Horizon crosing mode: k = aH = " << k_horizoncrossing << " [MPl], " << UC::kMPl_to_kMpc(k_horizoncrossing)  << " [Mpc^-1], knum = " << UC::kMPl_to_knum(k_horizoncrossing) << std::endl;
                 std::cout << "k = 0.7aH = " << k_horizoncrossing2 << " [MPl], " << UC::kMPl_to_kMpc(k_horizoncrossing2)  << " [Mpc^-1], knum = " << UC::kMPl_to_knum(k_horizoncrossing2) << std::endl << std::endl;
                Logout("-----------------------------------------------------\n\n");
                
            }
            
        }
    //output data
    zeromode_output << std::setw(6) << (la-BEGIN_EFOLD)*msigma/H << " " //t
    << std::setw(10) << la << " " //log(a)
    << std::setw(10) << tr[0] << " "  //buffer for yp_p[i][0] sigma
    << std::setw(20) << std::setprecision(20) << tr[1] << " " //buffer for yp_p[i][1] psi
    << std::setw(10) << tr[2] << " " //buffer for yp_p[i][1] phi
    << std::setw(10) << tr[3] << " "  //buffer for yp_p[i][0] sigma_dot
    << std::setw(20) << std::setprecision(20) << tr[4] << " " //buffer for yp_p[i][1] psi_dot
    << std::setw(10) << tr[5] << " " //buffer for yp_p[i][1] phi_dot
    << std::setw(10) << w << " "     //log10(H)
    << std::setw(10) << rho << " "
    << std::setw(10) << Kinetic << " "
    << std::setw(10) << V(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << rho_rad << " "
    << std::setw(10) << dda << " "
    << std::setw(10) << epsilon << " "
    << std::setw(10) << V_11(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << V_22(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << V_33(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << H << " "
    << std::setw(10) << V_1(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << V_2(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << V_3(tr[0],tr[1],tr[2]) << " "
    << std::setw(10) << pressure/rho << " "
    << std::setw(10) << k_comoving/a << " ";
        for(size_t loop = 0; loop < knum_zero.size(); loop++ ){
            if(loop ==  knum_zero.size()-1){
    zeromode_output << std::setw(10) << log10(UC::knum_to_kMPl(knum_zero[loop])/(a*H)) << "\n\n";
            }else{
                zeromode_output << std::setw(10) << log10(UC::knum_to_kMPl(knum_zero[loop])/(a*H)) << " ";
            };
        };
        
    };
}

//subroutine for k-analyze output
void kanalyze_output(const std::string dir, std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount, int knum, DP k_comoving ){
    //output results
    static int output_timecount = 0;
    static int latticerange_secondhalf_calc_count = 0;
    static int knum_static = knum;
    if(knum_static == knum){
        ++output_timecount;
    }else{

        if(knum > knum_static && latticerange_secondhalf_calc_count==0)
        {
            knum_static = knum;
            output_timecount = 1;
        }
        else
        {//This corresponds to the kanalyze output for latticerange_secondhalf_calc
            knum_static = knum;
            latticerange_secondhalf_calc_count++;
            if(latticerange_secondhalf_calc_count == N/2)
            {
                //This enables the function to output the upper range
                latticerange_secondhalf_calc_count=0;
            }
        }
        
    };
    
//    std::cout << knum_static << "\n";
//    std::cout << knum << "\n";
//    std::cout << output_timecount << "\n";
    
    std::stringstream ss;
    std::ofstream k_output;
    
    if(lattice_kmodes_switch ||
       (
        latticerange_switch &&  kfrom_knum_lattice <= knum &&  knum <= kto_knum_lattice)
        ){
        
        double kMpc, kMpc_int;
        kMpc = UC::kMPl_to_kMpc(k_comoving);
        kMpc_int = (int)round(100*kMpc); //Round to the second decimal place and make it an integer by multiplying by 100
        ss << "../" << par_set_name << "/" << dir << "/" << file << "_kMpc_" << std::setw(6) << std::setfill('0') << kMpc_int <<".txt";
    }else{
        ss << "../" << par_set_name << "/" << dir << "/" << file << "_knum_" << std::setw(4) << std::setfill('0') << knum <<".txt";
    }
    
    
    if( output_timecount == 1 )
    {
        k_output.open(ss.str().c_str(),std::ios::out);
    }else{
        k_output.open(ss.str().c_str(),std::ios::app);
    }
    
    DP H,la,rho,rhop,w,a,Pzeta,Pzeta_raw,PPot,P_sigma,P_psi,P_phi, P_sigma_raw, P_psi_raw, P_phi_raw, rho_rad, Kinetic, dda, pressure, epsilon;
    
    double k_horizoncrossing;
    double k_horizoncrossing2;
    static int epsilon_count=0;
    static int OSCSTART_count=0;
    double kMpc_max;
    static double OSCSTART_MAX = OSCSTART;
    static int kMpc_static = UC::kMPl_to_kMpc(k_comoving);

    Vec_DP zeta(6);
    Vec_DP tr(N_pert);
    int i,j;
    for (j=0;j<timecount;j++) {
        for (i=0;i<N_pert;i++) tr[i]=yp[i][j];
        la=xx[j];
//        if(j==0){std::cout <<  "la = " << la << "\n";
//            std::cout  <<  "tr[0] = " << tr[0] << "\n";
//            std::cout  << "tr[1] = " << tr[1] << "\n";
//            std::cout  << "tr[2] = " << tr[2] << "\n";
//        };
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
        P_sigma_raw = P_sigma;
        P_sigma = P_sigma/(2*M_PI*M_PI);
        P_sigma = log10(P_sigma) + 3*log10(k_comoving);
        
        P_psi = 0;
        for (i=0;i<3;i++) P_psi = P_psi + tr[i+10]*tr[i+10];
        for (i=0;i<3;i++) P_psi = P_psi + tr[i+34]*tr[i+34];
        P_psi_raw = P_psi;
        P_psi = P_psi/(2*M_PI*M_PI);
        P_psi = log10(P_psi) + 3*log10(k_comoving);
        
        P_phi = 0;
        for (i=0;i<3;i++) P_phi = P_phi + tr[i+13]*tr[i+13];
        for (i=0;i<3;i++) P_phi = P_phi + tr[i+37]*tr[i+37];
        P_phi_raw = P_phi;
        P_phi = P_phi/(2*M_PI*M_PI);
        P_phi = log10(P_phi) + 3*log10(k_comoving);
        
        rho_rad = tr[6];
        Kinetic = rho - V(tr[0],tr[1],tr[2]) - rho_rad;
        pressure = p_tot(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
        dda = -a*(rho+3*pressure)/6;
        epsilon = 1 - dda/(pw2(H)*a);
        
        
        if(latticerange_switch || lattice_kmodes_switch ){
            kMpc_max =  kto_Mpc_lattice;
        }else{
            kMpc_max = kto_Mpc;
        }
        

            //Set oscillation start to the ln(a) when epsilon=1
        if (epsilon > 1 && OSCSTART_count==0){

            if(kMpc_static == (int)UC::kMPl_to_kMpc(k_comoving)){
                ++epsilon_count;
            }else{
                kMpc_static = (int)UC::kMPl_to_kMpc(k_comoving);
                epsilon_count = 1;
            };
            
            if(epsilon_count == 1 )
            {
                
            OSCSTART = la;
            
            if (OSCSTART_MAX < OSCSTART)
            {
                OSCSTART_MAX = OSCSTART;
            }
                
//                std::cout << "OSCSTART_MAX = "  << OSCSTART_MAX  << std::endl;
//                std::cout << "OSCSTART = "  << OSCSTART  << std::endl;
            
            if( (int)UC::kMPl_to_kMpc(k_comoving) == (int)kMpc_max  ){
                
                
                OSCSTART = OSCSTART_MAX;
                k_horizoncrossing = a*H;
                k_horizoncrossing2 = 0.7*a*H;
                
                Logout("\n\n-----------------------------------------------------\n\n");
                std::cout << "OSCSTART (Perturbation) = " << OSCSTART  << std::endl ;
                std::cout << "Variable values at OSCSTART (Perturbation class)" << std::endl ;
                std::cout << "sigma = " << tr[0]<< ", psi = " << tr[1] << ", phi = " << tr[2]  << std::endl ;
                std::cout << "sigma_dot = " << tr[3]<< ", psi_dot = " << tr[4] << ", phi_dot = " << tr[5]  << std::endl ;
                std::cout << "rho_rad = " << tr[6] << std::endl ;
                std::cout << "dda = " << dda << std::endl ;
                std::cout << "Horizon crosing mode: k = aH = " << k_horizoncrossing << " [MPl], " << UC::kMPl_to_kMpc(k_horizoncrossing)  << " [Mpc^-1], knum = " << UC::kMPl_to_knum(k_horizoncrossing) << std::endl;
                std::cout << "k = 0.7aH = " << k_horizoncrossing2 << " [MPl], " << UC::kMPl_to_kMpc(k_horizoncrossing2)  << " [Mpc^-1], knum = " << UC::kMPl_to_knum(k_horizoncrossing2) << std::endl << std::endl;
                Logout("-----------------------------------------------------\n\n");
                
                OSCSTART_count++;
                
            }
            
        }
            
        }
        
        
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
        << std::setw(10) << log10(sqrt(P_sigma_raw)) << " "               //log10(dsigma)
        << std::setw(10) << log10(sqrt(P_psi_raw)) << " "                 //log10(dpsi)
        << std::setw(10) << log10(sqrt(P_phi_raw)) << " "                 //log10(dphi)
        << std::setw(10) << zeta[0]*zeta[0]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[1]*zeta[1]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[2]*zeta[2]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[3]*zeta[3]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[4]*zeta[4]/Pzeta_raw << " "                 //
        << std::setw(10) << zeta[5]*zeta[5]/Pzeta_raw << " "                 //
        << std::setw(10) << tr[0] << " "  //buffer for yp_p[i][0] sigma
        << std::setw(20) << std::setprecision(20) << tr[1] << " " //buffer for yp_p[i][1] psi
        << std::setw(10) << tr[2] << " " //buffer for yp_p[i][1] phi
        << std::setw(10) << tr[3] << " "  //buffer for yp_p[i][0] sigma_dot
        << std::setw(20) << std::setprecision(20) << tr[4] << " " //buffer for yp_p[i][1] psi_dot
        << std::setw(10) << tr[5] << " " //buffer for yp_p[i][1] phi_dot
        << std::setw(10) << rho << " "
        << std::setw(10) << Kinetic << " "
        << std::setw(10) << V(tr[0],tr[1],tr[2]) << " "
        << std::setw(10) << rho_rad << " "
        << std::setw(10) << dda << " "
        << std::setw(10) << epsilon << " "
        << "\n\n";
       
    };

}

//subroutine for spectrum output
void spectrum_output(const std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount, int knum, const DP k_comoving){
    
    static int output_timecount = 0;
    
//    std::cout << "1: output_timecount = " << output_timecount <<
//    " File Name: " << file << std::endl << std::endl;
    
    if(output_timecount == 0){
    sp_file_vec.push_back(file);
        ++output_timecount;
    }else
    {
        for (int fnum = 0; fnum < (int)sp_file_vec.size(); fnum++ )
                {
                    if(sp_file_vec[fnum] == file)//Already encountered this spectrum file
                    {
                        ++output_timecount;
                        break;
                    }else
                    {
                        if(fnum == (int)sp_file_vec.size() - 1 )//First time encountering this spectrum file
                        {
                            sp_file_vec.push_back(file);
                            output_timecount = 1;
                            break;
                        }
                        
                        
                    }
                        
                        
                    
                    
                }
    }
    
    std::stringstream ss;
    std::ofstream sp_output;
    ss << "../" << par_set_name << "/" << file;

//    std::cout << "2: output_timecount = " << output_timecount << " File Name: " << file << std::endl << std::endl;
    if( output_timecount == 1 )
    {
        sp_output.open(ss.str().c_str(),std::ios::out);
    }else{
        sp_output.open(ss.str().c_str(),std::ios::app);
    }
    
    DP H,a,la,rho,rhop,Pzeta,PPot,PSTR; //Pzeta_raw,,P_sigma,P_psi,P_phi,;

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
    << std::setw(10) << UC::kMPl_to_kMpc(k_comoving) << " "
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


//--------------------------------
// Lattice Simulation Subroutines
//--------------------------------


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
            {// zero-frequency (DC)
                in[idx][0] = f[idx]  ;
                in[idx][1] = 0  ;
                
            }else if(idx==N/2){//Nyquist frequency
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

void DFT_r2cD1( std::vector<double>& f, std::vector<double>& f_k)
{
    fftw_plan p;
    double* in;
    fftw_complex* out;
    size_t out_size;
    
    in = new double [N]();
    out_size = sizeof(fftw_complex) * (N/2+1);
    out  = (fftw_complex*)fftw_malloc( out_size );
    p = fftw_plan_dft_r2c_1d( N, in, out, FFTW_ESTIMATE );
    
    for( int j = 0; j < N; ++j ){
        int idx = j;
        in[idx] = f[idx];
    }
    
    fftw_execute(p);
    
    // Set output data
    for( int j = 0; j < N/2+1; ++j ){
        int idx = j;
        //f[][0] corresponds to Re(f_{i=0})
        //f[][2],f[][3] corresponds to the Re and Im of f_{i=1}
        // ...
        //f[][1] corresponds to Re(f_{i=N/2})
        if(idx==0)
        {// zero-frequency (DC)
             f_k[idx] = out[idx][0];
           // std::cout << " out[idx][0] = " << out[idx][0] << std::endl;
          //  std::cout << " out[idx][1] = " << out[idx][1] << std::endl;
        }else if(idx==N/2){//Nyquist frequency
            f_k[1] = out[idx][0];
        }else{
            f_k[2*idx] = out[idx][0];
            f_k[2*idx+1] = out[idx][1];
        }
        
    }
    
    if( p ) fftw_destroy_plan(p);
    if( out ) fftw_free(out);
    delete[] in;
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
void DFT_c2rD2( double* f, double* fnyquist )
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

void DFT_r2cD2(std::vector<double>& f, std::vector<double>& f_k, std::vector<double>& f_k_nyquist )
{
    fftw_plan p;
    double* in;
    fftw_complex *out;
    size_t out_size;
    
    in = new double [N*N]();
    out_size = sizeof(fftw_complex) * N * (N/2+1);
    out  = (fftw_complex*)fftw_malloc( out_size );
    p = fftw_plan_dft_r2c_2d( N, N, in, out, FFTW_ESTIMATE );
    
    for( int j = 0; j < N; ++j ){
        for( int k = 0; k < N; ++k ){
            int idx = j*N + k;
             in[idx] = f[idx];
        }
    }
    
    fftw_execute(p);
    
    // Set output data
    
    for( int j = 0; j < N; ++j ){
        for( int k = 0; k < N/2+1; ++k ){
            int idx = (N/2+1)*j + k;
            //  int idx = j*N + k;
            if(k == N/2){
                 f_k_nyquist[2*j] = out[idx][0];
                 f_k_nyquist[2*j+1] = out[idx][1];
            }else {
                f_k[j*N + 2*k] = out[idx][0];
                f_k[j*N + 2*k+1] = out[idx][1];
            }
        }
    }
    
    
    if( p ) fftw_destroy_plan(p);
    if( out ) fftw_free(out);
    delete[] in;
}




void DFT_c2rD3( double* f, double** fnyquist )
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

void DFT_r2cD3( std::vector<double>& f, std::vector<double>& f_k, std::vector<std::vector<double>>& f_k_nyquist )
{
    
    fftw_plan p;
    double* in;
    fftw_complex *out;
    size_t out_size;
    
    in = new double [N*N*N]();
    out_size = sizeof(fftw_complex) * N * N * (N/2+1);
    out  = (fftw_complex*)fftw_malloc( out_size );
    p = fftw_plan_dft_r2c_3d( N, N, N, in, out, FFTW_ESTIMATE );

  
    for( int j = 0; j < N; ++j ){
        for( int k = 0; k < N; ++k ){
            for( int l = 0; l < N; ++l ){
                int idx = (j*N + k)*N + l;
                in[idx] = f[idx];
            }
        }
    }

    fftw_execute(p);
    
    // Set output data
    
    
//    #pragma omp parallel for simd collapse(3) schedule( static ) num_threads( num_threads )
    for( int j = 0; j < N; ++j ){
        for( int k = 0; k < N; ++k ){
            for( int l = 0; l < N/2+1; ++l ){
                int idx = (j*N + k)*(N/2+1) + l;
                if(l == N/2){
                   f_k_nyquist[j][2*k] = out[idx][0];
                   f_k_nyquist[j][2*k+1] = out[idx][1];
                }
                else{
                    f_k[(j*N + k)*N + 2*l]  = out[idx][0];
                    f_k[(j*N + k)*N + 2*l+1] = out[idx][1];
                }
            }
        }
    }
    
    
    
    
    if( p ) fftw_destroy_plan(p);
    if( out ) fftw_free(out);
    delete[] in;
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
    
    double dx_Kpc;
    
    dx_Kpc = 1000*UC::xMPl_to_xMpc(dx/rescale_B);
	
    if( dim == 1 ){

		ss << "../" << par_set_name << "/" << dir_f << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".txt";
    	fout.open( ss.str().c_str() );	
 
    	for( int j = 0; j < N; j++ ){
			int idx = j;
			fout << idx*dx_Kpc << " " << f[idx] << std::endl;
		}
    }else{
    	ss << "../" << par_set_name << "/" << dir_f << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".vti";
    	fout.open( ss.str().c_str() );
    
  	  	fout << "<?xml version=\"1.0\"?>" << std::endl;
  		fout << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl <<std::endl;
        if(dim == 2){
            size = sizeof(double) * pow(N, 2);//8byte*pow(N,2)
            fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\" Origin=\"0 0 0\" Spacing=\"" << dx_Kpc << " " << dx_Kpc << " " << dx_Kpc << "\">" << std::endl;
            fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\">" << std::endl;
        }else if(dim == 3){
            size = sizeof(double) * pow(N, 3);
            fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 <<" 0 " << N-1 << "\" Origin=\"0 0 0\" Spacing=\"" << dx_Kpc << " " << dx_Kpc << " " << dx_Kpc << "\">" << std::endl;
            fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << N-1 << "\">" << std::endl;
        }else{
            std::cout << "dim must be 1~3" << std::endl;
            exit(1);
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
    double dx_Kpc;
    
    dx_Kpc = 1000*UC::xMPl_to_xMpc(dx/rescale_B);
    
    if( dim == 1 ){
        ss << "../" << par_set_name << "/" << dir_ed << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".txt";
        fout.open( ss.str().c_str() );
        
        for( int j = 0; j < N; j++ ){
            int idx = j;
            fout << idx*dx_Kpc << " " << f[idx] << std::endl;
        }
    }else{
        ss << "../" << par_set_name << "/" << dir_ed << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".vti";
        fout.open( ss.str().c_str() );
        
        fout << "<?xml version=\"1.0\"?>" << std::endl;
        fout << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl <<std::endl;
        if(dim == 2){
            size = sizeof(double) * pow(N, 2);//8byte*pow(N,2)
            fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\" Origin=\"0 0 0\" Spacing=\"" << dx_Kpc << " " << dx_Kpc << " " << dx_Kpc << "\">" << std::endl;
            fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\">" << std::endl;
        }else if(dim == 3){
            size = sizeof(double) * pow(N, 3);
            fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 <<" 0 " << N-1 << "\" Origin=\"0 0 0\" Spacing=\"" << dx_Kpc << " " << dx_Kpc << " " << dx_Kpc << "\">" << std::endl;
            fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << N-1 << "\">" << std::endl;
        }else{
            std::cout << "dim must be 1~3" << std::endl;
            exit(1);
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


void write_status( const std::string status_file, Field* field, LeapFrog* leapfrog, Energy* energy, double** f, double** df, double t )
{
	double a = leapfrog->a();
    double da = leapfrog->da();
    double efolds = leapfrog->efolds();
    
	std::ofstream ofs;
	
	if( t == t0 )
	{
		ofs.open( "../" + par_set_name +  "/" + status_file, std::ios::trunc );
		ofs << std::setw(3) << std::right << "  t_pr ";
        if( expansion ) {
            ofs << "  a ";
            ofs << "  efolds ";
        }
		for( int i = 0; i < num_fields; ++i ) ofs << "field_ave["  << i << "] ";
		for( int i = 0; i < num_fields; ++i ) ofs << "field_var["  << i << "] ";
        for( int i = 0; i < num_fields; ++i ) ofs << "field_deriv_ave["  << i << "] ";
        for( int i = 0; i < num_fields; ++i ) ofs << "field_deriv_var["  << i << "] ";
//        for( int i = 0; i < num_fields; ++i ) ofs << "energy_ave[" << i << "] ";
//        for( int i = 0; i < num_fields; ++i ) ofs << "energy_var[" << i << "] ";
        ofs << "total_energy_ave ";
        ofs << "total_energy_var ";
         ofs << "time_deriv_ave ";
         ofs << "gradient_ave ";
        ofs << "potential_ave ";
        ofs << "radiation ";
        ofs << "hubble ";
        ofs << "adotdot ";
        ofs << "energy_max ";
        ofs << "ddV_sigma ";
        ofs << "ddV_psi ";
        ofs << "ddV_phi ";
        ofs << "dV_sigma ";
        ofs << "dV_psi ";
        ofs << "dV_phi "<< std::endl;
	}
	else ofs.open( "../" + par_set_name + "/" + status_file, std::ios::app );
	 
	ofs << std::setw(3) << std::right << t << " ";
	if( expansion )
	{
        
		ofs << std::setw(3) << std::right << a << " ";
        ofs << std::setw(3) << std::right << efolds << " ";
        
        //f
        for( int i = 0; i < num_fields; ++i ){
            if(i == num_fields-1){
                ofs << std::showpos << std::scientific << std::setprecision(4) << field->average(f[i], i)/(a*a) << " "; //Gravitational Potential in Reduced Plank units
            }else{
                
//                if( i == 1){
//                    ofs << std::showpos << std::scientific << std::setprecision(4) <<  FIXPSI - field->average(f[i], i)/(rescale_A*a) << " "; //Scalar Field psi in Reduced Plank units
//                }else{
//                    ofs << std::showpos << std::scientific << std::setprecision(4) << field->average(f[i], i)/(rescale_A*a) << " "; //Scalar Fields other than psi in Reduced Plank units
//                      }
                
                 ofs << std::showpos << std::scientific << std::setprecision(4) << field->average(f[i], i)/(rescale_A*a) << " "; //Scalar Fields other than psi in Reduced Plank units
            }
        }
   
        for( int i = 0; i < num_fields; ++i ){
             if(i == num_fields-1){
            ofs << std::showpos << std::scientific << std::setprecision(4) << field->variance(f[i], i)/(a*a) << " ";//Variance of Gravitational Potential in Reduced Plank units
             }else{
                ofs << std::showpos << std::scientific << std::setprecision(4) << field->variance(f[i], i)/(rescale_A*a) << " ";//Variance of scalar fields in Reduced Plank units
             }
        }
        
        //df
      
        for( int i = 0; i < num_fields; ++i ){
            if(i == num_fields-1){
                ofs << std::showpos << std::scientific << std::setprecision(4) <<
                rescale_B*( field->average(df[i], i)  - 2*(da/a)*field->average(f[i], i) )/pow(a,3) << " ";//Derivative of Gravitational Potential in Reduced Plank units
            }else{
                
//                if( i == 1){
//                    ofs << std::showpos << std::scientific << std::setprecision(4) <<
//                    -(rescale_B/rescale_A)*( field->average(df[i], i)  - (da/a)*field->average(f[i], i) )/pow(a,2) << " "; //Derivative of Scalar Field psi in Reduced Plank units
//                   std::cout << "( field->average(df[i], i)  - (da/a)*field->average(f[i], i) ) = " << ( field->average(df[i], i)  - (da/a)*field->average(f[i], i) ) << std::endl;
//                   std::cout << "-(rescale_B/rescale_A)/pow(a,2) = " << -(rescale_B/rescale_A)/pow(a,2) << std::endl;
//
//                }else{
//                    ofs << std::showpos << std::scientific << std::setprecision(4) <<
//                    (rescale_B/rescale_A)*( field->average(df[i], i)  - (da/a)*field->average(f[i], i) )/pow(a,2) << " "; //Derivative of Scalar Fields other than psi in Reduced Plank units
                
                ofs << std::showpos << std::scientific << std::setprecision(4) <<
                                    (rescale_B/rescale_A)*( field->average(df[i], i)  - (da/a)*field->average(f[i], i) )/pow(a,2) << " "; //Derivative of Scalar Fields in Reduced Plank units
                //}
            }
        }
        
        
        for( int i = 0; i < num_fields; ++i ){ ofs << std::showpos << std::scientific << std::setprecision(4) << field->variance(df[i], i) << " ";//Variance of Derivative of fields all in Programming variables
        }
	}
	else
	{
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->average(f[i], i) << " ";//Reduced Plank units
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->variance(f[i], i) << " ";//Reduced Plank units
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->average(df[i], i) << " ";//Programming variable
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->variance(df[i], i) << " ";//Programming variable
	}
//    for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << energy->average(i) << " ";
//    for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << energy->variance(i) << " ";
	ofs << std::showpos << std::scientific << std::setprecision(4) << energy->total_average() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->total_variance() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->timederiv_average () << " ";
     ofs << std::showpos << std::scientific << std::setprecision(4) << energy->grad_average ()  << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->potential_average () << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->radiation () << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << leapfrog->hubble() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << leapfrog->adotdot() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->energy_max() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << field->ddV_lattice( f, 0, 0, a )*pow(rescale_B/a,2) << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << field->ddV_lattice( f, 1, 0, a )*pow(rescale_B/a,2) << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << field->ddV_lattice( f, 2, 0, a )*pow(rescale_B/a,2) << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << field->dV_lattice( f, 0, 0, a )*pow(rescale_B,2)/(pow(a,3)*rescale_A) << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << field->dV_lattice( f, 1, 0, a )*pow(rescale_B,2)/(pow(a,3)*rescale_A) << " ";
    //-field->dV_lattice( f, 1, 0, a )*pow(rescale_B,2)/(pow(a,3)*rescale_A) << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << field->dV_lattice( f, 2, 0, a )*pow(rescale_B,2)/(pow(a,3)*rescale_A) << std::endl;
    ofs << std::endl;
    
   
}

//const vector<vector<vector<double>>>&PS,
void kanalyze_output_lattice(const std::string dir, std::string file, Field* field, LeapFrog* leapfrog, double** f){
    //output results
    int  kMpc_int;//knum,
    double k_comoving, kMpc;

    DP H, la, w, a, a_bar;
    //rho,rhop,w,a,Pzeta,Pzeta_raw,PPot,P_sigma,P_psi,P_phi, P_sigma_raw, P_psi_raw, P_phi_raw, rho_rad, Kinetic, dda, pressure, epsilon, a_bar;
    
     a_bar = leapfrog->a();
     a = a_bar*exp(OSCSTART);
    
     H = leapfrog->hubble();
     la = log(a);
     w = log10(H);
    
    int j_start, j_end;
    static int time_count = 0;
    
    
    if (k_lattice_grid_min_MPl < kfrom_MPl_lattice)
    {
        j_start = outrange_num + 1;
    }else
    {
        j_start = 1;
    }
    
    j_end = N/2;
     
     for(int j=j_start; j < j_end + 1 ; j++) // j = 1,2...,N/2-1,N/2
        {
            
            k_comoving = (2*M_PI/(L/rescale_B))*j;// [MPl]
            
           // std::cout << "knum = " << knum << std::endl;
            kMpc = UC::kMPl_to_kMpc(k_comoving);
            kMpc_int = (int)round(100*kMpc); //Round to the second decimal place and make it an integer by multiplying by 100
            std::stringstream ss;
            std::ofstream k_output;
            
            
                 ss << "../" << par_set_name << "/" << dir << "/" << file << "_kMpc_" << std::setw(6) << std::setfill('0') << kMpc_int <<".txt";
                 
            
            
            if(k_lattice_startfromlattice_switch && time_count == 0){
                 k_output.open(ss.str().c_str(),std::ios::out);
            }else{
                k_output.open(ss.str().c_str(),std::ios::app);
            }
            
                 k_output <<  std::setw(20) << std::setprecision(20) << la << " "               //log(a)
                 << std::setw(10) << w << " "                        //log(H)
                 << std::setw(10) << 0 << " "//Pzeta << " "                    //log(Zeta)
            << std::setw(10) << log10(field->power_spectrum(f, 3, j)
                                      /pow(a_bar,4) ) << " "//PPot << " "                    //log(Gravitational Potential)
                 << std::setw(10) << 0 << " "//log10(rhop/rho) << " "        //log(p/rho)
                 << std::setw(10) << log10(H/a) << " "            //log(H/a)
                 << std::setw(10) << log10(k_comoving/(a*H)) << " "        //log(k/(a*H))
                                << std::setw(10) << log10(field->power_spectrum(f, 0, j)
                    /pow(rescale_A*a_bar,2) )<< " "               //log(P_sigma)
                                << std::setw(10) << log10(field->power_spectrum(f, 1, j)
                    /pow(rescale_A*a_bar,2)) << " "                 //log(P_psi)
                                << std::setw(10) << log10(field->power_spectrum(f, 2, j)
                    /pow(rescale_A*a_bar,2)) << " "
                         << std::setw(10) << 0 << " "//log10(sqrt(P_sigma_raw)) << " "               //log10(dsigma)
                         << std::setw(10) << 0 << " "//log10(sqrt(P_psi_raw)) << " "                 //log10(dpsi)
                         << std::setw(10) << 0 << " "//log10(sqrt(P_phi_raw)) << " "                 //log10(dphi)
                         << std::setw(10) << 0 << " "//zeta[0]*zeta[0]/Pzeta_raw << " "                 //
                         << std::setw(10) << 0 << " "//zeta[1]*zeta[1]/Pzeta_raw << " "                 //
                         << std::setw(10) << 0 << " "//zeta[2]*zeta[2]/Pzeta_raw << " "                 //
                         << std::setw(10) << 0 << " "//zeta[3]*zeta[3]/Pzeta_raw << " "                 //
                         << std::setw(10) << 0 << " "//zeta[4]*zeta[4]/Pzeta_raw << " "                 //
                         << std::setw(10) << 0 << " "//zeta[5]*zeta[5]/Pzeta_raw << " "                 //
                         << std::setw(10) << 0 << " "//tr[0] << " "  //buffer for yp_p[i][0] sigma
                         << std::setw(20) << 0 << " "//std::setprecision(20) << tr[1] << " " //buffer for yp_p[i][1] psi
                         << std::setw(10) << 0 << " "//tr[2] << " " //buffer for yp_p[i][1] phi
                         << std::setw(10) << 0 << " "//tr[3] << " "  //buffer for yp_p[i][0] sigma_dot
                         << std::setw(20) << 0 << " "//std::setprecision(20) << tr[4] << " " //buffer for yp_p[i][1] psi_dot
                         << std::setw(10) << 0 << " "//tr[5] << " " //buffer for yp_p[i][1] phi_dot
                         << std::setw(10) << 0 << " "//rho << " "
                         << std::setw(10) << 0 << " "//Kinetic << " "
                         << std::setw(10) << 0 << " "//V(tr[0],tr[1],tr[2]) << " "
                         << std::setw(10) << 0 << " "//rho_rad << " "
                         << std::setw(10) << 0 << " "//dda << " "
                         << std::setw(10) << 0 << " "//epsilon << " "
                        << "\n\n";
            
        }
    
    time_count++;
//    exit(1);
}

