#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <omp.h>
//#include <math.h>
#include "nr.h"
#include "equations.hpp"
#include "uc.hpp"

////Libraries for using JSON
#include <fstream>
#include <json.hpp>
using json = nlohmann::json;

//====================================
// Parameter Set defined in JSON file
//====================================

extern std::string par_set_name;




//=================
//Non-lattice Range
//=================

extern std::string exist_filename_zero, new_filename_zero, exist_dirname_k, new_dirname_k, filename_k, exist_filename_sp, new_filename_sp, new_filename_spbfosc, exist_filename_spbfosc;

extern bool k_switch;

extern bool zeromode_switch;
extern bool perturbation_switch;

//Array elements
extern const int N_zero, N_pert;
extern int k_target;

extern DP sigma_c;
extern DP FIXPHI;
extern DP FIXPSI;
extern DP CNT;
extern DP Gamma1,Gamma2,Gamma3;
extern DP OSCSTART;
extern double msigma;

extern std::vector<int> knum_zero;

extern bool kanalyze_switch;
extern bool spectrum_switch;
extern bool spectrum_bfosc_switch;

extern DP kfrom_Mpc;
extern DP kto_Mpc;
extern int kinterval_knum;
extern bool lattice_kmodes_switch;


//Choose Parameter Set
#define par_set_num 1
//Takayama's Master Thesis: 0
//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (a): 1
//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (b): 2
//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG2: 3
//Primordial seeds of SMBHs (peak at 2kMpc-1): 4
//Primordial seeds of SMBHs trying to shift the peak: 5

#if  par_set_num == 0

    //Takayama's Master Thesis
    #define GNORMAL 1.0E-11                //decay rate 1
    #define GLARGE 1.0E-11          //decay rate 2 (if decay rate changes during the calculation)
    #define GLARGE2 1.0E-11             //decay rate 3 (if decay rate changes during the calculation)
    #define CN_par 0.1                    //Potential paramater CN
    #define mu_par 2.7E-3                //Potential paramater mu
    #define Cv_par 6.4E-4                //Potential parameter Cv mu/4.
    #define M_par 1.6                    //Potential paramater M
    #define m_par 2.0                        //Potential paramater m
    #define n_par 10.0                        //Potential paramater n
    #define g_par 1.0                //Potential parameter g
    #define BEGIN_EFOLD -135.6            //initial ln(a)
    #define UNPERT_EFOLD -110.815                //ln(a) at which sigma and psi are fixed to the minimum
    #define NEWINF_END_EFOLD -90.0                //ln(a) at the beginning of oscillation of phi.-61.5
    #define END_EFOLD -70.0                //ln(a) at the end of calculation
    #define dla 1.0E-1                //stepsize for fixed step RQ-method
    #define itvl 1000.                //interval for output in fixed RQ-method
    #define sigma_init 0.8                //initial value of sigma

#elif par_set_num == 1

        //Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (a)
        #define GNORMAL 0.0                //decay rate 1
        #define GLARGE 1.0E-10          //decay rate 2 (if decay rate changes during the calculation)
        #define GLARGE2 1.0E-10            //decay rate 3 (if decay rate changes during the calculation)
        #define CN_par 0.04                    //Potential paramater CN
        #define mu_par 2.04E-3                //Potential paramater mu
        #define Cv_par 4.7E-4                //Potential parameter Cv mu/4.
        #define M_par 1.17                    //Potential paramater M
        #define m_par 2.0                        //Potential paramater m
        #define n_par 4.0                        //Potential paramater n
        #define g_par 2.0E-5                //Potential parameter g
        #define BEGIN_EFOLD -136.75            //initial ln(a)
        #define UNPERT_EFOLD -108.                //ln(a) at which sigma and psi are fixed to the minimum
        #define NEWINF_END_EFOLD -61.5                //ln(a) at the beginning of oscillation of phi.-61.5
        #define END_EFOLD -57.                //ln(a) at the end of calculation
        #define dla 1.0E-5                //stepsize for fixed step RQ-method
        #define itvl 1000.                //interval for output in fixed RQ-method
        #define sigma_init 0.4                //initial value of sigma

#elif par_set_num == 2
////Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (b)
#define GNOMAL 0.0                //decay rate 1
#define GLARGE 6.0E-8          //decay rate 2 (if decay rate changes during the calculation)
#define GLARGE2 6.0E-8            //decay rate 3 (if decay rate changes during the calculation)
#define CN_par 0.04                    //Potential paramater CN
#define mu_par 2.04E-3                //Potential paramater mu
#define Cv_par 4.7E-4                //Potential parameter Cv mu/4.
#define M_par 1.17                    //Potential paramater M
#define m_par 2.0                        //Potential paramater m
#define n_par 4.0                        //Potential paramater n
#define g_par 2.0E-5                //Potential parameter g
#define BEGIN_EFOLD -131.8            //initial ln(a)
#define UNPERT_EFOLD -111.0                //ln(a) at which sigma and psi are fixed to the minimum
#define NEWINF_END_EFOLD -61.5                //ln(a) at the beginning of oscillation of phi.-61.5
#define END_EFOLD -58.                //ln(a) at the end of calculation
#define dla 1.0E-5                //stepsize for fixed step RQ-method
#define itvl 1000.                //interval for output in fixed RQ-method
#define sigma_init 0.3                //initial value of sigma
        
#endif





 //Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG2
//#define GNOMAL 0.0                //decay rate 1
//#define GLARGE 1.0E-10          //decay rate 2 (if decay rate changes during the calculation)
//#define GLARGE2 1.0E-10            //decay rate 3 (if decay rate changes during the calculation)
//#define CN_par 0.04                    //Potential paramater CN
//#define mu_par 2.04E-3                //Potential paramater mu
//#define Cv_par 5.0E-4                //Potential parameter Cv mu/4.
//#define M_par 1.12                    //Potential paramater M
//#define m_par 2.0                        //Potential paramater m
//#define n_par 4.0                        //Potential paramater n
//#define g_par 2.0E-5                //Potential parameter g
////#define FIXPHI -0.473        //minimum value of Phi (set by hand according to zero-mode calculation) -0.472871 -0.475 -0.459
//#define SH -58.                //ln(a) at the end of calculation
//#define THRUNP -110.0                //ln(a) at which sigma and psi are fixed to the minimum
//#define THRLAST -62.0                //ln(a) at the beginning of oscillation of phi.-61.5
//#define OSCSTART -114.5        //ln(a) at the beginning of oscillation/
//#define dla 1.0E-5                //stepsize for fixed step RQ-method
//#define itvl 1000.0                //interval for output in fixed RQ-method
//#define Ini -131.8            //initial ln(a)  I1=0.3,Ini=-131.8
//#define sigma_init 0.3                //initial value of sigma
////////////////////////////////////////////
//////The following only holds for m=2//////////////
////////////////////////////////////////////
//#define msigma sqrt(8*pow(mu_par,3)/M_par)      //effective mass of sigma

//Primordial seeds of SMBHs (peak at 2kMpc-1)
//#define CAPTION "ver.06-04-05-02:00"
//#define Ck 2.626E-61            //k[MPl] for k=10000Mpc^{-1} // I think he means k[MPl] for k=10^-4 Mpc^{-1}
//#define Pi 3.14159265358979        //pi
//#define r2 1.41421356237310        //sqrt(2)
//#define GNOMAL 0                //decay rate 1
//#define GLARGE 2.5E-7            //decay rate 2 (if decay rate changes during the calculation)
//#define GLARGE2 2.5E-7            //decay rate 3 (if decay rate changes during the calculation)
//#define CN 0.05                    //Potential paramater CN
//#define mu 2.0E-3                //Potential paramater mu
//#define Cv mu/4                //Potential parameter Cv mu/4.
//#define M 0.59                    //Potential paramater M
//#define m 2                        //Potential paramater m
//#define n 4                        //Potential paramater n
//#define g 2.0E-5                //Potential parameter g
//#define FIXPHI -0.473        //minimum value of Phi (set by hand according to zero-mode calculation) -0.472871 -0.475 -0.459
//#define SH -67                //ln(a) at the end of calculation
//#define THRUNP -115                //ln(a) at which sigma and psi are fixed to the minimum
//#define THRLAST -72.5                //ln(a) at the beginning of oscillation of phi.-61.5
//#define OSCSTART -122        //ln(a) at the beginning of oscillation/
//#define dla 1.0E-5                //stepsize for fixed step RQ-method
//#define itvl 1000                //interval for output in fixed RQ-method
//#define Ini -150           //initial ln(a)  I1=0.3,Ini=-131.8
//#define I1 0.28               //initial value of sigma
//
//#define msigma sqrt(8*pow(mu,3)/M)     //effective mass of sigma

//Primordial seeds of SMBHs trying to shift the peak

//#define CAPTION "ver.06-04-05-02:00"
//#define Ck 2.626E-61            //k[MPl] for k=10000Mpc^{-1} // I think he means k[MPl] for k=10^-4 Mpc^{-1}
//#define Pi 3.14159265358979        //pi
//#define r2 1.41421356237310        //sqrt(2)
//#define GNOMAL 0                //decay rate 1
//#define GLARGE 5.0E-7            //decay rate 2 (if decay rate changes during the calculation)
//#define GLARGE2 5.0E-7            //decay rate 3 (if decay rate changes during the calculation)
//#define CN 0.05                    //Potential paramater CN
//#define mu 2.0E-3                //Potential paramater mu
//#define Cv mu/4                //Potential parameter Cv mu/4.
//#define M 0.59                    //Potential paramater M
//#define m 2                        //Potential paramater m
//#define n 4                        //Potential paramater n
//#define g 2.0E-5                //Potential parameter g
//#define FIXPHI -0.473        //minimum value of Phi (set by hand according to zero-mode calculation) -0.472871 -0.475 -0.459
//#define SH -62                //ln(a) at the end of calculation
//#define THRUNP -111                //ln(a) at which sigma and psi are fixed to the minimum
//#define THRLAST -80                //ln(a) at the beginning of oscillation of phi.-61.5
//#define OSCSTART -117        //ln(a) at the beginning of oscillation/
//#define dla 1.0E-5                //stepsize for fixed step RQ-method
//#define itvl 1000                //interval for output in fixed RQ-method
//#define Ini -150           //initial ln(a)  I1=0.3,Ini=-131.8
//#define I1 0.38               //initial value of sigma
//
//#define msigma sqrt(8*pow(mu,3)/M)     //effective mass of sigma


//=================
//Lattice Range
//=================


#define dim 1

extern std::string exist_dirname_ed, new_dirname_ed, exist_dirname_f, new_dirname_f, exist_filename_status, new_filename_status, exist_dirname_k_lattice, new_dirname_k_lattice, filename_k_lattice;


inline double pw2(double x) { return (x*x);}

extern int N;

extern bool k_lattice_switch;
extern bool k_lattice_startfromlattice_switch;

extern bool latticerange_switch;
extern bool initialize_perturb_switch;
extern int fluc_calc_switch;

extern double kfrom_Mpc_lattice;//[Mpc] Calculate from this k for lattice range
extern double kto_Mpc_lattice;//[Mpc] Calculate to this k for lattice range

extern double kfrom_MPl_lattice; //convert to MPl units
extern double kto_MPl_lattice; //convert to MPl units

//extern double sigma_initial;
//extern double rescale_A;
//extern double rescale_B;
//extern double L;
//
//extern double k_lattice_grid_min_pr;
extern double k_lattice_grid_min_MPl;
//
//extern double k_lattice_grid_max_pr;
extern double k_lattice_grid_max_MPl;

extern int outrange_num;
extern int latticerange_num;


extern int rnd;
extern int num_fields;
extern int num_threads;
extern const double initfield[];
extern const double initderivs[];
extern double Hinitial_pr;

extern int output_step;
extern int total_step;
extern int max_loop;
extern int st_output_step;
extern int st_max_loop;

extern double t0;
extern double dt;
//extern double dx;

extern const int expansion;
extern const int precision;

extern const int screen_latticeloop_number;

#endif
