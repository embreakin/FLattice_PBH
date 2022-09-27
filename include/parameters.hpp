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

//=================
//Non-lattice Range
//=================

extern std::string exist_filename_zero, new_filename_zero, exist_dirname_k, new_dirname_k, filename_k, exist_filename_sp, new_filename_sp, new_filename_spbfosc, exist_filename_spbfosc;

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

extern std::vector<int> knum_zero;

extern bool kanalyze_switch;
extern bool spectrum_switch;
extern bool spectrum_bfosc_switch;

extern DP kfrom_Mpc;
extern DP kto_Mpc;
extern int kinterval_knum;

//extern DP CN_par;                    //Potential paramater CN
//extern DP mu_par;                //Potential paramater mu
//extern DP Cv_par;                //Potential parameter Cv mu/4.
//extern DP M_par;                    //Potential paramater M
//extern DP m_par;                        //Potential paramater m
//extern DP n_par;                        //Potential paramater n
//extern DP g_par;                //Potential parameter g

//#define CAPTION "ver.06-04-05-02:00"
//#define Ck 2.626E-61            //k[MPl] for k=10000Mpc^{-1} // I think he means k[MPl] for k=10^-4 Mpc^{-1}
//#define Pi 3.14159265358979        //pi
//#define r2 1.41421356237310        //sqrt(2)
//#define GNOMAL 0                //decay rate 1
//#define GLARGE 6.0E-8            //decay rate 2 (if decay rate changes during the calculation)
//#define GLARGE2 6.0E-8            //decay rate 3 (if decay rate changes during the calculation)
//#define CN 0.04                    //Potential paramater CN
//#define mu 2.04E-3                //Potential paramater mu
//#define Cv 4.7E-4                //Potential parameter Cv mu/4.
//#define M 1.17                    //Potential paramater M
//#define m 2                        //Potential paramater m
//#define n 4                        //Potential paramater n
//#define g 2.0E-5                //Potential parameter g
//#define FIXPHI -0.459        //minimum value of Phi (set by hand according to zero-mode calculation) -0.472871 -0.475 -0.459
//#define SH -58.                //ln(a) at the end of calculation
//#define THRUNP -111                //ln(a) at which sigma and psi are fixed to the minimum
//#define THRLAST -61.5                //ln(a) at the beginning of oscillation of phi.-61.5
//#define OSCSTART -114.5        //ln(a) at the beginning of oscillation/
//#define dla 1.0E-5                //stepsize for fixed step RQ-method
//#define itvl 1000                //interval for output in fixed RQ-method
//#define Ini -131.8            //initial ln(a)  I1=0.3,Ini=-131.8
//#define I1 0.3                //initial value of sigma
//
//#define msigma sqrt(8*pow(mu,3)/M)      //effective mass of sigma

// mu 2.04E-3 Cv 4.7E-4  Ini -131.8  THRLAST -61.5 FIXPHI -0.458465


//#define CAPTION "ver.06-04-05-02:00"
//#define Ck 2.626E-61            //k[MPl] for k=10000Mpc^{-1} // I think he means k[MPl] for k=10^-4 Mpc^{-1}
//#define Pi 3.14159265358979        //pi
//#define r2 1.41421356237310        //sqrt(2)
//#define GNOMAL 0                //decay rate 1
//#define GLARGE 1.0E-10            //decay rate 2 (if decay rate changes during the calculation)
//#define GLARGE2 1.0E-10            //decay rate 3 (if decay rate changes during the calculation)
//#define CN 0.04                    //Potential paramater CN
//#define mu 2.04E-3                //Potential paramater mu
//#define Cv 4.7E-4                //Potential parameter Cv mu/4.
//#define M 1.17                    //Potential paramater M
//#define m 2                        //Potential paramater m
//#define n 4                        //Potential paramater n
//#define g 2.0E-5                //Potential parameter g
//#define FIXPHI -0.459        //minimum value of Phi (set by hand according to zero-mode calculation) -0.472871 -0.475 -0.459
//#define SH -57.                //ln(a) at the end of calculation
//#define THRUNP -111                //ln(a) at which sigma and psi are fixed to the minimum
//#define THRLAST -61.5                //ln(a) at the beginning of oscillation of phi.-61.5
//#define OSCSTART -114.5        //ln(a) at the beginning of oscillation/
//#define dla 1.0E-5                //stepsize for fixed step RQ-method
//#define itvl 1000                //interval for output in fixed RQ-method
//#define Ini -136.6           //initial ln(a)  I1=0.3,Ini=-131.8
//#define I1 0.4                //initial value of sigma
//
//#define msigma sqrt(8*pow(mu,3)/M)      //effective mass of sigma


//Takayama's Master Thesis
#define GNOMAL 0;//1.0E-11                 //decay rate 1
#define GLARGE 0;//1.0E-11            //decay rate 2 (if decay rate changes during the calculation)
#define GLARGE2 0;//1.0E-11            //decay rate 3 (if decay rate changes during the calculation)
#define CN_par 0.1                    //Potential paramater CN
#define Cv_par 6.4E-4                //Potential parameter Cv mu/4.
#define mu_par 2.7E-3                //Potential paramater mu
#define M_par 1.6                    //Potential paramater M
#define m_par 2.                        //Potential paramater m
#define n_par 10.                        //Potential paramater n
#define g_par 1.                //Potential parameter g
#define SH -70.                //ln(a) at the end of calculation
#define THRUNP -112                //ln(a) at which sigma and psi are fixed to the minimum
#define THRLAST -90                //ln(a) at the beginning of oscillation of phi.-61.5

#define dla 1.0E-1                //stepsize for fixed step RK-method
#define itvl 1000                //interval for output in fixed RK-method
#define Ini -135.6           //initial ln(a)  I1=0.3,Ini=-131.8
#define Init_sigma 0.8                //initial value of sigma
//10^-4[Mpc] corresponds to 0, 10^4[Mpc] to 800 in knum units
//////////////////////////////////////////
////The following only holds when m=2//////////////
//////////////////////////////////////////
#define msigma sqrt(8*pow(mu_par,3)/M_par)     //effective mass of sigma
//#define FIXPHI -2.31593e-03    //minimum value of Phi (set by hand according to zero-mode calculation) -0.472871 -0.475 -0.459
//#define FIXPSI  2*pow(mu_par*pow(M_par,m_par-1),1/m_par)  //PSI at the minimum



// FIG2
//#define CAPTION "ver.06-04-05-02:00"
//#define Ck 2.626E-61            //k[MPl] for k=10000Mpc^{-1} // I think he means Ck[MPl] for 10^-4 Mpc^{-1}
//#define Pi 3.14159265358979        //pi
//#define r2 1.41421356237310        //sqrt(2)
//#define GNOMAL 0                //decay rate 1
//#define GLARGE 1.0E-10          //decay rate 2 (if decay rate changes during the calculation)
//#define GLARGE2 1.0E-10            //decay rate 3 (if decay rate changes during the calculation)
//#define CN 0.04                    //Potential paramater CN
//#define mu 2.04E-3                //Potential paramater mu
//#define Cv 5.0E-4                //Potential parameter Cv mu/4.
//#define M 1.12                    //Potential paramater M
//#define m 2                        //Potential paramater m
//#define n 4                        //Potential paramater n
//#define g 2.0E-5                //Potential parameter g
//#define FIXPHI -0.473        //minimum value of Phi (set by hand according to zero-mode calculation) -0.472871 -0.475 -0.459
//#define SH -58.                //ln(a) at the end of calculation
//#define THRUNP -110                //ln(a) at which sigma and psi are fixed to the minimum
//#define THRLAST -62                //ln(a) at the beginning of oscillation of phi.-61.5
//#define OSCSTART -114.5        //ln(a) at the beginning of oscillation/
//#define dla 1.0E-5                //stepsize for fixed step RQ-method
//#define itvl 1000                //interval for output in fixed RQ-method
//#define Ini -131.8            //initial ln(a)  I1=0.3,Ini=-131.8
//#define I1 0.3                //initial value of sigma
//
//#define msigma sqrt(8*pow(mu,3)/M)      //effective mass of sigma

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

extern std::string exist_dirname_ed, new_dirname_ed, exist_dirname_f, new_dirname_f, exist_filename_status, new_filename_status;

inline double pw2(double x) { return (x*x);}

extern int N;

extern bool latticerange_switch;

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

#endif
