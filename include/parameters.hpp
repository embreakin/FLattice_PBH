#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <string>
//#include <math.h>
#include "nr.h"
#include "equations.hpp"

extern std::string exist_dirname_k, new_dirname_k, filename_k, exist_filename_sp, new_filename_sp;
extern DP k_comoving, CNT;
extern DP Gamma1,Gamma2,Gamma3;
extern DP dxsav;
extern int kmax,timecount;
extern Vec_DP *xp_p;
extern Mat_DP *yp_p;
extern Mat_DP *delp_p;
extern const int timecount_max;

extern bool kanalyze_switch;
extern bool spectrum_switch;

extern DP kfromMpc;
extern DP ktoMpc;
extern int kinterval;

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
#define CAPTION "ver.06-04-05-02:00"
#define Ck 2.626E-61            //k[MPl] for k=10^-4 Mpc^{-1}
#define GNOMAL 0                //decay rate 1
#define GLARGE 1.0E-11            //decay rate 2 (if decay rate changes during the calculation)
#define GLARGE2 1.0E-11            //decay rate 3 (if decay rate changes during the calculation)
#define CN_par 0.1                    //Potential paramater CN
#define mu_par 2.7E-3                //Potential paramater mu
#define Cv_par 6.4E-4                //Potential parameter Cv mu/4.
#define M_par 1.6                    //Potential paramater M
#define m_par 2                        //Potential paramater m
#define n_par 10                        //Potential paramater n
#define g_par 1                //Potential parameter g
#define SH -71.                //ln(a) at the end of calculation
#define THRUNP -102                //ln(a) at which sigma and psi are fixed to the minimum
#define THRLAST -75                //ln(a) at the beginning of oscillation of phi.-61.5
#define OSCSTART -117.75        //ln(a) at the beginning of oscillation/
#define dla 1.0E-5                //stepsize for fixed step RQ-method
#define itvl 1000                //interval for output in fixed RQ-method
#define Ini -135.6           //initial ln(a)  I1=0.3,Ini=-131.8
#define Init_sigma 0.8                //initial value of sigma
//10^-4[Mpc] corresponds to 0, 10^4[Mpc] to 800 in knum units
#define msigma sqrt(8*pow(mu_par,3)/M_par)     //effective mass of sigma
#define FIXPHI -0.325    //minimum value of Phi (set by hand according to zero-mode calculation) -0.472871 -0.475 -0.459
#define FIXPSI  2*sqrt(mu_par*M_par)  //PSI at the minimum



// FIG2
//#define CAPTION "ver.06-04-05-02:00"
//#define Ck 2.626E-61            //k[MPl] for k=10000Mpc^{-1} // I think he means k[MPl] for k=10^-4 Mpc^{-1}
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

#endif
