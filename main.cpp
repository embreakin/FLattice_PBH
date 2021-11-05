#include <chrono>// Measuring elapsed time
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "lattice.hpp"
#include "nr.h"
#include "equations.hpp"
#include "utilities.hpp"
#include "parameters.hpp"

//-------------------
// Declare variables





Vec_DP *xp_p;
Mat_DP *yp_p;
Mat_DP *delp_p;


//xp[i] stores integration variable (log(a)) for each steps[i].
//delp_p[i][j] stores variables for each steps[i]. [j] specifies variables as follows:
//    j=0-2  : zero modes of inflaton sigma, psi, and phi
//    j=3-5  : log(a) derivatives of inflaton sigma', psi', phi'
//    j=6    : energy density of radiation
//    j=7-15 : mode functions of field perturbation: delta_{sigma,sigma}, delta_{sigma,psi}, delta{sigma,phi}, delta_{psi,sigma}, etc...
//    j=16-24: log(a) derivatives of mode functions
//    j=25-27: mode functions of gravitational potential perturbation: delta Phi_{sigma}, delta Phi_{psi}, delta Phi_{psi}
//    j=28-30: log(a) derivatives of gravitational potential perturbation
//    j=31-54: complex conjugate of [7]-[30]
//unp[i] stores zero mode variables(delp_p[i][0-6]
//tr[j]: buffer for delp_p[i][j]

int main(int argc, char *argv[])//comand line arguments: #1: knum
{
    
    std::chrono::system_clock::time_point  time_start, time_end;
    time_start = std::chrono::system_clock::now(); // Start measuring elapsed time
    
    //Initial Data Directory Management
    dir_manage(exist_dirname_k, new_dirname_k);
    file_manage(exist_filename_sp);
   
    const int N1=55, N2=7;
    int timecount, mid,i,j,nbad,nok,p,kil,kfrom,kto,knum;
    DP a,xmid,H,rhop,theta,la, dxsav;
    DP epsSHI=1.0e-9,epsosc=1.0e-10,epsnew=1.0e-10,epslast=1.0e-10, eps=1.0E-6;
    //eps parameters set the allowed error in adapive step size RQ-method.
    DP h1=1.0E-4,h2=1.0E-6,hmin=0.0,xbegin=Ini,xend=SH,xint=-110;
    //h parameters set the initial trial step size in adaptive step size RQ-method.
    Vec_DP ystart(N2),delstart(N1);
    Vec_DP tr(N1);
    DP w,rho,rho1,rho2,rho3;
    
    Vec_DP zeta(6),dens(4),unp(N2),yout(N2),dydx(N2);
    Vec_DP adia(3),iso(3),fields(3),numdens(3),term(6);
    xp_p=new Vec_DP(timecount_max);
    delp_p=new Mat_DP(N1,timecount_max);
    Vec_DP &xp=*xp_p;
    Mat_DP &delp=*delp_p;
    kfrom = (int)( 100*(log(kfromMpc)/log(10.0)+4) ); // convert to original knum units
    kto = (int)( 100*(log(ktoMpc)/log(10.0)+4) ); // convert to original knum units
    int calcpercentage = 0;
    
    Logout( "kfrom =  %d, kto = %d, kinterval = %d \n", kfrom, kto, kinterval);
    Logout( "Calculation 0%% Complete\n");

    //Start for loop of k
    for (knum=kfrom;knum<=kto;knum=knum+kinterval){
        
    calcpercentage =  round( ( (knum - kfrom)/kinterval + 1 )*100 / ( (kto - kfrom)/kinterval + 1 ) );
    
//    k_comoving=Ck*exp(0.023025851*knum);            //knum -> k
        k_comoving=Ck*pow(10,knum/100);
    a=exp(xbegin);
    p=itvl;
    //set initial conditions
    unp[0] = Init_sigma;
    unp[1] = exp(log(Pow(4,m_par)/(2*m_par*(2*m_par-1)))/(2*(m_par-1)))*M_par*exp(log(mu_par/unp[0])/(m_par-1));
    unp[2] = -Cv_par*Cv_par*unp[0]/(mu_par*mu_par);
    unp[3] = 0;
    unp[4] = 0;
    unp[5] = 0;
    unp[6] = 0;
    //set decay rates
    Gamma1=GLARGE;
    Gamma2=GLARGE;
    Gamma3=GNOMAL;
    
    H=Fri(unp[0],unp[1],unp[2],unp[3],unp[4],unp[5],unp[6]);
    unpert(la,unp,dydx);
    //First, solve the evolution of zero mode, until the specified mode crosses the horizon.
    //Calculation for first order perturbation begins slightly before the horizon crossing.
    //Runge-Qutta method with constant step size is used.
    for (la=xbegin;k_comoving/(a*H)>500;la=la) {    //Here, zero-mode is solved until k=a*H*5000.
        if (p==itvl){
            rho=rho_tot(unp[0],unp[1],unp[2],unp[3],unp[4],unp[5],unp[6]);
            w=log10(H);
            p=0;
        };
        p++;
        H=Fri(unp[0],unp[1],unp[2],unp[3],unp[4],unp[5],unp[6]);
        unpert(la,unp,dydx);
        NR::rk4(unp,dydx,la,dla,yout,unpert);
        for (j=0;j<N2;j++){
            unp[j]=yout[j];
        };
        la=la+dla;
        a=exp(la);
    };
    //Hereafter, evolution of perturbation is solved.
    //Setting initial conditions for perturbation
    xmid=la;
    for (i=0;i<N2;i++) tr[i]=unp[i];
    
    a=exp(xmid);
    H=Fri(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
    rhop=rhoandp(tr[3],tr[4],tr[5],tr[6]);
    theta=M_PI/4;
    dxsav=(xend-xmid)/5000.0;
    for (i=N2;i<N1;i++) tr[i]=0;
    for (i=0;i<3;i++) tr[(i*4)+7]=cos(theta)/(sqrt(2*k_comoving)*a);
    for (i=0;i<3;i++) tr[(i*4)+31]=sin(theta)/(sqrt(2*k_comoving)*a);
    for (i=0;i<3;i++) tr[(i*4)+16]=(-1)*H*tr[(i*4)+7] + k_comoving*tr[(i*4)+31]/a;
    for (i=0;i<3;i++) tr[(i*4)+40]=(-1)*H*tr[(i*4)+31] - k_comoving*tr[(i*4)+7]/a;
    for (i=0;i<3;i++) tr[i+25]=0.5*a*a*(tr[3]*tr[i+16] + tr[4]*tr[i+19] + tr[5]*tr[i+22] + 3*H*(tr[3]*tr[i+7] + tr[4]*tr[i+10] + tr[5]*tr[i+13]) + V_1(tr[0],tr[1],tr[2])*tr[i+7] + V_2(tr[0],tr[1],tr[2])*tr[i+10] + V_3(tr[0],tr[1],tr[2])*tr[i+13])/(k_comoving*k_comoving);
    for (i=0;i<3;i++) tr[i+49]=0.5*a*a*(tr[3]*tr[i+40] + tr[4]*tr[i+43] + tr[5]*tr[i+46] + 3*H*(tr[3]*tr[i+31] + tr[4]*tr[i+34] + tr[5]*tr[i+37]) + V_1(tr[0],tr[1],tr[2])*tr[i+31] + V_2(tr[0],tr[1],tr[2])*tr[i+34] + V_3(tr[0],tr[1],tr[2])*tr[i+37])/(k_comoving*k_comoving);
    for (i=0;i<3;i++) tr[i+25]=tr[i+25]/(1 - (a*a*(tr[3]*tr[3] + tr[4]*tr[4] + tr[5]*tr[5])/(2*k_comoving*k_comoving)));
    for (i=0;i<3;i++) tr[i+49]=tr[i+49]/(1 - (a*a*(tr[3]*tr[3] + tr[4]*tr[4] + tr[5]*tr[5])/(2*k_comoving*k_comoving)));
    for (i=0;i<3;i++) tr[i+28]=(-1)*H*tr[i+25] - 0.5*(tr[3]*tr[i+7] + tr[4]*tr[i+10] + tr[5]*tr[i+13]);
    for (i=0;i<3;i++) tr[i+52]=(-1)*H*tr[i+49] - 0.5*(tr[3]*tr[i+31] + tr[4]*tr[i+34] + tr[5]*tr[i+37]);
    for (i=0;i<N1;i++) delstart[i]=tr[i];
    //Full evolution equations including first order perturbation of three inflaton fields are solved until oscillatory phase begins.
    //If horizon entering is later, this step is skipped.
    if(xmid < OSCSTART){
       NR::odeintpert(delstart,xmid,OSCSTART,epsSHI,h2,hmin,nok,nbad,timecount,dxsav,full,NR::rkqs);
        if(kanalyze_switch){
        kanalyze_output(new_dirname_k, filename_k, xp, delp, timecount, knum);
        }
        xmid=xp[timecount-1];
        a=exp(xmid);
        for (i=0;i<N1;i++) delstart[i]=delp[i][timecount-1];
    };
    Gamma1 = GLARGE2;
    Gamma2 = GLARGE2;
    //changing decay rates (if neccesary)
    //Full evolution equations including first order perturbation of three inflaton fields are solved until the contribution
    //from two inflatons sigma and psi become negligible at ln(a)=THRUNP.
    //THRUNP is set by hand, estimeted from the result of zero-mode evolution.
    //If horizon entering is later, this step is skipped.
    if(xmid < THRUNP){
        NR::odeintpert(delstart,xmid,THRUNP,epsosc,h2,hmin,nok,nbad,timecount,dxsav,full,NR::rkqs);
        if(kanalyze_switch){
        kanalyze_output(new_dirname_k, filename_k, xp, delp, timecount, knum);
        }
        xmid=xp[timecount-1];
        a=exp(xmid);
        for (i=0;i<N1;i++) delstart[i]=delp[i][timecount-1];
    };
    //Fixing sigma = psi = 0, in order to avoid solving oscillation of these two firlds which are negligible.
    delstart[0]=0;
    delstart[3]=0;
    delstart[1]=FIXPSI;
    delstart[4]=0;
    for (i=0;i<6;i++) delstart[7+i]=0;
    for (i=0;i<6;i++) delstart[16+i]=0;
    for (i=0;i<6;i++) delstart[31+i]=0;
    for (i=0;i<6;i++) delstart[40+i]=0;
    //Evolution equations for phi and its perturbations are solved, with zero-modes of sigma and psi are gixed to minimum
    //until phi begins oscillation at ln)a)=THRLAST.
    //THRLAST is set by hand according to the result for zero-mode.
    NR::odeintpert(delstart,xmid,THRLAST,epsnew,h2,hmin,nok,nbad,timecount,dxsav,newinf,NR::rkqs);
        if(kanalyze_switch){
    kanalyze_output(new_dirname_k, filename_k, xp, delp, timecount, knum);
        }
    xmid=xp[timecount-1];
    a=exp(xmid);
    for (i=0;i<N1;i++) delstart[i]=delp[i][timecount-1];
    delstart[0]=0;
    delstart[3]=0;
    delstart[1]=FIXPSI;
    delstart[4]=0;
    for (i=0;i<6;i++) delstart[7+i]=0;
    for (i=0;i<6;i++) delstart[16+i]=0;
    for (i=0;i<6;i++) delstart[31+i]=0;
    for (i=0;i<6;i++) delstart[40+i]=0;
    //Evolution equations for phi and all perturbations are solved, with zero-modes of sigma and psi are fixed to minimum.
    //until the amplitude of oscillation of phi becomed sufficiently small at ln(a)=xend.
    //xend is set by hand according to the result for zero-mode.
    //(if perturbations of of sigma and psi are not solved, superhorizon parturbations begin to decrease)
    NR::odeintpert(delstart,xmid,xend,epslast,h2,hmin,nok,nbad,timecount,dxsav,fixfix,NR::rkqs);
        if(kanalyze_switch){
    kanalyze_output(new_dirname_k, filename_k, xp, delp, timecount, knum);
        }
   
        
        //     cout << "Final Result delstart[2] = " << tr[2] << " delp[2][" << timecount-1 << "] = " << delp[2][timecount-1] << endl;
        //   for (i=0;i<N1;i++) cout << "Final Result tr[" << i << "]" << tr[i] << endl;
        if(spectrum_switch){
        spectrum_output(new_filename_sp, xp, delp, timecount, knum, k_comoving);
        }
        
        Logout( "Calculation %d%% Complete\n",calcpercentage);
    }; //End for loop of k
    
    time_end = std::chrono::system_clock::now();       // End measuring elapsed time

    int time_hours = std::chrono::duration_cast<std::chrono::hours>(time_end - time_start).count(); // Casting time
    int time_minutes = std::chrono::duration_cast<std::chrono::minutes>(time_end - time_start).count() - time_hours*60; // Casting time
    int time_seconds = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start).count() - time_hours*60*60 - time_minutes*60; // Casting time
    int time_days = time_hours / 24;
    time_hours = time_hours % 24;
    
    Logout( "Computation Time: %d d %d h %d m %d s\n",time_days,time_hours,time_minutes,time_seconds);

    delete delp_p;
    delete yp_p;
    delete xp_p;
    return 0;
}
