#include <chrono>
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "lattice.hpp"
#include "nr.h"
#include "equations.hpp"
#include "parameters.hpp"

using namespace std;

#define FILE "../unpWMAP5.txt"
#define TGT 634
//#define TGT 130

//-------------------

DP dxsav;
DP k,FIXPSI,CNT;
DP Gamma1,Gamma2,Gamma3;
int kmax,kount;
Vec_DP *xp_p;
Mat_DP *yp_p;

//xp[i] stores integration variable (log(a)) for each steps[i].
//yp_p[i][j] stores variables for each steps[i]. [j] specifies variables as follows:
//    j=0-2  : zero modes of inflaton sigma, psi, and phi
//    j=3-5  : log(a) derivatives of inflaton sigma', psi', phi'
//    j=6    : energy density of radiation
//unp[i] stores initial conditions for variables(delp_p[i][0-6]
//tr[j]: buffer for yp_p[i][j]

int nrhs;   // counts function evaluations

int main(void)
{
    const int N2=7, KMAX=100000, BINWIDTH=41,MIDPOINT=20;
    int mid,i,j,knum,nbad,nok,ncrs,p,q;
    DP a,H,rhop,la,hr,ave;
    DP eps=1.0E-6,eps2=1.0E-7,eps3=1.0E-8,eps4=1.0E-6,h1=1.0E-4,h2=1.0E-6,hmin=0.0,x1=Ini,x2=SH;
    //eps parameters set the allowed error in adapive step size RQ-method.
    //h parameters set the initial trial step size in adaptive step size RQ-method.
    Vec_DP tr(N2),scale(BINWIDTH),cross(5);
    DP w,rho;
    Vec_DP unp(N2);
    Vec_DP *bin_p;
    Vec_DP *binx_p;
    ofstream store;
    store.open(FILE);
    xp_p=new Vec_DP(KMAX);
    yp_p=new Mat_DP(N2,KMAX);
    bin_p=new Vec_DP(KMAX/20);
    binx_p=new Vec_DP(KMAX/20);
    Vec_DP &xp=*xp_p;
    Mat_DP &yp=*yp_p;
    Vec_DP &bin=*bin_p;
    Vec_DP &binx=*binx_p;
    FIXPSI = 2*sqrt(mu*M);
    CNT = (-1)*Vbare(0,FIXPSI,FIXPHI);
    knum=TGT;
    k=Ck*exp(0.023025851*knum);
    a=exp(x1);
    //set initial conditions
    unp[0] = I1;  //initial value of sigma
    unp[1] = exp(log(Pow(4,m)/(2*m*(2*m-1)))/(2*(m-1)))*M*exp(log(mu/unp[0])/(m-1));
    unp[2] = -Cv*Cv*unp[0]/(mu*mu);
    unp[3] = 0;
    unp[4] = 0;
    unp[5] = 0;
    unp[6] = 0;
    nrhs=0;
    //set decay rates
    Gamma1=GLARGE; //sigma
    Gamma2=GLARGE; // psi
    Gamma3=GNOMAL; //phi
    dxsav=(x2-x1)/50000.0;
    kmax=KMAX;
    //    int AA, BB, H0;
    H=Fri(unp[0],unp[1],unp[2],unp[3],unp[4],unp[5],unp[6]);
    //    H0=Fri(unp[0],unp[1],unp[2],unp[3],unp[4],unp[5],unp[6]);
    //    AA=(2/3)*Pow(a, -3/2)*(1/H0);
    //    BB=-exp(log(AA)+log(a)*3/2);
    //    cout << "AA = " << AA << "\n";
    //    cout << "BB = " << BB << "\n";
    cout << "H = " << H << "\n";
    
    
    //Full evolution equations of zero-mode variables are solved until oscillatory phase begins.
    //THRUNP is set by hand, estimeted from the result of zero-mode evolution.
    NR::odeint(unp,x1,THRUNP,eps4,h2,hmin,nok,nbad,unpert,NR::rkqs);
    for (j=0;j<kount;j++) {
        for (i=0;i<N2;i++) tr[i]=yp[i][j];
        la=xp[j];
        rho=rho_tot(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
        rhop=rhoandp(tr[3],tr[4],tr[5],tr[6]);
        H=Fri(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
        w=log10(H);
        a=exp(la);
        //output data
        store << setw(6) << (la-Ini)*msigma/H << " " //t
        << setw(10) << la << " " //log(a)
        << setw(10) << tr[0] << " "  //buffer for yp_p[i][0] sigma
        << setw(20) << setprecision(20) << tr[1] << " " //buffer for yp_p[i][1] psi
        << setw(10) << tr[2] << " " //buffer for yp_p[i][1] phi
        << setw(10) << w << " "     //log10(H)
        << setw(10) << rho << " "
        << setw(10) << V_11(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_12(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_13(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_22(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_23(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_33(tr[0],tr[1],tr[2]) << " "
        << setw(10) << log10(k/(a*H)) << "\n";
    };
    tr[0]=0;
    tr[1]=FIXPSI;
    tr[3]=0;
    tr[4]=0;
    Gamma1 = GLARGE2;
    Gamma2 = GLARGE2;
    //changing decay rates (if neccesary)
    //Evolution equations for phi is solved, with zero-modes of sigma and psi are fixed to minimum,
    //until the amplitude of oscillation of phi becomed sufficiently small at ln(a)=x2.
    //x2 is set by hand according to the result for zero-mode.
    NR::odeint(tr,la,x2,eps,h2,hmin,nok,nbad,unpertfixfix,NR::rkqs);
    for (j=0;j<kount;j++) {
        for (i=0;i<N2;i++) tr[i]=yp[i][j];
        la=xp[j];
        rho=rho_tot(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
        rhop=rhoandp(tr[3],tr[4],tr[5],tr[6]);
        H=Fri(tr[0],tr[1],tr[2],tr[3],tr[4],tr[5],tr[6]);
        w=log10(H);
        a=exp(la);
        //output data
        store << setw(6) << (la-Ini)*msigma/H << " "
        << setw(10) << la << " " //log(a)
        << setw(10) << tr[0] << " "
        << setw(10) << tr[1] << " "
        << setw(10) << tr[2] << " "
        << setw(10) << w << " "
        << setw(10) << rho << " "
        << setw(10) << V_11(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_12(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_13(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_22(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_23(tr[0],tr[1],tr[2]) << " "
        << setw(10) << V_33(tr[0],tr[1],tr[2]) << " "
        << setw(10) << log10(k/(a*H)) << "\n";
    };
    delete yp_p;
    delete xp_p;
    delete bin_p;
    delete binx_p;
    
}

