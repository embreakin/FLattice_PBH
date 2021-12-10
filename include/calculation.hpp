#ifndef _CALCULATION_H_
#define _CALCULATION_H_

#include "nr.h"
#include "parameters.hpp"
#include "equations.hpp"
#include "uc.hpp"

class Zeromode
{
    private:
    //eps parameters set the allowed error in adapive step size RK-method for zeromode calcualtion.
    DP eps=1.0E-6,eps2=1.0E-7,eps3=1.0E-8,eps4=1.0E-6;
    
    const int timecount_max_zero = 100000;
    
    Vec_DP *xp_p;
    Mat_DP *yp_p;
    
    Vec_DP &xp=*xp_p;
    Mat_DP &yp=*yp_p;
    
    DP k_comoving;
    
    int i,j,nbad,nok,timecount;
    DP a,H,w,rho,rhop,la,dxsav;
    
    //h parameters set the initial trial step size in adaptive step size RK-method for zeromode calcualtion and w/ perturbation calcualtion.
    DP h1=1.0E-4,h2=1.0E-6,hmin=0.0,xbegin=Ini,xend=SH;
    
   
    
    public:
    
    Zeromode (): xp_p(new Vec_DP(timecount_max_zero)), yp_p(new Mat_DP(N_zero,timecount_max_zero)){
        
    }
    
    ~Zeromode () {
        delete xp_p;
        delete yp_p;
    }
    
     void zeromode_calc();
     void zeromode_initial(Vec_DP &unp, DP &a, DP &H, DP &xbegin);
    
    
};



class Perturbation
{

    private:
    
    DP xmid, la2;
    DP theta = M_PI/4;;
    DP p = itvl;
     //eps parameters set the allowed error in adapive step size RK-method for w/-perturbation calcualtion.
    DP epsSHI=1.0e-9,epsosc=1.0e-10,epsnew=1.0e-10,epslast=1.0e-10;
    
    const int timecount_max_pert = 10000;
    
    Vec_DP *xp2_p;
    Mat_DP *delp_p;
    
    Vec_DP &xp2=*xp2_p;
    Mat_DP &delp=*delp_p;
    
    DP k_comoving;
    
    int i,j,nbad,nok,timecount;
    DP a,H,w,rho,rhop,la,dxsav;
    
    //h parameters set the initial trial step size in adaptive step size RQ-method for zeromode calcualtion and w/ perturbation calcualtion.
    DP h1=1.0E-4,h2=1.0E-6,hmin=0.0,xbegin=Ini,xend=SH;
    
    
    
    
    int knum;
    int percentage;

    
    public:
    
    Perturbation(): xp2_p(new Vec_DP(timecount_max_pert)), delp_p(new Mat_DP(N_pert,timecount_max_pert))
    {}
    
    ~Perturbation () {
        delete xp2_p;
        delete delp_p;
       
    }

    void perturbation_initial(Vec_DP &tr2, DP &k_comoving, DP &a, DP &H);
    void nonlatticerange_calc(int &k_begin, int &k_end, Zeromode &zeromode);
    void lattice_initialize( double** &latticep );
    void lattice_finalize( double** &latticep );
    void latticerange_firsthalf_calc( double** latticep, Zeromode &zeromode);
    void latticerange_secondhalf_calc(double** latticep);

};

#endif
