#include "nr.h"//Classic fourth-order Runge–Kutta methodvoid NR::rk4(Vec_I_DP &y, Vec_I_DP &dydx, const DP x, const DP h,	Vec_O_DP &yout,void derivs(const DP,Vec_I_DP &, Vec_O_DP &, const DP), const DP k_comoving){	int i;	DP xh,hh,h6;	int nn=y.size();	Vec_DP dym(nn),dyt(nn),yt(nn);	hh=h*0.5;	h6=h/6.0;	xh=x+hh;	for (i=0;i<nn;i++) yt[i]=y[i]+hh*dydx[i];	derivs(xh,yt,dyt,k_comoving);	for (i=0;i<nn;i++) yt[i]=y[i]+hh*dyt[i];	derivs(xh,yt,dym,k_comoving);	for (i=0;i<nn;i++) {		yt[i]=y[i]+h*dym[i];		dym[i] += dyt[i];	}	derivs(x+h,yt,dyt,k_comoving);	for (i=0;i<nn;i++)		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);    }// Runge–Kutta Cash–Karp methodvoid NR::rkck(Vec_I_DP &y, Vec_I_DP &dydx, const DP x,              const DP h, Vec_O_DP &yout, Vec_O_DP &yerr,              void derivs(const DP, Vec_I_DP &, Vec_O_DP &, const DP), const DP k_comoving){    static const DP a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,    b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9,    b43=1.2, b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0,    b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,    b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,    c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,    dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,    dc4=c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;    int i;        int n=y.size();    Vec_DP ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),ytemp(n);        //  cout << " Before 1 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;    for (i=0;i<n;i++) ytemp[i]=y[i]+b21*h*dydx[i];    //  cout << " After 1 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;        derivs(x+a2*h,ytemp,ak2,k_comoving);    //  cout << "derivs(x+a2*h,ytemp,ak2) = " << ak2[2] << endl;        //  cout << " Before 2 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;    for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);    // cout << " After 2 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;        derivs(x+a3*h,ytemp,ak3,k_comoving);    //  cout << "derivs(x+a3*h,ytemp,ak3) = " << ak3[2] << endl;        // cout << " Before 3 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;    for (i=0;i<n;i++)  ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);    // cout << " After 3 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;        derivs(x+a4*h,ytemp,ak4,k_comoving);    // cout << "derivs(x+a4*h,ytemp,ak4) = " << ak4[2] << endl;        //   cout << " Before 4 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;    for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);    //   cout << " After 4 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;        derivs(x+a5*h,ytemp,ak5,k_comoving);    //  cout << "derivs(x+a5*h,ytemp,ak5) = " << ak5[2] << endl;        //  cout << " Before 5 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;    for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);    //  cout << " After 5 y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;        derivs(x+a6*h,ytemp,ak6,k_comoving);    //  cout << "derivs(x+a6*h,ytemp,ak6) = " << ak6[2] << endl;        for (i=0;i<n;i++)        yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);    for (i=0;i<n;i++)        yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);}//extern ofstream STEP;// “quality-controlled” Runge Kuttavoid NR::rkqs(Vec_IO_DP &y, Vec_IO_DP &dydx, DP &x, const DP htry,              const DP eps, Vec_I_DP &yscal, DP &hdid, DP &hnext,              void derivs(const DP, Vec_I_DP &, Vec_O_DP &, const DP), const DP k_comoving){    const DP SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;    int i;    DP errmax,h,htemp,xnew;        int nn=y.size();    h=htry;    Vec_DP yerr(nn),ytemp(nn);    for (;;) {        //  cout << " Before rkck y[2] = " << y[2] << ", ytemp[2] = " << ytemp[2] << endl;        rkck(y,dydx,x,h,ytemp,yerr,derivs,k_comoving);        //   cout << " After rkck y[2] = " << y[2]  << ", ytemp[2] = " << ytemp[2] << endl;        errmax=0.0;        for (i=0;i<nn;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));        errmax /= eps;        if (errmax <= 1.0) break;        htemp=SAFETY*h*pow(errmax,PSHRNK);//*exp(-1.5*h);        // cout << "Before h = " <<  h << " htemp = " << htemp << "\n";        h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));        //cout << "After h = " <<  h << " htemp = " << htemp << "\n";                xnew=x+h;        //cout << setw(10) << errmax << "\n";        //cin >> i;        //STEP << " " << h;        if (h == 0) nrerror("stepsize underflow in rkqs");    }    //  cout << "lna = " << x << " errmax = " << errmax << "\n";    if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);//*exp(-1.5*h);    else hnext=5.0*h;    x += (hdid=h);    for (i=0;i<nn;i++) y[i]=ytemp[i];}//NR::odeint(unp,x1,THRUNP,eps4,h2,hmin,nok,nbad,unpert,NR::rkqs);void NR::odeint(Vec_IO_DP &ystart, const DP x1, const DP x2, const DP eps,                const DP h1, const DP hmin, int &nok, int &nbad, int &timecount, DP &dxsav,                void derivs(const DP, Vec_I_DP &, Vec_O_DP &, const DP),                void rkqs(Vec_IO_DP &, Vec_IO_DP &, DP &, const DP, const DP,                          Vec_I_DP &, DP &, DP &, void (*)(const DP, Vec_I_DP &, Vec_O_DP &,const DP),const DP),const DP k_comoving, Vec_DP *xp_p, Mat_DP *yp_p ,const int timecount_max_zero){    const int MAXSTP=10000000;    const DP TINY=1.0e-30;    int i,nstp;    DP xsav,x,hnext,hdid,h;        int nvar=ystart.size();    Vec_DP yscal(nvar),y(nvar),dydx(nvar);    Vec_DP &xp=*xp_p;    Mat_DP &yp=*yp_p;    x=x1;    h=SIGN(h1,x2-x1);    nok = nbad = timecount = 0;    for (i=0;i<nvar;i++) y[i]=ystart[i];    if (timecount_max_zero > 0) xsav=x-dxsav*2.0;    for (nstp=0;nstp<MAXSTP;) {        derivs(x,y,dydx,k_comoving);        for (i=0;i<nvar;i++)            yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;        if (timecount_max_zero > 0 && timecount < timecount_max_zero-1 && fabs(x-xsav) > fabs(dxsav)) {            for (i=0;i<nvar;i++) yp[i][timecount]=y[i];            xp[timecount++]=x;            xsav=x;        }        if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;        rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs,k_comoving);//       std::cout << x << "\n";        if (hdid == h) ++nok; else ++nbad;        if ((x-x2)*(x2-x1) >= 0.0) {            for (i=0;i<nvar;i++) ystart[i]=y[i];            if (timecount_max_zero != 0) {                for (i=0;i<nvar;i++) yp[i][timecount]=y[i];                xp[timecount++]=x;            }            return;        }        if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");        h=hnext;    }    nrerror("Too many steps in routine odeint");}void NR::odeintpert(Vec_IO_DP &ystart, const DP x1, const DP x2, const DP eps,                    const DP h1, const DP hmin, int &nok, int &nbad, int &timecount, DP &dxsav,                    void derivs(const DP, Vec_I_DP &, Vec_O_DP &, const DP),                    void rkqs(Vec_IO_DP &, Vec_IO_DP &, DP &, const DP, const DP,                              Vec_I_DP &, DP &, DP &, void (*)(const DP, Vec_I_DP &, Vec_O_DP &, const DP),const DP), const DP k_comoving, Vec_DP *xp2_p, Mat_DP *delp_p  ,const int timecount_max_pert){    const int MAXSTP=100000000;    const DP TINY=1.0e-30;    int i,j,nstp;    static int keepcount = 0; ++keepcount;    DP xsav,x,hnext,hdid,h;    DP H;        int nvar=ystart.size();    Vec_DP yscal(nvar),y(nvar),dydx(nvar);    Vec_DP &xp2=*xp2_p;    Mat_DP &delp=*delp_p;    x=x1;    h=SIGN(h1,x2-x1);    nok = nbad = timecount = 0;    for (i=0;i<nvar;i++) y[i]=ystart[i];//     if (keepcount % 4 == 2)   std::cout << "After h=SIGN(h1,x2-x1) h = " <<  h << "\n";    //  cout << "y[2] = " << y[2] << endl;        if (timecount_max_pert > 0) xsav=x-dxsav*2.0; //dxsav=(x2-x1)/5000.0;    for (nstp=0;nstp<MAXSTP;nstp++) {//        if (keepcount % 4 == 2)   std::cout << "nstp = " << nstp << ", check all delp[2][" << timecount << "] phi = " << delp[2][timecount]  << ", y[2] = " << y[2] << std::endl;        derivs(x,y,dydx,k_comoving);        for (i=0;i<nvar;i++)            yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;        if (timecount_max_pert > 0 && timecount < timecount_max_pert-1 && fabs(x-xsav) > fabs(dxsav)) {            for (i=0;i<nvar;i++) delp[i][timecount]=y[i];//            if (keepcount % 4 == 2)       std::cout << "nstp = " << nstp << ", delp[2][" << timecount << "] phi = " << delp[2][timecount] << ", y[2] = " << y[2] << std::endl;                        //cout << timecount << "\n";        //Writing time steps            xp2[timecount++]=x;            xsav=x;        }//        std::cout << "timecount_max_pert" << timecount_max_pert << "dxsav" << dxsav << "nvar" << nvar << "\n";        //cout << setw(10) << nstp << setw(10) << x << "\n";//        if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;        //cout <<setw(10) << hdid << "\n";        //cin >> i;//              if (keepcount % 4 == 2)   std::cout << "nstp = " << nstp << ",Before rkqs delp[2][" << timecount << "] phi = " << delp[2][timecount] << ", y[2] = " << y[2] << std::endl;        rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs,k_comoving);       // std::cout << x << "\n";//              if (keepcount % 4 == 2)      std::cout << h << "\n";        //cin >> i;        //        if (keepcount % 4 == 2)  cout << "nstp = " << nstp << ",After rkqs delp[2][" << timecount << "] phi = " << delp[2][timecount] << ", y[2] = " << y[2] << endl;                if (hdid == h) ++nok; else ++nbad;                if ((x-x2)*(x2-x1) >= 0.0 ) {            for (i=0;i<nvar;i++) ystart[i]=y[i];            if (timecount_max_pert != 0) {                for (i=0;i<nvar;i++) delp[i][timecount]=y[i];                xp2[timecount++]=x;//                           if (keepcount % 4 == 2)     std::cout << "nstp = " << nstp << ", timecount_max_pert != 0  delp[2][" << timecount << "] phi = " << delp[2][timecount] << ", y[2] = " << y[2] << std::endl;//                            }//                   if (keepcount % 4 == 2)    std::cout << "nstp = " << nstp << ", (x-x2)*(x2-x1) >= 0.0  delp[2][" << timecount << "] phi = " << delp[2][timecount] << ", y[2] = " << y[2] << std::endl;            return;        }        if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");        h=hnext;                //   cout << "nstp = " << nstp << ",(nstp=0;nstp<MAXSTP;nstp++) right before }  delp[2][" << timecount << "] phi = " << delp[2][timecount] << ", y[2] = " << y[2] << endl;        //cout << " delp[2][timecount] phi " << delp[2][timecount] << endl;//              if (keepcount % 4 == 2) std::cout << "keepcount = " << keepcount << " delp[2][timecount] phi " << delp[2][timecount] << std::endl;            }        std::cout << "calculation stopped at la=" << x << "with hdid=" << hdid << "\n";    nrerror("Too many steps in routine odeint");}