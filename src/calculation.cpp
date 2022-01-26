#include "calculation.hpp"
#include "utilities.hpp"

//------------------------------------
//Subroutines for zeromode calculation
//------------------------------------

//Initial conditions for zeromode calculation
void Zeromode::zeromode_initial(Vec_DP &unp, DP &a, DP &H, DP &xbegin){
        
        //set initial conditions
        a = exp(xbegin);
    
        unp[0] = Init_sigma;  //initial value of sigma
        unp[1] = exp(log(Pow(4,m_par)/(2*m_par*(2*m_par-1)))/(2*(m_par-1)))*M_par*exp(log(mu_par/unp[0])/(m_par-1));
        unp[2] = -Cv_par*Cv_par*unp[0]/(mu_par*mu_par);
        unp[3] = 0;
        unp[4] = 0;
        unp[5] = 0;
        unp[6] = 0;
        H=Fri(unp[0],unp[1],unp[2],unp[3],unp[4],unp[5],unp[6]);
        
        //set decay rates
        Gamma1=GLARGE; //sigma
        Gamma2=GLARGE; // psi
        Gamma3=GNOMAL; //phi
   
    }

//Zeromode calculation
void Zeromode::zeromode_calc(){
    
    Vec_DP unp(N_zero), tr(N_zero);
    
    for(i = 0; i < knum_zero.size(); i++ ){
        Logout("knum = %d, kMpc = %2.5e, kMPl = %2.5e \n\n",knum_zero[i], UC::knum_to_kMpc(knum_zero[i]), UC::knum_to_kMPl(knum_zero[i]));
    }
    
    k_comoving = UC::knum_to_kMPl(knum_zero[0]);
    
    Zeromode::zeromode_initial(unp, a, H, xbegin);
    
    dxsav=(xend-xbegin)/50000.0;
   
    //Full evolution equations of zero-mode variables are solved until oscillatory phase begins.
    //THRUNP is set by hand, estimeted from the result of zero-mode evolution.
    NR::odeint(unp,xbegin,THRUNP,eps4,h2,hmin,nok,nbad,timecount,dxsav,unpert,NR::rkqs,k_comoving,&xp,&yp, timecount_max_zero);
    zeromode_output(new_filename_zero, xp, yp, timecount);
    
    for (j=0;j<N_zero;j++) tr[j]=yp[j][timecount-1];
    la=xp[timecount-1];
    
    tr[0]=0;
    tr[1]=FIXPSI;
    tr[3]=0;
    tr[4]=0;
    
    Gamma1 = GLARGE2;
    Gamma2 = GLARGE2;
    //changing decay rates (if neccesary)
    //Evolution equations for phi is solved, with zero-modes of sigma and psi are fixed to minimum,
    //until the amplitude of oscillation of phi becomed sufficiently small at ln(a)=xend.
    //xend is set by hand according to the result for zero-mode.
    
    NR::odeint(tr,la,xend,eps,h2,hmin,nok,nbad,timecount,dxsav,unpertfixfix,NR::rkqs,k_comoving,&xp,&yp, timecount_max_zero);
    
    zeromode_output(new_filename_zero, xp, yp, timecount);
    
    
    
}

//------------------------------------------------------
//Subroutines for zeromode with perturbation calculation
//------------------------------------------------------

//Initial conditions for zeromode with perturbation calculation
void Perturbation::perturbation_initial(Vec_DP &tr2, DP &k_comoving, DP &a, DP &H)
{
    for (i=N_zero;i<N_pert;i++) tr2[i]=0;
    for (i=0;i<3;i++) tr2[(i*4)+7]=cos(theta)/(sqrt(2*k_comoving)*a);
    for (i=0;i<3;i++) tr2[(i*4)+31]=sin(theta)/(sqrt(2*k_comoving)*a);
    for (i=0;i<3;i++) tr2[(i*4)+16]=(-1)*H*tr2[(i*4)+7] + k_comoving*tr2[(i*4)+31]/a;
    for (i=0;i<3;i++) tr2[(i*4)+40]=(-1)*H*tr2[(i*4)+31] - k_comoving*tr2[(i*4)+7]/a;
    for (i=0;i<3;i++) tr2[i+25]=0.5*a*a*(tr2[3]*tr2[i+16] + tr2[4]*tr2[i+19] + tr2[5]*tr2[i+22] + 3*H*(tr2[3]*tr2[i+7] + tr2[4]*tr2[i+10] + tr2[5]*tr2[i+13]) + V_1(tr2[0],tr2[1],tr2[2])*tr2[i+7] + V_2(tr2[0],tr2[1],tr2[2])*tr2[i+10] + V_3(tr2[0],tr2[1],tr2[2])*tr2[i+13])/(k_comoving*k_comoving);
    for (i=0;i<3;i++) tr2[i+49]=0.5*a*a*(tr2[3]*tr2[i+40] + tr2[4]*tr2[i+43] + tr2[5]*tr2[i+46] + 3*H*(tr2[3]*tr2[i+31] + tr2[4]*tr2[i+34] + tr2[5]*tr2[i+37]) + V_1(tr2[0],tr2[1],tr2[2])*tr2[i+31] + V_2(tr2[0],tr2[1],tr2[2])*tr2[i+34] + V_3(tr2[0],tr2[1],tr2[2])*tr2[i+37])/(k_comoving*k_comoving);
    for (i=0;i<3;i++) tr2[i+25]=tr2[i+25]/(1 - (a*a*(tr2[3]*tr2[3] + tr2[4]*tr2[4] + tr2[5]*tr2[5])/(2*k_comoving*k_comoving)));
    for (i=0;i<3;i++) tr2[i+49]=tr2[i+49]/(1 - (a*a*(tr2[3]*tr2[3] + tr2[4]*tr2[4] + tr2[5]*tr2[5])/(2*k_comoving*k_comoving)));
    for (i=0;i<3;i++) tr2[i+28]=(-1)*H*tr2[i+25] - 0.5*(tr2[3]*tr2[i+7] + tr2[4]*tr2[i+10] + tr2[5]*tr2[i+13]);
    for (i=0;i<3;i++) tr2[i+52]=(-1)*H*tr2[i+49] - 0.5*(tr2[3]*tr2[i+31] + tr2[4]*tr2[i+34] + tr2[5]*tr2[i+37]);
}

//Zeromode with perturbation calculation (no lattice range)
void Perturbation::nonlatticerange_calc(int &k_begin, int &k_end, Zeromode &zeromode){

    Vec_DP delstart(N_pert);
    Vec_DP unp2(N_zero), tr2(N_pert);
    Vec_DP adia(3),iso(3),fields(3),numdens(3),term(6);
    Vec_DP zeta(6),dens(4),yout(N_zero),dydx(N_zero);
    
    static int nonlatticerange_count = 0;
    int k_loopend;
    
    if(latticerange_switch && nonlatticerange_count == 0){k_loopend = k_end;}
    else{ k_loopend = k_end + kinterval_knum;};
    
    for (knum = k_begin; knum < k_loopend; knum = knum + kinterval_knum)
    {
        percentage =  round( ( (knum - k_begin)/kinterval_knum + 1 )*100 / ( (k_loopend - k_begin)/kinterval_knum ) );
        
        k_comoving = UC::knum_to_kMPl(knum);
        
        Logout("knum = %d, kMpc = %2.5e, kMPl = %2.5e: ",knum, UC::knum_to_kMpc(knum), k_comoving);
        
         p=itvl;
        
            zeromode.zeromode_initial(unp2, a, H, xbegin);
    

            unpert(la2,unp2,dydx,k_comoving);
            //First, solve the evolution of zero mode, until the specified mode crosses the horizon.
            //Calculation for first order perturbation begins slightly before the horizon crossing.
            //Runge-Kutta method with constant step size is used.
     
            for (la2=xbegin;k_comoving/(a*H)>5000;la2=la2) {    //Here, zero-mode is solved until k=a*H*5000.
                if (p==itvl){
                    rho=rho_tot(unp2[0],unp2[1],unp2[2],unp2[3],unp2[4],unp2[5],unp2[6]);
                    w=log10(H);
                    p=0;
                };
                p++;
                H=Fri(unp2[0],unp2[1],unp2[2],unp2[3],unp2[4],unp2[5],unp2[6]);
                unpert(la2,unp2,dydx,k_comoving);
                NR::rk4(unp2,dydx,la2,dla,yout,unpert,k_comoving);
                for (j=0;j<N_zero;j++){
                    unp2[j]=yout[j];
                };
                la2=la2+dla;
                a=exp(la2);
            };
         std::cout << "a = " << a << "\n";

//            std::cout << "la2 = " << la2 << "\n";
            //Hereafter, evolution of perturbation is solved.
            //Setting initial conditions for perturbation
            xmid=la2;
            for (i=0;i<N_zero;i++) tr2[i]=unp2[i];

            a=exp(xmid);
            H=Fri(tr2[0],tr2[1],tr2[2],tr2[3],tr2[4],tr2[5],tr2[6]);
            rhop=rhoandp(tr2[3],tr2[4],tr2[5],tr2[6]);

            dxsav=(xend-xmid)/5000.0;
      
        perturbation_initial(tr2, k_comoving, a, H);
        
        for (i=0;i<N_pert;i++) delstart[i]=tr2[i];
            //Full evolution equations including first order perturbation of three inflaton fields are solved until oscillatory phase begins.
            //If horizon entering is later, this step is skipped.
        //        for (i=0;i<N_pert;i++) std::cout << "delstart[" << i << "] = " << delstart[i] << std::endl;
            if(xmid < OSCSTART){
               NR::odeintpert(delstart,xmid,OSCSTART,epsSHI,h2,hmin,nok,nbad,timecount,dxsav,full,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
        //         std::cout << "timecount = " << timecount << std::endl;
                if(kanalyze_switch){
                kanalyze_output(new_dirname_k, filename_k, xp2, delp, timecount, knum,k_comoving);
                }
                xmid=xp2[timecount-1];
                a=exp(xmid);
                for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
            };
            Gamma1 = GLARGE2;
            Gamma2 = GLARGE2;
        //        for (i=0;i<N_pert;i++) std::cout << "delstart[" << i << "] = " << delstart[i] << std::endl;
        //        std::cout << "xmid = " << xmid << std::endl;
        //        std::cout << "THRUNP = " << THRUNP << std::endl;
        //        std::cout << "epsosc = " << epsosc << std::endl;
        //        std::cout << "h2 = " << h2 << std::endl;
        //        std::cout << "hmin = " << hmin << std::endl;
        //        std::cout << "nok = " << nok << std::endl;
        //        std::cout << "nbad = " << nbad << std::endl;
        //        std::cout << "timecount = " << timecount << std::endl;
        //        std::cout << "dxsav = " << dxsav << std::endl;
            //changing decay rates (if neccesary)
            //Full evolution equations including first order perturbation of three inflaton fields are solved until the contribution
            //from two inflatons sigma and psi become negligible at ln(a)=THRUNP.
            //THRUNP is set by hand, estimeted from the result of zero-mode evolution.
            //If horizon entering is later, this step is skipped.
        
//        Logout("delstart[0] = %2.5e \n",delstart[i] );
        
        if(spectrum_bfosc_switch){
            spectrum_bfosc_output(new_filename_spbfosc, xp2, delp, timecount, knum, k_comoving);
        };


            if(xmid < THRUNP){
                NR::odeintpert(delstart,xmid,THRUNP,epsosc,h2,hmin,nok,nbad,timecount,dxsav,full,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
        //         std::cout << "timecount = " << timecount << std::endl;
                if(kanalyze_switch){
                kanalyze_output(new_dirname_k, filename_k, xp2, delp, timecount, knum,k_comoving);
                }
                xmid=xp2[timecount-1];
                a=exp(xmid);
                for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
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
            NR::odeintpert(delstart,xmid,THRLAST,epsnew,h2,hmin,nok,nbad,timecount,dxsav,newinf,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
        //         std::cout << "timecount = " << timecount << std::endl;
                if(kanalyze_switch){
            kanalyze_output(new_dirname_k, filename_k, xp2, delp, timecount, knum, k_comoving);
                }
            xmid=xp2[timecount-1];
            a=exp(xmid);
            for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
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
            NR::odeintpert(delstart,xmid,xend,epslast,h2,hmin,nok,nbad,timecount,dxsav,fixfix,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
        //         std::cout << "timecount = " << timecount << std::endl;
                if(kanalyze_switch){
            kanalyze_output(new_dirname_k, filename_k, xp2, delp, timecount, knum, k_comoving);
                }


                //     cout << "Final Result delstart[2] = " << tr2[2] << " delp[2][" << timecount-1 << "] = " << delp[2][timecount-1] << endl;
                //   for (i=0;i<N_pert;i++) cout << "Final Result tr2[" << i << "]" << tr2[i] << endl;
                if(spectrum_switch){
                spectrum_output(new_filename_sp, xp2, delp, timecount, knum, k_comoving);
                }

        Logout( "Calculation %d%% Complete\n\n",percentage);
            };
    
   nonlatticerange_count++;

    
}

//--------------------------
//Lattice Range Subroutines
//--------------------------

void Perturbation::lattice_initialize( double** &latticep )
{
    latticep = new double* [N/2];
    latticep[0] = new double [(N/2)*N_pert];
    for( int loop = 0; loop < N/2; ++loop )
    {
        latticep[loop] = latticep[0] + loop*N_pert;
        
    }
}

void Perturbation::lattice_finalize( double** &latticep )
{
    delete [] latticep[0];
    delete [] latticep;
}

void Perturbation::latticerange_firsthalf_calc( double** latticep, Zeromode &zeromode ){
    
    Vec_DP delstart(N_pert);
    Vec_DP unp2(N_zero), tr2(N_pert);
    Vec_DP adia(3),iso(3),fields(3),numdens(3),term(6);
    Vec_DP zeta(6),dens(4),yout(N_zero),dydx(N_zero);
    
    Logout("kfrom_MPl_lattice =  %2.5e, kto_MPl_lattice =  %2.5e ", kfrom_MPl_lattice, kto_MPl_lattice);
    
    for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){

           percentage =  round( ( lattice_loop + 1 )*100 / ( N/2 ) );

           k_comoving = kfrom_MPl_lattice + k_lattice_grid_min_MPl*lattice_loop;

           Logout("k_comoving =  %2.5e \n", k_comoving);
            
            p=itvl;
            
            zeromode.zeromode_initial(unp2, a, H, xbegin);
            
            unpert(la2,unp2,dydx,k_comoving);
            
            
            //First, solve the evolution of zero mode, until the specified mode crosses the horizon.
            //Calculation for first order perturbation begins slightly before the horizon crossing.
            //Runge-Kutta method with constant step size is used.
            
            for (la2=xbegin;k_comoving/(a*H)>5000;la2=la2) {    //Here, zero-mode is solved until k=a*H*5000.
                if (p==itvl){
                    rho=rho_tot(unp2[0],unp2[1],unp2[2],unp2[3],unp2[4],unp2[5],unp2[6]);
                    w=log10(H);
                    p=0;
                };
                p++;
                H=Fri(unp2[0],unp2[1],unp2[2],unp2[3],unp2[4],unp2[5],unp2[6]);
                unpert(la2,unp2,dydx,k_comoving);
                NR::rk4(unp2,dydx,la2,dla,yout,unpert,k_comoving);
                for (j=0;j<N_zero;j++){
                    unp2[j]=yout[j];
                };
                la2=la2+dla;
                a=exp(la2);
                
            };
         //    Logout("1:unp2[0] = %2.5e \n", unp2[0] );
              //  std::cout << "a = " << a << "\n";
            //Hereafter, evolution of perturbation is solved.
            //Setting initial conditions for perturbation
            xmid=la2;
           
            for (i=0;i<N_zero;i++) tr2[i]=unp2[i];

            a=exp(xmid);
            H=Fri(tr2[0],tr2[1],tr2[2],tr2[3],tr2[4],tr2[5],tr2[6]);
            rhop=rhoandp(tr2[3],tr2[4],tr2[5],tr2[6]);
            theta=M_PI/4;
            dxsav=(xend-xmid)/5000.0;

            perturbation_initial(tr2, k_comoving, a, H);
            
         //   for ( i=0;i<N_zero;i++) Logout("2:tr2[%d] = %2.5e \n", i , tr2[i] );
            
              for (i=0;i<N_pert;i++) delstart[i]=tr2[i];
            
//            Logout("delstart[0] = %2.5e \n",delstart[0] );
            //Full evolution equations including first order perturbation of three inflaton fields are solved until oscillatory phase begins.
            //If horizon entering is later, this step is skipped.
            //        for (i=0;i<N1;i++) std::cout << "delstart[" << i << "] = " << delstart[i] << std::endl;
            
//            std::cout << "xmid = " << xmid << "\n";
//            std::cout << "OSCSTART = " << OSCSTART << "\n";
            
            if(xmid < OSCSTART){
                NR::odeintpert(delstart,xmid,OSCSTART,epsSHI,h2,hmin,nok,nbad,timecount,dxsav,full,NR::rkqs, k_comoving, &xp2, &delp, timecount_max_pert);
                //         std::cout << "timecount = " << timecount << std::endl;
                if(kanalyze_switch){
                    kanalyze_output(new_dirname_k, filename_k, xp2, delp, timecount, knum, k_comoving);
                }
                xmid=xp2[timecount-1];
                a=exp(xmid);
                for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
            };
            Gamma1 = GLARGE2;
            Gamma2 = GLARGE2;
            //        for (i=0;i<N1;i++) std::cout << "delstart[" << i << "] = " << delstart[i] << std::endl;
            //        std::cout << "xmid = " << xmid << std::endl;
            //        std::cout << "THRUNP = " << THRUNP << std::endl;
            //        std::cout << "epsosc = " << epsosc << std::endl;
            //        std::cout << "h2 = " << h2 << std::endl;
            //        std::cout << "hmin = " << hmin << std::endl;
            //        std::cout << "nok = " << nok << std::endl;
            //        std::cout << "nbad = " << nbad << std::endl;
            //        std::cout << "timecount = " << timecount << std::endl;
            //        std::cout << "dxsav = " << dxsav << std::endl;
            //changing decay rates (if neccesary)
            //Full evolution equations including first order perturbation of three inflaton fields are solved until the contribution
            //from two inflatons sigma and psi become negligible at ln(a)=THRUNP.
            //THRUNP is set by hand, estimeted from the result of zero-mode evolution.
            //If horizon entering is later, this step is skipped.

           for (i=0;i<N_pert;i++) latticep[lattice_loop][i] = delstart[i];
            
//  for (i=0;i<N_pert;i++) Logout("latticep[%d][%d] = %2.5e \n",lattice_loop,i,latticep[lattice_loop][i] );
        Logout("latticep[%d][0] = %2.5e \n",lattice_loop,latticep[lattice_loop][0] );
        Logout( "Calculation %d%% Complete\n\n",percentage);

        }
    
    
}

void Perturbation::latticerange_secondhalf_calc( double** latticep ){
    
    Vec_DP delstart(N_pert);
    
     Logout("kfrom_MPl_lattice =  %2.5e, kto_MPl_lattice =  %2.5e ", kfrom_MPl_lattice, kto_MPl_lattice);
    
         for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){
    
             percentage =  round( ( lattice_loop + 1 )*100 / ( N/2 ) );
    
             k_comoving = kfrom_MPl_lattice + k_lattice_grid_min_MPl*lattice_loop;
             
            Logout("k_comoving =  %2.5e ", k_comoving);
           // for (i=0;i<N_pert;i++) Logout("latticep[%d][%d] = %2.5e \n",lattice_loop,i,latticep[lattice_loop][i] );
            for (i=0;i<N_pert;i++) delstart[i]  = latticep[lattice_loop][i];
    
            //Fixing sigma = psi = 0, in order to avoid solving oscillation of these two fields which are negligible.
            delstart[0]=0;
            delstart[3]=0;
            delstart[1]=FIXPSI;
            delstart[4]=0;
            for (i=0;i<6;i++) delstart[7+i]=0;
            for (i=0;i<6;i++) delstart[16+i]=0;
            for (i=0;i<6;i++) delstart[31+i]=0;
            for (i=0;i<6;i++) delstart[40+i]=0;
            //Evolution equations for phi and its perturbations are solved, with zero-modes of sigma and psi are fixed to minimum
            //until phi begins oscillation at ln(a)=THRLAST.
            //THRLAST is set by hand according to the result for zero-mode.
            NR::odeintpert(delstart,THRUNP,THRLAST,epsnew,h2,hmin,nok,nbad,timecount,dxsav,newinf,NR::rkqs, k_comoving, &xp2, &delp, timecount_max_pert);
            //         std::cout << "timecount = " << timecount << std::endl;
            if(kanalyze_switch){
                kanalyze_output(new_dirname_k, filename_k, xp2, delp, timecount, knum,k_comoving);
            }
            xmid=xp2[timecount-1];
            a=exp(xmid);
            for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
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
            NR::odeintpert(delstart,xmid,xend,epslast,h2,hmin,nok,nbad,timecount,dxsav,fixfix,NR::rkqs, k_comoving, &xp2, &delp, timecount_max_pert);
            //         std::cout << "timecount = " << timecount << std::endl;
            if(kanalyze_switch){
                kanalyze_output(new_dirname_k, filename_k, xp2, delp, timecount, knum, k_comoving);
            }
    
    
            //     cout << "Final Result delstart[2] = " << tr2[2] << " delp[2][" << timecount-1 << "] = " << delp[2][timecount-1] << endl;
            //   for (i=0;i<N1;i++) cout << "Final Result tr2[" << i << "]" << tr2[i] << endl;
            if(spectrum_switch){
                spectrum_output(new_filename_sp, xp2, delp, timecount, knum, k_comoving);
            }
    
            Logout( "Calculation %d%% Complete\n\n",percentage);
        }; //End 2nd for loop of k
    
    
}

