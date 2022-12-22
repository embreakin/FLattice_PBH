//Doxygen
/**
* @file    calculation.cpp
* @brief    Calculation source file
* @author   Francis Otani
* @date
* @details
*/


#include "calculation.hpp"
#include "utilities.hpp"
#include <chrono>// Measuring elapsed time




//------------------------------------
//Subroutines for zeromode calculation
//------------------------------------

//Initial conditions for zeromode calculation
void Zeromode::zeromode_initial(Vec_DP &unp, DP &a, DP &H, DP &xbegin){
        
        //set initial conditions
        a = exp(xbegin);
    
        unp[0] = sigma_init;  //initial value of sigma
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
        Gamma3=GNORMAL; //phi
    
    static int function_count = 0; ++function_count;
    
    if(function_count == 1){
    //Analytically estimated value of sigma at the end of hybrid inflation
     
       Logout("m_par: %2.5e \n\n", m_par );
        
        Logout("Potential_bare for initial value of fields: %2.5e \n\n", Vbare(unp[0], unp[1], unp[2]) );
        Logout("Potential for initial value of fields: %2.5e \n\n", V(unp[0], unp[1], unp[2]) );
    Logout("Hybrid Inflation initial values:\n sigma = %2.5e, psi = %2.5e, phi = %2.5e, \n sigma_dot = %2.5e,  psi_dot = %2.5e, phi_dot = %2.5e, radiation = %2.5e, Hubble =  %2.5e \n\n",unp[0], unp[1], unp[2], unp[3], unp[4], unp[5], unp[6], H );
        
    Logout("Field psi is set using the approx eq which uses the initial field value of sigma. sigma =  %2.5e >> %2.5e must be satisfied for it to hold. \n\n",unp[0], pow(mu_par*pow(M_par,m_par-1),1/m_par) );
        
    Logout("Estimated value of fields sigma and phi at the end of hybrid inflation (OSCSTART): sigma = %2.5e, phi = %2.5e \n\n", sigma_c, -pow(Cv_par/mu_par,2.0)*sigma_c);
        Logout("Estimated value of fields after oscillation period between hybrid and new inflation: sigma = %2.5e, psi = %2.5e, phi = %2.5e \n\n", 0.0, FIXPSI, -pow(Cv_par/mu_par,3.0)*sigma_c );
    Logout("Estimated value of field phi after oscillation that takes place after new inflation ends: phi = %2.5e \n\n", FIXPHI );
    Logout("Estimated value of CNT: %2.5e, Actual value of CNT: %2.5e \n\n",
               3*pow(n_par/(n_par+1),2.0)*pow(Cv_par,4.0)*pow(FIXPHI/sqrt(2),2.0), CNT );
        Logout("Potential_bare for final value of fields: %2.5e \n\n", Vbare(0, FIXPSI, FIXPHI) );
        Logout("Potential for final value of fields: %2.5e \n\n", V(0, FIXPSI, FIXPHI) );
    }
    
    }

//Zeromode calculation
void Zeromode::zeromode_calc(){
    
    std::chrono::system_clock::time_point  time_BEGIN_zero, time_UNPERT_zero, time_END_zero;
    time_BEGIN_zero = std::chrono::system_clock::now();
    
    Vec_DP unp(N_zero), tr(N_zero);
    
    for(size_t i = 0; i < knum_zero.size(); i++ ){
        Logout("knum = %d, kMpc = %2.5e, kMPl = %2.5e \n\n",knum_zero[i], UC::knum_to_kMpc(knum_zero[i]), UC::knum_to_kMPl(knum_zero[i]));
    }
    
    k_comoving = UC::knum_to_kMPl(k_target); 
    
    Logout("target wave mode actually used for calculation: knum =  %d\n\n",k_target);
    
    Zeromode::zeromode_initial(unp, a, H, xbegin);
    
    dxsav=(xend-xbegin)/50000.0;
   
    //Full evolution equations of zero-mode variables are solved until oscillatory phase begins.
    //THRUNP is set by hand, estimeted from the result of zero-mode evolution.
    NR::odeint(unp,xbegin,UNPERT_EFOLD,eps4,h2,hmin,nok,nbad,timecount,dxsav,unpert,NR::rkqs,k_comoving,&xp,&yp, timecount_max_zero);
    
    time_UNPERT_zero = std::chrono::system_clock::now();
    time_calc(time_BEGIN_zero,time_UNPERT_zero,"Zeromode BEGIN~UNPERT");
    
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
    
    
    time_END_zero = std::chrono::system_clock::now();
    time_calc(time_UNPERT_zero,time_END_zero,"Zeromode UNPERT~END");
    time_calc(time_BEGIN_zero,time_END_zero,"Zeromode BEGIN~END");
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

//Zero mode with perturbation calculation (non lattice range)
void Perturbation::nonlatticerange_calc(int &k_begin, int &k_end, Zeromode &zeromode){
        
    Vec_DP delstart(N_pert);
    Vec_DP unp2(N_zero), tr2(N_pert);
    Vec_DP adia(3),iso(3),fields(3),numdens(3),term(6);
    Vec_DP zeta(6),dens(4),yout(N_zero),dydx(N_zero);
    
    static int nonlatticerange_count = 0;
    int k_loopend = 0;
    double m_end;
    double k_begin_lattice = 0;
    
    if(latticerange_switch)
    {
        if(nonlatticerange_count == 0 && kfrom_Mpc < kfrom_Mpc_lattice )
        {
            //When there is a lattice range and we have a non-lattice lower range, we don't want kend of the lower range to be included in the non-lattice range, but rather included in the lattice range
            k_loopend = k_end;
            
            m_end = (double)(k_loopend - k_begin)/kinterval_knum;
          //  Logout("m_end = %e \n\n", m_end);
            
        }
        else
        {
            //When there is a lattice range, we want kend of the upper range to be included in the non-lattice range.
            k_loopend = k_end + kinterval_knum;
                
            m_end = (double)(k_loopend - k_begin)/kinterval_knum;
        }
    }
    else//no lattice range
    {
        //This corresponds to the case when there is no lattice range and lattice_kmodes_switch is on.
        if(lattice_kmodes_switch){
            
            //In this case, we start calculating from the mode that is just one mode above kfrom_MPl_lattice
            if(k_lattice_grid_min_MPl < kfrom_MPl_lattice)
            {
                
                outrange_num = floor(kfrom_MPl_lattice/k_lattice_grid_min_MPl);
                
                latticerange_num = (N/2) - outrange_num;
                
                k_begin_lattice = (outrange_num+1)*k_lattice_grid_min_MPl;
                
                m_end = latticerange_num;

                Logout("k_lattice_grid_min_MPl < kfrom_MPl_lattice\n\n");
                Logout("outrange_num = %d, latticerange_num = %d \n\n",outrange_num, latticerange_num);
                
            }else{
                //In this case, we simply start calculating from k_lattice_grid_min_MPl
                k_begin_lattice = k_lattice_grid_min_MPl;
                m_end = N/2;
                latticerange_num = N/2;
                Logout("kfrom_MPl_lattice < k_lattice_grid_min_MPl\n\n");
                Logout("latticerange_num = %d \n\n", latticerange_num);
                
            }
            
        }else
        {
            //This corresponds to the case when there is no lattice range and lattice_kmodes_switch is off.
                //In this case, kend is included in the range.
            k_loopend = k_end + kinterval_knum;
                
            m_end = (double)(k_loopend - k_begin)/kinterval_knum;
        }
    }
    
    
   std::chrono::system_clock::time_point  time_BEGIN_pert, time_OSCSTART_pert, time_OSCEND_pert, time_UNPERT_pert,time_NEWINF_END_pert, time_END_pert;
    

    for (int m = 0; m < m_end; m++)
    //for (knum = k_begin; knum < k_loopend; knum = knum + kinterval_knum)
    {
        if(m==0){
        time_BEGIN_pert = std::chrono::system_clock::now();
        }
       
        if(!latticerange_switch && lattice_kmodes_switch){
            
            k_comoving = k_begin_lattice + m*k_lattice_grid_min_MPl;
            
            knum = UC::kMPl_to_knum(k_comoving);
            
            percentage =   ( (k_comoving - k_begin_lattice)/k_lattice_grid_min_MPl + 1) / ( ceil(m_end) )*100;

           
        }
        else{
        knum = k_begin + m*kinterval_knum;
        
        percentage =   ( (knum - k_begin)/kinterval_knum + 1) / ( ceil(m_end)  )*100;
        
        k_comoving = UC::knum_to_kMPl(knum);

        }
        
         Logout("%d/%d: knum = %d, kMpc = %2.5e, kMPl = %2.5e: \n", m+1, (int)ceil(m_end) ,knum, UC::kMPl_to_kMpc(k_comoving), k_comoving);
        
         p=itvl;
        
            zeromode.zeromode_initial(unp2, a, H, xbegin);
    

            unpert(la2,unp2,dydx,k_comoving);
            //First, solve the evolution of zero mode, until the specified mode crosses the horizon.
            //Calculation for first order perturbation begins slightly before the horizon crossing.
            //Runge-Kutta method with constant step size is used.
//     std::cout << "la2 = " << la2 << "\n";
        la2=xbegin;
        while (k_comoving/(a*H)>5000)
        {
            //Here, zero-mode is solved until k=a*H*5000.
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
//                std::cout << "2 la2 = " << la2 << "\n";
        }
//         std::cout << "a = " << a << "\n";
//        std::cout << "2' la2 = " << la2 << "\n";
//            std::cout << "xmid = " << xmid << "\n";
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
               
                if(kanalyze_switch){
                kanalyze_output(new_dirname_k, xp2, delp, timecount, knum,k_comoving);
                }
                xmid=xp2[timecount-1];
                a=exp(xmid);
//                std::cout << "2: xmid = " << xmid << "\n";
                for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
            };
        
        if(m==0 && nonlatticerange_count == 0){
        time_OSCSTART_pert = std::chrono::system_clock::now();
        time_calc(time_BEGIN_pert,time_OSCSTART_pert,"Perturbation BEGIN~OSCSTART (Non-lattice Range)");
        }
        
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
            spectrum_output(new_filename_sp_bfosc, xp2, delp, timecount, knum, k_comoving);
        };
        
        if(spectrum_afosc_switch)
        {
            if(m==0 && nonlatticerange_count == 0){
            Logout("OSCEND_EFOLD = %2.5e \n\n", OSCEND_EFOLD);
            }
            
        if(xmid < OSCEND_EFOLD){
            NR::odeintpert(delstart,xmid,OSCEND_EFOLD,epsosc,h2,hmin,nok,nbad,timecount,dxsav,full,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
    //         std::cout << "timecount = " << timecount << std::endl;
            
            if(kanalyze_switch){
            kanalyze_output(new_dirname_k, xp2, delp, timecount, knum,k_comoving);
            }
            xmid=xp2[timecount-1];
            a=exp(xmid);
            for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
        };
        
        if(m==0 && nonlatticerange_count == 0){
        time_OSCEND_pert = std::chrono::system_clock::now();
        time_calc(time_OSCSTART_pert, time_OSCEND_pert,"Perturbation OSCSTART~OSCEND (Non-lattice Range)");
        }
            
            spectrum_output(new_filename_sp_afosc, xp2, delp, timecount, knum, k_comoving);
            
            if(xmid < UNPERT_EFOLD){
                NR::odeintpert(delstart,xmid,UNPERT_EFOLD,epsosc,h2,hmin,nok,nbad,timecount,dxsav,full,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
        //         std::cout << "timecount = " << timecount << std::endl;
                
                if(kanalyze_switch){
                kanalyze_output(new_dirname_k, xp2, delp, timecount, knum,k_comoving);
                }
                xmid=xp2[timecount-1];
                a=exp(xmid);
                for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
            };
            
            if(m==0 && nonlatticerange_count == 0){
            time_UNPERT_pert = std::chrono::system_clock::now();
            time_calc(time_OSCEND_pert, time_UNPERT_pert,"Perturbation OSCEND~UNPERT (Non-lattice Range)");
            }
            
        }else
        {
            if(xmid < UNPERT_EFOLD){
                NR::odeintpert(delstart,xmid,UNPERT_EFOLD,epsosc,h2,hmin,nok,nbad,timecount,dxsav,full,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
        //         std::cout << "timecount = " << timecount << std::endl;
                
                if(kanalyze_switch){
                kanalyze_output(new_dirname_k, xp2, delp, timecount, knum,k_comoving);
                }
                xmid=xp2[timecount-1];
                a=exp(xmid);
                for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
            };
            
            if(m==0 && nonlatticerange_count == 0){
            time_UNPERT_pert = std::chrono::system_clock::now();
            time_calc(time_OSCSTART_pert, time_UNPERT_pert,"Perturbation OSCSTART~UNPERT (Non-lattice Range)");
            }
        
        }
        //std::cout << "3: xmid = " << xmid << "\n";
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
            //until phi begins oscillation at ln(a)=THRLAST.
            //THRLAST is set by hand according to the result for zero-mode.
            NR::odeintpert(delstart,xmid,NEWINF_END_EFOLD,epsnew,h2,hmin,nok,nbad,timecount,dxsav,newinf,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
        //         std::cout << "timecount = " << timecount << std::endl;
        
                if(kanalyze_switch){
            kanalyze_output(new_dirname_k, xp2, delp, timecount, knum, k_comoving);
                }
            xmid=xp2[timecount-1];
            a=exp(xmid);
         //std::cout << "4: xmid = " << xmid << "\n";
            for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
            delstart[0]=0;
            delstart[3]=0;
            delstart[1]=FIXPSI;
            delstart[4]=0;
            for (i=0;i<6;i++) delstart[7+i]=0;
            for (i=0;i<6;i++) delstart[16+i]=0;
            for (i=0;i<6;i++) delstart[31+i]=0;
            for (i=0;i<6;i++) delstart[40+i]=0;
        
        if(m==0 && nonlatticerange_count == 0){
        time_NEWINF_END_pert = std::chrono::system_clock::now();
        time_calc(time_UNPERT_pert, time_NEWINF_END_pert,"Perturbation UNPERT~NEWINF_END (Non-lattice Range)");
        }
        
            //Evolution equations for phi and all perturbations are solved, with zero-modes of sigma and psi are fixed to minimum.
            //until the amplitude of oscillation of phi becomed sufficiently small at ln(a)=xend.
            //xend is set by hand according to the result for zero-mode.
            //(if perturbations of of sigma and psi are not solved, superhorizon parturbations begin to decrease)
            NR::odeintpert(delstart,xmid,xend,epslast,h2,hmin,nok,nbad,timecount,dxsav,fixfix,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
        //         std::cout << "timecount = " << timecount << std::endl;
        
                if(kanalyze_switch){
            kanalyze_output(new_dirname_k, xp2, delp, timecount, knum, k_comoving);
                }


                //     cout << "Final Result delstart[2] = " << tr2[2] << " delp[2][" << timecount-1 << "] = " << delp[2][timecount-1] << endl;
                //   for (i=0;i<N_pert;i++) cout << "Final Result tr2[" << i << "]" << tr2[i] << endl;
                if(spectrum_switch){
                spectrum_output(new_filename_sp_final, xp2, delp, timecount, knum, k_comoving);
                }
        
        if(m==0 && nonlatticerange_count == 0){
        time_END_pert = std::chrono::system_clock::now();
        time_calc(time_NEWINF_END_pert,time_END_pert,"Perturbation NEWINF_END~END (Non-lattice Range)");
        time_calc(time_BEGIN_pert,time_END_pert,"Perturbation BEGIN~END (Non-lattice Range)");
            Logout( "\n" );
        }
        
        Logout( "Calculation %d%% Complete\n\n",percentage);
        
            };
    
   nonlatticerange_count++;

    
}

//--------------------------
//Lattice Range Subroutines
//--------------------------

void Perturbation::lattice_initialize( double**& latticep )
{
    latticep = new double* [N/2];
    latticep[0] = new double [(N/2)*N_pert];
    for( int loop = 0; loop < N/2; ++loop )
    {
        latticep[loop] = latticep[0] + loop*N_pert;
        
    }
}

void Perturbation::lattice_finalize( double** latticep )
{
    delete [] latticep[0];
    delete [] latticep;
}

void Perturbation::latticerange_firsthalf_calc( double** latticep, Zeromode &zeromode ){
    
    Vec_DP delstart(N_pert);
    Vec_DP unp2(N_zero), tr2(N_pert);
    Vec_DP adia(3),iso(3),fields(3),numdens(3),term(6);
    Vec_DP zeta(6),dens(4),yout(N_zero),dydx(N_zero);
    
    Logout("MPl Units: k_lattice_grid_min_MPl = %2.5e, kfrom_MPl_lattice = %2.5e, kto_MPl_lattice = %2.5e\n\n",k_lattice_grid_min_MPl, kfrom_MPl_lattice, kto_MPl_lattice);
    Logout("Mpc^-1 Units: k_lattice_grid_min_MPl = %2.5e, kfrom_MPl_lattice = %2.5e, kto_MPl_lattice = %2.5e\n\n",UC::kMPl_to_kMpc(k_lattice_grid_min_MPl), UC::kMPl_to_kMpc(kfrom_MPl_lattice), UC::kMPl_to_kMpc(kto_MPl_lattice));
    Logout("knum Units: k_lattice_grid_min_MPl = %d, kfrom_MPl_lattice = %d, kto_MPl_lattice = %d\n\n",UC::kMPl_to_knum(k_lattice_grid_min_MPl), UC::kMPl_to_knum(kfrom_MPl_lattice), UC::kMPl_to_knum(kto_MPl_lattice));
    
    double k_comoving_start;
    
    if (k_lattice_grid_min_MPl < kfrom_MPl_lattice)
    {
        outrange_num = floor(kfrom_MPl_lattice/k_lattice_grid_min_MPl);
        
        latticerange_num = (N/2) - outrange_num;
        
        k_comoving_start = (outrange_num+1)*k_lattice_grid_min_MPl;
        
        Logout("k_lattice_grid_min_MPl < kfrom_MPl_lattice\n\n");
        Logout("outrange_num = %d, latticerange_num = %d \n\n",outrange_num, latticerange_num);

    }
    else{
        latticerange_num = N/2;
        k_comoving_start = k_lattice_grid_min_MPl;
        
        Logout("kfrom_MPl_lattice < k_lattice_grid_min_MPl\n\n");
         Logout("latticerange_num = %d \n\n", latticerange_num);
         
    }
    
     Logout("-----------------------------------------------------\n\n");
    
    std::chrono::system_clock::time_point  time_BEGIN_pert, time_OSCSTART_pert;
    
    for (int latticerange_loop = 0; latticerange_loop < latticerange_num; latticerange_loop++){

        if(latticerange_loop==0){
        time_BEGIN_pert = std::chrono::system_clock::now();
        }
        
           percentage =  round( ( latticerange_loop + 1 )*100 / ( latticerange_num ) );

           k_comoving = k_comoving_start + k_lattice_grid_min_MPl*latticerange_loop;
            knum = UC::kMPl_to_knum(k_comoving);
            Logout("%d/%d: knum = %d, kMpc = %2.5e, kMPl = %2.5e: \n", latticerange_loop+1, latticerange_num, knum, UC::kMPl_to_kMpc(k_comoving), k_comoving);
        
           p=itvl;
            
            zeromode.zeromode_initial(unp2, a, H, xbegin);
            
            unpert(la2,unp2,dydx,k_comoving);
            
            
            //First, solve the evolution of zero mode, until the specified mode crosses the horizon.
            //Calculation for first order perturbation begins slightly before the horizon crossing.
            //Runge-Kutta method with constant step size is used.
            
            la2=xbegin;
            while(k_comoving/(a*H)>5000)
            {
                //Here, zero-mode is solved until k=a*H*5000.
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
                
            }
            
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
                    kanalyze_output(new_dirname_k_lattice, xp2, delp, timecount, knum, k_comoving);
                }
                xmid=xp2[timecount-1];
                a=exp(xmid);
                for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
            };
        
        
            if(latticerange_loop==0){
            time_OSCSTART_pert = std::chrono::system_clock::now();
            time_calc(time_BEGIN_pert,time_OSCSTART_pert,"Perturbation BEGIN~OSCSTART (Lattice Range)");
            }
        
        
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
        
        if(spectrum_bfosc_switch){
            spectrum_output(new_filename_sp_bfosc, xp2, delp, timecount, knum, k_comoving);
        };
        
        if (k_lattice_grid_min_MPl < kfrom_MPl_lattice)
        {
            for (i=0;i<N_pert;i++){ latticep[outrange_num+latticerange_loop][i] = delstart[i];
            }
    
        }
        else{
            for (i=0;i<N_pert;i++) latticep[latticerange_loop][i] = delstart[i];
             
        }
           
            
//  for (i=0;i<N_pert;i++) Logout("latticep[%d][%d] = %2.5e \n",lattice_loop,i,latticep[lattice_loop][i] );
//        Logout("latticep[%d][0] = %2.5e \n",lattice_loop,latticep[lattice_loop][0] );
        Logout( "Calculation %d%% Complete\n\n",percentage);

        }
    
    if (k_lattice_grid_min_MPl < kfrom_MPl_lattice)
    {
        for(int outrange = 0; outrange < outrange_num; outrange++ )
        {
            for (i=0;i<N_pert;i++)
            {
                if(i < N_zero)
                {   //For zeromodes of the outrange wave modes, we use the values that are one step larger than kfrom_MPl_lattice
                    latticep[outrange][i] =
                    latticep[outrange_num][i];
                }
                else
                {
//                    For perturbations of the outrange wave modes, we temporarily set them to zero.
                    latticep[outrange][i] = 0;
                }
            }
        }
    

    }
    
    
    
}

void Perturbation::latticerange_secondhalf_calc( double** latticep ){
    
    
    
    Vec_DP delstart(N_pert);
    
     Logout("kfrom_MPl_lattice =  %2.5e, kto_MPl_lattice =  %2.5e, k_lattice_grid_min_MPl = %2.5e, floor(kfrom_MPl_lattice/k_lattice_grid_min_MPl) =  %2.5e\n", kfrom_MPl_lattice, kto_MPl_lattice, k_lattice_grid_min_MPl, floor(kfrom_MPl_lattice/k_lattice_grid_min_MPl));
    
    double k_comoving_start;
    
    if (k_lattice_grid_min_MPl < kfrom_MPl_lattice)
    {
        outrange_num = floor(kfrom_MPl_lattice/k_lattice_grid_min_MPl);
        
        latticerange_num = (N/2) - outrange_num;
        
        k_comoving_start = (outrange_num+1)*k_lattice_grid_min_MPl;
        
        Logout("k_lattice_grid_min_MPl < kfrom_MPl_lattice\n");
        Logout("outrange_num = %d, latticerange_num = %d \n\n",outrange_num, latticerange_num);

    }
    else{
        latticerange_num = N/2;
        k_comoving_start = k_lattice_grid_min_MPl;
        
        Logout("kfrom_MPl_lattice < k_lattice_grid_min_MPl\n\n");
         
    }
    
    Logout("-----------------------------------------------------\n\n");
    std::chrono::system_clock::time_point time_OSCEND_pert, time_UNPERT_pert, time_NEWINF_END_pert, time_END_pert;

    //LOOP
     for (int latticerange_loop = 0; latticerange_loop < latticerange_num; latticerange_loop++){
         
         
         
    
             percentage =  round( ( latticerange_loop + 1 )*100 / ( latticerange_num ) );
    
             k_comoving = k_comoving_start + k_lattice_grid_min_MPl*latticerange_loop;
             
         if(latticerange_loop==27){
            for (i=0;i<N_pert;i++) Logout("latticep[%d][%d] = %2.5e \n",latticerange_loop,i,latticep[latticerange_loop][i] );
         }
             
             knum = UC::kMPl_to_knum(k_comoving);
             Logout("%d/%d: knum = %d, kMpc = %2.5e, kMPl = %2.5e: \n", latticerange_loop+1, latticerange_num, knum, UC::kMPl_to_kMpc(k_comoving), k_comoving);
             
             //set delstart
             if (k_lattice_grid_min_MPl < kfrom_MPl_lattice)
             {
                 for (i=0;i<N_pert;i++) delstart[i]  =  latticep[outrange_num+latticerange_loop][i];
             }else{
            for (i=0;i<N_pert;i++) delstart[i]  = latticep[latticerange_loop][i];
                
             }
             
             //set efold
             xmid = log(a_lattice_end);
             
                 
             if(latticerange_loop==0)
             {
             Logout("OSCEND_EFOLD = %2.5e \n", xmid);
             Logout("UNPERT_EFOLD = %2.5e \n", UNPERT_EFOLD);
             }
                 if(latticerange_loop==27){
                 for (i=0;i<N_pert;i++)     Logout("delstart[%d] = %2.5e \n", i, delstart[i]);
                 }
         
         timecount = 1;
         xp2[timecount-1] = xmid;
         for (i=0;i<N_pert;i++) delp[i][timecount-1] = delstart[i];
         
         if(spectrum_afosc_switch){
             spectrum_output(new_filename_sp_afosc, xp2, delp, timecount, knum, k_comoving);
         };
         
         
         if(latticerange_loop==0){
            time_OSCEND_pert = std::chrono::system_clock::now();
         }
         
                 if(latticerange_loop==27){
                 for (i=0;i<N_pert;i++)     Logout("2: delstart[%d] = %2.5e \n", i, delstart[i]);
                 }
             if(xmid < UNPERT_EFOLD){
                 NR::odeintpert(delstart,xmid,UNPERT_EFOLD,epsosc,h2,hmin,nok,nbad,timecount,dxsav,full,NR::rkqs,k_comoving, &xp2, &delp, timecount_max_pert);
         //         std::cout << "timecount = " << timecount << std::endl;
                 
                 if(latticerange_loop==27){
                 for (i=0;i<N_pert;i++)     Logout("3: delstart[%d] = %2.5e \n", i, delstart[i]);
                 }
                 
                 if(kanalyze_switch){
                 kanalyze_output(new_dirname_k_lattice, xp2, delp, timecount, knum,k_comoving);
                 }
                 xmid=xp2[timecount-1];
                 a=exp(xmid);
                 for (i=0;i<N_pert;i++) delstart[i]=delp[i][timecount-1];
             };
             
             
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
            //UNPERT_EFOLD is set by hand according to the result for zero-mode.
         
         
         if(latticerange_loop==0){
         time_UNPERT_pert = std::chrono::system_clock::now();
         time_calc(time_OSCEND_pert, time_UNPERT_pert, "Perturbation OSCEND~UNPERT (Lattice Range)");
         }
         
            NR::odeintpert(delstart,UNPERT_EFOLD,NEWINF_END_EFOLD,epsnew,h2,hmin,nok,nbad,timecount,dxsav,newinf,NR::rkqs, k_comoving, &xp2, &delp, timecount_max_pert);
            //         std::cout << "timecount = " << timecount << std::endl;
             
            if(kanalyze_switch){
                kanalyze_output(new_dirname_k_lattice, xp2, delp, timecount, knum,k_comoving);
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
            //until the amplitude of oscillation of phi becomes sufficiently small at ln(a)=xend.
            //xend is set by hand according to the result for zero-mode.
            //(if perturbations of sigma and psi are not solved, superhorizon parturbations begin to decrease)
         
         if(latticerange_loop==0 ){
         time_NEWINF_END_pert = std::chrono::system_clock::now();
         time_calc(time_UNPERT_pert, time_NEWINF_END_pert,"Perturbation UNPERT~NEWINF_END (Lattice Range)");
         }
         
             NR::odeintpert(delstart,xmid,xend,epslast,h2,hmin,nok,nbad,timecount,dxsav,fixfix,NR::rkqs, k_comoving, &xp2, &delp, timecount_max_pert);
           
            //         std::cout << "timecount = " << timecount << std::endl;
            if(kanalyze_switch){
                kanalyze_output(new_dirname_k_lattice, xp2, delp, timecount, knum, k_comoving);
            }
    
    
            //     cout << "Final Result delstart[2] = " << tr2[2] << " delp[2][" << timecount-1 << "] = " << delp[2][timecount-1] << endl;
            //   for (i=0;i<N1;i++) cout << "Final Result tr2[" << i << "]" << tr2[i] << endl;
            if(spectrum_switch){
                spectrum_output(new_filename_sp_final, xp2, delp, timecount, knum, k_comoving);
            }
         
         if(latticerange_loop==0 ){
         time_END_pert = std::chrono::system_clock::now();
         time_calc(time_NEWINF_END_pert,time_END_pert,"Perturbation NEWINF_END~END (Lattice Range)");
             Logout( "\n" );
         }
    
            Logout( "Calculation %d%% Complete\n\n",percentage);
        }; //End 2nd for loop of k
    
    
}

