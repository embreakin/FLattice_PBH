//Doxygen
/**
* @file    main.cpp
* @brief    Main source file
* @author   Francis Otani
* @date
* @details
*/

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
#include "calculation.hpp"





//xp[i] stores integration variable (log(a)) for each steps[i].
//delp_p[i][j] stores variables for each steps[i]. [j] specifies variables as follows:
//    j=0-2  : zero modes of inflaton sigma, psi, and phi
//    j=3-5  : log(a) derivatives of (<-this part is probably wrong) inflaton sigma', psi', phi'
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
    
    Logout("=====================================================\n");
    
   
    //Output Data File/Directory Management
    dir_manage(new_dirname_k);
    dir_manage(new_dirname_k_lattice);
    
    
    //------------------------
    //Name of Parameter Set
    //------------------------
    
    Logout("=====================================================\n\n");
    
    Logout("Parameter Set: %s \n\n", par_set_name.c_str() );
    
    Logout("=====================================================\n\n");
    //------------------------
    //calculation of zeromode
    //------------------------
    if(zeromode_switch){
    Logout("=====================================================\n\n");
     
    Logout( "Calculating Zero mode...\n\n");
    
    //Instantiate zeromode
    Zeromode Zero;
    
    //Calculate zeromode
    Zero.zeromode_calc();
    

    Logout( "\nZero mode Calculation Complete\n\n");
   
    }
    
    
    
    
    //---------------------------------------
    //calculation of zeromode w/ perturbation
    //---------------------------------------
     if(perturbation_switch){
         
    Logout("=====================================================\n");
    Logout("=====================================================\n\n");
    Logout( "Calculating Zero mode with Perturbation...\n\n");
    Logout("=====================================================\n");
    Logout("=====================================================\n\n");

     //Instantiate zeromode
     Zeromode Zero2;
    //Instantiate perturbabtion
     Perturbation Perturb;

    
    //Calculate with perturbation
  if(latticerange_switch){ // When there is lattice range

         Logout("Start knum lower range\n\n");
         Logout("-----------------------------------------------------\n\n");

          // convert to knum units

          Logout("Range: kfrom_knum = %d, kfrom_knum_lattice = %d, kinterval_knum = %d \n\n",kfrom_knum, kfrom_knum_lattice,kinterval_knum);
    //
          Perturb.nonlatticerange_calc(kfrom_knum, kfrom_knum_lattice, Zero2);
    //
          Logout("=====================================================\n\n");
          Logout( "Start Lattice Range\n\n");

          double **latticep;
      
          //Initialize latticep
          Perturb.lattice_initialize(latticep);
      
      
          Logout("=====================================================\n\n");
          Logout( "Start loop calculation up to OSCSTART\n\n");
          Logout("-----------------------------------------------------\n\n");
      
      
          Perturb.latticerange_firsthalf_calc(latticep, Zero2);

          Logout("-----------------------------------------------------\n\n");
          Logout("Start Lattice Simulation\n\n");
          Logout("-----------------------------------------------------\n\n");
      
      
      
          lattice(latticep);
      
          Logout("-----------------------------------------------------\n\n");
          Logout("Start loop calculation after Lattice Simulation\n\n");
          Logout("-----------------------------------------------------\n\n");
      
//      for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){
//
//          for (int i=0;i<N_zero;i++) Logout("lattice_var[%d][%d] = %2.5e \n",lattice_loop, i , latticep[lattice_loop][i] );
//
//      }

          Perturb.latticerange_secondhalf_calc(latticep);

          Logout("=====================================================\n\n");
          Logout("Start knum upper range\n\n");
          Logout("-----------------------------------------------------\n\n");


          Logout("Range: kstart_knum = %d, kto_knum = %d, kinterval_knum = %d \n\n",kstart_knum, kto_knum, kinterval_knum);

          Perturb.nonlatticerange_calc(kstart_knum, kto_knum, Zero2);

          Perturb.lattice_finalize(latticep);



        }else{ // When there is no lattice range
            
            if(lattice_kmodes_switch)
            {
                 Perturb.nonlatticerange_calc(kfrom_knum, kto_knum, Zero2);
            }else
            {
                Logout("Range: kfrom_knum = %d, kto_knum = %d, kinterval_knum = %d \n\n",kfrom_knum, kto_knum,kinterval_knum);
          Perturb.nonlatticerange_calc(kfrom_knum, kto_knum, Zero2);
            }
       }
         
     }
    
    time_end = std::chrono::system_clock::now();       // End measuring elapsed time
//    int time_hours = std::chrono::duration_cast<std::chrono::hours>(time_end - time_start).count(); // Casting time
//    int time_minutes = std::chrono::duration_cast<std::chrono::minutes>(time_end - time_start).count() - time_hours*60; // Casting time
//    int time_seconds = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start).count() - time_hours*60*60 - time_minutes*60; // Casting time
//    int time_days = time_hours / 24;
//    time_hours = time_hours % 24;
//     Logout("=====================================================\n\n");
//     Logout( "Total Computation Time: %d d %d h %d m %d s\n\n",time_days,time_hours,time_minutes,time_seconds);
//     Logout("=====================================================\n");
   
    Logout("=====================================================\n");
    time_calc(time_start, time_end, "Total Computation Time");
    Logout("=====================================================\n");
//
//
//
//    delete latticep_p;
//    delete delp_p;
//    delete xp2_p;
    return 0;
}
