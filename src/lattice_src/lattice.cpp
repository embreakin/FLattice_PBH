//Doxygen
/**
* @file    lattice_initialize.cpp
* @brief    Lattice source file
* @author   Francis Otani
* @date
* @details
*/


#include <chrono>
#include "parameters.hpp"
#include <sys/stat.h>
#include "lattice_initialize.hpp"
#include "utilities.hpp"
#include "lattice.hpp"



double rescale_A;
double rescale_B;
double L;//L_pr
double dx; //dx_pr

void lattice(double** lattice_var)
{
    //Start Counting Time
    std::chrono::high_resolution_clock::time_point start, loop_start, current;
    std::chrono::milliseconds init_elapsed, elapsed;

    start = std::chrono::high_resolution_clock::now();


    //Output Data File/Directory Management
    dir_manage(new_dirname_ed);
    dir_manage(new_dirname_f);
    
    Logout( "\n----------------------------------------------\n" );
     
    double sigma_initial = 0;
    
    for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){
        
        sigma_initial += lattice_var[lattice_loop][0];
        //          Logout("latticep[%d][0] = %2.5e",lattice_loop, latticep[lattice_loop][0]);
        
    }
    sigma_initial /= (N/2);
    
     rescale_A = 1/sigma_initial;
    rescale_B = sqrt(V_11(0,FIXPSI,0))*sigma_initial;
     L = N*M_PI*rescale_B/(kto_MPl_lattice);
    
    double k_lattice_grid_min_pr = 2*M_PI/L;
    
    double k_lattice_grid_max_pr = N*M_PI/L;
    
    dx = 1.* L/N;
    
    
    Logout("kfrom_lattice =  %2.5e [MPl], kto_lattice =  %2.5e [MPl], k_lattice_grid_min [MPl] =  %2.5e \n\n",kfrom_MPl_lattice, kto_MPl_lattice, k_lattice_grid_min_MPl);
    
    
    Logout("kfrom_lattice =  %2.5e [Mpc^-1], kto_lattice =  %2.5e [Mpc^-1], k_lattice_grid_min =  %2.5e [Mpc^-1] \n\n",kfrom_Mpc_lattice, kto_Mpc_lattice, kto_Mpc_lattice/(N/2));
    
    Logout("kfrom_lattice_pr =  %2.5e, kto_lattice_pr =  %2.5e, k_lattice_grid_min_pr =  %2.5e \n\n",kfrom_MPl_lattice/rescale_B, kto_MPl_lattice/rescale_B, k_lattice_grid_min_pr);
    
    
    Logout("sigma_initial = %2.5e \n\n", sigma_initial);
    Logout("L = %2.5e [Mpc], dx = %2.5e [Mpc] \n\n", UC::xMPl_to_xMpc(L/rescale_B), UC::xMPl_to_xMpc(dx/rescale_B));
    Logout("Range of k_pr in lattice: %2.5e <= |k_pr| <= %2.5e \n\n", k_lattice_grid_min_pr, k_lattice_grid_max_pr);
    

    Logout( "\n----------------------------------------------\n" );
    Logout( "            SIMULATION PARAMETERS            \n\n" );

    Logout( " Dimension    =  %d\n", dim );
    Logout( " Box Size   =  %2.2e\n", L );
    Logout( " Grid Number   =  %d\n", N );
    Logout( " Initial Time =  %4.1lf\n", t0 );
    Logout( " Final Time   =  %4.2e\n", t0 + total_step*dt );
    Logout( " dt_pr     =  %2.2e\n", dt );
    Logout( " dx_pr     =  %2.2e\n", dx );
    if(dt/dx <  1/sqrt(dim)){
    Logout( " dt_pr/dx_pr     =  %2.2e < 1/sqrt(%d) =  %2.2e \n", dt/dx ,dim, 1/sqrt(dim));}
    Logout( " Number of fields   =  %d\n", num_fields );
    Logout( " Number of threads  =  %d\n", num_threads );
    Logout("  rescale_A = %2.2e \n", rescale_A);
    Logout("  rescale_B = %2.2e \n\n", rescale_B);

    //--------------------------------------------------
    //       SETTING INITIAL CONDITIONS
    //--------------------------------------------------
    Logout("----------------------------------------------\n");
    Logout("            STARTING INITIALIZATION                \n\n");

    //lattice_var[knum_lattice][j] [j] specifies variables as follows:
    //    j=0-2  : zero modes of inflaton sigma, psi, and phi
    //    j=3-5  : t derivatives of inflaton sigma', psi', phi'
    //    j=6    : energy density of radiation
    //    j=7-15 : mode functions of field perturbation: delta_{sigma,sigma}, delta_{sigma,psi}, delta{sigma,phi}, delta_{psi,sigma}, etc...
    //    j=16-24: t derivatives of mode functions
    //    j=25-27: mode functions of gravitational potential perturbation: delta Phi_{sigma}, delta Phi_{psi}, delta Phi_{phi}
    //    j=28-30: t derivatives of gravitational potential perturbation
    //    j=31-54: complex conjugate of [7]-[30]
//    for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){
//
//     for (int i=0;i<N_zero;i++) Logout("lattice_var[%d][%d] = %2.5e \n",lattice_loop, i , lattice_var[lattice_loop][i] );
//
//    }
    //Declare fields and their derivatives
    double **f;
    double **df;
    //Declare radiation
    double radiation = 0;
//
    //Instantiate fields
    Field field;
//
    Logout("radiation1 = %2.5e \n", radiation);
//    Initialize fields and their derivatives
   
    initialize( f, df, &field, radiation, lattice_var);
    
    
    Logout("radiation2 = %2.5e \n", radiation);
    
    Logout( "f[0][4] = %2.5e \n", f[0][4]);
    //Instantiate leapfrog and initialize
    LeapFrog leapfrog(&field, f, df, radiation);
   
////    //Instantiate energy and initialize
    Energy energy;
////
////    // 1: self-consistent, 2: radiation dominant, 3: matter dominant -> No Output of vti files for the field
        //  0: no expansion -> Output vti files for the field
//    if(expansion){
//    }else{
//        write_VTK_f(new_dirname_f, f[0], "field", -1 );
//    }
    Logout("radiation3 = %2.5e \n", radiation);
    // Calculate all necessary initial data regarding energy density
    energy.energy_calc( &field, &leapfrog, f, df ,radiation);

    Logout("radiation4 = %2.5e \n", radiation);
    // Output vti files for the energy density of the field
    write_VTK_ed( new_dirname_ed, energy.value, "energy", -1  );
    Logout("write_status \n");
    // Write data to status.txt
    write_status( new_filename_lattice, &field, &leapfrog, &energy, f, df, t0 );
    // add data to kanalyze files (power spectrum)
    kanalyze_output_lattice(new_dirname_k_lattice, &field, &leapfrog, &energy, f, df);
    
    //Calculate Initialization Time
    current = std::chrono::high_resolution_clock::now();
    init_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( current - start );
    double hourresidue = fmod(init_elapsed.count()*1.e-3, 3600);
    Logout( " Initialization Time: %d h %d m %2.3f s \n", int(init_elapsed.count()*1.e-3)/3600,int(hourresidue)/60, fmod(hourresidue, 60.0));
////    //--------------------------------------------------
////    //       TIME ITERATION LOOP
////    //--------------------------------------------------
////
    Logout("----------------------------------------------\n");
    Logout("            STARTING TIME ITERATION LOOP                \n\n");

    int latticeloop_count = 1;
    int latticeloop_interval = max_loop/screen_latticeloop_number;

    for( int loop = 0; loop < max_loop; ++loop ) // This many times vti files will be created
    {
        double t = t0 + loop*output_step*dt;

        if((loop+1) - 1 == latticeloop_interval*(latticeloop_count-1))
        {
            loop_start = std::chrono::high_resolution_clock::now();
        }
        


        for( int st_loop = 0; st_loop < st_max_loop; ++st_loop ) // This many times data will be added to status.txt between the output of vti files
        {
        // 1: self-consistent, 2: radiation dominant, 3: matter dominant -> Compute using the function that takes expansion into account
        //  0: no expansion -> Compute using the function that doesn't take expansion into account
            if(expansion){

            leapfrog.evolution_expansion( &field, f, df, radiation );

            }else{
                //leapfrog.evolution( &field, f, df );
            }

            // Calculate all necessary data regarding energy density.
            energy.energy_calc( &field, &leapfrog, f, df , radiation );
            //  std::cout << "t1 = " << t << std::endl;

            // Add data to status.txt
            write_status( new_filename_lattice, &field, &leapfrog, &energy, f, df,  t+st_output_step*dt );
            //  std::cout << "t2 = " << t << std::endl;

            // Evolve time
            t = t + st_output_step*dt;
            // std::cout << "t3 = " << t << std::endl;
        }


        // 1: self-consistent, 2: radiation dominant, 3: matter dominant -> No Output of vti files for the field
        //  0: no expansion -> Output vti files for the field
//        if(expansion){
//        }else{
//            write_VTK_f( new_dirname_f, f[0], "field", loop );
//        }

        // Output vti files for the energy density of the field
        write_VTK_ed( new_dirname_ed, energy.value, "energy", loop );
        
        // add data to kanalyze files (power spectrum)
        kanalyze_output_lattice(new_dirname_k_lattice, &field, &leapfrog, &energy, f, df);
        
        current = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( current - loop_start );

        if(loop)
        {
                
                if(loop+1 == latticeloop_interval*latticeloop_count)
                {
                    Logout( " Loop %d/%d: %2.3f s \n", latticeloop_interval*latticeloop_count, max_loop, elapsed.count()*1.e-3 );
                    
                    latticeloop_count++;
                }

        }else{ // This branch is run only for the first loop (loop = 0) to calculate the total estimate time of the simulation judging from the first loop and the initialization time.
            double estimate_time = (init_elapsed.count()+elapsed.count()*max_loop)*1.e-3;

            hourresidue = fmod( estimate_time , 3600);

            Logout( " Loop %d/%d: %2.3f s \t Estimated Total Time: %d h %d m %2.3f s \n", loop+1, max_loop, elapsed.count()*1.e-3, int(estimate_time)/3600,int(hourresidue)/60, fmod(hourresidue, 60.0) );
        }

    }

    
    //Calculate total elapsed time
    current = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( current - start );
    hourresidue = fmod(elapsed.count()*1.e-3, 3600);
    Logout( " Total Time: %d h %d m %2.3f s \n", int(elapsed.count()*1.e-3)/3600,int(hourresidue)/60, fmod(hourresidue, 60.0));


    //Release all memory of fields and their derivatives
    field.finalize( f, df, &leapfrog, radiation, lattice_var );

    
//    for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){
//
//
//        for (int i=0;i<N_pert;i++) Logout("lattice_var[%d][%d] = %2.5e \n",lattice_loop, i , lattice_var[lattice_loop][i] );
//
//    }

    Logout( "\n----------------------------------------------\n" );
    Logout( "            SIMULATION PARAMETERS            \n\n" );
    
    Logout( " Dimension    =  %d\n", dim );
    Logout( " Box Size  =  %2.2e\n", L );
    Logout( " Grid Number   =  %d\n", N );
    Logout( " Initial Time =  %4.1lf\n", t0 );
    Logout( " Final Time   =  %4.2e\n", t0 + total_step*dt );
    Logout( " dt_pr     =  %2.2e\n", dt );
    Logout( " dx_pr     =  %2.2e\n", dx );
    if(dt/dx <  1/sqrt(dim)){
        Logout( " dt/dx     =  %2.2e < 1/sqrt(%d) =  %2.2e \n", dt/dx ,dim, 1/sqrt(dim));}
    Logout( " Number of fields   =  %d\n", num_fields );
    Logout( " Number of threads  =  %d\n\n", num_threads );
    Logout("----------------------------------------------\n");
    Logout( "              PROGRAM FINISHED            \n\n" );
    Logout( "----------------------------------------------\n\n" );
//
}
