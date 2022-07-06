#include <chrono>
#include "parameters.hpp"
#include <sys/stat.h>
#include "utilities.hpp"
#include "lattice.hpp"

double rescale_A;
double rescale_B;
double L;//L_pr
double dx; //dx_pr

void lattice(double**& lattice_var)
{
    //Start Counting Time
    std::chrono::high_resolution_clock::time_point start, loop_start, current;
    std::chrono::milliseconds init_elapsed, elapsed;

    start = std::chrono::high_resolution_clock::now();


    //Output Data File/Directory Management
    dir_manage(exist_dirname_ed, new_dirname_ed);
    dir_manage(exist_dirname_f, new_dirname_f);
    file_manage(exist_filename_status);
    
    
    
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
    
    
    Logout("sigma_initial = %2.5e \n", sigma_initial);
    Logout("rescale_A = %2.5e \n", rescale_A);
    Logout("kfrom_MPl_lattice =  %2.5e, kto_MPl_lattice =  %2.5e, k_lattice_grid_min_MPl =  %2.5e \n\n",kfrom_MPl_lattice, kto_MPl_lattice, k_lattice_grid_min_MPl);
    Logout("kfrom_pr_lattice =  %2.5e, kto_pr_lattice =  %2.5e, k_lattice_grid_min_pr =  %2.5e \n\n",kfrom_MPl_lattice/rescale_B, kto_MPl_lattice/rescale_B, k_lattice_grid_min_pr);
    Logout("rescale_B = %2.5e, L_pr = %2.5e, N = %d \n\n",rescale_B, L, N);
    Logout("Range of k_pr in lattice: %2.5e <= |k_pr| <= %2.5e \n\n", k_lattice_grid_min_pr, k_lattice_grid_max_pr);
    

    Logout( "\n----------------------------------------------\n" );
    Logout( "            SIMULATION PARAMETERS            \n\n" );

    Logout( " Dimension    =  %d\n", dim );
    Logout( " Box Size     =  %d\n", L );
    Logout( " Grid Number   =  %d\n", N );
    Logout( " Initial Time =  %4.1lf\n", t0 );
    Logout( " Final Time   =  %4.2e\n", t0 + total_step*dt );
    Logout( " dt     =  %2.2e\n", dt );
    Logout( " dx     =  %2.2e\n", dx );
    if(dt/dx <  1/sqrt(dim)){
    Logout( " dt/dx     =  %2.2e < 1/sqrt(%d) =  %2.2e \n", dt/dx ,dim, 1/sqrt(dim));}
    Logout( " Number of fields   =  %d\n", num_fields );
    Logout( " Number of threads  =  %d\n\n", num_threads );


    //--------------------------------------------------
    //       SETTING INITIAL CONDITIONS
    //--------------------------------------------------
    Logout("----------------------------------------------\n");
    Logout("            STARTING INITIALIZATION                \n\n");

    //lattice_var[knum_lattice][j] [j] specifies variables as follows:
    //    j=0-2  : zero modes of inflaton sigma, psi, and phi
    //    j=3-5  : log(a) derivatives of inflaton sigma', psi', phi'
    //    j=6    : energy density of radiation
    //    j=7-15 : mode functions of field perturbation: delta_{sigma,sigma}, delta_{sigma,psi}, delta{sigma,phi}, delta_{psi,sigma}, etc...
    //    j=16-24: log(a) derivatives of mode functions
    //    j=25-27: mode functions of gravitational potential perturbation: delta Phi_{sigma}, delta Phi_{psi}, delta Phi_{phi}
    //    j=28-30: log(a) derivatives of gravitational potential perturbation
    //    j=31-54: complex conjugate of [7]-[30]
    for (int lattice_loop = 0; lattice_loop < N/2; lattice_loop++){

     for (int i=0;i<N_zero;i++) Logout("lattice_var[%d][%d] = %2.5e \n",lattice_loop, i , lattice_var[lattice_loop][i] );

    }
    //Declare fields and their derivatives
    double **f, **df;
    //Declare radiation
    double radiation;
//
    //Instantiate fields
    Field field;
//
//    Initialize fields and their derivatives
    initialize( f, df, &field, radiation, lattice_var);
//
    //Instantiate leapfrog and initialize
    LeapFrog leapfrog(&field, f, df, radiation);
   
////    //Instantiate energy and initialize
////    Energy energy;
////
////    // 1: self-consistent, 2: radiation dominant, 3: matter dominant -> No Output of vti files for the field
////    //  0: no expansion -> Output vti files for the field
////    if(expansion){
////    }else{
////       write_VTK_f(new_dirname_f, f[0], "field", -1 );
////    }
////
////    // Calculate all necessary initial data regarding energy density
////    energy.energy_calc( &field, &leapfrog, f, df );
////
////    // Output vti files for the energy density of the field
////    write_VTK_ed( new_dirname_ed, energy.value[0], "energy", -1  );
////
////    // Write data to status.txt
////    write_status( new_filename_status, &field, &leapfrog, &energy, f, t0 );
////
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



    for( int loop = 0; loop < max_loop; ++loop ) // This many times vti files will be created
    {
////        double t = t0 + loop*output_step*dt;
////
////        loop_start = std::chrono::high_resolution_clock::now();
////
////
////        for( int st_loop = 0; st_loop < st_max_loop; ++st_loop ) // This many times data will be added to status.txt between the output of vti files
////        {
////        // 1: self-consistent, 2: radiation dominant, 3: matter dominant -> Compute using the function that takes expansion into account
////        //  0: no expansion -> Compute using the function that doesn't take expansion into account
////            if(expansion){
////
////            leapfrog.evolution_expansion( &field, f, df, t );
////
////            }else{
////                leapfrog.evolution( &field, f, df );
////            }
////
////            // Calculate all necessary data regarding energy density.
////            energy.energy_calc( &field, &leapfrog, f, df );
////            //  std::cout << "t1 = " << t << std::endl;
////
////            // Add data to status.txt
////            write_status( new_filename_status, &field, &leapfrog, &energy, f, t+st_output_step*dt );
////            //  std::cout << "t2 = " << t << std::endl;
////
////            // Evolve time
////            t = t + st_output_step*dt;
////            // std::cout << "t3 = " << t << std::endl;
////        }
////
////
////        // 1: self-consistent, 2: radiation dominant, 3: matter dominant -> No Output of vti files for the field
////        //  0: no expansion -> Output vti files for the field
////        if(expansion){
////        }else{
////            write_VTK_f( new_dirname_f, f[0], "field", loop );
////        }
////
////        // Output vti files for the energy density of the field
////        write_VTK_ed( new_dirname_ed, energy.value[0], "energy", loop );
////
////
////        //Calculate the elapsed time between the output of vti files
////        current = std::chrono::high_resolution_clock::now();
////        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( current - loop_start );
////
        if(loop)
        {
            Logout( " Loop %d/%d: %2.3f s \n", loop+1, max_loop, elapsed.count()*1.e-3 );

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
    finalize( f, df );

    Logout( "\n----------------------------------------------\n" );
    Logout( "            SIMULATION PARAMETERS            \n\n" );

    Logout( " Dimension    =  %d\n", dim );
    Logout( " Box Size     =  %d\n", L );
    Logout( " Grid Number   =  %d\n", N );
    Logout( " Initial Time =  %4.1lf\n", t0 );
    Logout( " Final Time   =  %4.2e\n", t0 + total_step*dt );
    Logout( " dt     =  %2.2e\n", dt );
    Logout( " dx     =  %2.2e\n", dx );
    if(dt/dx <  1/sqrt(dim)){
        Logout( " dt/dx     =  %2.2e < 1/sqrt(%d) =  %2.2e \n", dt/dx ,dim, 1/sqrt(dim));}
    Logout( " Number of fields   =  %d\n", num_fields );
    Logout( " Number of threads  =  %d\n\n", num_threads );
    Logout("----------------------------------------------\n");
    Logout( "              PROGRAM FINISHED            \n\n" );
    Logout( "----------------------------------------------\n\n" );
//
}
