#include <chrono>
#include "parameter.hpp"
//#include <filesystem>
#include <boost/filesystem.hpp>
//#include <experimental/filesystem>
#include <sys/stat.h>
#include "utilities.hpp"

namespace fs = boost::filesystem;
//namespace fs = std::filesystem;

int main( int argc, char** argv )
{   
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
	
	double **f, **df;
    Field field;
    LeapFrog leapfrog;
    
	initialize( f, df, &field);
   // fs::remove_all("data");
  // std::experimental::filesystem::path::remove("../data");
    
    const fs::path path("../data");
    try {
        fs::remove_all(path);
    }
    catch (fs::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }
   
    boost::system::error_code error;
    const bool result = fs::create_directory(path, error);
    if (!result || error) {
        std::cout << "failed to create directory" << std::endl;
    }
     /*
    if(mkdir("../data",0755)==0){
        printf("successfully created directory\n");
    }else{
        printf("failed to create directory\n");
    }*/
    // fs::remove_all("data");
    //fs::create_directory("data");
     /*
    const fs::path path("data2");
    try {
        fs::remove_all(path);
    }
    catch (fs::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }
    
    boost::system::error_code error;
    const bool result = fs::create_directory(path, error);
    if (!result || error) {
        std::cout << "failed to create directory" << std::endl;
    }*/
    //std::filesystem::remove_all("data");
   // std::filesystem::create_directory("data");
    if(expansion){
    }else{
       write_VTK_f( f[0], "field", -1 );
    }
    
        
	Energy energy( &field, &leapfrog, f, df );
   write_VTK_ed( energy.value[0], "energy", -1  );
    
	write_status( &field, &leapfrog, &energy, f, t0 );

	//--------------------------------------------------
	//       TIME ITERATION LOOP
	//--------------------------------------------------

	Logout("----------------------------------------------\n");
	Logout("            STARTING COMPUTATION                \n\n");

	std::chrono::high_resolution_clock::time_point start, loop_start, current;
	std::chrono::milliseconds elapsed;

	start = std::chrono::high_resolution_clock::now();
    
    for( int loop = 0; loop < max_loop; ++loop )
	{
        //std::cout << leapfrog.da() << std::endl;
	    double t = t0 + loop*output_step*dt;
	    loop_start = std::chrono::high_resolution_clock::now();

            for( int st_loop = 0; st_loop < st_max_loop; ++st_loop ){
            
        if(expansion){
            
        leapfrog.evolution_expansion( &field, f, df, t );
            
        }else{
            leapfrog.evolution( &field, f, df );
        }
		
            Energy energy( &field, &leapfrog, f, df );
        //  std::cout << "t1 = " << t << std::endl;
           write_status( &field, &leapfrog, &energy, f, t+st_output_step*dt );
          //  std::cout << "t2 = " << t << std::endl;
             t = t + st_output_step*dt;
           // std::cout << "t3 = " << t << std::endl;
        }
		  Energy energy( &field, &leapfrog, f, df );
        
        if(expansion){
        }else{
            write_VTK_f( f[0], "field", loop );
        }
        
		write_VTK_ed( energy.value[0], "energy", loop );
       
		//write_status( &field, &leapfrog, &energy, f, t+output_step*dt );

        current = std::chrono::high_resolution_clock::now();
	    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( current - loop_start );
        Logout( " Loop %d/%d: %2.3f s \n", loop+1, max_loop, elapsed.count()*1.e-3 );
    }

    current = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( current - start );
    double hourresidue = fmod(elapsed.count()*1.e-3, 3600);
    Logout( " Total time: %d h %d m %2.3f s \n", int(elapsed.count()*1.e-3)/3600,int(hourresidue)/60, fmod(hourresidue, 60.0));
    
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
}
