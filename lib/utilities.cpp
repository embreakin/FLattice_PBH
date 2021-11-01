#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <sstream>
#include <fftw3.h>
#include "parameter.hpp"
#include "utilities.hpp"
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;


void dir_manage(const std::string exist_dir, const std::string new_dir )
{
    
    //Path of the existing data directory
    const fs::path existpath("../" + exist_dir );
    
    //Tries to remove all files in there. If it fails, it throws an error.
    try {
        fs::remove_all(existpath);
    }
    catch (fs::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }
    
    //Path of the newly set directory
    const fs::path newpath("../" + new_dir );
    
    //Tries to create a directory. If it fails, it throws an error.
    boost::system::error_code error;
    const bool result = fs::create_directory(newpath, error);
    if (!result || error) {
        std::cout << "failed to create directory" << std::endl;
    }
    
}

void file_manage(const std::string exist_file)
{
    //Path of the existing status file
    const fs::path statuspath("../" + exist_file);
    
    //Tries to remove the file. If it fails, it throws an error.
    try {
        fs::remove(statuspath);
    }
    catch (fs::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }
    
}


double rand_uniform(void)
{
    std::random_device rnd;
    std::mt19937 mt( rnd() );
    std::uniform_real_distribution<> rand( 0, 1 );
 //   std::cout << rand(mt) << "\n";
    return (rand(mt));
}

void set_mode(double p2, double m2, double *field, double *deriv, int real)
{
    double phase,phase2, amplitude, rms_amplitude, omega;
    double re_f_left, im_f_left, re_f_right, im_f_right;
#if  dim==1
    static double norm = m*pow(L/(dx*dx),.5)/sqrt(4*M_PI);
#elif  dim==2
    static double norm =  m*pow(L/(dx*dx),1)/(sqrt(2*M_PI));
#elif  dim==3
    static double norm =  m*pow(L/(dx*dx),1.5)/sqrt(2);
#endif
    static int tachyonic = 0; //Avoid printing the same error repeatedly
    
    if(p2+m2>0)
        omega=sqrt(p2+m2);
    else
    {
        if(tachyonic==0)
            std::cout <<"Warning: Tachyonic mode(s) may be initialized inaccurately"<< std::endl;
        omega=sqrt(p2);
        tachyonic=1;
    }
    
    if(omega>0.)
        rms_amplitude=norm/sqrt(omega)*pow(p2,.75-(double)dim/4.);
    else
        rms_amplitude=0.;
        
        //Amplitude = RMS amplitude x Rayleigh distributed random number
        // The same amplitude is used for left and right moving waves to generate standing waves. The extra 1/sqrt(2) normalizes the initial occupation number correctly.
        
    amplitude = rms_amplitude/sqrt(2.)*sqrt(log(1./rand_uniform()));
    phase = 2*M_PI*rand_uniform();
  // std::cout << "phase1 " << phase/(2*M_PI) << std::endl;
        //Left moving component
        re_f_left = amplitude * cos( phase );
        im_f_left = amplitude * sin( phase );
        //Right moving component
    phase2 = 2*M_PI*rand_uniform();
  //  std::cout << "phase2 " << phase/(2*M_PI) << std::endl;
        re_f_right = amplitude * cos( phase2 );
        im_f_right = amplitude * sin( phase2 );
    
    field[0] = re_f_left + re_f_right;
    field[1] = im_f_left + im_f_right;
    deriv[0] = omega*(im_f_left - im_f_right);
    deriv[1] = -omega*(re_f_left - re_f_right);
    if(real==1)
    {
        field[1]=0;
        deriv[1]=0;
    }
    return;
   }

void DFT_c2rD1( double* f)
{
    fftw_plan p;
    double* out;
    fftw_complex *in;
    size_t in_size;
    
    in_size = sizeof(fftw_complex) * (N/2+1);
    in  = (fftw_complex*)fftw_malloc( in_size );
    out = new double [N]();
    p = fftw_plan_dft_c2r_1d( N, in, out, FFTW_ESTIMATE );
    
        for( int j = 0; j < N/2+1; ++j ){
            int idx = j;
            if(idx==0)
            {
                in[idx][0] = f[idx]  ;
                in[idx][1] = 0  ;
                
            }else if(idx==N/2){
                in[idx][0] = f[1]  ;
                in[idx][1] = 0  ;
            }else{
                in[idx][0] = f[2*idx]  ;
                in[idx][1] = f[2*idx+1]  ;
            }
            
        }
    
    fftw_execute(p);
    
    // Set output data
    for( int j = 0; j < N; ++j ){
        int idx = j;
        f[idx] = out[idx]/N;
    }
    
    if( p ) fftw_destroy_plan(p);
    if( in ) fftw_free(in);
    delete[] out;
}
/*
void DFT_c2rD2d( double* df,double* fdnyquist )
{
    fftw_plan p;
    double* out;
    fftw_complex *in;
    size_t in_size;
    
    in_size = sizeof(fftw_complex) * N * (N/2+1);
    in  = (fftw_complex*)fftw_malloc( in_size );
    out = new double [N*N]();
    p = fftw_plan_dft_c2r_2d( N, N, in, out, FFTW_ESTIMATE );
    
    
    
    for( int j = 0; j < N; ++j ){
        //#pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int k = 0; k < N/2+1; ++k ){
            int idx = (N/2+1)*j + k;
            // int idx = j*N + k;
            if(k == N/2){
                in[idx][0] = (double)fdnyquist[2*j] ;
                in[idx][1] = (double)fdnyquist[2*j+1] ;
               // std::cout << "in[" << idx << "][0] = " << in[idx][0] << std::endl;
               // std::cout << "in[" << idx << "][1] = " << in[idx][1] << std::endl;
            }else {
                in[idx][0] =  (double)df[j*N + 2*k] ;
                in[idx][1] = (double)df[j*N + 2*k+1] ;
               // std::cout << "in[" << idx << "][0] = " << in[idx][0] << std::endl;
               // std::cout << "in[" << idx << "][1] = " << in[idx][1] << std::endl;
            }
        }
    }
    fftw_execute(p);
    
    // Set output data
    
    for( int j = 0; j < N; ++j ){
        //#pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int k = 0; k < N; ++k ){
            int idx = j*N + k;
            std::cout << "out[" << idx << "] = " << out[idx] << std::endl;
            df[idx] = out[idx]/(N*N);
        }
    }
    
    if( p ) fftw_destroy_plan(p);
    if( in ) fftw_free(in);
    delete[] out;
}
*/
void DFT_c2rD2( double* f,double* fnyquist )
{
    fftw_plan p;
    double* out;
    fftw_complex *in;
    size_t in_size;
    
    in_size = sizeof(fftw_complex) * N * (N/2+1);
    in  = (fftw_complex*)fftw_malloc( in_size );
    out = new double [N*N]();
    p = fftw_plan_dft_c2r_2d( N, N, in, out, FFTW_ESTIMATE );
    
    
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N/2+1; ++k ){
                int idx = (N/2+1)*j + k;
              //  int idx = j*N + k;
                if(k == N/2){
                    in[idx][0] =  fnyquist[2*j] ;
                    in[idx][1] =  fnyquist[2*j+1] ;
                }else {
                    in[idx][0] =  f[j*N + 2*k]; 
                    in[idx][1] =  f[j*N + 2*k+1] ;
                }
            }
        }
        fftw_execute(p);
        
        // Set output data

     for( int j = 0; j < N; ++j ){
        for( int k = 0; k < N; ++k ){
            int idx = j*N + k;
            f[idx] = out[idx]/(N*N);
        }
     }
    
    if( p ) fftw_destroy_plan(p);
    if( in ) fftw_free(in);
    delete[] out;
}




void DFT_c2rD3( double* f,double** fnyquist )
{
    fftw_plan p;
    double* out;
    fftw_complex *in;
    size_t in_size;
    
    in_size = sizeof(fftw_complex) * N * N * (N/2+1);
    in  = (fftw_complex*)fftw_malloc( in_size );
    out = new double [N*N*N]();
    p = fftw_plan_dft_c2r_3d( N, N, N, in, out, FFTW_ESTIMATE );
    
   
//#pragma omp parallel for simd collapse(3) schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            for( int k = 0; k < N; ++k ){
                for( int l = 0; l < N/2+1; ++l ){
                    int idx = (j*N + k)*(N/2+1) + l;
                    if(l == N/2){
                        in[idx][0] = fnyquist[j][2*k] ;
                        in[idx][1] = fnyquist[j][2*k+1] ;
                    }
                    else{
                        in[idx][0] =  f[(j*N + k)*N + 2*l] ;
                        in[idx][1] =  f[(j*N + k)*N + 2*l+1] ;
                    }
                }
            }
        }
        fftw_execute(p);
        
        // Set output data


         for( int j = 0; j < N; ++j ){
        for( int k = 0; k < N; ++k ){
            for( int l = 0; l < N; ++l ){
                int idx = (j*N + k)*N + l;
                f[idx] = out[idx]/(N*N*N);
            }
        }
    }
    
    if( p ) fftw_destroy_plan(p);
    if( in ) fftw_free(in);
    delete[] out;
}

/*
void DFT_c2r( double** f,double* fnyquist )
{
	  fftw_plan p;
		double* out;
    fftw_complex *in;
		size_t in_size;
	  
		switch( dim )
		{
			case 1:
				in_size = sizeof(fftw_complex) * (N/2+1);
				in  = (fftw_complex*)fftw_malloc( in_size );
				out = new double [N]();
				p = fftw_plan_dft_c2r_1d( N, in, out, FFTW_ESTIMATE );
				break;
			case 2:
				in_size = sizeof(fftw_complex) * N * (N/2+1);
				in  = (fftw_complex*)fftw_malloc( in_size );
				out = new double [N*N]();
				p = fftw_plan_dft_c2r_2d( N, N, in, out, FFTW_ESTIMATE );
				break;
			case 3:
				in_size = sizeof(fftw_complex) * N * N * (N/2+1);
				in  = (fftw_complex*)fftw_malloc( in_size );
				out = new double [N*N*N]();
				p = fftw_plan_dft_c2r_3d( N, N, N, in, out, FFTW_ESTIMATE );
				break;
		}

		for( int i = 0; i < num_fields; ++i )
		{
				
				// Create input data
				switch( dim ){
					case 1:
						#pragma omp parallel for schedule( static ) num_threads( num_threads )
						for( int j = 0; j < N/2+1; ++j ){
							int idx = j;
							if(idx==0)
                            {
                                in[idx][0] = f[i][idx]  ;
                                in[idx][1] = 0  ;
                                
                            }else if(idx==N/2){
                                in[idx][0] = f[i][1]  ;
                                in[idx][1] = 0  ;
                            }else{
                                in[idx][0] = f[i][2*idx]  ;
                                in[idx][1] = f[i][2*idx+1]  ;
                            }
							
						}
						break;
					case 2:
						#pragma omp parallel for schedule( static ) num_threads( num_threads )
						for( int j = 0; j < N; ++j ){
							for( int k = 0; k < N/2+1; k=k+2 ){
									int idx = j*N + k;
                                if(k == N/2){
                                    in[idx][0] = fnyquist[2*j] ;
                                    in[idx][1] = fnyquist[2*j+1] ;
                                }else {
                                    in[idx][0] =  f[i][idx] ;
                                    in[idx][1] =  f[i][idx+1] ;
                                }
							}
						}
						break;
					case 3:
						#pragma omp parallel for schedule( static ) num_threads( num_threads )
						for( int j = 0; j < N; ++j ){
							for( int k = 0; k < N; ++k ){
									for( int l = 0; l < N/2+1; l=l+2 ){
											int idx = (j*N + k)*N + l;
                                        if(l == N/2){
                                        in[idx][0] = fnyquist[j][2*k] ;
                                        in[idx][1] = fnyquist[j][2*k+1] ;
                                        }
                                        else{
                                            in[idx][0] =  f[i][idx] ;
                                            in[idx][1] =  f[i][idx+1] ;
                                        }
                                        }
							}
						}
						break;
				}
		
				fftw_execute(p);
	
				// Set output data
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            switch( dim ){
							case 1:
                int idx = j;
                f[i][idx] = out[idx]/N;
								break;
            	case 2:
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f[i][idx] = out[idx]/(N*N);
                }
								break;
            	case 3:
                for( int k = 0; k < N; ++k ){
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f[i][idx] = out[idx]/(N*N*N);
                    }
                }
								break;
						}
        }
    }

		if( p ) fftw_destroy_plan(p);
		if( in ) fftw_free(in);
		delete[] out;
}
*/

void write_VTK_f( const std::string dir_f, double* f, std::string str, int loop )
{
   // double a = leapfrog->a();
	unsigned int size;
	std::stringstream ss;
	std::ofstream fout;
    
    
	
	if( dim == 1 ){
		ss << "../" << dir_f << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".txt";
    	fout.open( ss.str().c_str() );	
 
    	for( int j = 0; j < N; j++ ){
			int idx = j;
			fout << idx*dx << " " << f[idx] << std::endl;
		}
    }else{
    	ss << "../" << dir_f << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".vti";
    	fout.open( ss.str().c_str() );
    
  	  	fout << "<?xml version=\"1.0\"?>" << std::endl;
  		fout << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl <<std::endl;
    	switch( dim ){
    		case 2:
				size = sizeof(double) * pow(N, 2);//8byte*pow(N,2)
				fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
				fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\">" << std::endl;
    			break;
			case 3:
				size = sizeof(double) * pow(N, 3);
				fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 <<" 0 " << N-1 << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
				fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << N-1 << "\">" << std::endl;
				break;
    	}
    	fout << "<PointData Scalars=\"field\">" << std::endl;
    	fout << "<DataArray type=\"Float64\" Name=\"field\" format=\"appended\" offset=\"0\" />" << std::endl;
    	fout << "</PointData>" << std::endl;
	    fout << "<CellData>" << std::endl;
	    fout << "</CellData>" << std::endl;
	    fout << "</Piece>" << std::endl << std::endl;
	    
	    fout << "</ImageData>" << std::endl << std::endl;
	    fout << "<AppendedData encoding=\"raw\">" << std::endl;
	    fout << "_" ;
	    fout.close();
	    
	    fout.open( ss.str().c_str(), std::ios::binary | std::ios::app);
	    fout.write( (char*) &size, sizeof(unsigned int) );//4byte
       
           fout.write( (char*) f, size );
	    fout.close();	
	    fout.open( ss.str().c_str(), std::ios::app);
	    fout << std::endl << "</AppendedData>" << std::endl;
		fout << "</VTKFile>" ;
	}
	fout.close();
}

void write_VTK_ed( const std::string dir_ed, double* f, std::string str, int loop )
{
  //  double a = leapfrog->a();
   // std::cout << "aaa = " <<  a << std::endl;
    unsigned int size;
    std::stringstream ss;
    std::ofstream fout;
    
    
    
    if( dim == 1 ){
        ss << "../" << dir_ed << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".txt";
        fout.open( ss.str().c_str() );
        
        for( int j = 0; j < N; j++ ){
            int idx = j;
            fout << idx*dx << " " << f[idx] << std::endl;
        }
    }else{
        ss << "../" << dir_ed << "/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".vti";
        fout.open( ss.str().c_str() );
        
        fout << "<?xml version=\"1.0\"?>" << std::endl;
        fout << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl <<std::endl;
        switch( dim ){
            case 2:
                size = sizeof(double) * pow(N, 2);//8byte*pow(N,2)
                fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
                fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\">" << std::endl;
                break;
            case 3:
                size = sizeof(double) * pow(N, 3);
                fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 <<" 0 " << N-1 << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
                fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << N-1 << "\">" << std::endl;
                break;
        }
        fout << "<PointData Scalars=\"energy\">" << std::endl;
        fout << "<DataArray type=\"Float64\" Name=\"energy density\" format=\"appended\" offset=\"0\" />" << std::endl;
        fout << "</PointData>" << std::endl;
        fout << "<CellData>" << std::endl;
        fout << "</CellData>" << std::endl;
        fout << "</Piece>" << std::endl << std::endl;
        
        fout << "</ImageData>" << std::endl << std::endl;
        fout << "<AppendedData encoding=\"raw\">" << std::endl;
        fout << "_" ;
        fout.close();
        
        fout.open( ss.str().c_str(), std::ios::binary | std::ios::app);
        fout.write( (char*) &size, sizeof(unsigned int) );//4byte
        
            fout.write( (char*) f, size );
        
        fout.close();
        fout.open( ss.str().c_str(), std::ios::app);
        fout << std::endl << "</AppendedData>" << std::endl;
        fout << "</VTKFile>" ;
    }
    fout.close();
}


void write_status( const std::string status_file, Field* field, LeapFrog* leapfrog, Energy* energy, double** f, double t )
{
	double a = leapfrog->a();
	std::ofstream ofs;
	
	if( t == t0 )
	{
		ofs.open( "../" + status_file, std::ios::trunc );

		ofs << std::setw(3) << std::right << "  t ";
		if( expansion ) ofs << "  a ";
		for( int i = 0; i < num_fields; ++i ) ofs << "field_ave["  << i << "] ";
		for( int i = 0; i < num_fields; ++i ) ofs << "field_var["  << i << "] ";
        for( int i = 0; i < num_fields; ++i ) ofs << "field_deriv_ave["  << i << "] ";
        for( int i = 0; i < num_fields; ++i ) ofs << "field_deriv_var["  << i << "] ";
		for( int i = 0; i < num_fields; ++i ) ofs << "energy_ave[" << i << "] ";
        for( int i = 0; i < num_fields; ++i ) ofs << "energy_var[" << i << "] ";
        ofs << "total_energy_ave ";
         ofs << "time_deriv_ave ";
         ofs << "gradient_ave ";
        ofs << "potential_ave ";
        ofs << "hubble ";
        ofs << "adotdot ";
        ofs << "energy_max" << std::endl;
	}
	else ofs.open( "../" + status_file, std::ios::app );
	
	ofs << std::setw(3) << std::right << t << " ";
	if( expansion )
	{
		ofs << std::setw(3) << std::right << a << " ";
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->f_average(f[i], i)/a << " "; //Reduced Plank units
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->f_variance(f[i], i)/a << " ";//Reduced Plank units
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->df_average(f[i], i) << " ";//Programming variable
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->df_variance(f[i], i) << " ";//Programming variable
	}
	else
	{
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->f_average(f[i], i) << " ";//Reduced Plank units
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->f_variance(f[i], i) << " ";//Reduced Plank units
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->df_average(f[i], i) << " ";//Programming variable
        for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->df_variance(f[i], i) << " ";//Programming variable
	}
	for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << energy->average(i) << " ";
    for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << energy->variance(i) << " ";
	ofs << std::showpos << std::scientific << std::setprecision(4) << energy->total_average() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->timederiv_average () << " ";
     ofs << std::showpos << std::scientific << std::setprecision(4) << energy->grad_average ()  << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->potential_average () << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << leapfrog->hubble() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << leapfrog->adotdot() << " ";
    ofs << std::showpos << std::scientific << std::setprecision(4) << energy->energy_max() << std::endl;
    
    
}


