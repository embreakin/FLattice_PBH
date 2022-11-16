#ifndef _LATTICEFIELD_H_
#define _LATTICEFIELD_H_

#include <cmath>
#include "parameters.hpp"
#include "lattice_evol.hpp"


class Field
{
	private:
		double* _average;
		double* _variance;
    
        double* f_MPl;
        int fld;
    
        int J,K;
        int m_start;
    
#if  dim==1
    
    
    //Use 2D vector for power spectrum (fields, wave modes). Initialize by 0
    std::vector<std::vector<double>> PS = std::vector<std::vector<double>>(num_fields, std::vector<double>(N/2+1, 0));
    
    //This 2D vector receives only the fluctuation of the input field
    std::vector<std::vector<double>> f_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N, 0));
    //This 2D vector receives the output of r2c Fourier transform
    std::vector<std::vector<double>> f_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N, 0));
    
    //This 2D vector receives only the fluctuation of the input field derivative
    std::vector<std::vector<double>> df_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N, 0));
    //This 2D vector receives the output of r2c Fourier transform
    std::vector<std::vector<double>> df_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N, 0));
    
    
#elif  dim==2
    
    //Use 2D vector for power spectrum (fields, wave modes). Initialize by 0
    std::vector<std::vector<double>> PS = std::vector<std::vector<double>>(num_fields, std::vector<double>(N/2+1, 0));
    
    //This 2D vector receives only the fluctuation of the input field
    std::vector<std::vector<double>> f_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N, 0));
    //This 2D vector receives the output of r2c Fourier transform
    std::vector<std::vector<double>> f_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N, 0));
    //This 2D vector receives the output of r2c Fourier transform (Nyquist frequency corresponding to k (z-axis) = N/2)
    std::vector<std::vector<double>> f_fluc_k_nyquist =  std::vector<std::vector<double>>(num_fields, std::vector<double>(2*N, 0));
    
    //This 2D vector receives only the fluctuation of the input field
    std::vector<std::vector<double>> df_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N, 0));
    //This 2D vector receives the output of r2c Fourier transform
    std::vector<std::vector<double>> df_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N, 0));
    //This 2D vector receives the output of r2c Fourier transform (Nyquist frequency corresponding to k (z-axis) = N/2)
    std::vector<std::vector<double>> df_fluc_k_nyquist =  std::vector<std::vector<double>>(num_fields, std::vector<double>(2*N, 0));
    
#elif  dim==3
    
    //Use 2D vector for power spectrum (fields, wave modes). Initialize by 0
    std::vector<std::vector<double>> PS = std::vector<std::vector<double>>(num_fields, std::vector<double>(N/2+1, 0));
    
    //This 2D vector receives only the fluctuation of the input field
    std::vector<std::vector<double>> f_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N*N, 0));
    //This 2D vector receives the output of r2c Fourier transform
    std::vector<std::vector<double>> f_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N*N, 0));
    //This 3D vector receives the output of r2c Fourier transform (Nyquist frequency corresponding to l (z-axis) = N/2)
    std::vector<std::vector<std::vector<double>>> f_fluc_k_nyquist =  std::vector<std::vector<std::vector<double>>>(num_fields, std::vector<std::vector<double>>(N, std::vector<double>(2*N, 0)));
    
    //This 2D vector receives only the fluctuation of the input field
    std::vector<std::vector<double>> df_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N*N, 0));
    //This 2D vector receives the output of r2c Fourier transform
    std::vector<std::vector<double>> df_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N*N, 0));
    //This 3D vector receives the output of r2c Fourier transform (Nyquist frequency corresponding to l (z-axis) = N/2)
    std::vector<std::vector<std::vector<double>>> df_fluc_k_nyquist =  std::vector<std::vector<std::vector<double>>>(num_fields, std::vector<std::vector<double>>(N, std::vector<double>(2*N, 0)));
    
    
#endif
    
    
	public:
    
   
    
        Field (): _average(new double [num_fields]()), _variance(new double [num_fields]()), f_MPl(new double [num_fields]())  {
            
           
//            PS = new double** [num_fields];
//
//            switch( dim )
//            {
//                case 1:
//                    PS[0] = new double* [num_fields*(N/2+1)];
//
//                    for( int i = 0; i < num_fields; ++i )
//                    {
//                        PS[i] = PS[0] + i*(N/2+1);
//
//                        //Initialize all elements to 0 first. If we don't do this, a very small (or in some cases very large) value will be assigned instead and it can cause errors. Also, for some reason std::fill() doesn't work for some cases, so I've avoided using it. //
//                        for( int j = 0; j < N; ++j ){
//                            int idx = j;
//
//                            PS[i][idx] = 0;
//                        }
//
//                    }
//
//                    break;
//                case 2:
//                    PS[0] = new double [num_fields*N*N];
//
//                    for( int i = 0; i < num_fields; ++i )
//                    {
//                        PS[i] = PS[0] + i*N*N;
//
//                        for( int j = 0; j < N; ++j ){
//                            for( int k = 0; k < N; ++k ){
//                                int idx = j*N + k;
//                                PS[i][idx] = 0;
//                            }
//                        }
//
//                    }
//                    break;
//                case 3:
//                    PS[0] = new double [num_fields*N*N*N];
//
//                    for( int i = 0; i < num_fields; ++i )
//                    {
//                        PS[i] = PS[0] + i*N*N*N;
//
//                        for( int j = 0; j < N; ++j ){
//                            for( int k = 0; k < N; ++k ){
//                                for( int l = 0; l < N; ++l ){
//                                    int idx = (j*N + k)*N + l;
//                                    PS[i][idx] = 0;
//                                }
//                            }
//                        }
//                    }
//                    break;
//                default:
//                    std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
//                    exit(1);
//            }
//
            
        }
    
        ~Field () {
            delete [] _average;
            delete [] _variance;
            delete [] f_MPl;
        }
    
    //void zeromode_initialize();
		double laplacian        ( double* f, int j, int k = 0, int l = 0 ); //You can omit k and l if they are zero
        double gradient_energy_eachpoint( double** f ,int i, int idx );
		double gradient_energy  ( double* f ) ;
		double potential_energy ( double** f, double a );
        double average  ( double* f, int i );
        double variance ( double* f, int i );
        double V_lattice  ( double** f, int idx, double a = 1 );
        double dV_lattice ( double** f, int i, int idx, double a = 1);
        double ddV_lattice ( double** f, int i, int idx, double a = 1 );
        void effective_mass(double mass_sq[], double *field_values);
        double power_spectrum( double** f, int i, int j);
    void finalize(double** f, double** df, LeapFrog* leapfrog, double radiation_pr, double** lattice_var );

};




#endif
