//Doxygen
/**
* @file   lattice_field.hpp
* @brief    Lattice field header file
* @author   Francis Otani
* @date
* @details
*/

#ifndef _LATTICEFIELD_H_
#define _LATTICEFIELD_H_

#include <cmath>
#include "parameters.hpp"

class LeapFrog;

class Field
{
	private:
		double* _average;
		double* _variance;
    
        double* f_MPl;
        int fld;
    
        int J,K;
        int m_start;
    
    std::vector<std::vector<double>> PS;
    std::vector<std::vector<double>> f_fluc;
    std::vector<std::vector<double>> f_fluc_k;
    std::vector<std::vector<double>> df_fluc;
    std::vector<std::vector<double>> df_fluc_k;
    std::vector<std::vector<double>> f_fluc_k_nyquist_2d;
    std::vector<std::vector<double>>
        df_fluc_k_nyquist_2d;
    std::vector<std::vector<std::vector<double>>>
    f_fluc_k_nyquist_3d;
    std::vector<std::vector<std::vector<double>>>
    df_fluc_k_nyquist_3d;
    
	public:
    
   
    
        Field (): _average(new double [num_fields]()), _variance(new double [num_fields]()), f_MPl(new double [num_fields]())  {
            
            switch (dim){
            
            case 1:
                {
            //Use 2D vector for power spectrum (fields, wave modes). Initialize by 0
             PS = std::vector<std::vector<double>>(num_fields+1, std::vector<double>(N/2+1, 0));
            
            //This 2D vector receives only the fluctuation of the input field
            f_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N, 0));
            //This 2D vector receives the output of r2c Fourier transform
             f_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N, 0));
            
            //This 2D vector receives only the fluctuation of the input field derivative
             df_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N, 0));
            //This 2D vector receives the output of r2c Fourier transform
             df_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N, 0));
                }
                    break;
        case 2:
                {
            //Use 2D vector for power spectrum (fields, wave modes). Initialize by 0
             PS = std::vector<std::vector<double>>(num_fields+1, std::vector<double>(N/2+1, 0));
            
            //This 2D vector receives only the fluctuation of the input field
             f_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N, 0));
            //This 2D vector receives the output of r2c Fourier transform
             f_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N, 0));
            //This 2D vector receives the output of r2c Fourier transform (Nyquist frequency corresponding to k (z-axis) = N/2)
             f_fluc_k_nyquist_2d =  std::vector<std::vector<double>>(num_fields, std::vector<double>(2*N, 0));
            
            //This 2D vector receives only the fluctuation of the input field
             df_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N, 0));
            //This 2D vector receives the output of r2c Fourier transform
            df_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N, 0));
            //This 2D vector receives the output of r2c Fourier transform (Nyquist frequency corresponding to k (z-axis) = N/2)
             df_fluc_k_nyquist_2d =  std::vector<std::vector<double>>(num_fields, std::vector<double>(2*N, 0));
                }
                    break;
        case 3:
                {
            //Use 2D vector for power spectrum (fields, wave modes). Initialize by 0
             PS = std::vector<std::vector<double>>(num_fields+1, std::vector<double>(N/2+1, 0));
            
            //This 2D vector receives only the fluctuation of the input field
             f_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N*N, 0));
            //This 2D vector receives the output of r2c Fourier transform
             f_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N*N, 0));
            //This 3D vector receives the output of r2c Fourier transform (Nyquist frequency corresponding to l (z-axis) = N/2)
             f_fluc_k_nyquist_3d =  std::vector<std::vector<std::vector<double>>>(num_fields, std::vector<std::vector<double>>(N, std::vector<double>(2*N, 0)));
            
            //This 2D vector receives only the fluctuation of the input field
             df_fluc =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N*N, 0));
            //This 2D vector receives the output of r2c Fourier transform
             df_fluc_k =  std::vector<std::vector<double>>(num_fields, std::vector<double>(N*N*N, 0));
            //This 3D vector receives the output of r2c Fourier transform (Nyquist frequency corresponding to l (z-axis) = N/2)
            df_fluc_k_nyquist_3d =  std::vector<std::vector<std::vector<double>>>(num_fields, std::vector<std::vector<double>>(N, std::vector<double>(2*N, 0)));
                }
                    break;
        } //switch
            
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
