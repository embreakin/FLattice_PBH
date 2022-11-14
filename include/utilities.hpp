#ifndef _UTILITIES_H_
#define _UTILITIES_H_


#include "lattice_calc.hpp"
#include "nr.h"
#include "uc.hpp"
#include <string>
#include <chrono>

#define Logout(...) 	do { printf(__VA_ARGS__); fflush(stdout); } while(0)

//--------------------------------
// Directory & File Management
//--------------------------------

void dir_manage(const std::string exist_dir, const std::string new_dir);

void file_manage(const std::string exist_file);

//--------------------------------
// Elapsed Time Calculation
//--------------------------------

void time_calc(std::chrono::system_clock::time_point time_start, std::chrono::system_clock::time_point time_end, std::string time_name);

//--------------------------------
// Double Inflation Subroutines
//--------------------------------
void zeromode_output(const std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount);

void kanalyze_output(const std::string dir, std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount, int knum, DP k_comoving);

void spectrum_bfosc_output(const std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount, int knum, const DP k_comoving);

void spectrum_output(const std::string file, Vec_I_DP &xx, Mat_I_DP &yp, int timecount, int knum, const DP k_comoving);

//--------------------------------
// Lattice Simulation Subroutines
//--------------------------------

void DFT_c2rD1( double* f);
void DFT_r2cD1( std::vector<double>& f, std::vector<double>& f_k);

void DFT_c2rD2( double* f,double* fnyquist );
void DFT_r2cD2(std::vector<double>& f, std::vector<double>& f_k, std::vector<double>& f_k_nyquist );

void DFT_c2rD3( double* f,double** fnyquist );
void DFT_r2cD3( std::vector<double>& f, std::vector<double>& f_k, std::vector<std::vector<double>>& f_k_nyquist );


void write_VTK_f  ( const std::string dir_f, double* f, std::string str, int loop );

void write_VTK_ed  (const std::string dir_ed,  double* f, std::string str, int loop);

void write_status ( const std::string status_file, Field* field, LeapFrog* leapfrog, Energy* energy, double** f, double** df, double t );

void kanalyze_output_lattice(const std::string dir, std::string file, Field* field, LeapFrog* leapfrog, double** f);

#endif
