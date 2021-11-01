#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <string>
#include "calculation.hpp"

#define Logout(...) 	do { printf(__VA_ARGS__); fflush(stdout); } while(0)


void dir_manage(const std::string exist_dir, const std::string new_dir);

void file_manage(const std::string exist_file);

void set_mode(double p2, double m2, double *field, double *deriv, int real);

void DFT_c2rD1( double* f);

void DFT_c2rD2( double* f,double* fnyquist );

void DFT_c2rD3( double* f,double** fnyquist );

void write_VTK_f  ( const std::string dir_f, double* f, std::string str, int loop );

void write_VTK_ed  (const std::string dir_ed,  double* f, std::string str, int loop);

void write_status ( const std::string status_file, Field* field, LeapFrog* leapfrog, Energy* energy, double** f, double t );


#endif
