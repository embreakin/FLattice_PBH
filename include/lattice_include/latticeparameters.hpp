#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <string>
#include <omp.h>

#define dim 3

extern std::string exist_dirname_ed, new_dirname_ed, exist_dirname_f, new_dirname_f, exist_filename_status, new_filename_status;

inline double pw2(double x) { return (x*x);}

extern int N;
extern int L;
extern int rnd;
extern int num_fields;
extern int num_threads;
extern const double initfield[];
extern const double initderivs[];
extern const double m;//rescale_B
extern const double ENGRESCALE;
extern const double Hinitial;

extern int output_step;
extern int total_step;
extern int max_loop;
extern int st_output_step;
extern int st_max_loop;

extern double t0;
extern double dt;
extern double dx;

extern const int expansion;
extern const int precision;
extern bool restart;

//potential parameters
extern const double aa;// = 2*M_PI;
extern const double AA;// = 10;
extern const double W0;// = -pow(10,-5);
extern const double D ;//= 4.6824231*pow(10,-12);

#endif
