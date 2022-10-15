#ifndef _LATTICEINITIALIZE_H_
#define _LATTICEINITIALIZE_H_

#include "lattice_field.hpp"


double Fk_log_int_calc(int k_int , double** lattice_var, int num_field);

double dFk_log_int_calc(int k_int , double** lattice_var, int num_field);

void fdf_calc(double distance, double** lattice_var, double *field, double *deriv, int num_field);

void set_mode(double p2, double m2, double *field, double *deriv, int i, int real);

void initialize_perturb(double** f, double** df, double** lattice_var, double mass_sq[]);

void initialize( double**& f, double**& df, Field* field, double &radiation_pr, double**& lattice_var);

void finalize( double** f, double** df );

#endif
