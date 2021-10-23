#include "parameter.hpp"




int N = 1024;
int L = 20;
int rnd = 1;
const int num_fields  = 1;
int num_threads = 8;

const double Hinitial = 5.98147171787852*pow(10,-6)/(sqrt(8*M_PI));
const double m = sqrt(5.67402574831172*pow(10,-10))/(sqrt(8*M_PI));//rescale_B, sqrt(V''(phi))
const double initfield[] = {2.2105/(sqrt(8*M_PI))};//{sqrt(2)*1.49652/(sqrt(8*M_PI))};
const double initderivs[] = {(1*Hinitial*initfield[0])/m};//{(1*Hinitial*initfield[0])/m}; //no expansion -> 0, expansion -> rescale_r*Hinitial*f_pr/rescale_B -> (1*Hinitial*initfield[0])/m

int output_step = 1.5e+1;
int total_step  = 1.2e+4; 
int max_loop    = total_step/output_step; // This many times vti files will be created
int st_output_step = 1;
int st_max_loop = output_step/st_output_step; // This many times data will be added to status.txt between the output of vti files 

double t0 = 0;
double dt = 5.e-3;
double dx = 1.* L/N;

int expansion = 3; // 0: no expansion, 1: self-consistent, 2: radiation dominant, 3: matter dominant
int precision = 2;

bool restart = false;  // cannot use yet

//potential parameters
const double aa = 2*M_PI;
const double AA = 10;
const double W0 = -pow(10,-5);//-pow(10,-5);
const double D = 3.169403285436*pow(10,-11);//4.6824231*pow(10,-12);

/*f/M_p=f_pr/a
 dx_pr=mdx
 dt_pr=mdt/a
 
 A = 1/M_p (Therefore 1 in Planck units.), B = m, r = 1, s = -1,
 
 */

