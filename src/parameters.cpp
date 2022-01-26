#include "parameters.hpp"

//Global Variables

//====================================================================
//                      Non-lattice Range
//====================================================================

//-----------
//File names
//-----------
std::string exist_filename_zero = "unpWMAP5.txt";
std::string new_filename_zero = "unpWMAP5.txt";

std::string exist_dirname_k = "dataDI"; //remove this existing directory for k-analyze txt files
std::string new_dirname_k = "dataDI"; //create a new directory for k-analyze txt files
std::string filename_k = "kAnalyze"; // Head of the file name for k-analyze txt files

std::string exist_filename_sp  = "spectrum.txt";// remove this existing spectrum file
std::string new_filename_sp  = "spectrum.txt"; // create this new spectrum file

std::string exist_filename_spbfosc  = "spectrum_bfosc.txt";// remove this existing spectrum file
std::string new_filename_spbfosc  = "spectrum_bfosc.txt"; // create this new spectrum file


//----------------------------------
//Variables for zeromode calculation
//----------------------------------
bool zeromode_switch = true;

//Array elements
const int N_zero=7;

DP CNT = (-1)*Vbare(0,FIXPSI,FIXPHI);
//Vbare(0,FIXPSI,FIXPHI);       //constant term in the potential (set V=0 at the minimum)
//calculate log(Gravitational Potential) and log(Zeta) for each knum specified.

DP Gamma1,Gamma2,Gamma3;

std::vector<int> knum_zero = {200, 400, 634, 700};//knum for calculating zeromode


//-------------------------------------------------
//Variables for zeromode w/ perturbation calculation
//-------------------------------------------------
bool perturbation_switch = true;

//Outputs
bool kanalyze_switch = true;// true:Calculate k-analyze, false:Don't calculate k-analyze
bool spectrum_switch = true;// true:Calculate final spectrum, false:Don't calculate final spectrum
bool spectrum_bfosc_switch = true; // true:Calculate spectrum before oscillation starts, false:Don't calculate

//Array elements
const int  N_pert=55;

DP kfrom_Mpc = 0.01;//[Mpc^-1] Calculate from this k
DP kto_Mpc = 1000;//[Mpc^-1] Calculate to this k
int kinterval_knum = 1;// [knum units] Calculate with this interval of knum



//====================================================================
//                             Lattice Range
//====================================================================

//-----------
//File names
//-----------
std::string exist_dirname_ed = "dataKKLT3DN64L20MDwparallel"; //remove this existing directory for energy vti files
std::string new_dirname_ed = "dataKKLT3DN64L20MDwparallel"; //create a new directory for energy vti files

std::string exist_dirname_f = "data"; //remove this existing directory for field vti files
std::string new_dirname_f = "data"; //create a new directory for field vti files


std::string exist_filename_status = "statusKKLT3DN64L20MDwparallel.txt"; // remove this existing status file
std::string new_filename_status = "statusKKLT3DN64L20MDwparallel.txt";// create this new status file

//-------------------------------------------------
//Variables for calculating lattice range
//-------------------------------------------------

bool latticerange_switch = false; // true:Set lattice range and calculate, false:Don't set lattice range and calculate

double kfrom_Mpc_lattice = 600;//[Mpc] Calculate from this k for lattice range
double kto_Mpc_lattice = 900;//[Mpc] Calculate to this k for lattice range

int N = 64; //Should be 2^n

double kfrom_MPl_lattice = UC::kMpc_to_kMPl(kfrom_Mpc_lattice); //convert to MPl units
double kto_MPl_lattice = UC::kMpc_to_kMPl(kto_Mpc_lattice); //convert to MPl units

double rescale_A = 1;
double rescale_B = sqrt(V_11(0,FIXPSI,0));
double L = N*M_PI*rescale_B/(kto_MPl_lattice);

double k_lattice_grid_min_pr = 2*M_PI/L;
double k_lattice_grid_min_MPl = k_lattice_grid_min_pr*rescale_B;

double k_lattice_grid_max_pr = N*M_PI/L;
double k_lattice_grid_max_MPl = k_lattice_grid_max_pr*rescale_B;

int rnd = 1;
int num_fields  = 3;
int num_threads = omp_get_num_procs()/2;

const double Hinitial = 5.98147171787852*pow(10,-6);
//const double m = sqrt(5.67402574831172*pow(10,-10));//rescale_B, sqrt(V''(phi))
const double ENGRESCALE = pw2(rescale_B); // Used to rescale the energies from reduced Planck units to program variables and vice versa
const double initfield[] = {2.2105};//{sqrt(2)*1.49652/(sqrt(8*M_PI))};
const double initderivs[] = {(1*Hinitial*initfield[0])/rescale_B};//{(1*Hinitial*initfield[0])/m}; //no expansion -> 0, expansion -> rescale_r*Hinitial*f_pr/rescale_B -> (1*Hinitial*initfield[0])/m

int output_step = 1.5e+1;
int total_step  = 1.2e+4;
int max_loop    = total_step/output_step; // This many times vti files will be created
int st_output_step = 1;
int st_max_loop = output_step/st_output_step; // This many times data will be added to status.txt between the output of vti files

double t0 = 0;
double dt = 5.e-3;
double dx = 1.* L/N;

const int expansion = 3; // 0: no expansion, 1: self-consistent, 2: radiation dominant, 3: matter dominant
const int precision = 2;

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
