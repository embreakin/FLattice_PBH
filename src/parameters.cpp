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

std::string exist_dirname_k = "kAnalyze_nonlattice"; //remove this existing directory for k-analyze txt files
std::string new_dirname_k = "kAnalyze_nonlattice"; //create a new directory for k-analyze txt files
std::string filename_k = "kAnalyze"; // Head of the file name for k-analyze txt files

bool k_switch = false; //If this is true, the existing dir and files for kAnalyze are deleted and new dir and files are created. If false, the dir and files remain as it is.

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

DP sigma_c = pow((m_par/(2*m_par-1))*((3*m_par-1)/(m_par-1)),(m_par-1)/(4*m_par-2))*2/pow(2*m_par*(2*m_par-1),m_par/(4*m_par-2))*pow(mu_par*pow(M_par,m_par-1),1/(2*m_par -1));

DP FIXPSI = 2*pow(mu_par*pow(M_par,m_par-1),1/m_par);
DP FIXPHI = -sqrt(2)*pow(pw2(Cv_par)/g_par,1/n_par);

DP CNT = (-1)*Vbare(0.0,FIXPSI,FIXPHI);

//Vbare(0,FIXPSI,FIXPHI);       //constant term in the potential (set V=0 at the minimum)
//calculate log(Gravitational Potential) and log(Zeta) for each knum specified.

DP Gamma1,Gamma2,Gamma3;
DP OSCSTART;

std::vector<int> knum_zero = {200, 400, 634, 700};//knum for calculating zeromode

int k_target = knum_zero[3]; //target wave mode actually used for zeromode calculation


//-------------------------------------------------
//Variables for zeromode w/ perturbation calculation
//-------------------------------------------------
bool perturbation_switch = true;//This needs to be true for perturbation calculation including lattice simulation calculation

bool lattice_kmodes_switch = true; //If this is true, it will use the modes calculated in lattice simulation for non lattice zeromode w/ perturb calculation

//Outputs
bool kanalyze_switch = true;// true:Calculate k-analyze, false:Don't calculate k-analyze
bool spectrum_switch = true;// true:Calculate final spectrum, false:Don't calculate final spectrum
bool spectrum_bfosc_switch = true; // true:Calculate spectrum before oscillation starts, false:Don't calculate

//Array elements
const int  N_pert=55;

//If lattice_kmodes_switch is false, we use the following range for calculation
DP kfrom_Mpc = 3;//[Mpc^-1] Calculate from this k (Must be greater than 0)
DP kto_Mpc = 1000;//[Mpc^-1] Calculate to this k
int kinterval_knum = 10;// [knum units] Calculate with this interval of knum





//====================================================================
//                             Lattice Range
//====================================================================

//-----------
//File names
//-----------
std::string exist_dirname_ed = "dataSHNI_energy"; //remove this existing directory for energy vti files
std::string new_dirname_ed = "dataSHNI_energy"; //create a new directory for energy vti files

std::string exist_dirname_f = "data"; //remove this existing directory for field vti files
std::string new_dirname_f = "data"; //create a new directory for field vti files


std::string exist_filename_status = "statusSHNI.txt"; // remove this existing status file
std::string new_filename_status = "statusSHNI.txt";// create this new status file

//lattice simulation version of kAnalyze
std::string exist_dirname_k_lattice = "kAnalyze_lattice"; //remove this existing directory for k-analyze txt files
std::string new_dirname_k_lattice = "kAnalyze_lattice"; //create a new directory for k-analyze txt files
std::string filename_k_lattice = "kAnalyze"; // Head of the file name for k-analyze txt files

bool k_lattice_switch = true; //If this is true, the existing dir and files for kAnalyze using lattice simulation are deleted and new dir and files are created. If false, the dir and files remain as it is.

bool k_lattice_startfromlattice_switch = false; //If this is true, data output of kAnalyze_lattice starts from the time when lattice simulation starts. If false, then data output starts from the beginning (the beginning of hybrid inflation)

//-------------------------------------------------
//Variables for calculating lattice range
//-------------------------------------------------

bool latticerange_switch = true; // true:Set lattice range and calculate, false:Don't set lattice range and calculate
bool initialize_perturb_switch = true; // true:Initialize fluctuation, false:Don't initialize fluctuation (only calculate zeromode for lattice)

int fluc_calc_switch  = 1;//Choose type of fluctuation initialization for scalar fields (for gravitational fluctuation, it is set to 1 regardless) 0:LatticeEasy case (when amplitudes of fluctuations are not predetermined) 1:when we use the amplitudes of predetermined fluctuations


double kfrom_Mpc_lattice = 3;//[Mpc^-1] Calculate from this k for lattice range
double kto_Mpc_lattice = 400;//[Mpc^-1] Calculate to this k for lattice range

int N = 32;//512; //Should be 2^n

double kfrom_MPl_lattice = UC::kMpc_to_kMPl(kfrom_Mpc_lattice); //convert to MPl units
double kto_MPl_lattice = UC::kMpc_to_kMPl(kto_Mpc_lattice); //convert to MPl units


double k_lattice_grid_min_MPl = kto_MPl_lattice/(N/2);
double k_lattice_grid_max_MPl = kto_MPl_lattice;

int outrange_num; //the number of wave modes that are not in the lattice range
int latticerange_num;//the number of wave modes that are in the lattice range


int rnd = 1;
int num_fields  = 4; //0:sigma 1:psi 2:phi 3:metric perturbation

//number of threads
int num_threads = omp_get_num_procs()/2;

double Hinitial_pr;
//const double m = sqrt(5.67402574831172*pow(10,-10));//rescale_B, sqrt(V''(phi))

//const double initfield[] = {2.2105};//{sqrt(2)*1.49652/(sqrt(8*M_PI))};
//const double initderivs[] = {(1*Hinitial*initfield[0])/rescale_B};//{(1*Hinitial*initfield[0])/m}; //no expansion -> 0, expansion -> rescale_r*Hinitial*f_pr/rescale_B -> (1*Hinitial*initfield[0])/m

int output_step = 2.0e+1;
int total_step  = 4.1e+3;//1.75e+4;//8.75e+3;
int max_loop    = total_step/output_step; // This many times vti files will be created
int st_output_step = 10;
int st_max_loop = output_step/st_output_step; // This many times data will be added to status.txt between the output of vti files

double t0 = 0;
double dt = 1.e-4;//5.0e-4;//1.e-3; //dt_pr


const int expansion = 1; // 0: no expansion, 1: self-consistent, 2: radiation dominant, 3: matter dominant
const int precision = 2;

const int screen_latticeloop_number = 100; //This many times loop will be displayed on the terminal for lattice simulation. If you want to show all loops, then set this number to max_loop.

/*f/M_p=f_pr/a
 dx_pr=mdx
 dt_pr=mdt/a
 
 A = 1/M_p (Therefore 1 in Planck units.), B = m, r = 1, s = -1,
 
 */
