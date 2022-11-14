#include "parameters.hpp"

//====================================
// Parameter Sets
//====================================

std::ifstream jsonfile("../include/parameters.json");
json par_set = json::parse(jsonfile);

//Choose Parameter Set defined in the JSON file
int par_set_num = 1;
//Takayama's Master Thesis: 0
//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (a): 1
//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (b): 2
//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG2: 3
//Primordial seeds of SMBHs (peak at 2kMpc-1): 4
//Primordial seeds of SMBHs trying to shift the peak: 5


double GNORMAL = par_set[par_set_num]["GNORMAL"].get<double>();//decay rate of phi set at the very beginning
double GLARGE = par_set[par_set_num]["GLARGE"].get<double>();//decay rate of sigma and psi set at the very beginning
double GLARGE2 = par_set[par_set_num]["GLARGE2"].get<double>(); //Decay rate of sigma and psi set at the end of oscillation period
double CN_par = par_set[par_set_num]["CN_par"].get<double>();//Potential paramater CN
double mu_par = par_set[par_set_num]["mu_par"].get<double>();//Potential paramater mu
double Cv_par = par_set[par_set_num]["Cv_par"].get<double>();//Potential parameter Cv
double M_par = par_set[par_set_num]["M_par"].get<double>();//Potential paramater M
//double m_par = par_set[par_set_num]["m_par"].get<double>();//Potential paramater m
//double n_par = par_set[par_set_num]["n_par"].get<double>();//Potential paramater n
//double g_par = par_set[par_set_num]["g_par"].get<double>();//Potential parameter g
//double BEGIN_EFOLD = par_set[par_set_num]["BEGIN_EFOLD"].get<double>();//initial ln(a)
//double UNPERT_EFOLD = par_set[par_set_num]["UNPERT_EFOLD"].get<double>();//ln(a) at which sigma and psi are fixed to the minimum after oscillation period
//double NEWINF_END_EFOLD = par_set[par_set_num]["NEWINF_END_EFOLD"].get<double>();//ln(a) at the beginning of oscillation of phi.
//double END_EFOLD = par_set[par_set_num]["END_EFOLD"].get<double>();//ln(a) at the very end
//double dla = par_set[par_set_num]["dla"].get<double>();//stepsize for fixed step RK-method
//double itvl = par_set[par_set_num]["itvl"].get<double>(); //interval for output in fixed RK-method
//double sigma_init = par_set[par_set_num]["sigma_init"].get<double>();//initial value of sigma
std::string par_set_name = par_set[par_set_num]["Name"].get<std::string>();//Parameter set name
//====================================
//The following only holds for m=2
//====================================
//double msigma = sqrt(8*pow(mu_par,3)/M_par);//effective mass of sigma





//====================================================================
//                      Non-lattice Range
//====================================================================
//-----------
//File names
//-----------
std::string exist_filename_zero = par_set_name + "_unpWMAP5.txt";
std::string new_filename_zero = par_set_name + "_unpWMAP5.txt";

std::string exist_dirname_k = par_set_name + "_kAnalyze_nonlattice"; //remove this existing directory for k-analyze txt files
std::string new_dirname_k = par_set_name + "_kAnalyze_nonlattice"; //create a new directory for k-analyze txt files
std::string filename_k = par_set_name + "_kAnalyze"; // Head of the file name for k-analyze txt files

bool k_switch = false; //If this is true, the existing dir and files for kAnalyze are deleted and new dir and files are created. If false, the dir and files remain as it is.

std::string exist_filename_sp  = par_set_name + "_spectrum.txt";// remove this existing spectrum file
std::string new_filename_sp  = par_set_name + "_spectrum.txt"; // create this new spectrum file

std::string exist_filename_spbfosc  = par_set_name + "_spectrum_bfosc.txt";// remove this existing spectrum file
std::string new_filename_spbfosc  = par_set_name + "_spectrum_bfosc.txt"; // create this new spectrum file


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

bool lattice_kmodes_switch = false; //If this is true, it will use the modes calculated in lattice simulation for non lattice zeromode w/ perturb calculation

//Outputs
bool kanalyze_switch = true;// true:Calculate k-analyze, false:Don't calculate k-analyze
bool spectrum_switch = true;// true:Calculate final spectrum, false:Don't calculate final spectrum
bool spectrum_bfosc_switch = true; // true:Calculate spectrum before oscillation starts, false:Don't calculate

//Array elements
const int  N_pert=55;

//If lattice_kmodes_switch is false, we use the following range for calculation

//Takayama Master Thesis
//DP kfrom_Mpc = 3;//[Mpc^-1] Calculate from this k (Must be greater than 0)
//DP kto_Mpc = 1000;//[Mpc^-1] Calculate to this k
//int kinterval_knum = 10;// [knum units] Calculate with this interval of knum

//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (a)
DP kfrom_Mpc = 1.0e-4;//[Mpc^-1] Calculate from this k (Must be greater than 0)
DP kto_Mpc = 1.0e+4;//[Mpc^-1] Calculate to this k
int kinterval_knum = 1;// [knum units] Calculate with this interval of knum







//====================================================================
//                             Lattice Range
//====================================================================

//-----------
//File names
//-----------


std::string exist_dirname_ed = par_set_name + "_energy"; //remove this existing directory for energy vti files
std::string new_dirname_ed = par_set_name + "_energy"; //create a new directory for energy vti files

std::string exist_dirname_f = par_set_name + "_field"; //remove this existing directory for field vti files
std::string new_dirname_f = par_set_name + "_field"; //create a new directory for field vti files


std::string exist_filename_status = par_set_name + "_status.txt"; // remove this existing status file
std::string new_filename_status = par_set_name + "_status.txt";// create this new status file

//lattice simulation version of kAnalyze
std::string exist_dirname_k_lattice = par_set_name + "_kAnalyze_lattice"; //remove this existing directory for k-analyze txt files
std::string new_dirname_k_lattice = par_set_name + "_kAnalyze_lattice"; //create a new directory for k-analyze txt files
std::string filename_k_lattice = par_set_name + "_kAnalyze"; // Head of the file name for k-analyze txt files

bool k_lattice_switch = true; //If this is true, the existing dir and files for kAnalyze using lattice simulation are deleted and new dir and files are created. If false, the dir and files remain as it is.

bool k_lattice_startfromlattice_switch = false; //If this is true, data output of kAnalyze_lattice starts from the time when lattice simulation starts. If false, then data output starts from the beginning (the beginning of hybrid inflation)

//-------------------------------------------------
//Variables for calculating lattice range
//-------------------------------------------------

bool latticerange_switch = false; // true:Set lattice range and calculate, false:Don't set lattice range and calculate
bool initialize_perturb_switch = true; // true:Initialize fluctuation, false:Don't initialize fluctuation (only calculate zeromode for lattice)

int fluc_calc_switch  = 1;//Choose type of fluctuation initialization for scalar fields (for gravitational fluctuation, it is set to 2 regardless) 0:LatticeEasy case (when amplitudes of fluctuations are not predetermined) 1:when we use the amplitudes of predetermined fluctuations 2:gravitational perturbation

//Takayama Master Thesis
//double kfrom_Mpc_lattice = 3;//[Mpc^-1] Calculate from this k for lattice range
//double kto_Mpc_lattice = 400;//[Mpc^-1] Calculate to this k for lattice range

//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (a)
double kfrom_Mpc_lattice = 1;//[Mpc^-1] Calculate from this k for lattice range
double kto_Mpc_lattice = 3000;//[Mpc^-1] Calculate to this k for lattice range

int N = 64;//512; //Should be 2^n

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
int total_step  = 8.2e+4;//1.75e+4;//8.75e+3;
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
