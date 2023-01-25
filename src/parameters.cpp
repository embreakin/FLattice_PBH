//Doxygen
/**
* @file    parameters.cpp
* @brief      Parameters source file
* @author   Francis Otani
* @date
* @details
*/


#include "parameters.hpp"



//====================================
//         General Parameters
//====================================
std::ifstream jsonfile("../include/parameters.json");
json par_set = json::parse(jsonfile);

bool exist_par_set_rmall_switch = false;//If this is true, the entire directory of the previously calculated parameter set will be deleted.

//Name of the parameter set

//This is the name of the parameter set that you are about to simulate
std::string par_set_name = par_set[par_set_num]["Name"].get<std::string>();

//This is the name of the parameter set that you have simulated previously and will be removed if exist_par_set_rmall_switch is set to true
std::string par_set_name_rm = par_set[par_set_num_rm]["Name"].get<std::string>();

//This is the name of the condition you are about to simulate
std::string condition_name;

void condition_name_set()
{
    
    condition_name = "zs_" + boost::lexical_cast<std::string>(par_set[par_set_num]["zeromode_switch"]);
    
    
    if(par_set[par_set_num]["perturbation_switch"])
    {
        condition_name += "-ps_" + boost::lexical_cast<std::string>(par_set[par_set_num]["perturbation_switch"])
        + "-ls_" + boost::lexical_cast<std::string>(par_set[par_set_num]["latticerange_switch"]);
        
        if(par_set[par_set_num]["latticerange_switch"])
        {
            
            condition_name +=
            "-kfrom_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kfrom_Mpc"])
           + "-kto_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kto_Mpc"])
            + "-iknum_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kinterval_knum"])
            +
            "-kfroml_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kfrom_Mpc_lattice"])
           + "-ktol_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kto_Mpc_lattice"])
           + "-N_" + boost::lexical_cast<std::string>(par_set[par_set_num]["N"])
           + "-dim_" + boost::lexical_cast<std::string>(par_set[par_set_num]["dim"]);
        }else
        {
            condition_name += "-lks_" +  boost::lexical_cast<std::string>(par_set[par_set_num]["lattice_kmodes_switch"]);
            
        if(par_set[par_set_num]["lattice_kmodes_switch"])
           {
            condition_name += "-kfroml_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kfrom_Mpc_lattice"])
            + "-ktol_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kto_Mpc_lattice"])
            + "-N_" + boost::lexical_cast<std::string>(par_set[par_set_num]["N"]);
            }
        else
           {
            condition_name +=
            "-kfrom_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kfrom_Mpc"])
           + "-kto_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kto_Mpc"])
            + "-iknum_" + boost::lexical_cast<std::string>(par_set[par_set_num]["kinterval_knum"]);
        
            }
        }
    }

}




//=========================================
//The following parameter only holds for m=2
//=========================================
double msigma = sqrt(8*pow(mu_par,3.0)/M_par);//effective mass of sigma





//====================================================================
//                      Non-lattice Range
//====================================================================
//-----------
//File names
//-----------
std::string new_filename_zero = "zeromode.txt";

std::string new_dirname_k = "k-analyze_nonlattice"; //create a new directory for k-analyze txt files
std::string new_filename_sp_bfosc  = "spectrum_bfosc.txt"; // create this new spectrum file
std::string new_filename_sp_afosc = "spectrum_afosc.txt"; // create this new spectrum file
std::string new_filename_sp_final  = "spectrum_final.txt"; // create this new spectrum file

std::string new_filename_sp_unpert = "spectrum_unpert.txt"; // create this new spectrum file
std::string new_filename_sp_newinfend = "spectrum_newinfend.txt"; // create this new spectrum file

//----------------------------------
//Variables for zeromode calculation
//----------------------------------
bool zeromode_switch = par_set[par_set_num]["zeromode_switch"];

//Array elements
const int N_zero=7;

DP sigma_c = pow((m_par/(2*m_par-1))*((3*m_par-1)/(m_par-1)),(m_par-1)/(4*m_par-2))*2/pow(2*m_par*(2*m_par-1),m_par/(4*m_par-2))*pow(mu_par*pow(M_par,m_par-1),1/(2*m_par -1));

DP FIXPSI = 2*pow(mu_par*pow(M_par,m_par-1),1/m_par);
DP FIXPHI = -sqrt(2)*pow(pw2(Cv_par)/g_par,1/n_par);

DP CNT = (-1)*Vbare(0.0,FIXPSI,FIXPHI);

//Vbare(0,FIXPSI,FIXPHI);       //constant term in the potential (set V=0 at the minimum)
//calculate log(Gravitational Potential) and log(Zeta) for each knum specified.

DP Gamma1,Gamma2,Gamma3;
DP OSCSTART,NEWINF_END_EFOLD;

std::vector<int> knum_zero = {200, 400, 634, 700};//knum for calculating k/aH

int k_target = knum_zero[3]; //target wave mode actually used for zeromode calculation

int OSCSTART_switch = 1;
//0: OSCSTART is set as the efold when sigma = sigma_c
//1: OSCSTART is set as the efold when epsilon = 1


//-------------------------------------------------
//Variables for zeromode w/ perturbation calculation
//-------------------------------------------------
bool perturbation_switch = par_set[par_set_num]["perturbation_switch"];//This needs to be true for perturbation calculation including lattice simulation calculation

bool lattice_kmodes_switch =  par_set[par_set_num]["lattice_kmodes_switch"]; //If this is true, it will use the modes calculated in lattice simulation for non lattice zeromode w/ perturb calculation

//Outputs
bool kanalyze_switch = true;// true:Calculate k-analyze, false:Don't calculate k-analyze
bool spectrum_switch = true;// true:Calculate final spectrum, false:Don't calculate final spectrum
bool spectrum_bfosc_switch = true; // true:Calculate spectrum right before oscillation period starts, false:Don't calculate
bool spectrum_afosc_switch = par_set[par_set_num]["spectrum_afosc_switch"];
// true:Calculate spectrum right after oscillation period ends, false:Don't calculate
bool spectrum_unpert_switch = true;
bool spectrum_newinfend_switch = true;

std::vector<std::string> sp_file_vec;

//Array elements
const int  N_pert=55;

//If lattice_kmodes_switch is false, we use the following range for calculation

//Takayama Master Thesis
//DP kfrom_Mpc = 3;//[Mpc^-1] Calculate from this k (Must be greater than 0)
//DP kto_Mpc = 1000;//[Mpc^-1] Calculate to this k
//int kinterval_knum = 10;// [knum units] Calculate with this interval of knum

//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (a)
DP kfrom_Mpc = par_set[par_set_num]["kfrom_Mpc"];//1.0e-4;//[Mpc^-1] Calculate from this k (Must be greater than 0)
DP kto_Mpc = par_set[par_set_num]["kto_Mpc"];//1.0e+4;//[Mpc^-1] Calculate to this k
int kfrom_knum = UC::kMpc_to_knum(kfrom_Mpc);
int kto_knum = UC::kMpc_to_knum(kto_Mpc);

int kinterval_knum = par_set[par_set_num]["kinterval_knum"];//100;// [knum units] Calculate with this interval of knum









//====================================================================
//                             Lattice Range
//====================================================================

//-----------
//File names
//-----------


std::string new_dirname_ed = "energy"; //create a new directory for energy vti files

std::string new_dirname_f = "field"; //create a new directory for field vti files

std::string new_filename_lattice = "lattice.txt";// create this new status file

//lattice simulation version of kAnalyze
std::string new_dirname_k_lattice = "k-analyze_lattice"; //create a new directory for k-analyze txt files

std::string new_dirname_sp_lattice = "spectrum_lattice"; //create a new directory for field vti files




//-------------------------------------------------
//Variables for calculating lattice range
//-------------------------------------------------

bool latticerange_switch = par_set[par_set_num]["latticerange_switch"]; // true:Set lattice range and calculate, false:Don't set lattice range and calculate
bool initialize_perturb_switch = true; // true:Initialize fluctuation, false:Don't initialize fluctuation (only calculate zeromode for lattice simulation)

int fluc_calc_switch  = 1;//Choose type of fluctuation initialization for scalar fields (for gravitational fluctuation, it is set to 2 regardless) 0:LatticeEasy case (when amplitudes of fluctuations are not predetermined) 1:when we use the amplitudes of predetermined fluctuations 2:gravitational perturbation

bool k_lattice_startfromlattice_switch = false; //If this is true, data output of kAnalyze_lattice starts from the time when lattice simulation starts. If false, then data output starts from the beginning (the beginning of hybrid inflation)

//Takayama Master Thesis
//double kfrom_Mpc_lattice = 3;//[Mpc^-1] Calculate from this k for lattice range
//double kto_Mpc_lattice = 400;//[Mpc^-1] Calculate to this k for lattice range

//Paper "Power spectrum of the density perturbations from smooth hybrid new inflation model" FIG1 (a)
double kfrom_Mpc_lattice = par_set[par_set_num]["kfrom_Mpc_lattice"];//1;//[Mpc^-1] Calculate from this k for lattice range
double kto_Mpc_lattice = par_set[par_set_num]["kto_Mpc_lattice"];//3000;//[Mpc^-1] Calculate to this k for lattice range
int kfrom_knum_lattice = UC::kMpc_to_knum(kfrom_Mpc_lattice);
int kto_knum_lattice = UC::kMpc_to_knum(kto_Mpc_lattice);
int kres_knum = (kto_knum_lattice - kfrom_knum) % kinterval_knum;
//This corresponds to the first mode in the upper range
int kstart_knum = kto_knum_lattice + ( kinterval_knum - kres_knum );


int dim = par_set[par_set_num]["dim"];
int N = par_set[par_set_num]["N"];//512; //Should be 2^n

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
double a_lattice_end;

int output_step = 2.0e+1;
int total_step  = 7.2e+4;//1.75e+4;//8.75e+3;
int max_loop    = total_step/output_step; // This many times vti files will be created
int st_output_step = 10;
int st_max_loop = output_step/st_output_step; // This many times data will be added to status.txt between the output of vti files
int spectrum_lattice_number = 30;// This many times spectrum will be created during lattice simulation

double t0 = 0;
double dt = 1.e-4;//5.0e-4;//1.e-3; //dt_pr


const int expansion = 1; // 0: no expansion, 1: self-consistent, 2: radiation dominant, 3: matter dominant
const int precision = 2;

const int screen_latticeloop_number = 100; //This many times loop will be displayed on the terminal for lattice simulation. If you want to show all loops, then set this number to max_loop.

const double metric_amp_rescale = 1000; // Right before lattice simulation ends, metric perturbation can be scaled down to a value metric_amp_rescale times lower if necessary. If this value is 1 then it means that the amplitude of the metric perturbation won't change.

/*f/M_p=f_pr/a
 dx_pr=mdx
 dt_pr=mdt/a
 
 A = 1/M_p (Therefore 1 in Planck units.), B = m, r = 1, s = -1,
 
 */



