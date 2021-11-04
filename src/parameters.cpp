#include "parameters.hpp"

std::string exist_dirname_k = "dataDI"; //remove this existing directory for k-analyze txt files
std::string new_dirname_k = "dataDI"; //create a new directory for k-analyze txt files
std::string filename_k = "kAnalyze"; // Head of the file name for k-analyze txt files

std::string exist_filename_sp  = "spectrum.txt";// remove this existing spectrum file
std::string new_filename_sp  = "spectrum.txt"; // create this new spectrum file

//Global Variables
DP CNT = (-1)*Vbare(0,FIXPSI,FIXPHI);
//Vbare(0,FIXPSI,FIXPHI);       //constant term in the potential (set V=0 at the minimum)
//calculate log(Gravitational Potential) and log(Zeta) for each knum specified.
DP k_comoving;
DP Gamma1,Gamma2,Gamma3;
const int timecount_max = 10000;

bool kanalyze_switch = true;// true:Calculate k-analyze, false:Don't calculate k-analyze
bool spectrum_switch = true;// true:Calculate spectrum, false:Don't calculate spectrum

DP kfromMpc = 220;//[Mpc] Calculate from this k
DP ktoMpc =  350;//[Mpc] Calculate to this k
int kinterval = 10;// [knum units] Calculate with this interval of knum
