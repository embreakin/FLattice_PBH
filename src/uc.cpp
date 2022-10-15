#include "uc.hpp"

//-----------------
//Unit Conversions
//-----------------

DP UC::knum_to_kMPl(int &k){
    
    return Ck*pow(10,k/100.);
    
}

int UC::kMPl_to_knum(DP &k){
    
    return 100*log10(k/Ck);
    
}

DP UC::knum_to_kMpc(int &k){
    
    return pow(10, k/100. - 4);
    
}

int UC::kMpc_to_knum(DP &k){
    
    return 100*(log10(k) + 4);
    
}

DP UC::kMpc_to_kMPl(DP &k){
    
    return Ck*pow(10, log10(k) + 4);
    
}

DP UC::kMPl_to_kMpc(DP &k){
    
    return pow(10, log10(k/Ck) - 4);
    
}

