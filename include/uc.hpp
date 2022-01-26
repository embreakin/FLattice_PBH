#ifndef _UC_H_
#define _UC_H_

#include "nr.h"

//Unit Conversion Functions
namespace UC {
    
    const DP Ck = 2.626E-61;  //Ck[MPl] for k = 10^-4 Mpc^{-1}
    
    DP knum_to_kMPl(int &k);
    int kMPl_to_knum(DP &k);
    
    DP knum_to_kMpc(int &k);
    int kMpc_to_knum(DP &k);
    
    DP kMpc_to_kMPl(DP &k);
    DP kMPl_to_kMpc(DP &k);
    
}
#endif
