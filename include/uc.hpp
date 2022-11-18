//Doxygen
/**
* @file   uc.hpp
* @brief    Unit conversions header file
* @author   Francis Otani
* @date
* @details
*/


#ifndef _UC_H_
#define _UC_H_

#include "nr.h"

//Unit Conversion Functions
namespace UC {
    
    const double Ck = 2.626E-61;  //Ck[MPl] for k = 10^-4 Mpc^{-1}
    
    double knum_to_kMPl(int k);
    int kMPl_to_knum(double k);
    
    double knum_to_kMpc(int k);
    int kMpc_to_knum(double k);
    
    double kMpc_to_kMPl(double k);
    double kMPl_to_kMpc(double k);
    
    double xMPl_to_xMpc(double x);
    
}
#endif
