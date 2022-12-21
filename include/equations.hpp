//Doxygen
/**
* @file    equations.hpp
* @brief    Equations header file
* @author   Francis Otani
* @date
* @details
*/

#ifndef _EQUATIONS_H_
#define _EQUATIONS_H_


#include "nr.h"
#include "parameters.hpp"

//-------------------
//Pow: Pow(x,i)
DP Pow(DP x,int i);

//Term inside the first parentheses of the first term in (24) of Takayama's published paper
DP f(DP x);

//Psi part of the second term in (24)
DP f1(DP x);

//Derivative of f1 with respect to x
DP f2(DP x);

//Derivative of f2 with respect to x
DP f3(DP x);


DP g0(DP x);

//First term in (25) of Takayama's thesis
DP g1(DP x);

//Derivative of g1 with respect to x
DP g2(DP x);

//Derivative of g2 with respect to x
DP g3(DP x);

DP Vbare(DP x,DP y,DP z);

DP V(DP x,DP y,DP z);

//1st potential derivative with respect to sigma
DP V_1(DP x,DP y,DP z);

//1st potential derivative with respect to psi
DP V_2(DP x,DP y,DP z);

//1st potential derivative with respect to phi
DP V_3(DP x,DP y,DP z);

//2nd potential derivative with respect to sigma-sigma
DP V_11(DP x,DP y,DP z);

//2nd potential derivative with respect to sigma-psi
DP V_12(DP x,DP y,DP z);

//2nd potential derivative with respect to sigma-phi
DP V_13(DP x,DP y,DP z);

//2nd potential derivative with respect to psi-psi
DP V_22(DP x,DP y,DP z);

//2nd potential derivative with respect to psi-phi
DP V_23(DP x,DP y,DP z);

//2nd potential derivative with respect to phi-phi
DP V_33(DP x,DP y,DP z);

//unperturbed scalar evolution1: dif_vphi1(phi1,phi2,phi3,vphi1,H)
DP dif_vphi1 (DP x,DP y,DP z,DP vx,DP HH);

//unperturbed scalar evolution2: dif_vphi2(phi1,phi2,phi3,vphi2,H)
DP dif_vphi2 (DP x,DP y,DP z,DP vy,DP HH);


//unperturbed scalar evolution3: dif_vphi3(phi1,phi2,phi3,vphi3,H)
DP dif_vphi3 (DP x,DP y,DP z,DP vz,DP HH);

//unperturbed total Kinetic energy:
DP K_tot (DP vx,DP vy,DP vz);

//unperturbed total rho: rho_tot(phi1,phi2,phi3,vphi1,vphi2,vphi3,rho_rad)
DP rho_tot (DP x,DP y,DP z,DP vx,DP vy,DP vz,DP rho_rad);

/*
 //unperturbed total rho: rho_tot(phi1,phi2,vphi1,vphi2,rho_rad)
 DP rho_12 (DP x,DP y,DP z,DP vx,DP vy);
 */

//unperturbed total p: p_tot(phi1,phi2,phi3,vphi1,vphi2,vphi3,rho_rad)
DP p_tot (DP x,DP y,DP z,DP vx,DP vy,DP vz,DP rho_rad);

//unperturbed rho+p: rhoandp(vphi1,vphi2,vphi3,rho_rad)
DP rhoandp(DP vphi1,DP vphi2,DP vphi3,DP rho_rad);

//unperturved Friedmann: Fri(phi1,phi2,phi3,v1,v2,v3,rho_rad)
DP Fri (DP x, DP y,DP z,DP vx, DP vy,DP vz, DP rr);

//unperturbed rad evolution: dif_rad(v1,v2,v3,rho_rad,H)
DP dif_rad (DP vx, DP vy,DP vz,DP rr,DP HH);

//perturbed vdelpsi_11 evolution: dif_DELv11(phi1,phi2,phi3,v1,psi11,psi21,psi31,DELv11,la,H,Phil,vPhil)
DP dif_DELv11 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot,DP k_comoving);

//perturbed vdelpsi_12 evolution: dif_DELv12(phi1,phi2,phi3,v1,psi12,psi22,psi32,DELv12,la,H,Phil,vPhil)
DP dif_DELv12 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot,DP k_comoving);

//perturbed vdelpsi_13 evolution: dif_DELv13(phi1,phi2,phi3,v1,psi13,psi23,psi33,DELv13,la,H,Phil,vPhil)
DP dif_DELv13 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot,DP k_comoving);

//perturbed vdelpsi_21 evolution: dif_DELv21(phi1,phi2,phi3,v2,psi11,psi21,psi31,DELv21,la,H,Phil,vPhil)
DP dif_DELv21 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot,DP k_comoving);

//perturbed vdelpsi_22 evolution: dif_DELv22(phi1,phi2,phi3,v2,psi12,psi22,psi32,DELv22,la,H,Phil,vPhil)
DP dif_DELv22 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot,DP k_comoving);

//perturbed vdelpsi_23 evolution: dif_DELv23(phi1,phi2,phi3,v2,psi12,psi22,psi32,DELv23,la,H,Phil,vPhil)
DP dif_DELv23 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot,DP k_comoving);

//perturbed vdelpsi_31 evolution: dif_DELv31(phi1,phi2,phi3,v2,psi11,psi21,psi31,DELv31,la,H,Phil,vPhil)
DP dif_DELv31 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot,DP k_comoving);

//perturbed vdelpsi_32 evolution: dif_DELv32(phi1,phi2,phi3,v2,psi12,psi22,psi32,DELv32,la,H,Phil,vPhil)
DP dif_DELv32 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot,DP k_comoving);

//perturbed vdelpsi_33 evolution: dif_DELv33(phi1,phi2,phi3,v2,psi13,psi23,psi33,DELv33,la,H,Phil,vPhil)
DP dif_DELv33 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot,DP k_comoving);

//perturbed vPot evolution: dif_vPot(phi1,phi2,phi3,vphi1,vphi2,vphi3,psi1l,psi2l,psi3l,vpsi1l,vpsi2l,vpsi3l,Phil,vPhil,la,H)
DP dif_vPot (DP x,DP y,DP z,DP vx,DP vy,DP vz,DP dx,DP dy,DP dz,DP vdx,DP vdy,DP vdz,DP Pot,DP vPot,DP la,DP H,DP k_comoving);

void fieldperturbation(Vec_I_DP &y, Vec_O_DP &fields);

void adiabaticity(Vec_I_DP &y, Vec_O_DP &adia);

void isocurvatureness(Vec_I_DP &y, Vec_O_DP &iso);

void density(Vec_I_DP &y, Vec_O_DP &dens);

void numberdens(const DP x, Vec_I_DP &y, Vec_O_DP &numdens, const DP k_comoving);

//    j=0-2  : zero modes of inflaton sigma, psi, and phi
//    j=3-5  : log(a) derivatives of inflaton sigma', psi', phi'
//    j=6    : energy density of radiation
//evolution equations for zero-mode
void unpert(const DP x, Vec_I_DP &y, Vec_O_DP &dydx, const DP k_comoving);

//evoluition equations for zero mode with fixed psi
void unpertfix(const DP x, Vec_I_DP &y, Vec_O_DP &dydx, const DP k_comoving);

//evolution equations for zero mode with fixed sigma and psi
void unpertfixfix(const DP x, Vec_I_DP &y, Vec_O_DP &dydx, const DP k_comoving);

//evolution equations for zero-mode and perturbation
void full(const DP x,Vec_I_DP &y, Vec_O_DP &dydx, const DP k_comoving);

//evolution equations for zero-mode with fixed sigma and psi perturbations
void sigma_psi_nopert(const DP x,Vec_I_DP &y, Vec_O_DP &dydx, const DP k_comoving);

//evolution equations for zero-mode and perturbation without gravitational potential perturbation
void nogravpert(const DP x,Vec_I_DP &y, Vec_O_DP &dydx, const DP k_comoving);

//xp[i] stores integration variable (log(a)) for each steps[i].
//delp_p[i][j] stores variables for each steps[i]. [j] specifies variables as follows:
//    j=0-2  : zero modes of inflaton sigma, psi, and phi
//    j=3-5  : log(a) derivatives of inflaton sigma', psi', phi'
//    j=6    : energy density of radiation
//    j=7-15 : mode functions of field perturbation: delta_{sigma,sigma}, delta_{sigma,psi}, delta{sigma,phi}, delta_{psi,sigma}, etc...
//    j=16-24: log(a) derivatives of mode functions
//    j=25-27: mode functions of gravitational potential perturbation: delta Phi_{sigma}, delta Phi_{psi}, delta Phi_{psi}
//    j=28-30: log(a) derivatives of gravitational potential perturbation
//    j=31-54: complex conjugate of [7]-[30]
//unp[i] stores zero mode variables(delp_p[i][0-6]
//tr[j]: buffer for delp_p[i][j]

//evolution equations for zero-mode and perturbation with fixed perturbations
void unpfull(const DP x,Vec_I_DP &y, Vec_O_DP &dydx, const DP k_comoving);

//evolution equations for zero-mode and perturbaton with fixed sigma and psi and their perturbaions
void newinf(const DP x,Vec_I_DP &y, Vec_O_DP &dydx, const DP k_comoving);

//evolution equations for zero-mode and perturbaton with fixed sigma and psi
void fixfix(const DP x,Vec_I_DP &y, Vec_O_DP &dydx, const DP k_comoving);


#endif
