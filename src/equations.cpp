#include <sstream>
#include "nr.h"
#include "parameters.hpp"

DP Pow(DP x,int i){
   DP val;
   int j;
   val=1.0;
   j=0;
   for (j=0;j<i;j++){
      val = val*x;
   };
   return val;
}

//Term inside the first parentheses of the first term in (24) of Takayama's published paper
DP f(DP x){
        DP val;
        val = (Pow(x,2*m_par)/(Pow(2,2*m_par)*Pow(M_par,2*(m_par-1)))) - mu_par*mu_par;
		if (x==FIXPSI) val=0;
        return val;
}

//Psi part of the second term in (24)
DP f1(DP x){
        DP val;
        val = m_par*Pow(x,2*m_par-1)/(Pow(2,2*m_par-1)*Pow(M_par,2*(m_par-1)));
        return val;
}

//Derivative of f1 with respect to x
DP f2(DP x){
        DP val;
        val = m_par*(2*m_par-1)*Pow(x,2*m_par-2)/(Pow(2,2*m_par-1)*Pow(M_par,2*(m_par-1)));
        return val;
}

//Derivative of f2 with respect to x
DP f3(DP x){
        DP val;
        val = m_par*(2*m_par-1)*(2*m_par-2)*Pow(x,2*m_par-3)/(Pow(2,2*m_par-1)*Pow(M_par,2*(m_par-1)));
        return val;
}


DP g0(DP x){
        DP val;
        val = Cv_par*Cv_par*x - g_par*Pow(x,n_par+1)/((n_par+1)*Pow(sqrt(2.0),n_par));
        return val;
}

//First term in (25) of Takayama's thesis
DP g1(DP x){
        DP val;
        val = Cv_par*Cv_par - g_par*Pow(x,n_par)/Pow(sqrt(2.0),n_par);
        return val;
}

//Derivative of g1 with respect to x
DP g2(DP x){
        DP val;
        val = (-1)*n_par*g_par*Pow(x,n_par-1)/Pow(sqrt(2.0),n_par);
        return val;
}

//Derivative of g2 with respect to x
DP g3(DP x){
        DP val;
        val = (-1)*n_par*(n_par-1)*g_par*Pow(x,n_par-2)/Pow(sqrt(2.0),n_par);
        return val;
}


DP Vbare(DP x,DP y,DP z){
        DP val;
        val = f(y)*f(y)*(1 + x*x*x*x/8 + y*y/2 + z*z/2)  // x -> sigma, y -> psi, z-> phi
             + g1(z)*g1(z) - CN_par*Cv_par*Cv_par*Cv_par*Cv_par*z*z/2
			 + x*x*f1(y)*f1(y) + f(y)*x*x*y*f1(y) - Cv_par*Cv_par*x*z*f(y);
        return val;
}

DP V(DP x,DP y,DP z){
        DP val;
        val = f(y)*f(y)*(1 + x*x*x*x/8 + y*y/2 + z*z/2)
             + g1(z)*g1(z) - CN_par*Cv_par*Cv_par*Cv_par*Cv_par*z*z/2
			 + x*x*f1(y)*f1(y) + f(y)*x*x*y*f1(y) - Cv_par*Cv_par*x*z*f(y)
             + CNT;
        return val;
}

//1st potential derivative with respect to sigma
DP V_1(DP x,DP y,DP z){
        DP val;
        val = x*x*x*f(y)*f(y)/2 + 2*x*f1(y)*f1(y) + 2*x*y*f(y)*f1(y) - Cv_par*Cv_par*z*f(y);
        return val;
}

//1st potential derivative with respect to psi
DP V_2(DP x,DP y,DP z){
        DP val;
        val = 2*f(y)*f1(y)*(1 + x*x*x*x/8 + y*y/2 + z*z/2)
		      + y*f(y)*f(y) + 2*x*x*f1(y)*f2(y) + x*x*f(y)*f1(y) + x*x*y*f(y)*f2(y) + x*x*y*f1(y)*f1(y) - Cv_par*Cv_par*x*z*f1(y);
		//if (y==0.13145341380001) val = 0;	  
		return val;
}

//1st potential derivative with respect to phi
DP V_3(DP x,DP y,DP z){
        DP val;
        val = z*f(y)*f(y) + 2*g1(z)*g2(z) - CN_par*Cv_par*Cv_par*Cv_par*Cv_par*z - Cv_par*Cv_par*x*f(y);
        return val;
}

//2nd potential derivative with respect to sigma-sigma
DP V_11(DP x,DP y,DP z){
        DP val;
        val = 3*x*x*f(y)*f(y)/2 + 2*f1(y)*f1(y) + 2*y*f(y)*f1(y);
		return val;
}

//2nd potential derivative with respect to sigma-psi
DP V_12(DP x,DP y,DP z){
        DP val;
        val = x*x*x*f(y)*f1(y) + 4*x*f1(y)*f2(y) + 2*x*f(y)*f1(y)
			 + 2*x*y*f(y)*f2(y) + 2*x*y*f1(y)*f1(y) - Cv_par*Cv_par*z*f1(y);
		return val;
}

//2nd potential derivative with respect to sigma-phi
DP V_13(DP x,DP y,DP z){
        DP val;
        val = (-1)*Cv_par*Cv_par*f(y);
		return val;
}

//2nd potential derivative with respect to psi-psi
DP V_22(DP x,DP y,DP z){
        DP val;
        val = (2*f1(y)*f1(y) + 2*f(y)*f2(y))*(1 + x*x*x*x/8 + y*y/2 + z*z/2)
			+ 4*y*f(y)*f1(y) + f(y)*f(y) + 2*x*x*f2(y)*f2(y) + 2*x*x*f1(y)*f3(y) + 2*x*x*f1(y)*f1(y)
	                + 3*x*x*y*f1(y)*f2(y) + 2*x*x*f(y)*f2(y) + x*x*y*f(y)*f3(y) - Cv_par*Cv_par*x*z*f2(y);
		return val;
}

//2nd potential derivative with respect to psi-phi
DP V_23(DP x,DP y,DP z){
        DP val;
        val = 2*z*f(y)*f1(y) - Cv_par*Cv_par*x*f1(y);
		return val;
}

//2nd potential derivative with respect to phi-phi
DP V_33(DP x,DP y,DP z){
        DP val;
        val = f(y)*f(y) + 2*g2(z)*g2(z) + 2*g1(z)*g3(z) - CN_par*Cv_par*Cv_par*Cv_par*Cv_par;
		return val;
}

//unperturbed scalar evolution1: dif_vphi1(phi1,phi2,phi3,vphi1,H)
DP dif_vphi1 (DP x,DP y,DP z,DP vx,DP HH){
	DP val;
	val = (-3)*vx - V_1(x,y,z)/HH - Gamma1*vx/HH;
	return val;
}

//unperturbed scalar evolution2: dif_vphi2(phi1,phi2,phi3,vphi2,H)
DP dif_vphi2 (DP x,DP y,DP z,DP vy,DP HH){
	DP val;
	val = (-3)*vy - V_2(x,y,z)/HH - Gamma2*vy/HH;
	return val;
}


//unperturbed scalar evolution3: dif_vphi3(phi1,phi2,phi3,vphi3,H)
DP dif_vphi3 (DP x,DP y,DP z,DP vz,DP HH){
	DP val;
	val = (-3)*vz - V_3(x,y,z)/HH - Gamma3*vz/HH;
	return val;
}


//unperturbed total rho: rho_tot(phi1,phi2,phi3,vphi1,vphi2,vphi3,rho_rad)
DP rho_tot (DP x,DP y,DP z,DP vx,DP vy,DP vz,DP rho_rad){
	DP val_rho;
	val_rho = vx*vx/2 + vy*vy/2 + vz*vz/2 + V(x,y,z) + rho_rad;
	return val_rho;
}

/*
//unperturbed total rho: rho_tot(phi1,phi2,vphi1,vphi2,rho_rad)
DP rho_12 (DP x,DP y,DP z,DP vx,DP vy){
	DP val_rho;
	val_rho = vx*vx/2 + vy*vy/2 + VSH(x,y,z);
	return val_rho;
}
*/

//unperturbed total p: p_tot(phi1,phi2,phi3,vphi1,vphi2,vphi3,rho_rad)
DP p_tot (DP x,DP y,DP z,DP vx,DP vy,DP vz,DP rho_rad){
	DP val_p;
	val_p = vx*vx/2 + vy*vy/2 + vz*vz/2 - V(x,y,z) + rho_rad/3;
	return val_p;
}

//unperturbed rho+p: rhoandp(vphi1,vphi2,vphi3,rho_rad)
DP rhoandp(DP vphi1,DP vphi2,DP vphi3,DP rho_rad){
	DP val_rhoandp;
	val_rhoandp = vphi1*vphi1 + vphi2*vphi2 + vphi3*vphi3 + 4*rho_rad/3;
	return val_rhoandp;
}

//unperturved Friedmann: Fri(phi1,phi2,phi3,v1,v2,v3,rho_rad)
DP Fri (DP x, DP y,DP z,DP vx, DP vy,DP vz, DP rr){
	DP val_H;
	val_H = sqrt(fabs(rho_tot(x,y,z,vx,vy,vz,rr)/3));
	return val_H;
}

//unperturbed rad evolution: dif_rad(v1,v2,v3,rho_rad,H)
DP dif_rad (DP vx, DP vy,DP vz,DP rr,DP HH){
	DP val_rad;
	val_rad = (-4)*rr + (Gamma1*vx*vx + Gamma2*vy*vy + Gamma3*vz*vz)/HH;
	return val_rad;
}

//perturbed vdelpsi_11 evolution: dif_DELv11(phi1,phi2,phi3,v1,psi11,psi21,psi31,DELv11,la,H,Phil,vPhil)
DP dif_DELv11 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot){
	DP val_DELv11;
	DP vala;
	vala=k_comoving/exp(la);
	val_DELv11 = (-1)*(3 + Gamma1/H)*vdx - vala*vala*dx/H - V_11(x,y,z)*dx/H - V_12(x,y,z)*dy/H -V_13(x,y,z)*dz/H + (2*V_1(x,y,z) + vx*Gamma1)*Pot/H - 4*vx*vPot/H;
	return val_DELv11;
}

//perturbed vdelpsi_12 evolution: dif_DELv12(phi1,phi2,phi3,v1,psi12,psi22,psi32,DELv12,la,H,Phil,vPhil)
DP dif_DELv12 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot){
	DP val_DELv12;
	DP vala;
	vala=k_comoving/exp(la);
	val_DELv12 = (-1)*(3 + Gamma1/H)*vdx - vala*vala*dx/H - V_11(x,y,z)*dx/H - V_12(x,y,z)*dy/H -V_13(x,y,z)*dz/H + (2*V_1(x,y,z) + vx*Gamma1)*Pot/H - 4*vx*vPot/H;
	return val_DELv12;
}

//perturbed vdelpsi_13 evolution: dif_DELv13(phi1,phi2,phi3,v1,psi13,psi23,psi33,DELv13,la,H,Phil,vPhil)
DP dif_DELv13 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot){
	DP val_DELv13;
	DP vala;
	vala=k_comoving/exp(la);
	val_DELv13 = (-1)*(3 + Gamma1/H)*vdx - vala*vala*dx/H - V_11(x,y,z)*dx/H - V_12(x,y,z)*dy/H -V_13(x,y,z)*dz/H + (2*V_1(x,y,z) + vx*Gamma1)*Pot/H - 4*vx*vPot/H;
	return val_DELv13;
}

//perturbed vdelpsi_21 evolution: dif_DELv21(phi1,phi2,phi3,v2,psi11,psi21,psi31,DELv21,la,H,Phil,vPhil)
DP dif_DELv21 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot){
	DP val_DELv21;
	DP vala;
	vala=k_comoving/exp(la);
	val_DELv21 = (-1)*(3 + Gamma2/H)*vdx - vala*vala*dy/H - V_12(x,y,z)*dx/H - V_22(x,y,z)*dy/H - V_23(x,y,z)*dz/H + (2*V_2(x,y,z) + vx*Gamma2)*Pot/H - 4*vx*vPot/H;
	return val_DELv21;
}

//perturbed vdelpsi_22 evolution: dif_DELv22(phi1,phi2,phi3,v2,psi12,psi22,psi32,DELv22,la,H,Phil,vPhil)
DP dif_DELv22 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot){
	DP val_DELv22;
	DP vala;
	vala=k_comoving/exp(la);
	val_DELv22 = (-1)*(3 + Gamma2/H)*vdx - vala*vala*dy/H - V_12(x,y,z)*dx/H - V_22(x,y,z)*dy/H - V_23(x,y,z)*dz/H + (2*V_2(x,y,z) + vx*Gamma2)*Pot/H - 4*vx*vPot/H;
	return val_DELv22;
}

//perturbed vdelpsi_23 evolution: dif_DELv23(phi1,phi2,phi3,v2,psi12,psi22,psi32,DELv23,la,H,Phil,vPhil)
DP dif_DELv23 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot){
	DP val_DELv23;
	DP vala;
	vala=k_comoving/exp(la);
	val_DELv23 = (-1)*(3 + Gamma2/H)*vdx - vala*vala*dy/H - V_12(x,y,z)*dx/H - V_22(x,y,z)*dy/H - V_23(x,y,z)*dz/H + (2*V_2(x,y,z) + vx*Gamma2)*Pot/H - 4*vx*vPot/H;
	return val_DELv23;
}

//perturbed vdelpsi_31 evolution: dif_DELv31(phi1,phi2,phi3,v2,psi11,psi21,psi31,DELv31,la,H,Phil,vPhil)
DP dif_DELv31 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot){
	DP val_DELv31;
	DP vala;
	vala=k_comoving/exp(la);
	val_DELv31 = (-1)*(3 + Gamma3/H)*vdx - vala*vala*dz/H - V_13(x,y,z)*dx/H - V_23(x,y,z)*dy/H - V_33(x,y,z)*dz/H + (2*V_3(x,y,z) + vx*Gamma3)*Pot/H - 4*vx*vPot/H;
	return val_DELv31;
}

//perturbed vdelpsi_32 evolution: dif_DELv32(phi1,phi2,phi3,v2,psi12,psi22,psi32,DELv32,la,H,Phil,vPhil)
DP dif_DELv32 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot){
	DP val_DELv32;
	DP vala;
	vala=k_comoving/exp(la);
	val_DELv32 = (-1)*(3 + Gamma3/H)*vdx - vala*vala*dz/H - V_13(x,y,z)*dx/H - V_23(x,y,z)*dy/H - V_33(x,y,z)*dz/H + (2*V_3(x,y,z) + vx*Gamma3)*Pot/H - 4*vx*vPot/H;
	return val_DELv32;
}

//perturbed vdelpsi_33 evolution: dif_DELv33(phi1,phi2,phi3,v2,psi13,psi23,psi33,DELv33,la,H,Phil,vPhil)
DP dif_DELv33 (DP x,DP y,DP z,DP vx,DP dx,DP dy,DP dz,DP vdx,DP la,DP H,DP Pot,DP vPot){
	DP val_DELv33;
	DP vala;
	vala=k_comoving/exp(la);
	val_DELv33 = (-1)*(3 + Gamma3/H)*vdx - vala*vala*dz/H - V_13(x,y,z)*dx/H - V_23(x,y,z)*dy/H - V_33(x,y,z)*dz/H + (2*V_3(x,y,z) + vx*Gamma3)*Pot/H - 4*vx*vPot/H;
	return val_DELv33;
}

//perturbed vPot evolution: dif_vPot(phi1,phi2,phi3,vphi1,vphi2,vphi3,psi1l,psi2l,psi3l,vpsi1l,vpsi2l,vpsi3l,Phil,vPhil,la,H)
DP dif_vPot (DP x,DP y,DP z,DP vx,DP vy,DP vz,DP dx,DP dy,DP dz,DP vdx,DP vdy,DP vdz,DP Pot,DP vPot,DP la,DP H){
	DP val_vPot;
	DP vala;
	vala=k_comoving/exp(la);
	val_vPot = (-5)*vPot - (vala*vala/3 + 4*V(x,y,z)/3)*Pot/H + (2*V_1(x,y,z)*dx + 2*V_2(x,y,z)*dy + 2*V_3(x,y,z)*dz - vx*vdx - vy*vdy - vz*vdz)/(3*H);
	return val_vPot;
}

void fieldperturbation(Vec_I_DP &y, Vec_O_DP &fields){
        int i,j;
        for(i=0;i<3;i++) fields[i]=0;
        for(j=0;j<3;j++){
                for(i=0;i<3;i++){
                        fields[j]=fields[j]+y[7+j*3+i]*y[7+j*3+i] + y[31+j*3+i]*y[31+j*3+i];
                };
		fields[j]=fields[j]/(2*M_PI*M_PI);		
        };
}

void adiabaticity(Vec_I_DP &y, Vec_O_DP &adia){
        int i,j;
        for(i=0;i<3;i++) adia[i]=0;
        for(j=0;j<3;j=j+2){
                for(i=0;i<3;i++){
                        adia[j]=adia[j]+y[7+j*3+i]*y[7+j*3+i] + y[31+j*3+i]*y[31+j*3+i];
                };
                adia[j]=adia[j]/(y[3+j]*y[3+j]);
        };
        if (y[4]!=0){
           j=1;
           for(i=0;i<3;i++){
                adia[j]=adia[j]+y[7+j*3+i]*y[7+j*3+i] + y[31+j*3+i]*y[31+j*3+i];
           };
           adia[j]=adia[j]/(y[3+j]*y[3+j]);
        };
        if (y[4]==0) adia[1]=1;
}

void isocurvatureness(Vec_I_DP &y, Vec_O_DP &iso){
        int i,j;
        for(i=0;i<3;i++) iso[i]=0;
        if (y[4]!=0){
                for(i=0;i<3;i++){
                        iso[0]=iso[0]+(y[7+i]/y[3] - y[10+i]/y[4])*(y[7+i]/y[3] - y[10+i]/y[4]) + (y[31+i]/y[3] - y[34+i]/y[4])*(y[31+i]/y[3] - y[34+i]/y[4]);
                };
        };
        for(i=0;i<3;i++){
                iso[1]=iso[1]+(y[7+i]/y[3] - y[13+i]/y[5])*(y[7+i]/y[3] - y[13+i]/y[5]) + (y[31+i]/y[3] - y[37+i]/y[5])*(y[31+i]/y[3] - y[37+i]/y[5]);
        };
        if (y[4]!=0){
                for(i=0;i<3;i++){
                        iso[2]=iso[2]+(y[10+i]/y[4] - y[13+i]/y[5])*(y[10+i]/y[4] - y[13+i]/y[5]) + (y[34+i]/y[4] - y[37+i]/y[5])*(y[34+i]/y[4] - y[37+i]/y[5]);
                };
        };
        if(y[4]==0){
                iso[0]=1;
                iso[2]=1;
        };
}

void density(Vec_I_DP &y, Vec_O_DP &dens){
        int i,j;
        DP H,med;
        H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
        for(i=0;i<4;i++) dens[i]=0;
        for(i=0;i<3;i++){
                med = y[3]*y[3]*y[25+i] + y[3]*y[16+i] + 3*H*y[3]*y[7+i] + V_1(y[0],y[1],y[2])*y[7+i];
                dens[0] = dens[0]+med*med;
                med = y[3]*y[3]*y[49+i] + y[3]*y[40+i] + 3*H*y[3]*y[31+i] + V_1(y[0],y[1],y[2])*y[31+i];
                dens[0] = dens[0]+med*med;
        };
        for(i=0;i<3;i++){
                med = y[4]*y[4]*y[25+i] + y[4]*y[19+i] + 3*H*y[4]*y[10+i] + V_2(y[0],y[1],y[2])*y[10+i];
                dens[1] = dens[1]+med*med;
                med = y[4]*y[4]*y[49+i] + y[4]*y[43+i] + 3*H*y[4]*y[34+i] + V_2(y[0],y[1],y[2])*y[34+i];
                dens[1] = dens[1]+med*med;
        };
        for(i=0;i<3;i++){
                med = y[5]*y[5]*y[25+i] + y[5]*y[22+i] + 3*H*y[5]*y[13+i] + V_3(y[0],y[1],y[2])*y[13+i];
                dens[2] = dens[2]+med*med;
                med = y[5]*y[5]*y[49+i] + y[5]*y[46+i] + 3*H*y[5]*y[37+i] + V_3(y[0],y[1],y[2])*y[37+i];
                dens[2] = dens[2]+med*med;
        };
        for(i=0;i<3;i++){
                med = y[3]*y[3]*y[25+i] + y[4]*y[4]*y[25+i] + y[5]*y[5]*y[25+i] + y[3]*y[16+i]+y[4]*y[19+i]+y[5]*y[22+i] + 3*H*(y[3]*y[7+i]+y[4]*y[10+i]+y[5]*y[13+i]) + V_1(y[0],y[1],y[2])*y[7+i]+V_2(y[0],y[1],y[2])*y[10+i]+V_3(y[0],y[1],y[2])*y[13+i];
                dens[3] = dens[3] + med*med;
                med = y[3]*y[3]*y[49+i] + y[4]*y[4]*y[49+i] + y[5]*y[5]*y[49+i] + y[3]*y[40+i]+y[4]*y[43+i]+y[5]*y[46+i] + 3*H*(y[3]*y[31+i]+y[4]*y[34+i]+y[5]*y[37+i]) + V_1(y[0],y[1],y[2])*y[31+i]+V_2(y[0],y[1],y[2])*y[34+i]+V_3(y[0],y[1],y[2])*y[37+i];
                dens[3] = dens[3] + med*med;
        };
}

void numberdens(DP x, Vec_I_DP &y, Vec_O_DP &numdens){
        int i,j;
		DP vala;
		vala = k_comoving/exp(x);
		vala = vala*vala;
        for(i=0;i<3;i++) numdens[i]=0;
        for(i=0;i<3;i++){
                numdens[0] = numdens[0] + y[7+i]*y[7+i] + y[31+i]*y[31+i];
                numdens[0] = numdens[0] + (y[16+i]*y[16+i] + y[40+i]*y[40+i])/(vala + V_11(y[0],y[1],y[2]));
                numdens[1] = numdens[1] + y[10+i]*y[10+i] + y[34+i]*y[34+i];
                numdens[1] = numdens[1] + (y[19+i]*y[19+i] + y[43+i]*y[43+i])/(vala + V_22(y[0],y[1],y[2]));
                numdens[2] = numdens[2] + y[13+i]*y[13+i] + y[37+i]*y[37+i];
                numdens[2] = numdens[2] + (y[22+i]*y[22+i] + y[46+i]*y[46+i])/(vala + V_33(y[0],y[1],y[2]));
        };
        numdens[0] = sqrt(fabs(vala + V_11(y[0],y[1],y[2])))*numdens[0]/2 - 0.5;
        numdens[1] = sqrt(fabs(vala + V_22(y[0],y[1],y[2])))*numdens[1]/2 - 0.5;
        numdens[2] = sqrt(fabs(vala + V_33(y[0],y[1],y[2])))*numdens[2]/2 - 0.5;
}

//    j=0-2  : zero modes of inflaton sigma, psi, and phi
//    j=3-5  : log(a) derivatives of inflaton sigma', psi', phi'
//    j=6    : energy density of radiation
//evolution equations for zero-mode
void unpert(const DP x, Vec_I_DP &y, Vec_O_DP &dydx)
{
    DP H;
    int ii;
	H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
	dydx[0] = y[3]/H;
	dydx[1] = y[4]/H;
	dydx[2] = y[5]/H;
	dydx[3] = dif_vphi1(y[0],y[1],y[2],y[3],H);
	dydx[4] = dif_vphi2(y[0],y[1],y[2],y[4],H);
	dydx[5] = dif_vphi3(y[0],y[1],y[2],y[5],H);
	dydx[6] = dif_rad(y[3],y[4],y[5],y[6],H);
}

//evoluition equations for zero mode with fixed psi
void unpertfix(const DP x, Vec_I_DP &y, Vec_O_DP &dydx)
{
    DP H;
    int ii;
	H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
	dydx[0] = y[3]/H;
	dydx[1] = 0;
	dydx[2] = y[5]/H;
	dydx[3] = dif_vphi1(y[0],y[1],y[2],y[3],H);
	dydx[4] = 0;
	dydx[5] = dif_vphi3(y[0],y[1],y[2],y[5],H);
	dydx[6] = dif_rad(y[3],y[4],y[5],y[6],H);
}

//evolution equations for zero mode with fixed sigma and psi
void unpertfixfix(const DP x, Vec_I_DP &y, Vec_O_DP &dydx)
{
    DP H;
    int ii;
	H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
	dydx[0] = 0;
	dydx[1] = 0;
	dydx[2] = y[5]/H;
	dydx[3] = 0;
	dydx[4] = 0;
	dydx[5] = dif_vphi3(y[0],y[1],y[2],y[5],H);
	dydx[6] = dif_rad(y[3],y[4],y[5],y[6],H);
}


//evolution equations for zero-mode and perturbation
void full(const DP x,Vec_I_DP &y, Vec_O_DP &dydx)
{
    DP H;
    int ii;
    H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
    dydx[0] = y[3]/H;
    dydx[1] = y[4]/H;
    dydx[2] = y[5]/H;
    dydx[3] = dif_vphi1(y[0],y[1],y[2],y[3],H);
    dydx[4] = dif_vphi2(y[0],y[1],y[2],y[4],H);
    dydx[5] = dif_vphi3(y[0],y[1],y[2],y[5],H);
    dydx[6] = dif_rad(y[3],y[4],y[5],y[6],H);
    for (ii=0;ii<9;ii++) dydx[ii+7] = y[ii+16]/H;
    dydx[16] = dif_DELv11(y[0],y[1],y[2],y[3],y[7],y[10],y[13],y[16],x,H,y[25],y[28]);
    dydx[17] = dif_DELv12(y[0],y[1],y[2],y[3],y[8],y[11],y[14],y[17],x,H,y[26],y[29]);
    dydx[18] = dif_DELv13(y[0],y[1],y[2],y[3],y[9],y[12],y[15],y[18],x,H,y[27],y[30]);
    dydx[19] = dif_DELv21(y[0],y[1],y[2],y[4],y[7],y[10],y[13],y[19],x,H,y[25],y[28]);
    dydx[20] = dif_DELv22(y[0],y[1],y[2],y[4],y[8],y[11],y[14],y[20],x,H,y[26],y[29]);
    dydx[21] = dif_DELv23(y[0],y[1],y[2],y[4],y[9],y[12],y[15],y[21],x,H,y[27],y[30]);
    dydx[22] = dif_DELv31(y[0],y[1],y[2],y[5],y[7],y[10],y[13],y[22],x,H,y[25],y[28]);
    dydx[23] = dif_DELv32(y[0],y[1],y[2],y[5],y[8],y[11],y[14],y[23],x,H,y[26],y[29]);
    dydx[24] = dif_DELv33(y[0],y[1],y[2],y[5],y[9],y[12],y[15],y[24],x,H,y[27],y[30]);
    for (ii=0;ii<3;ii++) dydx[ii+25] = y[ii+28]/H;
    for (ii=0;ii<3;ii++) dydx[ii+28] = dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[7+ii],y[10+ii],y[13+ii],y[16+ii],y[19+ii],y[22+ii],y[25+ii],y[28+ii],x,H);
    for (ii=0;ii<9;ii++) dydx[ii+31] = y[ii+40]/H;
    dydx[40] = dif_DELv11(y[0],y[1],y[2],y[3],y[31],y[34],y[37],y[40],x,H,y[49],y[52]);
    dydx[41] = dif_DELv12(y[0],y[1],y[2],y[3],y[32],y[35],y[38],y[41],x,H,y[50],y[53]);
    dydx[42] = dif_DELv13(y[0],y[1],y[2],y[3],y[33],y[36],y[39],y[42],x,H,y[51],y[54]);
    dydx[43] = dif_DELv21(y[0],y[1],y[2],y[4],y[31],y[34],y[37],y[43],x,H,y[49],y[52]);
    dydx[44] = dif_DELv22(y[0],y[1],y[2],y[4],y[32],y[35],y[38],y[44],x,H,y[50],y[53]);
    dydx[45] = dif_DELv23(y[0],y[1],y[2],y[4],y[33],y[36],y[39],y[45],x,H,y[51],y[54]);
    dydx[46] = dif_DELv31(y[0],y[1],y[2],y[5],y[31],y[34],y[37],y[46],x,H,y[49],y[52]);
    dydx[47] = dif_DELv32(y[0],y[1],y[2],y[5],y[32],y[35],y[38],y[47],x,H,y[50],y[53]);
    dydx[48] = dif_DELv33(y[0],y[1],y[2],y[5],y[33],y[36],y[39],y[48],x,H,y[51],y[54]);
    for (ii=0;ii<3;ii++) dydx[ii+49] = y[ii+52]/H;
    for (ii=0;ii<3;ii++) dydx[ii+52] = dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[31+ii],y[34+ii],y[37+ii],y[40+ii],y[43+ii],y[46+ii],y[49+ii],y[52+ii],x,H);
}

//evolution equations for zero-mode with fixed sigma and psi perturbations
void sigma_psi_nopert(const DP x,Vec_I_DP &y, Vec_O_DP &dydx)
{
    DP H;
    int ii;
    H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
    dydx[0] = y[3]/H;
    dydx[1] = y[4]/H;
    dydx[2] = y[5]/H;
    dydx[3] = dif_vphi1(y[0],y[1],y[2],y[3],H);
    dydx[4] = dif_vphi2(y[0],y[1],y[2],y[4],H);
    dydx[5] = dif_vphi3(y[0],y[1],y[2],y[5],H);
    dydx[6] = dif_rad(y[3],y[4],y[5],y[6],H);
    for (ii=0;ii<6;ii++) dydx[ii+7] = 0;
    for (ii=0;ii<3;ii++) dydx[ii+13] = y[ii+22]/H;
    for (ii=0;ii<6;ii++) dydx[ii+16] = 0;
    dydx[22] = dif_DELv31(y[0],y[1],y[2],y[5],0,0,y[13],y[22],x,H,y[25],y[28]);
    dydx[23] = dif_DELv32(y[0],y[1],y[2],y[5],0,0,y[14],y[23],x,H,y[26],y[29]);
    dydx[24] = dif_DELv33(y[0],y[1],y[2],y[5],0,0,y[15],y[24],x,H,y[27],y[30]);
    for (ii=0;ii<3;ii++) dydx[ii+25] = y[ii+28]/H;
    for (ii=0;ii<3;ii++) dydx[ii+28] = dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],0,0,y[13+ii],0, 0,y[22+ii],y[25+ii],y[28+ii],x,H);//dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[7+ii],y[10+ii],y[13+ii],y[16+ii],y[19+ii],y[22+ii],y[25+ii],y[28+ii],x,H);
    for (ii=0;ii<6;ii++) dydx[ii+31] = 0;
    for (ii=0;ii<3;ii++) dydx[ii+37] = y[ii+46]/H;
    for (ii=0;ii<6;ii++) dydx[ii+40] = 0;
    dydx[46] = dif_DELv31(y[0],y[1],y[2],y[5],0,0,y[37],y[46],x,H,y[49],y[52]);
    dydx[47] = dif_DELv32(y[0],y[1],y[2],y[5],0,0,y[38],y[47],x,H,y[50],y[53]);
    dydx[48] = dif_DELv33(y[0],y[1],y[2],y[5],0,0,y[39],y[48],x,H,y[51],y[54]);
    for (ii=0;ii<3;ii++) dydx[ii+49] = y[ii+52]/H;
    for (ii=0;ii<3;ii++) dydx[ii+52] = dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],0,0,y[37+ii],0,0,y[46+ii],y[49+ii],y[52+ii],x,H);//dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[31+ii],y[34+ii],y[37+ii],y[40+ii],y[43+ii],y[46+ii],y[49+ii],y[52+ii],x,H);
}

//evolution equations for zero-mode and perturbation without gravitational potential perturbation
void nogravpert(const DP x,Vec_I_DP &y, Vec_O_DP &dydx)
{
    DP H;
    int ii;
    H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
    dydx[0] = y[3]/H;
    dydx[1] = y[4]/H;
    dydx[2] = y[5]/H;
    dydx[3] = dif_vphi1(y[0],y[1],y[2],y[3],H);
    dydx[4] = dif_vphi2(y[0],y[1],y[2],y[4],H);
    dydx[5] = dif_vphi3(y[0],y[1],y[2],y[5],H);
    dydx[6] = dif_rad(y[3],y[4],y[5],y[6],H);
    for (ii=0;ii<9;ii++) dydx[ii+7] = y[ii+16]/H;
    dydx[16] = dif_DELv11(y[0],y[1],y[2],y[3],y[7],y[10],y[13],y[16],x,H,0,0);
    dydx[17] = dif_DELv12(y[0],y[1],y[2],y[3],y[8],y[11],y[14],y[17],x,H,0,0);
    dydx[18] = dif_DELv13(y[0],y[1],y[2],y[3],y[9],y[12],y[15],y[18],x,H,0,0);
    dydx[19] = dif_DELv21(y[0],y[1],y[2],y[4],y[7],y[10],y[13],y[19],x,H,0,0);
    dydx[20] = dif_DELv22(y[0],y[1],y[2],y[4],y[8],y[11],y[14],y[20],x,H,0,0);
    dydx[21] = dif_DELv23(y[0],y[1],y[2],y[4],y[9],y[12],y[15],y[21],x,H,0,0);
    dydx[22] = dif_DELv31(y[0],y[1],y[2],y[5],y[7],y[10],y[13],y[22],x,H,0,0);
    dydx[23] = dif_DELv32(y[0],y[1],y[2],y[5],y[8],y[11],y[14],y[23],x,H,0,0);
    dydx[24] = dif_DELv33(y[0],y[1],y[2],y[5],y[9],y[12],y[15],y[24],x,H,0,0);
    for (ii=0;ii<3;ii++) dydx[ii+25] = 0;
    for (ii=0;ii<3;ii++) dydx[ii+28] = dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[7+ii],y[10+ii],y[13+ii],y[16+ii],y[19+ii],y[22+ii],0,0,x,H);
    for (ii=0;ii<9;ii++) dydx[ii+31] = y[ii+40]/H;
    dydx[40] = dif_DELv11(y[0],y[1],y[2],y[3],y[31],y[34],y[37],y[40],x,H,0,0);
    dydx[41] = dif_DELv12(y[0],y[1],y[2],y[3],y[32],y[35],y[38],y[41],x,H,0,0);
    dydx[42] = dif_DELv13(y[0],y[1],y[2],y[3],y[33],y[36],y[39],y[42],x,H,0,0);
    dydx[43] = dif_DELv21(y[0],y[1],y[2],y[4],y[31],y[34],y[37],y[43],x,H,0,0);
    dydx[44] = dif_DELv22(y[0],y[1],y[2],y[4],y[32],y[35],y[38],y[44],x,H,0,0);
    dydx[45] = dif_DELv23(y[0],y[1],y[2],y[4],y[33],y[36],y[39],y[45],x,H,0,0);
    dydx[46] = dif_DELv31(y[0],y[1],y[2],y[5],y[31],y[34],y[37],y[46],x,H,0,0);
    dydx[47] = dif_DELv32(y[0],y[1],y[2],y[5],y[32],y[35],y[38],y[47],x,H,0,0);
    dydx[48] = dif_DELv33(y[0],y[1],y[2],y[5],y[33],y[36],y[39],y[48],x,H,0,0);
    for (ii=0;ii<3;ii++) dydx[ii+49] = 0;
    for (ii=0;ii<3;ii++) dydx[ii+52] = dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[31+ii],y[34+ii],y[37+ii],y[40+ii],y[43+ii],y[46+ii],0,0,x,H);
}

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
void unpfull(const DP x,Vec_I_DP &y, Vec_O_DP &dydx)
{
    DP H;
	int ii;
	H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
	dydx[0] = y[3]/H;
	dydx[1] = y[4]/H;
	dydx[2] = y[5]/H;
	dydx[3] = dif_vphi1(y[0],y[1],y[2],y[3],H);
	dydx[4] = dif_vphi2(y[0],y[1],y[2],y[4],H);
	dydx[5] = dif_vphi3(y[0],y[1],y[2],y[5],H);
	dydx[6] = dif_rad(y[3],y[4],y[5],y[6],H);
	for (ii=7;ii<55;ii++) dydx[ii] = 0;
}

//evolution equations for zero-mode and perturbaton with fixed sigma and psi and their perturbaions
void newinf(const DP x,Vec_I_DP &y, Vec_O_DP &dydx)
{
    DP H;
	int ii;
	H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
	dydx[0] = 0;
	dydx[1] = 0;
	dydx[2] = y[5]/H;
	dydx[3] = 0;
	dydx[4] = 0;
	dydx[5] = dif_vphi3(y[0],y[1],y[2],y[5],H);
	dydx[6] = dif_rad(y[3],y[4],y[5],y[6],H);
	for (ii=0;ii<6;ii++) dydx[ii+7] = 0;
	for (ii=0;ii<3;ii++) dydx[ii+13] = y[ii+22]/H;
	for (ii=0;ii<6;ii++) dydx[ii+16] = 0;
	dydx[22] = dif_DELv31(y[0],y[1],y[2],y[5],y[7],y[10],y[13],y[22],x,H,y[25],y[28]);
	dydx[23] = dif_DELv32(y[0],y[1],y[2],y[5],y[8],y[11],y[14],y[23],x,H,y[26],y[29]);
	dydx[24] = dif_DELv33(y[0],y[1],y[2],y[5],y[9],y[12],y[15],y[24],x,H,y[27],y[30]);
	for (ii=0;ii<3;ii++) dydx[ii+25] = y[ii+28]/H;
    for (ii=0;ii<3;ii++) dydx[ii+28] =  dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[7+ii],y[10+ii],y[13+ii],y[16+ii],y[19+ii],y[22+ii],y[25+ii],y[28+ii],x,H);
	for (ii=0;ii<6;ii++) dydx[ii+31] = 0;
	for (ii=0;ii<3;ii++) dydx[ii+37] = y[ii+46]/H;
	for (ii=0;ii<6;ii++) dydx[ii+40] = 0;
	dydx[46] = dif_DELv31(y[0],y[1],y[2],y[5],y[31],y[34],y[37],y[46],x,H,y[49],y[52]);
	dydx[47] = dif_DELv32(y[0],y[1],y[2],y[5],y[32],y[35],y[38],y[47],x,H,y[50],y[53]);
	dydx[48] = dif_DELv33(y[0],y[1],y[2],y[5],y[33],y[36],y[39],y[48],x,H,y[51],y[54]);
	for (ii=0;ii<3;ii++) dydx[ii+49] = y[ii+52]/H;
    for (ii=0;ii<3;ii++) dydx[ii+52] =  dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[31+ii],y[34+ii],y[37+ii],y[40+ii],y[43+ii],y[46+ii],y[49+ii],y[52+ii],x,H);
}
//evolution equations for zero-mode and perturbaton with fixed sigma and psi
void fixfix(const DP x,Vec_I_DP &y, Vec_O_DP &dydx)
{
    DP H;
	int ii;
	H = Fri(y[0],y[1],y[2],y[3],y[4],y[5],y[6]);
	dydx[0] = 0;
	dydx[1] = 0;
	dydx[2] = y[5]/H;
	dydx[3] = 0;
	dydx[4] = 0;
	dydx[5] = dif_vphi3(y[0],y[1],y[2],y[5],H);
	dydx[6] = dif_rad(y[3],y[4],y[5],y[6],H);
	for (ii=0;ii<9;ii++) dydx[ii+7] = y[ii+16]/H;
	dydx[16] = dif_DELv11(y[0],y[1],y[2],y[3],y[7],y[10],y[13],y[16],x,H,y[25],y[28]);
	dydx[17] = dif_DELv12(y[0],y[1],y[2],y[3],y[8],y[11],y[14],y[17],x,H,y[26],y[29]);
	dydx[18] = dif_DELv13(y[0],y[1],y[2],y[3],y[9],y[12],y[15],y[18],x,H,y[27],y[30]);
	dydx[19] = dif_DELv21(y[0],y[1],y[2],y[4],y[7],y[10],y[13],y[19],x,H,y[25],y[28]);
	dydx[20] = dif_DELv22(y[0],y[1],y[2],y[4],y[8],y[11],y[14],y[20],x,H,y[26],y[29]);
	dydx[21] = dif_DELv23(y[0],y[1],y[2],y[4],y[9],y[12],y[15],y[21],x,H,y[27],y[30]);
	dydx[22] = dif_DELv31(y[0],y[1],y[2],y[5],y[7],y[10],y[13],y[22],x,H,y[25],y[28]);
	dydx[23] = dif_DELv32(y[0],y[1],y[2],y[5],y[8],y[11],y[14],y[23],x,H,y[26],y[29]);
	dydx[24] = dif_DELv33(y[0],y[1],y[2],y[5],y[9],y[12],y[15],y[24],x,H,y[27],y[30]);
	for (ii=0;ii<3;ii++) dydx[ii+25] = y[ii+28]/H;
	for (ii=0;ii<3;ii++) dydx[ii+28] = dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[7+ii],y[10+ii],y[13+ii],y[16+ii],y[19+ii],y[22+ii],y[25+ii],y[28+ii],x,H);
	for (ii=0;ii<9;ii++) dydx[ii+31] = y[ii+40]/H;
	dydx[40] = dif_DELv11(y[0],y[1],y[2],y[3],y[31],y[34],y[37],y[40],x,H,y[49],y[52]);
	dydx[41] = dif_DELv12(y[0],y[1],y[2],y[3],y[32],y[35],y[38],y[41],x,H,y[50],y[53]);
	dydx[42] = dif_DELv13(y[0],y[1],y[2],y[3],y[33],y[36],y[39],y[42],x,H,y[51],y[54]);
	dydx[43] = dif_DELv21(y[0],y[1],y[2],y[4],y[31],y[34],y[37],y[43],x,H,y[49],y[52]);
	dydx[44] = dif_DELv22(y[0],y[1],y[2],y[4],y[32],y[35],y[38],y[44],x,H,y[50],y[53]);
	dydx[45] = dif_DELv23(y[0],y[1],y[2],y[4],y[33],y[36],y[39],y[45],x,H,y[51],y[54]);
	dydx[46] = dif_DELv31(y[0],y[1],y[2],y[5],y[31],y[34],y[37],y[46],x,H,y[49],y[52]);
	dydx[47] = dif_DELv32(y[0],y[1],y[2],y[5],y[32],y[35],y[38],y[47],x,H,y[50],y[53]);
	dydx[48] = dif_DELv33(y[0],y[1],y[2],y[5],y[33],y[36],y[39],y[48],x,H,y[51],y[54]);
	for (ii=0;ii<3;ii++) dydx[ii+49] = y[ii+52]/H;
	for (ii=0;ii<3;ii++) dydx[ii+52] = dif_vPot(y[0],y[1],y[2],y[3],y[4],y[5],y[31+ii],y[34+ii],y[37+ii],y[40+ii],y[43+ii],y[46+ii],y[49+ii],y[52+ii],x,H);
}


