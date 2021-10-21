#include "smpar.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

smpar::smpar(){pi=M_PI;
    scale_mu=4.8;
    m_d=0.00467,m_c=1.27,m_s=0.093,m_b=4.18;
    m_e=0.511*pow(10.0,-3.0),m_mu=105.658*pow(10.0,-3.0),m_tau=1.77686;
    
    m_Bd=5.27965,m_Bs=5.36688,m_K=0.497614,m_Kst=0.896;
    
    tau_Bd=1.519*1.52*pow(10.0,12.0),tau_Bs=1.527*1.52*pow(10.0,12.0),D_Gamma_dbar=0.001,D_Gamma_sbar=0.124,f_Bd=0.1905,f_Bs=0.2303;
    
    abs_VtbVtdStr=0.00853871,abs_VtbVtsStr=0.0408726,G_F=1.16637*pow(10.0,-5.0),alphae=1.0/127.944;
    
    C_1=-0.257,C_2=1.009,C_3=-0.005,C_4=-0.078,C_5=0.000,C_6=0.001,C_7effRe=-0.304,C_7effIm=0.0,C_8eff=-0.167,C_9Re=4.211,C_9Im=0.0,C_10Re=-4.103,C_10Im=0.0;
    
    C_7RHRe=0.0,C_7RHIm=0.0,C_9RHRe=0.0,C_9RHIm=0.0,C_10RHRe=0.0,C_10RHIm=0.0,C_SRe=0.0,C_SIm=0.0,C_SRHRe=0.0,C_SRHIm=0.0,C_PRe=0.0,C_PIm=0.0,C_PRHRe=0.0,C_PRHIm=0.0,C_TRe=0.0,C_TIm=0.0,C_T5Re=0.0,C_T5Im=0.0;
}
//Quark_Mass
double smpar::mu(){return scale_mu;}
double smpar::md(){return m_d;}
double smpar::mc(){return m_c;}
double smpar::ms(){return m_s;}
double smpar::mb(){return m_b;}
//Lepton_Mass
double smpar::me(){return m_e;}
double smpar::mmu(){return m_mu;}
double smpar::mtau(){return m_tau;}
//Meson_Mass
double smpar::mBd(){return m_Bd;}
double smpar::mBs(){return m_Bs;}
double smpar::mK(){return m_K;}
double smpar::mKst(){return m_Kst;}
//Width_Lifetime
double smpar::tauBd(){return tau_Bd;}
double smpar::tauBs(){return tau_Bs;}
double smpar::DGamma_dbar(){return D_Gamma_dbar;}
double smpar::DGamma_sbar(){return D_Gamma_sbar;}
double smpar::fBd(){return f_Bd;}
double smpar::fBs(){return f_Bs;}
//Coupling_Constants
double smpar::absVtbVtdStr(){return abs_VtbVtdStr;}
double smpar::absVtbVtsStr(){return abs_VtbVtsStr;}
double smpar::GF(){return G_F;}
double smpar::alpha_e(){return alphae;}
//SM_Wilson_Coefficients
double smpar::C1(){return C_1;}
double smpar::C2(){return C_2;}
double smpar::C3(){return C_3;}
double smpar::C4(){return C_4;}
double smpar::C5(){return C_5;}
double smpar::C6(){return C_6;}
double smpar::C7effRe(){return C_7effRe;}
double smpar::C7effIm(){return C_7effIm;}
double smpar::C8eff(){return C_8eff;}
double smpar::C9Re(){return C_9Re;}
double smpar::C9Im(){return C_9Im;}
double smpar::C10Re(){return C_10Re;}
double smpar::C10Im(){return C_10Im;}
//NP_Wilson_Coefficients
double smpar::C7RHRe(){return C_7RHRe;}
double smpar::C7RHIm(){return C_7RHIm;}
double smpar::C9RHRe(){return C_9RHRe;}
double smpar::C9RHIm(){return C_9RHIm;}
double smpar::C10RHRe(){return C_10RHRe;}
double smpar::C10RHIm(){return C_10RHIm;}
double smpar::CSRe(){return C_SRe;}
double smpar::CSIm(){return C_SIm;}
double smpar::CSRHRe(){return C_SRHRe;}
double smpar::CSRHIm(){return C_SRHIm;}
double smpar::CPRe(){return C_PRe;}
double smpar::CPIm(){return C_PIm;}
double smpar::CPRHRe(){return C_PRHRe;}
double smpar::CPRHIm(){return C_PRHIm;}
double smpar::CTRe(){return C_TRe;}
double smpar::CTIm(){return C_TIm;}
double smpar::CT5Re(){return C_T5Re;}
double smpar::CT5Im(){return C_T5Im;}
//
/*double smpar::tpd(){return pow(m_Bd+m_Kst,2.0);}
 double smpar::tps(){return pow(m_Bs+m_Kst,2.0);}
 double smpar::tmd(){return pow(m_Bd-m_Kst,2.0);}
 double smpar::tms(){return pow(m_Bs-m_Kst,2.0);}
 double smpar::tzd(){return tpd()*(1.0-sqrt(1.0-tmd()/tpd()));}
 double smpar::tzs(){return tps()*(1.0-sqrt(1.0-tms()/tps()));}*/

