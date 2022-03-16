#include "SM_parameters_mc.h"

smpar_mc::smpar_mc()
{
    pi = M_PI;
    scale_mu = 4.8;
    m_d = 0.00467, m_c = 1.27, m_s = 0.093, m_b = 4.18;

    m_Bd = 5.27965, m_Bs = 5.36688, m_K = 0.497614, m_Kst = 0.896, m_phi = 1.019461;

    tau_Bd = 1.519*1.52*pow(10.0,12.0), tau_Bs = 1.527*1.52*pow(10.0,12.0), D_Gamma_dbar = 0.001, D_Gamma_sbar = 0.124, f_Bd = 0.1905, f_Bs = 0.2303;

    G_F = 1.1663787*pow(10.0,-5.0), alphae = 1.0/127.944;

    C_1 = -0.257, C_2 = 1.009, C_3 = -0.005, C_4 = -0.078, C_5 = 0.000, C_6 = 0.001, C_7effRe = -0.304, C_7effIm = 0.0, C_8eff = -0.167,
    C_9Re = 4.211, C_9Im = 0.0, C_10Re = -4.103, C_10Im = 0.0;

    C_7RHRe = -0.006, C_7RHIm = 0.0, C_9RHRe = 0.0, C_9RHIm = 0.0, C_10RHRe = 0.0, C_10RHIm = 0.0, C_SRe = 0.0, C_SIm = 0.0, C_SRHRe = 0.0, C_SRHIm = 0.0,
    C_PRe = 0.0, C_PIm = 0.0, C_PRHRe = 0.0, C_PRHIm = 0.0, C_TRe = 0.0, C_TIm = 0.0, C_T5Re = 0.0, C_T5Im = 0.0;
}
double smpar_mc::cen_smpar[] = {/*m_e*/0.5109989461*pow(10.0,-3.0), /*m_mu*/105.6583745*pow(10.0,-3.0), /*m_tau*/1.77686,
                                /*m_Bd*/5.27965, /*m_Bs*/5.36688, /*m_K*/0.497614, /*m_Kst*/0.896, /*m_phi*/1.019461};

double smpar_mc::cen_smpar[] = {/*m_e*/0.0000000031*pow(10.0,-3.0), /*m_mu*/000.0000024*pow(10.0,-3.0), /*m_tau*/0.00012,
                                /*m_Bd*/0.00012, /*m_Bs*/5.36688, /*m_K*/0.497614, /*m_Kst*/0.896, /*m_phi*/1.019461};
//Quark_Mass
double smpar_mc::mu(){return scale_mu;}
double smpar_mc::md(){return m_d;}
double smpar_mc::mc(){return m_c;}
double smpar_mc::ms(){return m_s;}
double smpar_mc::mb(){return m_b;}
//Lepton_Mass
double smpar_mc::me(){return m_e;}
double smpar_mc::mmu(){return m_mu;}
double smpar_mc::mtau(){return m_tau;}
//Meson_Mass
double smpar_mc::mBd(){return m_Bd;}
double smpar_mc::mBs(){return m_Bs;}
double smpar_mc::mK(){return m_K;}
double smpar_mc::mKst(){return m_Kst;}
double smpar_mc::mphi(){return m_phi;}
//Width_Lifetime
double smpar_mc::tauBd(){return tau_Bd;}
double smpar_mc::tauBs(){return tau_Bs;}
double smpar_mc::DGamma_dbar(){return D_Gamma_dbar;}
double smpar_mc::DGamma_sbar(){return D_Gamma_sbar;}
double smpar_mc::fBd(){return f_Bd;}
double smpar_mc::fBs(){return f_Bs;}
//Coupling_Constants
double smpar_mc::GF(){return G_F;}
double smpar_mc::alpha_e(){return alphae;}
//SM_Wilson_Coefficients
double smpar_mc::C1(){return C_1;}
double smpar_mc::C2(){return C_2;}
double smpar_mc::C3(){return C_3;}
double smpar_mc::C4(){return C_4;}
double smpar_mc::C5(){return C_5;}
double smpar_mc::C6(){return C_6;}
double smpar_mc::C7effRe(){return C_7effRe;}
double smpar_mc::C7effIm(){return C_7effIm;}
double smpar_mc::C8eff(){return C_8eff;}
double smpar_mc::C9Re(){return C_9Re;}
double smpar_mc::C9Im(){return C_9Im;}
double smpar_mc::C10Re(){return C_10Re;}
double smpar_mc::C10Im(){return C_10Im;}
//NP_Wilson_Coefficients
double smpar_mc::C7RHRe(){return C_7RHRe;}
double smpar_mc::C7RHIm(){return C_7RHIm;}
double smpar_mc::C9RHRe(){return C_9RHRe;}
double smpar_mc::C9RHIm(){return C_9RHIm;}
double smpar_mc::C10RHRe(){return C_10RHRe;}
double smpar_mc::C10RHIm(){return C_10RHIm;}
double smpar_mc::CSRe(){return C_SRe;}
double smpar_mc::CSIm(){return C_SIm;}
double smpar_mc::CSRHRe(){return C_SRHRe;}
double smpar_mc::CSRHIm(){return C_SRHIm;}
double smpar_mc::CPRe(){return C_PRe;}
double smpar_mc::CPIm(){return C_PIm;}
double smpar_mc::CPRHRe(){return C_PRHRe;}
double smpar_mc::CPRHIm(){return C_PRHIm;}
double smpar_mc::CTRe(){return C_TRe;}
double smpar_mc::CTIm(){return C_TIm;}
double smpar_mc::CT5Re(){return C_T5Re;}
double smpar_mc::CT5Im(){return C_T5Im;}


