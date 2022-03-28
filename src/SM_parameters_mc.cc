#include "SM_parameters_mc.h"

smpar_mc::smpar_mc()
{
    pi = M_PI;
    scale_mu = 4.8;//the renormalization scale
    G_F = 1.1663787*pow(10.0,-5.0), alphae = 1.0/128.541;//errors on the Fermi constant and fine-structure constant are too small to consider
}

double smpar_mc::cen_smpar[] = {/*m_d*/0.00467, /*m_c*/1.27, /*m_s*/0.093, /*m_b*/4.18,
                                /*m_e*/0.5109989461*pow(10.0,-3.0), /*m_mu*/105.6583745*pow(10.0,-3.0), /*m_tau*/1.77686,
                                /*m_Bd*/5.27965, /*m_Bs*/5.36688, /*m_K*/0.497611, /*m_Kst*/0.89555, /*m_phi*/1.019461,
                                /*tau_Bd*/1.519*1.52*pow(10.0,12.0), /*tau_Bs*/1.527*1.52*pow(10.0,12.0),
                                /*D_Gamma_dbar*/0.001, /*D_Gamma_sbar*/0.124, /*f_Bd*/0.1905, /*f_Bs*/0.2303};

double smpar_mc::unc_smpar[] = {/*m_d*/0.00048, /*m_c*/0.02, /*m_s*/0.011, /*m_b*/0.03,
                                /*m_e*/0.0000000031*pow(10.0,-3.0), /*m_mu*/000.0000024*pow(10.0,-3.0), /*m_tau*/0.00012,
                                /*m_Bd*/0.00012, /*m_Bs*/0.00014, /*m_K*/0.000013, /*m_Kst*/0.00020, /*m_phi*/0.000016,
                                /*tau_Bd*/0.004*1.52*pow(10.0,12.0), /*tau_Bs*/0.011*1.52*pow(10.0,12.0),
                                /*D_Gamma_dbar*/0.010, /*D_Gamma_sbar*/0.008, /*f_Bd*/0.0013, /*f_Bs*/0.0013};

double smpar_mc::chd_smpar[][MXdm] =  /*m_d*/     {{0.00048, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_c*/     {0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_s*/     {0, 0, 0.011, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_b*/     {0, 0, 0, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_e*/     {0, 0, 0, 0, 0.0000000031*pow(10.0,-3.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_mu*/    {0, 0, 0, 0, 0, 0.0000024*pow(10.0,-3.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_tau*/   {0, 0, 0, 0, 0, 0, 0.00012, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_Bd*/    {0, 0, 0, 0, 0, 0, 0, 0.00012, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_Bs*/    {0, 0, 0, 0, 0, 0, 0, 0, 0.00014, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_K*/     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000013, 0, 0, 0, 0, 0, 0, 0, 0},
                                /*m_Kst*/   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00020, 0, 0, 0, 0, 0, 0, 0},
                                /*m_phi*/   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000016, 0, 0, 0, 0, 0, 0},
                                /*tau_Bd*/  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.004*1.52*pow(10.0,12.0), 0, 0, 0, 0, 0},
                                /*tau_Bs*/  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.011*1.52*pow(10.0,12.0), 0, 0, 0, 0},
                        /*D_Gamma_dbar*/    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.010, 0, 0, 0},
                        /*D_Gamma_sbar*/    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0, 0},
                                /*f_Bd*/    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0013, 0},
                                /*f_Bs*/    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0013}};

//double smpar_mc::cen_wilsonSM[] = {/*C_1*/-0.257, /*C_2*/1.009, /*C_3*/-0.005, /*C_4*/-0.078, /*C_5*/0.000, /*C_6*/0.001, /*C_7effRe*/-0.304, /*C_7effIm*/0.0,
//                                    /*C_8eff*/-0.167, /*C_9Re*/4.211, /*C_9Im*/0.0, /*C_10Re*/-4.103, /*C_10Im*/0.0};

//double smpar_mc::unc_wilsonSM[] = {/*C_1*/0.0, /*C_2*/0.0, /*C_3*/0.0, /*C_4*/0.0, /*C_5*/0.0, /*C_6*/0.0, /*C_7effRe*/0.0, /*C_7effIm*/0.0,
//                                    /*C_8eff*/0.0, /*C_9Re*/0.0, /*C_9Im*/0.0, /*C_10Re*/0.0, /*C_10Im*/0.0};

double smpar_mc::chd_wilsonSM[][MXdm] = /*C_1*/   {{0,0,0,0,0,0,0,0,0,0,0,0,0},
                                  /*C_2*/   {0,0,0,0,0,0,0,0,0,0,0,0,0},
                                  /*C_3*/   {0,0,0,0,0,0,0,0,0,0,0,0,0},
                                  /*C_4*/   {0,0,0,0,0,0,0,0,0,0,0,0,0},
                                  /*C_5*/   {0,0,0,0,0,0,0,0,0,0,0,0,0},
                                  /*C_6*/   {0,0,0,0,0,0,0,0,0,0,0,0,0},
                            /*C_7effRe*/    {0,0,0,0,0,0,0,0,0,0,0,0,0},
                            /*C_7effIm*/    {0,0,0,0,0,0,0,0,0,0,0,0,0},
                                /*C_8eff*/  {0,0,0,0,0,0,0,0,0,0,0,0,0},
                                /*C_9Re*/   {0,0,0,0,0,0,0,0,0,0,0,0,0},
                                /*C_9Im*/   {0,0,0,0,0,0,0,0,0,0,0,0,0},
                                /*C_10Re*/  {0,0,0,0,0,0,0,0,0,0,0,0,0},
                                /*C_10Im*/  {0,0,0,0,0,0,0,0,0,0,0,0,0}};

//double smpar_mc::cen_wilsonNP[] = {/*C_7RHRe*/-0.006, /*C_7RHIm*/0.0, /*C_9RHRe*/0.0, /*C_9RHIm*/0.0, /*C_10RHRe*/0.0, /*C_10RHIm*/0.0, /*C_SRe*/0.0,
//                                    /*C_SIm*/0.0, /*C_SRHRe*/0.0, /*C_SRHIm*/0.0, /*C_PRe*/0.0, /*C_PIm*/0.0, /*C_PRHRe*/0.0, /*C_PRHIm*/0.0,
//                                    /*C_TRe*/0.0, /*C_TIm*/0.0, /*C_T5Re*/0.0, /*C_T5Im*/0.0};

//double smpar_mc::unc_wilsonNP[] = {/*C_7RHRe*/0.0, /*C_7RHIm*/0.0, /*C_9RHRe*/0.0, /*C_9RHIm*/0.0, /*C_10RHRe*/0.0, /*C_10RHIm*/0.0, /*C_SRe*/0.0,
//                                    /*C_SIm*/0.0, /*C_SRHRe*/0.0, /*C_SRHIm*/0.0, /*C_PRe*/0.0, /*C_PIm*/0.0, /*C_PRHRe*/0.0, /*C_PRHIm*/0.0,
//                                    /*C_TRe*/0.0, /*C_TIm*/0.0, /*C_T5Re*/0.0, /*C_T5Im*/0.0};

double smpar_mc::chd_wilsonNP[][MXdm] =   /*C_7RHRe*/     {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_7RHIm*/     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_9RHRe*/     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_9RHIm*/     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_10RHRe*/    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_10RHIm*/    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_SRe*/       {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_SIm*/       {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_SRHRe*/     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_SRHIm*/     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_PRe*/       {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_PIm*/       {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_PRHRe*/     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_PRHIm*/     {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_TRe*/       {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_TIm*/       {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_T5Re*/      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    /*C_T5Im*/      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};


//Renormalization Scale
double smpar_mc::mu(){return scale_mu;}
//Coupling_Constants
double smpar_mc::GF(){return G_F;}
double smpar_mc::alpha_e(){return alphae;}
//Quark_Mass
double smpar_mc::md(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 0);}
double smpar_mc::mc(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 1);}
double smpar_mc::ms(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 2);}
double smpar_mc::mb(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 3);}
//Lepton_Mass
double smpar_mc::me(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 4);}
double smpar_mc::mmu(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 5);}
double smpar_mc::mtau(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 6);}
//Meson_Mass
double smpar_mc::mBd(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 7);}
double smpar_mc::mBs(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 8);}
double smpar_mc::mK(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 9);}
double smpar_mc::mKst(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 10);}
double smpar_mc::mphi(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 11);}
//Width_Lifetime
double smpar_mc::tauBd(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 12);}
double smpar_mc::tauBs(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 13);}
double smpar_mc::DGamma_dbar(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 14);}
double smpar_mc::DGamma_sbar(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 15);}
double smpar_mc::fBd(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 16);}
double smpar_mc::fBs(double unv[]){return mnd_cov_smpar(cen_smpar, chd_smpar, unv, 17);}
//SM_Wilson_Coefficients
double smpar_mc::C1(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 0);}
double smpar_mc::C2(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 1);}
double smpar_mc::C3(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 2);}
double smpar_mc::C4(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 3);}
double smpar_mc::C5(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 4);}
double smpar_mc::C6(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 5);}
double smpar_mc::C7effRe(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 6);}
double smpar_mc::C7effIm(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 7);}
double smpar_mc::C8eff(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 8);}
double smpar_mc::C9Re(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 9);}
double smpar_mc::C9Im(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 10);}
double smpar_mc::C10Re(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 11);}
double smpar_mc::C10Im(double smwc[], double unv[]){return mnd_cov_wcSM(smwc, chd_wilsonSM, unv, 12);}
//NP_Wilson_Coefficients
double smpar_mc::C7RHRe(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 0);}
double smpar_mc::C7RHIm(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 1);}
double smpar_mc::C9RHRe(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 2);}
double smpar_mc::C9RHIm(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 3);}
double smpar_mc::C10RHRe(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 4);}
double smpar_mc::C10RHIm(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 5);}
double smpar_mc::CSRe(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 6);}
double smpar_mc::CSIm(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 7);}
double smpar_mc::CSRHRe(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 8);}
double smpar_mc::CSRHIm(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 9);}
double smpar_mc::CPRe(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 10);}
double smpar_mc::CPIm(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 11);}
double smpar_mc::CPRHRe(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 12);}
double smpar_mc::CPRHIm(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 13);}
double smpar_mc::CTRe(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 14);}
double smpar_mc::CTIm(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 15);}
double smpar_mc::CT5Re(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 16);}
double smpar_mc::CT5Im(double npwc[], double unv[]){return mnd_cov_wcNP(npwc, chd_wilsonNP, unv, 17);}


/*References/////////////////////////////////////////////////
masses and lifetimes are taken from pdgLive interactive listings: https://pdglive.lbl.gov/Viewer.action
decay constants(i.e. fBd, fBs): Eq-84.20 of https://pdg.lbl.gov/2019/reviews/rpp2019-rev-pseudoscalar-meson-decay-cons.pdf
Fermi constant and fine-structure constant: https://pdg.lbl.gov/2019/reviews/rpp2019-rev-phys-constants.pdf
*/




