#ifndef smpar_h
#define smpar_h

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <string>
#include <bits/stdc++.h>
#define mxdm 5
#define MXdm 21
using namespace std;

class smpar {
public:
    smpar ();
    double pi;
    double mu();
    double md(),mc(),ms(),mb();
    double me(),mmu(),mtau();
    double mBd(),mBs(),mK(),mKst();
    double tauBd(),tauBs(),DGamma_dbar(),DGamma_sbar(),fBd(),fBs();
    double absVtbVtdStr(),absVtbVtsStr(),GF(),alpha_e();
    double C1(),C2(),C3(),C4(),C5(),C6(),C7effRe(),C7effIm(),C8eff(),C9Re(),C9Im(),C10Re(),C10Im();
    double C7RHRe(),C7RHIm(),C9RHRe(),C9RHIm(),C10RHRe(),C10RHIm(),CSRe(),CSIm(),CSRHRe(),CSRHIm(),CPRe(),CPIm(),CPRHRe(),CPRHIm(),CTRe(),CTIm(),CT5Re(),CT5Im();

private:
    double scale_mu;
    double m_d,m_c,m_s,m_b;
    double m_e,m_mu,m_tau;
    double m_Bd,m_Bs,m_K,m_Kst;
    double tau_Bd,tau_Bs,D_Gamma_dbar,D_Gamma_sbar,f_Bd,f_Bs;
    double abs_VtbVtdStr,abs_VtbVtsStr,G_F,alphae;
    double C_1,C_2,C_3,C_4,C_5,C_6,C_7effRe,C_7effIm,C_8eff,C_9Re,C_9Im,C_10Re,C_10Im;
    double C_7RHRe,C_7RHIm,C_9RHRe,C_9RHIm,C_10RHRe,C_10RHIm,C_SRe,C_SIm,C_SRHRe,C_SRHIm,C_PRe,C_PIm,C_PRHRe,C_PRHIm,C_TRe,C_TIm,C_T5Re,C_T5Im;
};

#endif
