#include <iostream>
#include <cmath>
#include <string>
using namespace std;

class smpar {
public:
    smpar ();
    double pi;
    double md(),mc(),ms(),mb();
    double me(),mmu(),mtau();
    double mBd(),mBs(),mK(),mKst();
    double tauBd(),tauBs(),DGamma_dbar(),DGamma_sbar(),fBd(),fBs();
    double absVtbVtdStr(),absVtbVtsStr(),GF(),alpha_e();
    double C1(),C2(),C3(),C4(),C5(),C6(),C7effRe(),C8eff(),C9Re(),C9Im(),C10Re(),C10Im();
    double C7RHRe(),C7RHIm(),C9RHRe(),C9RHIm(),C10RHRe(),C10RHIm(),CSRe(),CSIm(),CSRHRe(),CSRHIm(),CPRe(),CPIm(),CPRHRe(),CPRHIm(),CTRe(),CTIm(),CT5Re(),CT5Im();
    
private:
    double m_d,m_c,m_s,m_b;
    double m_e,m_mu,m_tau;
    double m_Bd,m_Bs,m_K,m_Kst;
    double tau_Bd,tau_Bs,D_Gamma_dbar,D_Gamma_sbar,f_Bd,f_Bs;
    double abs_VtbVtdStr,abs_VtbVtsStr,G_F,alphae;
    double C_1,C_2,C_3,C_4,C_5,C_6,C_7effRe,C_8eff,C_9Re,C_9Im,C_10Re,C_10Im;
    double C_7RHRe,C_7RHIm,C_9RHRe,C_9RHIm,C_10RHRe,C_10RHIm,C_SRe,C_SIm,C_SRHRe,C_SRHIm,C_PRe,C_PIm,C_PRHRe,C_PRHIm,C_TRe,C_TIm,C_T5Re,C_T5Im;
};
///////////////////////
smpar::smpar(){pi=M_PI;
    m_d=0.00467,m_c=1.27,m_s=0.093,m_b=4.18;
    m_e=0.511*pow(10.0,-3.0),m_mu=105.658*pow(10.0,-3.0),m_tau=1.77686;
    
    m_Bd=5.27965,m_Bs=5.36688,m_K=0.497614,m_Kst=0.896;
    
    tau_Bd=1.519*1.52*pow(10.0,12.0),tau_Bs=1.527*1.52*pow(10.0,12.0),D_Gamma_dbar=0.001,D_Gamma_sbar=0.124,f_Bd=0.1905,f_Bs=0.2303;
    
    abs_VtbVtdStr=0.00853871,abs_VtbVtsStr=0.0408726,G_F=1.16637*pow(10.0,-5.0),alphae=1.0/127.944;
    
    C_1=-0.257,C_2=1.009,C_3=-0.005,C_4=-0.078,C_5=0.000,C_6=0.001,C_7effRe=-0.304,C_8eff=-0.167,C_9Re=4.211,C_9Im=0.0,C_10Re=-4.103,C_10Im=0.0;
    
    C_7RHRe=0.0,C_7RHIm=0.0,C_9RHRe=0.0,C_9RHIm=0.0,C_10RHRe=0.0,C_10RHIm=0.0,C_SRe=0.0,C_SIm=0.0,C_SRHRe=0.0,C_SRHIm=0.0,C_PRe=0.0,C_PIm=0.0,C_PRHRe=0.0,C_PRHIm=0.0,C_TRe=0.0,C_TIm=0.0,C_T5Re=0.0,C_T5Im=0.0;
}
//Quark_Mass
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
///////////////////////
int main(){
    smpar p;
    cout << pow(p.mBd()-p.mKst()+p.C9Re()-p.CT5Im(),2.0) << endl;
    return 0;
}

