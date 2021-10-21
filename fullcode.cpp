#include <iostream>
#include <cmath>
#include <string>
using namespace std;
/////smpar.h
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
    /*double tpd();
    double tps();
    double tmd();
    double tms();
    double tzd(),tzs();*/
    
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
/////smpar.cc
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

/////ffpar.h
class ffpar : public smpar {
public:
    ffpar();
    double tpd,tps,tmd,tms,tzd,tzs;
    double zd(double qsq),zs(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(),Va1(),Va2(),A0a0(),A0a1(),A0a2(),A1a0(),A1a1(),A1a2(),A12a0(),A12a1(),A12a2(),T1a0(),T1a1(),T1a2(),T2a0(),T2a1(),T2a2(),T23a0(),T23a1(),T23a2();
private:
    //Bs->phill
    double m_PSbs,m_Vbs,m_Abs;
    double  V_a0,V_a1,V_a2,A0_a0,A0_a1,A0_a2,A1_a0,A1_a1,A1_a2,A12_a0,A12_a1,A12_a2,T1_a0,T1_a1,T1_a2,T2_a0,T2_a1,T2_a2,T23_a0,T23_a1,T23_a2;
};

//////ffpar.cc
ffpar::ffpar(){
    tpd=pow(mBd()+mKst(),2.0);tps=pow(mBs()+mKst(),2.0);
    tmd=pow(mBd()-mKst(),2.0);tms=pow(mBs()-mKst(),2.0);
    tzd=tpd*(1.0-sqrt(1.0-tmd/tpd));
    tzs=tps*(1.0-sqrt(1.0-tms/tps));
    
    m_PSbs=5.366,m_Vbs=5.415,m_Abs=5.829;
    
    V_a0=0.376313,V_a1=-1.16597,V_a2=2.42443,
    A0_a0=0.369196,A0_a1=-1.36584,A0_a2=0.128191,
    A1_a0=0.29725,A1_a1=0.392378,A1_a2=1.18916,
    A12_a0=0.265375,A12_a1=0.533638,A12_a2=0.483166,
    T1_a0=0.312055,T1_a1=-1.00893,T1_a2=1.5272,
    T2_a0=0.312055,T2_a1=0.496846,T2_a2=1.61431,
    T23_a0=0.667412,T23_a1=1.31812,T23_a2=3.82334;}

double ffpar::zd(double qsq){
    return (sqrt(tpd-qsq)-sqrt(tpd-tzd))/(sqrt(tpd-qsq) +sqrt(tpd-tzd));}
double ffpar::zs(double qsq){
    return (sqrt(tps-qsq)-sqrt(tps-tzs))/(sqrt(tps-qsq) +sqrt(tps-tzs));}

double ffpar::mPSbs(){return m_PSbs;}
double ffpar::mVbs(){return m_Vbs;}
double ffpar::mAbs(){return m_Abs;}
double ffpar::Va0(){return V_a0;}
double ffpar::Va1(){return V_a1;}
double ffpar::Va2(){return V_a2;}
double ffpar::A0a0(){return A0_a0;}
double ffpar::A0a1(){return A0_a1;}
double ffpar::A0a2(){return A0_a2;}
double ffpar::A1a0(){return A1_a0;}
double ffpar::A1a1(){return A1_a1;}
double ffpar::A1a2(){return A1_a2;}
double ffpar::A12a0(){return A12_a0;}
double ffpar::A12a1(){return A12_a1;}
double ffpar::A12a2(){return A12_a2;}
double ffpar::T1a0(){return T1_a0;}
double ffpar::T1a1(){return T1_a1;}
double ffpar::T1a2(){return T1_a2;}
double ffpar::T2a0(){return T2_a0;}
double ffpar::T2a1(){return T2_a1;}
double ffpar::T2a2(){return T2_a2;}
double ffpar::T23a0(){return T23_a0;}
double ffpar::T23a1(){return T23_a1;}
double ffpar::T23a2(){return T23_a2;}


/////fffun.h
class fffun : public ffpar {
public:
    //Bs->phill
    double V(double qsq),A0(double qsq),A1(double qsq),A12(double qsq),T1(double qsq),T2(double qsq),T23(double qsq);
};

/////fffun.cc
double fffun::V(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( Va0()*pow(zs(qsq)-zs(0.0),0.0) + Va1()*pow(zs(qsq)-zs(0.0),1.0) + Va2()*pow(zs(qsq)-zs(0.0),2.0) );}

double fffun::A0(double qsq){
    return 1.0/(1.0-qsq/pow(mPSbs(),2.0))*( A0a0()*pow(zs(qsq)-zs(0.0),0.0) + A0a1()*pow(zs(qsq)-zs(0.0),1.0) + A0a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double fffun::A1(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A1a0()*pow(zs(qsq)-zs(0.0),0.0) + A1a1()*pow(zs(qsq)-zs(0.0),1.0) + A1a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double fffun::A12(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A12a0()*pow(zs(qsq)-zs(0.0),0.0) + A12a1()*pow(zs(qsq)-zs(0.0),1.0) + A12a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double fffun::T1(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( T1a0()*pow(zs(qsq)-zs(0.0),0.0) + T1a1()*pow(zs(qsq)-zs(0.0),1.0) + T1a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double fffun::T2(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T2a0()*pow(zs(qsq)-zs(0.0),0.0) + T2a1()*pow(zs(qsq)-zs(0.0),1.0) + T2a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double fffun::T23(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T23a0()*pow(zs(qsq)-zs(0.0),0.0) + T23a1()*pow(zs(qsq)-zs(0.0),1.0) + T23a2()*pow(zs(qsq)-zs(0.0),2.0) );}



//////amp.h
class amp : public fffun {
public:
    double hRe(double qsq, double mq);
    double hIm(double qsq, double mq);
    double hReZ(double qsq);
    double hImZ(double qsq);
    double YRe(double qsq);
    double YIm(double qsq);
    double C9effRe(double qsq);
    double C9effIm(double qsq);
    double betal(double qsq, double ml);
    double lambda(double qsq);
    double nf(double qsq, double ml);
    
    double ApLRe(double qsq, double ml);
    double ApLIm(double qsq, double ml);
    double ApRRe(double qsq, double ml);
    double ApRIm(double qsq, double ml);
    double AaLRe(double qsq, double ml);
    double AaLIm(double qsq, double ml);
    double AaRRe(double qsq, double ml);
    double AaRIm(double qsq, double ml);
    double AzLRe(double qsq, double ml);
    double AzLIm(double qsq, double ml);
    double AzRRe(double qsq, double ml);
    double AzRIm(double qsq, double ml);
    double AtRe(double qsq, double ml);
    double AtIm(double qsq, double ml);
    double ASRe(double qsq, double ml);
    double ASIm(double qsq, double ml);
};

//////amp.cc
double amp::hRe(double qsq, double mq){
    if (qsq < 4.0*pow(mq,2.0)){
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(4.0*pow(mq,2.0)/qsq-1.0)* atan(1.0/sqrt(4.0*pow(mq,2.0)/qsq-1.0)); }
    else {
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq)* log((1.0+sqrt(1.0-4.0*pow(mq,2.0)/qsq))/sqrt(4.0*pow(mq,2.0)/qsq)); } }

double amp::hIm(double qsq, double mq){
    if (qsq >= 4.0*pow(mq,2.0)){
        return 2.0*pi/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq); }
    else {
        return 0.0; } }

double amp::hReZ(double qsq){
    return 8.0/27.0 - 4.0/9.0*log(qsq/pow(mu(),2.0)); }

double amp::hImZ(double qsq){
    return 4.0/9.0*pi; }

double amp::YRe(double qsq){
    return hRe(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hRe(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hReZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()) + 4.0/3.0*C3() + 64.0/9.0*C5() + 64.0/27.0*C6(); }

double amp::YIm(double qsq){
    return hIm(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hIm(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hImZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()); }

double amp::C9effRe(double qsq){return 4.211 + YRe(qsq);}
double amp::C9effIm(double qsq){return YIm(qsq);}

double amp::betal(double qsq, double ml){return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double amp::lambda(double qsq){
    return pow(qsq,2.0)+pow(mBs(),4.0)+pow(mKst(),4.0)-2.0*(pow(mBs()*mKst(),2.0)+qsq*pow(mBs(),2.0)+qsq*pow(mKst(),2.0));}

double amp:: nf(double qsq, double ml){
    return absVtbVtsStr()*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBs(),3.0)) );}

double amp::ApLRe(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}
double amp::ApLIm(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}

double amp::ApRRe(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}
double amp::ApRIm(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}

double amp::AaLRe(double qsq, double ml){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq));}
double amp::AaLIm(double qsq, double ml){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq));}

double amp::AaRRe(double qsq, double ml){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq));}
double amp::AaRIm(double qsq, double ml){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq));}

double amp::AzLRe(double qsq, double ml){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq) + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) );}
double amp::AzLIm(double qsq, double ml){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq) + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) );}

double amp::AzRRe(double qsq, double ml){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq) + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) );}
double amp::AzRIm(double qsq, double ml){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq) + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) );}

double amp::AtRe(double qsq, double ml){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq);}
double amp::AtIm(double qsq, double ml){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq);}

double amp::ASRe(double qsq, double ml){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq);}
double amp::ASIm(double qsq, double ml){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq);}



/////obs.h
class obs : public amp {
public:
    double J1s(double qsq, double ml);
    double J1c(double qsq, double ml);
    double J2s(double qsq, double ml);
    double J2c(double qsq, double ml);
    double J3(double qsq, double ml);
    double J4(double qsq, double ml);
    double J5(double qsq, double ml);
    double J6s(double qsq, double ml);
    double J6c(double qsq, double ml);
    double J7(double qsq, double ml);
    double J8(double qsq, double ml);
    double J9(double qsq, double ml);
    double diffWidth(double qsq, double ml);
    double FL(double qsq, double ml);
    double AFB(double qsq, double ml);
};

//////obs.cc
double obs::J1s(double qsq, double ml){
    return (4.0*pow(ml,2.0)*(AaLIm(qsq,ml)*AaRIm(qsq,ml) + AaLRe(qsq,ml)*AaRRe(qsq,ml) + ApLIm(qsq,ml)*ApRIm(qsq,ml) +
                             ApLRe(qsq,ml)*ApRRe(qsq,ml)))/qsq +
    ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
      pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double obs::J1c(double qsq, double ml){
    return pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtIm(qsq,ml),2.0) + pow(AtRe(qsq,ml),2.0) +
                      2.0*(AzLIm(qsq,ml)*AzRIm(qsq,ml) + AzLRe(qsq,ml)*AzRRe(qsq,ml))))/qsq +
    (pow(ASIm(qsq,ml),2.0) + pow(ASRe(qsq,ml),2.0))*pow(betal(qsq,ml),2.0) ;}

double obs::J2s(double qsq, double ml){
    return ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double obs::J2c(double qsq, double ml){
    return -((pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double obs::J3(double qsq, double ml){
    return ((-pow(AaLIm(qsq,ml),2.0) - pow(AaLRe(qsq,ml),2.0) - pow(AaRIm(qsq,ml),2.0) - pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double obs::J4(double qsq, double ml){
    return ((AaLIm(qsq,ml)*AzLIm(qsq,ml) + AaLRe(qsq,ml)*AzLRe(qsq,ml) + AaRIm(qsq,ml)*AzRIm(qsq,ml) +
             AaRRe(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double obs::J5(double qsq, double ml){
    return sqrt(2.0)*(-((ml*(AaLIm(qsq,ml)*ASIm(qsq,ml) + AaRIm(qsq,ml)*ASIm(qsq,ml) + AaLRe(qsq,ml)*ASRe(qsq,ml) +
                             AaRRe(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq)) + ApLIm(qsq,ml)*AzLIm(qsq,ml) +
                      ApLRe(qsq,ml)*AzLRe(qsq,ml) - ApRIm(qsq,ml)*AzRIm(qsq,ml) - ApRRe(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ;}

double obs::J6s(double qsq, double ml){
    return 2.0*(AaLIm(qsq,ml)*ApLIm(qsq,ml) + AaLRe(qsq,ml)*ApLRe(qsq,ml) - AaRIm(qsq,ml)*ApRIm(qsq,ml) -
                AaRRe(qsq,ml)*ApRRe(qsq,ml))*betal(qsq,ml) ;}

double obs::J6c(double qsq, double ml){
    return (4.0*ml*(ASIm(qsq,ml)*AzLIm(qsq,ml) + ASRe(qsq,ml)*AzLRe(qsq,ml) + ASIm(qsq,ml)*AzRIm(qsq,ml) +
                    ASRe(qsq,ml)*AzRRe(qsq,ml))*betal(qsq,ml))/sqrt(qsq) ; }

double obs::J7(double qsq, double ml){
    return sqrt(2.0)*((ml*(-(ApLRe(qsq,ml)*ASIm(qsq,ml)) - ApRRe(qsq,ml)*ASIm(qsq,ml) + ApLIm(qsq,ml)*ASRe(qsq,ml) +
                           ApRIm(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq) + AaLRe(qsq,ml)*AzLIm(qsq,ml) -
                      AaLIm(qsq,ml)*AzLRe(qsq,ml) - AaRRe(qsq,ml)*AzRIm(qsq,ml) + AaRIm(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ;}

double obs::J8(double qsq, double ml){
    return ((ApLRe(qsq,ml)*AzLIm(qsq,ml) - ApLIm(qsq,ml)*AzLRe(qsq,ml) + ApRRe(qsq,ml)*AzRIm(qsq,ml) -
             ApRIm(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double obs::J9(double qsq, double ml){
    return (AaLRe(qsq,ml)*ApLIm(qsq,ml) - AaLIm(qsq,ml)*ApLRe(qsq,ml) + AaRRe(qsq,ml)*ApRIm(qsq,ml) -
            AaRIm(qsq,ml)*ApRRe(qsq,ml))*pow(betal(qsq,ml),2.0) ;}

double obs::diffWidth(double qsq, double ml){
    return (3.0*(J1c(qsq,ml) + 2.0*J1s(qsq,ml) + (-J2c(qsq,ml) - 2.0*J2s(qsq,ml))/3.0))/4.0;}

double obs::FL(double qsq, double ml){
    return (pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))/
    (pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
     pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0) +
     pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0));}

double obs::AFB(double qsq, double ml){
    return (3.0*(J6c(qsq,ml) + 2.0*J6s(qsq,ml)))/(8.0*diffWidth(qsq,ml));}

/////main.cc

int main(){
    smpar p;
    obs o;
    double qsq;
    cout << "Type a value for qsq:";
    cin >> qsq;
    cout << "Observable values at this point are:"
    << o.diffWidth(qsq,p.mmu()) << "," << o.FL(qsq,p.mmu()) << "," << o.AFB(qsq,p.mmu()) << endl;
    return 0;
}
                                                
                                                



















