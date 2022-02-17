#include "ffpar.h"

//Bs->Kstr,ll
BdtoKstrll_ffpar::BdtoKstrll_ffpar(){
    tpd=pow(mBd()+mKst(),2.0);
    tmd=pow(mBd()-mKst(),2.0);
    tzd=tpd*(1.0-sqrt(1.0-tmd/tpd));

    m_PSbs = 5.336, m_Vbs = 5.412, m_Abs = 5.829;

    A0_a0 = 0.369196, A0_a1 = -1.36584, A0_a2 = 0.128191,
    A1_a0 = 0.29725, A1_a1 = 0.392378, A1_a2 = 1.18916,
    A12_a0 = 0.265375, A12_a1 = 0.533638, A12_a2 = 0.483166,
    V_a0 = 0.376313, V_a1 = -1.16597, V_a2 = 2.42443,
    T1_a0 = 0.312055, T1_a1 = -1.00893, T1_a2 = 1.5272,
    T2_a0 = 0.312055, T2_a1 = 0.496846, T2_a2 = 1.61431,
    T23_a0 = 0.667412, T23_a1 = 1.31812, T23_a2 = 3.82334;}

double BdtoKstrll_ffpar::cov_A0[][mxdm] = {{0.000837635, 0.00471065, 0.00272131}, {0.00471065, 0.0659712, 0.204148}, {0.00272131, 0.204148, 2.67119}};
double BdtoKstrll_ffpar::cov_A1[][mxdm] = {{0.000694641, 0.00347885, 0.00415084}, {0.00347885, 0.0353042, 0.144041}, {0.00415084, 0.144041, 1.05126}};
double BdtoKstrll_ffpar::cov_A12[][mxdm] = {{0.000432775, 0.00165482, 0.000605487}, {0.00165482, 0.0165823, 0.059165}, {0.000605487, 0.059165, 0.430695}};
double BdtoKstrll_ffpar::cov_V[][mxdm] = {{0.00110852, 0.00658825, -0.0202371}, {0.00658825, 0.0682608, 0.0138459}, {-0.0202371, 0.0138459, 2.34402}};
double BdtoKstrll_ffpar::cov_T1[][mxdm] = {{0.000756032, 0.00372113, -0.0213639}, {0.00372113, 0.0359388, 0.00028928}, {-0.0213639, 0.00028928, 2.68844}};
double BdtoKstrll_ffpar::cov_T2[][mxdm] = {{0.000756032, 0.00316639, -0.00503577}, {0.00316639, 0.0276287, 0.0511444}, {-0.00503577, 0.0511444, 0.646067}};
double BdtoKstrll_ffpar::cov_T23[][mxdm] = {{0.0040111, 0.00460369, -0.0993374}, {0.00460369, 0.0493411, 0.140689}, {-0.0993374, 0.140689, 4.85487}};

double BdtoKstrll_ffpar::zd(double qsq){
    return (sqrt(tpd-qsq)-sqrt(tpd-tzd))/(sqrt(tpd-qsq) +sqrt(tpd-tzd));}

double BdtoKstrll_ffpar::mPSbs(){return m_PSbs;}
double BdtoKstrll_ffpar::mVbs(){return m_Vbs;}
double BdtoKstrll_ffpar::mAbs(){return m_Abs;}
double BdtoKstrll_ffpar::Va0(){return V_a0;}
double BdtoKstrll_ffpar::Va1(){return V_a1;}
double BdtoKstrll_ffpar::Va2(){return V_a2;}
double BdtoKstrll_ffpar::A0a0(){return A0_a0;}
double BdtoKstrll_ffpar::A0a1(){return A0_a1;}
double BdtoKstrll_ffpar::A0a2(){return A0_a2;}
double BdtoKstrll_ffpar::A1a0(){return A1_a0;}
double BdtoKstrll_ffpar::A1a1(){return A1_a1;}
double BdtoKstrll_ffpar::A1a2(){return A1_a2;}
double BdtoKstrll_ffpar::A12a0(){return A12_a0;}
double BdtoKstrll_ffpar::A12a1(){return A12_a1;}
double BdtoKstrll_ffpar::A12a2(){return A12_a2;}
double BdtoKstrll_ffpar::T1a0(){return T1_a0;}
double BdtoKstrll_ffpar::T1a1(){return T1_a1;}
double BdtoKstrll_ffpar::T1a2(){return T1_a2;}
double BdtoKstrll_ffpar::T2a0(){return T2_a0;}
double BdtoKstrll_ffpar::T2a1(){return T2_a1;}
double BdtoKstrll_ffpar::T2a2(){return T2_a2;}
double BdtoKstrll_ffpar::T23a0(){return T23_a0;}
double BdtoKstrll_ffpar::T23a1(){return T23_a1;}
double BdtoKstrll_ffpar::T23a2(){return T23_a2;}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll
Bstophill_ffpar::Bstophill_ffpar(){
    tps=pow(mBs()+mphi(),2.0);
    tms=pow(mBs()-mphi(),2.0);
    tzs=tps*(1.0-sqrt(1.0-tms/tps));

    m_PSbs=5.366,m_Vbs=5.415,m_Abs=5.829;

    V_a0=0.376313,V_a1=-1.16597,V_a2=2.42443,
    A0_a0=0.369196,A0_a1=-1.36584,A0_a2=0.128191,
    A1_a0=0.29725,A1_a1=0.392378,A1_a2=1.18916,
    A12_a0=0.265375,A12_a1=0.533638,A12_a2=0.483166,
    T1_a0=0.312055,T1_a1=-1.00893,T1_a2=1.5272,
    T2_a0=0.312055,T2_a1=0.496846,T2_a2=1.61431,
    T23_a0=0.667412,T23_a1=1.31812,T23_a2=3.82334;}

double Bstophill_ffpar::cov_A0[][mxdm] = {{0.00057847, 0.00392743, -0.0114296}, {0.00392743, 0.0564395, 0.0390906}, {-0.0114296, 0.0390906, 1.84713}};
double Bstophill_ffpar::cov_A1[][mxdm] = {{0.00011185, 0.000275488, -0.000881354}, {0.000275488, 0.0107771, 0.0643267}, {-0.000881354, 0.0643267, 0.624569}};
double Bstophill_ffpar::cov_A12[][mxdm] = {{0.000232401, 0.0014885, 0.00309781}, {0.0014885, 0.0158455, 0.048424}, {0.00309781, 0.048424, 0.229259}};
double Bstophill_ffpar::cov_V[][mxdm] = {{0.000199807, 0.000574973, -0.00300535}, {0.000574973, 0.0269538, 0.209736}, {-0.00300535, 0.209736, 2.98191}};
double Bstophill_ffpar::cov_T1[][mxdm] = {{0.000145656, 0.00045251, -0.00679762}, {0.00045251, 0.00697628, -0.0142247}, {-0.00679762, -0.0142247, 1.00649}};
double Bstophill_ffpar::cov_T2[][mxdm] = {{0.000145656, 0.000206198, -0.00363963}, {0.000206198, 0.00645713, 0.0200771}, {-0.00363963, 0.0200771, 0.370725}};
double Bstophill_ffpar::cov_T23[][mxdm] = {{0.00127784, 0.00917743, 0.0367908}, {0.00917743, 0.108954, 0.532838}, {0.0367908, 0.532838, 3.22296}};

double Bstophill_ffpar::zs(double qsq){
    return (sqrt(tps-qsq)-sqrt(tps-tzs))/(sqrt(tps-qsq) +sqrt(tps-tzs));}

double Bstophill_ffpar::mPSbs(){return m_PSbs;}
double Bstophill_ffpar::mVbs(){return m_Vbs;}
double Bstophill_ffpar::mAbs(){return m_Abs;}
double Bstophill_ffpar::Va0(){return V_a0;}
double Bstophill_ffpar::Va1(){return V_a1;}
double Bstophill_ffpar::Va2(){return V_a2;}
double Bstophill_ffpar::A0a0(){return A0_a0;}
double Bstophill_ffpar::A0a1(){return A0_a1;}
double Bstophill_ffpar::A0a2(){return A0_a2;}
double Bstophill_ffpar::A1a0(){return A1_a0;}
double Bstophill_ffpar::A1a1(){return A1_a1;}
double Bstophill_ffpar::A1a2(){return A1_a2;}
double Bstophill_ffpar::A12a0(){return A12_a0;}
double Bstophill_ffpar::A12a1(){return A12_a1;}
double Bstophill_ffpar::A12a2(){return A12_a2;}
double Bstophill_ffpar::T1a0(){return T1_a0;}
double Bstophill_ffpar::T1a1(){return T1_a1;}
double Bstophill_ffpar::T1a2(){return T1_a2;}
double Bstophill_ffpar::T2a0(){return T2_a0;}
double Bstophill_ffpar::T2a1(){return T2_a1;}
double Bstophill_ffpar::T2a2(){return T2_a2;}
double Bstophill_ffpar::T23a0(){return T23_a0;}
double Bstophill_ffpar::T23a1(){return T23_a1;}
double Bstophill_ffpar::T23a2(){return T23_a2;}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
BdtoKll_ffpar::BdtoKll_ffpar(){
    tp=pow(mBd()+mK(),2.0);
    tm=pow(mBd()-mK(),2.0);
    tz=tp*(1.0-sqrt(1.0-tm/tp));

    m_Vbs=5.412,m_Sbs=5.630;

    /*fz_a0=0.0,*/fz_a1=0.195117,fz_a2=-0.446126,
    fT_a0=0.299383,fT_a1=-0.773546,fT_a2=0.00955438,
    fp_a0=0.32909,fp_a1=-0.866947,fp_a2=0.00609567;}

double BdtoKll_ffpar::cov_fz[][mxdm] = {{0.0282565, 0.057892}, {0.057892, 0.167237}};
double BdtoKll_ffpar::cov_fT[][mxdm] = {{0.000669362, -0.00119402, -0.0131777}, {-0.00119402, 0.0225385, 0.0781834}, {-0.0131777, 0.0781834, 0.759645}};
double BdtoKll_ffpar::cov_fp[][mxdm] = {{0.000768854, -0.00063249, -0.0119012}, {-0.00063249, 0.0190837, 0.0795904}, {-0.0119012, 0.0795904, 0.56374}};

double BdtoKll_ffpar::z(double qsq){return (sqrt(tp-qsq)-sqrt(tp-tz))/(sqrt(tp-qsq)+sqrt(tp-tz));}
double BdtoKll_ffpar::mVbs(){return m_Vbs;}
double BdtoKll_ffpar::mSbs(){return m_Sbs;}
double BdtoKll_ffpar::fza1(){return fz_a1;}
double BdtoKll_ffpar::fza2(){return fz_a2;}
double BdtoKll_ffpar::fTa0(){return fT_a0;}
double BdtoKll_ffpar::fTa1(){return fT_a1;}
double BdtoKll_ffpar::fTa2(){return fT_a2;}
double BdtoKll_ffpar::fpa0(){return fp_a0;}
double BdtoKll_ffpar::fpa1(){return fp_a1;}
double BdtoKll_ffpar::fpa2(){return fp_a2;}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////










