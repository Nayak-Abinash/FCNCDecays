#include "smpar.h"
#include "myfun.h"
#include "fferrpar.h"
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <string>
#include <iostream>

using namespace std;

//Bs->Kstr,ll
BdtoKstrll_fferrpar::BdtoKstrll_fferrpar(){
    tpd=pow(mBd()+mKst(),2.0);
    tmd=pow(mBd()-mKst(),2.0);
    tzd=tpd*(1.0-sqrt(1.0-tmd/tpd));

    m_PSbs=5.366,m_Vbs=5.415,m_Abs=5.829;

    V_a0=0.376313,V_a1=-1.16597,V_a2=2.42443,
    A0_a0=0.369196,A0_a1=-1.36584,A0_a2=0.128191,
    A1_a0=0.29725,A1_a1=0.392378,A1_a2=1.18916,
    A12_a0=0.265375,A12_a1=0.533638,A12_a2=0.483166,
    T1_a0=0.312055,T1_a1=-1.00893,T1_a2=1.5272,
    T2_a0=0.312055,T2_a1=0.496846,T2_a2=1.61431,
    T23_a0=0.667412,T23_a1=1.31812,T23_a2=3.82334;

    eV_a0=0.0376313,eV_a1=0.116597,eV_a2=0.242443,
    eA0_a0=0.0369196,eA0_a1=0.136584,eA0_a2=0.0128191,
    eA1_a0=0.029725,eA1_a1=0.0392378,eA1_a2=0.118916,
    eA12_a0=0.0265375,eA12_a1=0.0533638,eA12_a2=0.0483166,
    eT1_a0=0.0312055,eT1_a1=0.100893,eT1_a2=0.15272,
    eT2_a0=0.0312055,eT2_a1=0.0496846,eT2_a2=0.161431,
    eT23_a0=0.0667412,eT23_a1=0.131812,eT23_a2=0.382334;}

double BdtoKstrll_fferrpar::zd(double qsq){
    return (sqrt(tpd-qsq)-sqrt(tpd-tzd))/(sqrt(tpd-qsq) +sqrt(tpd-tzd));}

double BdtoKstrll_fferrpar::mPSbs(){return m_PSbs;}
double BdtoKstrll_fferrpar::mVbs(){return m_Vbs;}
double BdtoKstrll_fferrpar::mAbs(){return m_Abs;}

double BdtoKstrll_fferrpar::Va0(){return mnd(V_a0,eV_a0);}
double BdtoKstrll_fferrpar::Va1(){return mnd(V_a1,eV_a1);}
double BdtoKstrll_fferrpar::Va2(){return mnd(V_a2,eV_a2);}
double BdtoKstrll_fferrpar::A0a0(){return mnd(A0_a0,eA0_a0);}
double BdtoKstrll_fferrpar::A0a1(){return mnd(A0_a1,eA0_a1);}
double BdtoKstrll_fferrpar::A0a2(){return mnd(A0_a2,eA0_a2);}
double BdtoKstrll_fferrpar::A1a0(){return mnd(A1_a0,eA1_a0);}
double BdtoKstrll_fferrpar::A1a1(){return mnd(A1_a1,eA1_a1);}
double BdtoKstrll_fferrpar::A1a2(){return mnd(A1_a2,eA1_a2);}
double BdtoKstrll_fferrpar::A12a0(){return mnd(A12_a0,eA12_a0);}
double BdtoKstrll_fferrpar::A12a1(){return mnd(A12_a1,eA12_a1);}
double BdtoKstrll_fferrpar::A12a2(){return mnd(A12_a2,eA12_a2);}
double BdtoKstrll_fferrpar::T1a0(){return mnd(T1_a0,eT1_a0);}
double BdtoKstrll_fferrpar::T1a1(){return mnd(T1_a1,eT1_a1);}
double BdtoKstrll_fferrpar::T1a2(){return mnd(T1_a2,eT1_a2);}
double BdtoKstrll_fferrpar::T2a0(){return mnd(T2_a0,eT2_a0);}
double BdtoKstrll_fferrpar::T2a1(){return mnd(T2_a1,eT2_a1);}
double BdtoKstrll_fferrpar::T2a2(){return mnd(T2_a2,eT2_a2);}
double BdtoKstrll_fferrpar::T23a0(){return mnd(T23_a0,eT23_a0);}
double BdtoKstrll_fferrpar::T23a1(){return mnd(T23_a1,eT23_a1);}
double BdtoKstrll_fferrpar::T23a2(){return mnd(T23_a2,eT23_a2);}


//Bs->phi,ll
Bstophill_fferrpar::Bstophill_fferrpar(){
    tps=pow(mBs()+mKst(),2.0);
    tms=pow(mBs()-mKst(),2.0);
    tzs=tps*(1.0-sqrt(1.0-tms/tps));

    m_PSbs=5.366,m_Vbs=5.415,m_Abs=5.829;

    V_a0=0.376313,V_a1=-1.16597,V_a2=2.42443,
    A0_a0=0.369196,A0_a1=-1.36584,A0_a2=0.128191,
    A1_a0=0.29725,A1_a1=0.392378,A1_a2=1.18916,
    A12_a0=0.265375,A12_a1=0.533638,A12_a2=0.483166,
    T1_a0=0.312055,T1_a1=-1.00893,T1_a2=1.5272,
    T2_a0=0.312055,T2_a1=0.496846,T2_a2=1.61431,
    T23_a0=0.667412,T23_a1=1.31812,T23_a2=3.82334;

    eV_a0=0.0376313,eV_a1=0.116597,eV_a2=0.242443,
    eA0_a0=0.0369196,eA0_a1=0.136584,eA0_a2=0.0128191,
    eA1_a0=0.029725,eA1_a1=0.0392378,eA1_a2=0.118916,
    eA12_a0=0.0265375,eA12_a1=0.0533638,eA12_a2=0.0483166,
    eT1_a0=0.0312055,eT1_a1=0.100893,eT1_a2=0.15272,
    eT2_a0=0.0312055,eT2_a1=0.0496846,eT2_a2=0.161431,
    eT23_a0=0.0667412,eT23_a1=0.131812,eT23_a2=0.382334;}

double Bstophill_fferrpar::zs(double qsq){
    return (sqrt(tps-qsq)-sqrt(tps-tzs))/(sqrt(tps-qsq) +sqrt(tps-tzs));}

double Bstophill_fferrpar::mPSbs(){return m_PSbs;}
double Bstophill_fferrpar::mVbs(){return m_Vbs;}
double Bstophill_fferrpar::mAbs(){return m_Abs;}

double Bstophill_fferrpar::Va0(){return mnd(V_a0,eV_a0);}
double Bstophill_fferrpar::Va1(){return mnd(V_a1,eV_a1);}
double Bstophill_fferrpar::Va2(){return mnd(V_a2,eV_a2);}
double Bstophill_fferrpar::A0a0(){return mnd(A0_a0,eA0_a0);}
double Bstophill_fferrpar::A0a1(){return mnd(A0_a1,eA0_a1);}
double Bstophill_fferrpar::A0a2(){return mnd(A0_a2,eA0_a2);}
double Bstophill_fferrpar::A1a0(){return mnd(A1_a0,eA1_a0);}
double Bstophill_fferrpar::A1a1(){return mnd(A1_a1,eA1_a1);}
double Bstophill_fferrpar::A1a2(){return mnd(A1_a2,eA1_a2);}
double Bstophill_fferrpar::A12a0(){return mnd(A12_a0,eA12_a0);}
double Bstophill_fferrpar::A12a1(){return mnd(A12_a1,eA12_a1);}
double Bstophill_fferrpar::A12a2(){return mnd(A12_a2,eA12_a2);}
double Bstophill_fferrpar::T1a0(){return mnd(T1_a0,eT1_a0);}
double Bstophill_fferrpar::T1a1(){return mnd(T1_a1,eT1_a1);}
double Bstophill_fferrpar::T1a2(){return mnd(T1_a2,eT1_a2);}
double Bstophill_fferrpar::T2a0(){return mnd(T2_a0,eT2_a0);}
double Bstophill_fferrpar::T2a1(){return mnd(T2_a1,eT2_a1);}
double Bstophill_fferrpar::T2a2(){return mnd(T2_a2,eT2_a2);}
double Bstophill_fferrpar::T23a0(){return mnd(T23_a0,eT23_a0);}
double Bstophill_fferrpar::T23a1(){return mnd(T23_a1,eT23_a1);}
double Bstophill_fferrpar::T23a2(){return mnd(T23_a2,eT23_a2);}


//Bd->K,ll
BdtoKll_fferrpar::BdtoKll_fferrpar(){
    tp=pow(mBd()+mK(),2.0);
    tm=pow(mBd()-mK(),2.0);
    tz=tp*(1.0-sqrt(1.0-tm/tp));

    lcsr_m_Vbs=5.412,lcsr_m_Sbs=5.630;

    c_fp0=0.32909,c_fp1=-0.866947,c_fp2=0.00609567,c_fz0=0.0,c_fz1=0.195117,
    c_fz2=-0.446126,c_ft0=0.299383,c_ft1=-0.773546,c_ft2=0.00955438;

    ec_fp0=0.032909,ec_fp1=0.0866947,ec_fp2=0.000609567,ec_fz0=0.00,ec_fz1=0.0195117,
    ec_fz2=0.0446126,ec_ft0=0.0299383,ec_ft1=0.0773546,ec_ft2=0.000955438;

    lat_m_Vbs=5.4154,lat_m_Sbs=5.711;

    b0_p=0.466,b1_p=-0.885,b2_p=-0.213,b0_z=0.292,b1_z=0.281,b2_z=0.150,
    b0_t=0.460,b1_t=-1.089,b2_t=-1.114;

    eb0_p=0.0466,eb1_p=0.0885,eb2_p=0.0213,eb0_z=0.0292,eb1_z=0.0281,eb2_z=0.0150,
    eb0_t=0.0460,eb1_t=0.1089,eb2_t=0.1114;}

double BdtoKll_fferrpar::z(double qsq){return (sqrt(tp-qsq)-sqrt(tp-tz))/(sqrt(tp-qsq)+sqrt(tp-tz));}
double BdtoKll_fferrpar::lcsr_mVbs(){return lcsr_m_Vbs;}
double BdtoKll_fferrpar::lcsr_mSbs(){return lcsr_m_Sbs;}

double BdtoKll_fferrpar::cfp0(){return mnd(c_fp0,ec_fp0);}
double BdtoKll_fferrpar::cfp1(){return mnd(c_fp1,ec_fp1);}
double BdtoKll_fferrpar::cfp2(){return mnd(c_fp2,ec_fp2);}
double BdtoKll_fferrpar::cfz0(){return mnd(c_fz0,ec_fz0);}
double BdtoKll_fferrpar::cfz1(){return mnd(c_fz1,ec_fz1);}
double BdtoKll_fferrpar::cfz2(){return mnd(c_fz2,ec_fz2);}
double BdtoKll_fferrpar::cft0(){return mnd(c_ft0,ec_ft0);}
double BdtoKll_fferrpar::cft1(){return mnd(c_ft1,ec_ft1);}
double BdtoKll_fferrpar::cft2(){return mnd(c_ft2,ec_ft2);}

double BdtoKll_fferrpar::lat_mVbs(){return lat_m_Vbs;}
double BdtoKll_fferrpar::lat_mSbs(){return lat_m_Sbs;}

double BdtoKll_fferrpar::b0p(){return mnd(b0_p,eb0_p);}
double BdtoKll_fferrpar::b1p(){return mnd(b1_p,eb1_p);}
double BdtoKll_fferrpar::b2p(){return mnd(b2_p,eb2_p);}
double BdtoKll_fferrpar::b0z(){return mnd(b0_z,eb0_z);}
double BdtoKll_fferrpar::b1z(){return mnd(b1_z,eb1_z);}
double BdtoKll_fferrpar::b2z(){return mnd(b2_z,eb2_z);}
double BdtoKll_fferrpar::b0t(){return mnd(b0_t,eb0_t);}
double BdtoKll_fferrpar::b1t(){return mnd(b1_t,eb1_t);}
double BdtoKll_fferrpar::b2t(){return mnd(b2_t,eb2_t);}












