#include "smpar.h"
#include "ffpar.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

//Bs->Kstr,ll
BdtoKstrll_ffpar::BdtoKstrll_ffpar(){
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
    T23_a0=0.667412,T23_a1=1.31812,T23_a2=3.82334;}

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

//Bs->phi,ll
Bstophill_ffpar::Bstophill_ffpar(){
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
    T23_a0=0.667412,T23_a1=1.31812,T23_a2=3.82334;}

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

//Bd->K,ll
BdtoKll_ffpar::BdtoKll_ffpar(){
    tp=pow(mBd()+mK(),2.0);
    tm=pow(mBd()-mK(),2.0);
    tz=tp*(1.0-sqrt(1.0-tm/tp));
    
    lcsr_m_Vbs=5.412,lcsr_m_Sbs=5.630;
    
    c_fp0=0.32909,c_fp1=-0.866947,c_fp2=0.00609567,c_fz0=0.0,c_fz1=0.195117,
    c_fz2=-0.446126,c_ft0=0.299383,c_ft1=-0.773546,c_ft2=0.00955438;
    
    lat_m_Vbs=5.4154,lat_m_Sbs=5.711;
    
    b0_p=0.466,b1_p=-0.885,b2_p=-0.213,b0_z=0.292,b1_z=0.281,b2_z=0.150,
    b0_t=0.460,b1_t=-1.089,b2_t=-1.114;}

double BdtoKll_ffpar::z(double qsq){return (sqrt(tp-qsq)-sqrt(tp-tz))/(sqrt(tp-qsq)+sqrt(tp-tz));}
double BdtoKll_ffpar::lcsr_mVbs(){return lcsr_m_Vbs;}
double BdtoKll_ffpar::lcsr_mSbs(){return lcsr_m_Sbs;}
double BdtoKll_ffpar::cfp0(){return c_fp0;}
double BdtoKll_ffpar::cfp1(){return c_fp1;}
double BdtoKll_ffpar::cfp2(){return c_fp2;}
double BdtoKll_ffpar::cfz0(){return c_fz0;}
double BdtoKll_ffpar::cfz1(){return c_fz1;}
double BdtoKll_ffpar::cfz2(){return c_fz2;}
double BdtoKll_ffpar::cft0(){return c_ft0;}
double BdtoKll_ffpar::cft1(){return c_ft1;}
double BdtoKll_ffpar::cft2(){return c_ft2;}
double BdtoKll_ffpar::lat_mVbs(){return lat_m_Vbs;}
double BdtoKll_ffpar::lat_mSbs(){return lat_m_Sbs;}
double BdtoKll_ffpar::b0p(){return b0_p;}
double BdtoKll_ffpar::b1p(){return b1_p;}
double BdtoKll_ffpar::b2p(){return b2_p;}
double BdtoKll_ffpar::b0z(){return b0_z;}
double BdtoKll_ffpar::b1z(){return b1_z;}
double BdtoKll_ffpar::b2z(){return b2_z;}
double BdtoKll_ffpar::b0t(){return b0_t;}
double BdtoKll_ffpar::b1t(){return b1_t;}
double BdtoKll_ffpar::b2t(){return b2_t;}
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
