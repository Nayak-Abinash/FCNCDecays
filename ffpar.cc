#include "smpar.h"
#include "ffpar.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

//Bs->phill
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
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
