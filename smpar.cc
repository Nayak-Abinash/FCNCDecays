#include "smpar.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

smpar::smpar(){pi=M_PI;}
//Quark_Mass
double smpar::md(){return 0.00467;}
double smapr::mc(){return 1.27;}
double smapr::ms(){return 0.093;}
double smapr::mb(){return 4.18;}
//Lepton_Mass
double smapr::me(){return 0.511*pow(10.0,-3.0);}
double smapr::mmu(){return 105.658*pow(10.0,-3.0);}
double smapr::mtau(){return 1.77686;}
//Meson_Mass
double smapr::mBd(){return 5.27965;}
double smapr::mBs(){return 5.36688;}
double smapr::mK(){return 0.497614;}
double smapr::mKst(){return 0.896;}
//Width_Lifetime
double smapr::tauBd(){return 1.519*1.52*pow(10.0,12.0);}
double smapr::tauBs(){return 1.527*1.52*pow(10.0,12.0);}
double smapr::DGamma_dbar(){return 0.001;}
double smapr::DGamma_sbar(){return 0.124;}
double smapr::fBd(){return 0.1905;}
double smapr::fBs(){return 0.2303;}
//Coupling_Constants
double smapr::absVtbVtdStr(){return 0.00853871;}
double smapr::absVtbVtsStr(){return 0.0408726;}
double smapr::GF(){return 1.16637*pow(10.0,-5.0);}
double smapr::alpha_e(){return 1.0/127.944;}
//SM_Wilson_Coefficients
double smapr::C1(){return -0.257;}
double smapr::C2(){return 1.009;}
double smapr::C3(){return -0.005;}
double smapr::C4(){return -0.078;}
double smapr::C5(){return 0.000;}
double smapr::C6(){return 0.001;}
double smapr::C7effRe(){return -0.304;}
double smapr::C8eff(){return -0.167;}
double smapr::C9Re(){return 4.211;}
double smapr::C9Im(){return 0.0;}
double smapr::C10Re(){return -4.103;}
double smapr::C10Im(){return 0.0;}
//NP_Wilson_Coefficients
double smapr::C7RHRe(){return 0.0;}
double smapr::C7RHIm(){return 0.0;}
double smapr::C9RHRe(){return 0.0;}
double smapr::C9RHIm(){return 0.0;}
double smapr::C10RHRe(){return 0.0;}
double smapr::C10RHIm(){return 0.0;}
double smapr::CSRe(){return 0.0;}
double smapr::CSIm(){return 0.0;}
double smapr::CSRHRe(){return 0.0;}
double smapr::CSRHIm(){return 0.0;}
double smapr::CPRe(){return 0.0;}
double smapr::CPIm(){return 0.0;}
double smapr::CPRHRe(){return 0.0;}
double smapr::CPRHIm(){return 0.0;}
double smapr::CTRe(){return 0.0;}
double smapr::CTIm(){return 0.0;}
double smapr::CT5Re(){return 0.0;}
double smapr::CT5Im(){return 0.0;}

