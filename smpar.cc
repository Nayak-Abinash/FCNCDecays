
#include "smpar.h"
#include <iostream>
#include <cmath>
#include <string>

using namespace std;


smpar :: smpar() {
    //SM parameter values:
    pi=M_PI;
    mB=5.27953,mKst=0.896,mu=4.8,mc=1.27,mb=4.8;
    C1=-0.257,C2=1.009,C3=-0.005,C4=-0.078,C5=0.000,C6=0.001,C7effRe=-0.304,C8eff=-0.167,C10effRe=-4.103;
    
    tauB=1.519*1.52*pow(10.0,12.0),alpha_e=1.0/133.0,GF1.16637*pow(10.0,-5.0),absVtbVtsStr=0.0409,me=0.511*pow(10.0,-3.0),mmu=105.658*pow(10.0,-3.0);
}

