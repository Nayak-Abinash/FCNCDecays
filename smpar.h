
#ifndef smpar_h
#define smpar_h

#include <iostream>
#include <cmath>
#include <string>

using namespace std;

//B->ll:
class smpar {
public:
    double pi,mB,mKst,mu,mc,mb;
    
    double C1,C2,C3,C4,C5,C6,C7effRe,C8eff,C10effRe;
    
    double tauB,alpha_e,GF,absVtbVtsStr,me,mmu;
    
    smpar ();
};

#endif
