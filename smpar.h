#ifndef smpar_h
#define smpar_h

#include <iostream>
#include <cmath>
#include <string>
using namespace std;

class smpar {
public:
    double pi;
    double mb(),md(),ms();
    double me(),mmu(),mtau();
    double mBd(),mBs(),mKst();
    double tauBd(),tauBs(),DGamma_dbar(),DGamma_sbar(),fBd(),fBs();
    double absVtbVtdStr(),absVtbVtsStr(),GF(),alpha_e();
    double C1(),C2(),C3(),C4(),C5(),C6(),C7effRe(),C8eff(),C10effRe(),C9Re(),C10Re(),C9Im();
    double C9RHRe(),C9RHIm(),C10Im(),C10RHRe(),C10RHIm(),CSRe(),CSIm(),CSRHRe(),CSRHIm(),CPRe(),CPIm(),CPRHRe(),CPRHIm();
    
    smpar ();
};

#endif
