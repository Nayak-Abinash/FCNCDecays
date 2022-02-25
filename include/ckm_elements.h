#ifndef ckm_elements_h
#define ckm_elements_h

#include "ref_fun.h"

class ckm_elements {
public:
    ckm_elements();
    double ckmA(), ckmlambda(), ckmrhobar(), ckmetabar();
    double s12(), s23(), s13pidRe(), s13pidIm();
    double c12(), c23(), c13();
    double ReVud(), ReVus(), ReVub(), ReVcd(), ReVcs(), ReVcb(), ReVtd(), ReVts(), ReVtb();
    double ImVud(), ImVus(), ImVub(), ImVcd(), ImVcs(), ImVcb(), ImVtd(), ImVts(), ImVtb();
private:
    double ckm_A, ckm_lambda, ckm_rhobar, ckm_etabar;
};

#endif // ckm_elements_h
