#ifndef ckm_elements_mc_h
#define ckm_elements_mc_h

#include "model_functions.h"

class ckm_elements_mc : public model_functions
{
public:
    //ckm_elements_mc();
    double static cen_ckm[mxdm], unc_ckm[mxdm], chd_ckm[mxdm][mxdm];
    double ckmA(double unv[]), ckmlambda(double unv[]), ckmrhobar(double unv[]), ckmetabar(double unv[]);
    double s12(double unv[]), s23(double unv[]), s13pidRe(double unv[]), s13pidIm(double unv[]);
    double c12(double unv[]), c23(double unv[]), c13(double unv[]);
    double ReVud(double unv[]), ReVus(double unv[]), ReVub(double unv[]), ReVcd(double unv[]), ReVcs(double unv[]),
        ReVcb(double unv[]), ReVtd(double unv[]), ReVts(double unv[]), ReVtb(double unv[]);
    double ImVud(double unv[]), ImVus(double unv[]), ImVub(double unv[]), ImVcd(double unv[]), ImVcs(double unv[]),
        ImVcb(double unv[]), ImVtd(double unv[]), ImVts(double unv[]), ImVtb(double unv[]);
    double absVtbVtdStr(double unv[]), absVtbVtsStr(double unv[]);
};

#endif
