//Here, we define the CKM matrix in Wolfenstein parameterization (https://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf)
//Though the matrix has a scale dependence above weak scale, below \mu = m_W the CKM elements can be treated as constants.
//Since, we are mostly working with effective theories much below the weak scale, they are assumed to be constants.
#include "ckm_elements_mc.h"

//the following values are taken from (http://ckmfitter.in2p3.fr/www/results/plots_spring21/num/ckmEval_results_spring21.pdf)
//the probability distribution for these fitted parameters are not Gaussian. The errors are asymmetric.
//However, when we do the MC simulation we shall assume a Gaussian form with the largest error being the standard deviation.
double ckm_elements_mc::cen_ckm[] = {0.8132, 0.22500, 0.1566, 0.3475};
double ckm_elements_mc::unc_ckm[] = {0.0119, 0.00024, 0.0085, 0.0118};
double ckm_elements_mc::chd_ckm[][mxdm] = {{0.0119, 0, 0, 0}, {0, 0.00024, 0, 0}, {0, 0, 0.0085, 0}, {0, 0, 0, 0.0118}};

double ckm_elements_mc::ckmA(double unv[]){return mnd_cov_ckm(cen_ckm, chd_ckm, unv, 0);}
double ckm_elements_mc::ckmlambda(double unv[]){return mnd_cov_ckm(cen_ckm, chd_ckm, unv, 1);}
double ckm_elements_mc::ckmrhobar(double unv[]){return mnd_cov_ckm(cen_ckm, chd_ckm, unv, 2);}
double ckm_elements_mc::ckmetabar(double unv[]){return mnd_cov_ckm(cen_ckm, chd_ckm, unv, 3);}

double ckm_elements_mc::s12(double unv[]){return ckmlambda(unv);}
double ckm_elements_mc::s23(double unv[]){return ckmA(unv)*pow(ckmlambda(unv),2);}

double ckm_elements_mc::s13pidRe(double unv[]){
    return (ckmA(unv)*pow(ckmlambda(unv),3)*sqrt(1 - pow(ckmA(unv),2)*pow(ckmlambda(unv),4))*(-(pow(ckmA(unv),2)*pow(ckmetabar(unv),2)*pow(ckmlambda(unv),4))
                + ckmrhobar(unv) - pow(ckmA(unv),2)*pow(ckmlambda(unv),4)*pow(ckmrhobar(unv),2)))/
            (sqrt(1 - pow(ckmlambda(unv),2))*(1 - 2*pow(ckmA(unv),2)*pow(ckmlambda(unv),4)*ckmrhobar(unv)
                + pow(ckmA(unv),4)*pow(ckmlambda(unv),8)*(pow(ckmetabar(unv),2) + pow(ckmrhobar(unv),2))));}

double ckm_elements_mc::s13pidIm(double unv[]){
    return (ckmA(unv)*ckmetabar(unv)*pow(ckmlambda(unv),3)*sqrt(1 - pow(ckmA(unv),2)*pow(ckmlambda(unv),4)))/
            (sqrt(1 - pow(ckmlambda(unv),2))*(1 - 2*pow(ckmA(unv),2)*pow(ckmlambda(unv),4)*ckmrhobar(unv)
                    + pow(ckmA(unv),4)*pow(ckmlambda(unv),8)*(pow(ckmetabar(unv),2) + pow(ckmrhobar(unv),2))));}

double ckm_elements_mc::c12(double unv[]){return sqrt(1 - pow(s12(unv),2));}
double ckm_elements_mc::c23(double unv[]){return sqrt(1 - pow(s23(unv),2));}
double ckm_elements_mc::c13(double unv[]){return sqrt(1 - pow(s13pidRe(unv),2) - pow(s13pidIm(unv),2));}

double ckm_elements_mc::ReVud(double unv[]){return c12(unv)*c13(unv);}
double ckm_elements_mc::ImVud(double unv[]){return 0.0;}

double ckm_elements_mc::ReVus(double unv[]){return c13(unv)*s12(unv);}
double ckm_elements_mc::ImVus(double unv[]){return 0.0;}

double ckm_elements_mc::ReVub(double unv[]){return s13pidRe(unv);}
double ckm_elements_mc::ImVub(double unv[]){return -s13pidIm(unv);}

double ckm_elements_mc::ReVcd(double unv[]){return -c23(unv)*s12(unv) - c12(unv)*s23(unv)*s13pidRe(unv);}
double ckm_elements_mc::ImVcd(double unv[]){return -c12(unv)*s23(unv)*s13pidIm(unv);}

double ckm_elements_mc::ReVcs(double unv[]){return c12(unv)*c23(unv) - s12(unv)*s23(unv)*s13pidRe(unv);}
double ckm_elements_mc::ImVcs(double unv[]){return -s12(unv)*s23(unv)*s13pidIm(unv);}

double ckm_elements_mc::ReVcb(double unv[]){return c13(unv)*s23(unv);}
double ckm_elements_mc::ImVcb(double unv[]){return 0.0;}

double ckm_elements_mc::ReVtd(double unv[]){return -c12(unv)*c23(unv)*s13pidRe(unv) + s12(unv)*s23(unv);}
double ckm_elements_mc::ImVtd(double unv[]){return -c12(unv)*c23(unv)*s13pidIm(unv);}

double ckm_elements_mc::ReVts(double unv[]){return -c23(unv)*s12(unv)*s13pidRe(unv) - c12(unv)*s23(unv);}
double ckm_elements_mc::ImVts(double unv[]){return -c23(unv)*s12(unv)*s13pidIm(unv);}

double ckm_elements_mc::ReVtb(double unv[]){return c13(unv)*c23(unv);}
double ckm_elements_mc::ImVtb(double unv[]){return 0.0;}

double ckm_elements_mc::absVtbVtdStr(double unv[]){return sqrt(pow(ReVtb(unv)*ReVtd(unv) + ImVtb(unv)*ImVtd(unv), 2)
                                                               + pow(ImVtb(unv)*ReVtd(unv) - ReVtb(unv)*ImVtd(unv), 2));}

double ckm_elements_mc::absVtbVtsStr(double unv[]){return sqrt(pow(ReVtb(unv)*ReVts(unv) + ImVtb(unv)*ImVts(unv), 2)
                                                               + pow(ImVtb(unv)*ReVts(unv) - ReVtb(unv)*ImVts(unv), 2));}

double ckm_elements_mc::ReVtbVtdStr(double unv[]){return ReVtb(unv)*ReVtd(unv) + ImVtb(unv)*ImVtd(unv);}

double ckm_elements_mc::ImVtbVtdStr(double unv[]){return ImVtb(unv)*ReVtd(unv) - ReVtb(unv)*ImVtd(unv);}

double ckm_elements_mc::ReVtbVtsStr(double unv[]){return ReVtb(unv)*ReVts(unv) + ImVtb(unv)*ImVts(unv);}

double ckm_elements_mc::ImVtbVtsStr(double unv[]){return ImVtb(unv)*ReVts(unv) - ReVtb(unv)*ImVts(unv);}
































