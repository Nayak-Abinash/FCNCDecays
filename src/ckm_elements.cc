//Here, we define the CKM matrix in Wolfenstein parameterization (https://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf)
//Though the matrix has a scale dependence above weak scale, below \mu = m_W the CKM elements can be treated as constants.
//Since, we are mostly working with effective theories much below the weak scale, they are assumed to be constants.
#include "ckm_elements.h"

ckm_elements::ckm_elements(){
    //the following values are taken from (http://ckmfitter.in2p3.fr/www/results/plots_spring21/num/ckmEval_results_spring21.pdf)
    //the probability distribution for these fitted parameters are not Gaussian. The errors are asymmetric.
    //However, when we do the MC simulation we shall assume a Gaussian form with the largest error being the standard deviation.
    ckm_A = 0.8132;
    ckm_lambda = 0.22500;
    ckm_rhobar = 0.1566;
    ckm_etabar = 0.3475;}

double ckm_elements::ckmA(){return ckm_A;}
double ckm_elements::ckmlambda(){return ckm_lambda;}
double ckm_elements::ckmrhobar(){return ckm_rhobar;}
double ckm_elements::ckmetabar(){return ckm_etabar;}

double ckm_elements::s12(){return ckmlambda();}
double ckm_elements::s23(){return ckmA()*pow(ckmlambda(),2);}

double ckm_elements::s13pidRe(){
    return (ckmA()*pow(ckmlambda(),3)*sqrt(1 - pow(ckmA(),2)*pow(ckmlambda(),4))*(-(pow(ckmA(),2)*pow(ckmetabar(),2)*pow(ckmlambda(),4))
                + ckmrhobar() - pow(ckmA(),2)*pow(ckmlambda(),4)*pow(ckmrhobar(),2)))/
            (sqrt(1 - pow(ckmlambda(),2))*(1 - 2*pow(ckmA(),2)*pow(ckmlambda(),4)*ckmrhobar()
                + pow(ckmA(),4)*pow(ckmlambda(),8)*(pow(ckmetabar(),2) + pow(ckmrhobar(),2))));}

double ckm_elements::s13pidIm(){
    return (ckmA()*ckmetabar()*pow(ckmlambda(),3)*sqrt(1 - pow(ckmA(),2)*pow(ckmlambda(),4)))/
            (sqrt(1 - pow(ckmlambda(),2))*(1 - 2*pow(ckmA(),2)*pow(ckmlambda(),4)*ckmrhobar()
                    + pow(ckmA(),4)*pow(ckmlambda(),8)*(pow(ckmetabar(),2) + pow(ckmrhobar(),2))));}

double ckm_elements::c12(){return sqrt(1 - pow(s12(),2));}
double ckm_elements::c23(){return sqrt(1 - pow(s23(),2));}
double ckm_elements::c13(){return sqrt(1 - pow(s13pidRe(),2) - pow(s13pidIm(),2));}

double ckm_elements::ReVud(){return c12()*c13();}
double ckm_elements::ImVud(){return 0.0;}

double ckm_elements::ReVus(){return c13()*s12();}
double ckm_elements::ImVus(){return 0.0;}

double ckm_elements::ReVub(){return s13pidRe();}
double ckm_elements::ImVub(){return -s13pidIm();}

double ckm_elements::ReVcd(){return -c23()*s12() - c12()*s23()*s13pidRe();}
double ckm_elements::ImVcd(){return -c12()*s23()*s13pidIm();}

double ckm_elements::ReVcs(){return c12()*c23() - s12()*s23()*s13pidRe();}
double ckm_elements::ImVcs(){return -s12()*s23()*s13pidIm();}

double ckm_elements::ReVcb(){return c13()*s23();}
double ckm_elements::ImVcb(){return 0.0;}

double ckm_elements::ReVtd(){return -c12()*c23()*s13pidRe() + s12()*s23();}
double ckm_elements::ImVtd(){return -c12()*c23()*s13pidIm();}

double ckm_elements::ReVts(){return -c23()*s12()*s13pidRe() - c12()*s23();}
double ckm_elements::ImVts(){return -c23()*s12()*s13pidIm();}

double ckm_elements::ReVtb(){return c13()*c23();}
double ckm_elements::ImVtb(){return 0.0;}

double ckm_elements::absVtbVtdStr(){return sqrt(pow(ReVtb()*ReVtd() + ImVtb()*ImVtd(), 2) + pow(ImVtb()*ReVtd() - ReVtb()*ImVtd(), 2));}

double ckm_elements::absVtbVtsStr(){return sqrt(pow(ReVtb()*ReVts() + ImVtb()*ImVts(), 2) + pow(ImVtb()*ReVts() - ReVtb()*ImVts(), 2));}

double ckm_elements::ReVtbVtdStr(){return ReVtb()*ReVtd() + ImVtb()*ImVtd();}

double ckm_elements::ImVtbVtdStr(){return ImVtb()*ReVtd() - ReVtb()*ImVtd();}

double ckm_elements::ReVtbVtsStr(){return ReVtb()*ReVts() + ImVtb()*ImVts();}

double ckm_elements::ImVtbVtsStr(){return ImVtb()*ReVts() - ReVtb()*ImVts();}

































