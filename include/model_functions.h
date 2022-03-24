#ifndef model_functions_h
#define model_functions_h

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFrame.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TLine.h"
#include "TH1.h"
#include "TRandom3.h"

#define mxdm 5
#define MXdm 25
using namespace std;

class model_functions : public TRandom3
{
public:
    model_functions();
    TRandom3* rndm;
    int mylength(double lst[], int iter, double a, double b);
    double mtrx_tp(double lpar[], double covm[][mxdm]);
    double mnd(double mu, double sigma);
    double mnd_default();
    double mnd_cov(double mu[], double cholesky[][MXdm], double untnormlv[], int i);
    double mnd_cov_ckm(double mu[], double cholesky[][mxdm], double untnormlv[], int i);
    double mnd_cov_wcSM(double mu[], double cholesky[][MXdm], double untnormlv[], int i);
    double mnd_cov_wcNP(double mu[], double cholesky[][MXdm], double untnormlv[], int i);
    double mnd_cov_smpar(double mu[], double cholesky[][MXdm], double untnormlv[], int i);
};


#endif



















