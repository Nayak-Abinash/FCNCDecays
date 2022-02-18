#ifndef myfun_h
#define myfun_h

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <bits/stdc++.h>
//ROOTFiles
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TLine.h"
#include "TH1.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"

#define mxdm 5
#define MXdm 21
using namespace std;

class myfun : public TRandom3
{
public:
    myfun();
    TRandom3* rndm;
    int mylength(double lst[], int iter, double a, double b);
    double mtrx_tp(double lpar[], double covm[][mxdm]);
    double mnd(double mu, double sigma);
    double mnd_default();
    double mnd_cov(double mu[], double cholesky[][MXdm], double untnormlv[], int i);
};


#endif
