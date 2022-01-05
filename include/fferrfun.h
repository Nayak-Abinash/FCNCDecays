#ifndef fferrfun_h
#define fferrfun_h

#include "smpar.h"
#include "myfun.h"
#include "fferrpar.h"
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <string>
#include <iostream>
using namespace std;

//Bd->Kstr,ll
class BdtoKstrll_fferrfun : public BdtoKstrll_fferrpar {
public:
    double V(double qsq),A0(double qsq),A1(double qsq),A12(double qsq),T1(double qsq),T2(double qsq),T23(double qsq);
};

//Bs->phi,ll
class Bstophill_fferrfun : public Bstophill_fferrpar {
public:
    double V(double qsq),A0(double qsq),A1(double qsq),A12(double qsq),T1(double qsq),T2(double qsq),T23(double qsq);
};

//Bd->K,ll
class BdtoKll_fferrfun : public BdtoKll_fferrpar {
public:
    double lcsrfp(double qsq),lcsrfz(double qsq),lcsrft(double qsq),latfp(double qsq),latfz(double qsq),latft(double qsq);
    double fp(double qsq),fz(double qsq),ft(double qsq);
};

#endif
