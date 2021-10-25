#ifndef fffun_h
#define fffun_h

#include "smpar.h"
#include "ffpar.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

//Bd->Kstr,ll
class BdtoKstrll_fffun : public BdtoKstrll_ffpar {
public:
    double V(double qsq),A0(double qsq),A1(double qsq),A12(double qsq),T1(double qsq),T2(double qsq),T23(double qsq);
};

//Bs->phi,ll
class Bstophill_fffun : public Bstophill_ffpar {
public:
    double V(double qsq),A0(double qsq),A1(double qsq),A12(double qsq),T1(double qsq),T2(double qsq),T23(double qsq);
};

//Bd->K,ll
class BdtoKll_fffun : public BdtoKll_ffpar {
public:
    double lcsrfp(double qsq),lcsrfz(double qsq),lcsrft(double qsq),latfp(double qsq),latfz(double qsq),latft(double qsq);
    double fp(double qsq),fz(double qsq),ft(double qsq);
};

#endif
