#ifndef fffun_h
#define fffun_h

#include "smpar.h"
#include "ffpar.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

class fffun : public ffpar {
public:
    //Bs->phill
    double V(double qsq),A0(double qsq),A1(double qsq),A12(double qsq),T1(double qsq),T2(double qsq),T23(double qsq);
};

#endif
