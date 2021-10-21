#ifndef obs_h
#define obs_h

#include "smpar.h"
#include "ffpar.h"
#include "fffun.h"
#include "amp.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

//Bs->phill:
class obs : public amp {
public:
    double J1s(double qsq, double ml);
    double J1c(double qsq, double ml);
    double J2s(double qsq, double ml);
    double J2c(double qsq, double ml);
    double J3(double qsq, double ml);
    double J4(double qsq, double ml);
    double J5(double qsq, double ml);
    double J6s(double qsq, double ml);
    double J6c(double qsq, double ml);
    double J7(double qsq, double ml);
    double J8(double qsq, double ml);
    double J9(double qsq, double ml);
    double diffWidth(double qsq, double ml);
    double FL(double qsq, double ml);
    double AFB(double qsq, double ml);
};

#endif

