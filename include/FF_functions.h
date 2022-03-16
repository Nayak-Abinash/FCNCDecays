#ifndef FF_functions_h
#define FF_functions_h

#include "FF_parameters.h"

//Bd->Kstr,ll
class BdtoKstrll_fffun : public BdtoKstrll_ffpar {
public:
    double V(double qsq),A0(double qsq),A1(double qsq),A12(double qsq),T1(double qsq),T2(double qsq),T23(double qsq);
    double ErV(double qsq),ErA0(double qsq),ErA1(double qsq),ErA12(double qsq),ErT1(double qsq),ErT2(double qsq),ErT23(double qsq);
};

//Bs->phi,ll
class Bstophill_fffun : public Bstophill_ffpar {
public:
    double V(double qsq),A0(double qsq),A1(double qsq),A12(double qsq),T1(double qsq),T2(double qsq),T23(double qsq);
    double ErV(double qsq),ErA0(double qsq),ErA1(double qsq),ErA12(double qsq),ErT1(double qsq),ErT2(double qsq),ErT23(double qsq);
};

//Bd->K,ll
class BdtoKll_fffun : public BdtoKll_ffpar {
public:
    double fz(double qsq),fT(double qsq),fp(double qsq);
    double Erfz(double qsq),ErfT(double qsq),Erfp(double qsq);
};

#endif
