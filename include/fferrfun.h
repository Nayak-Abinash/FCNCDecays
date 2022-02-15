#ifndef fferrfun_h
#define fferrfun_h

#include "myfun.h"
#include "smpar.h"
#include "fferrpar.h"

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
    double fz(double qsq),fT(double qsq),fp(double qsq);
};

#endif
