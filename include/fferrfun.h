#ifndef fferrfun_h
#define fferrfun_h

#include "fferrpar.h"

//Bd->Kstr,ll
class BdtoKstrll_fferrfun : public BdtoKstrll_fferrpar {
public:
    double V(double qsq, double unv[]), A0(double qsq, double unv[]), A1(double qsq, double unv[]), A12(double qsq, double unv[]),
        T1(double qsq, double unv[]), T2(double qsq, double unv[]), T23(double qsq, double unv[]);
};

//Bs->phi,ll
class Bstophill_fferrfun : public Bstophill_fferrpar {
public:
    double V(double qsq, double unv[]), A0(double qsq, double unv[]), A1(double qsq, double unv[]), A12(double qsq, double unv[]),
        T1(double qsq, double unv[]), T2(double qsq, double unv[]), T23(double qsq, double unv[]);
};

//Bd->K,ll
class BdtoKll_fferrfun : public BdtoKll_fferrpar {
public:
    double fz(double qsq, double unv[]), fT(double qsq, double unv[]), fp(double qsq, double unv[]);
};

#endif
