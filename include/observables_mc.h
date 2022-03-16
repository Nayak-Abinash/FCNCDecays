#ifndef observables_mc_h
#define observables_mc_h

#include "amplitudes_mc.h"


//Bd->Kstr,ll:
class BdtoKstrll_obserr : public BdtoKstrll_amperr {
public:
    double J1s(double qsq, double ml, double unv[]);
    double J1c(double qsq, double ml, double unv[]);
    double J2s(double qsq, double ml, double unv[]);
    double J2c(double qsq, double ml, double unv[]);
    double J3(double qsq, double ml, double unv[]);
    double J4(double qsq, double ml, double unv[]);
    double J5(double qsq, double ml, double unv[]);
    double J6s(double qsq, double ml, double unv[]);
    double J6c(double qsq, double ml, double unv[]);
    double J7(double qsq, double ml, double unv[]);
    double J8(double qsq, double ml, double unv[]);
    double J9(double qsq, double ml, double unv[]);
    double diffWidth(double qsq, double ml, double unv[]);
    double FL(double qsq, double ml, double unv[]);
    double AFB(double qsq, double ml, double unv[]);
};

//Bs->phi,ll:
class Bstophill_obserr : public Bstophill_amperr {
public:
    double J1s(double qsq, double ml, double unv[]);
    double J1c(double qsq, double ml, double unv[]);
    double J2s(double qsq, double ml, double unv[]);
    double J2c(double qsq, double ml, double unv[]);
    double J3(double qsq, double ml, double unv[]);
    double J4(double qsq, double ml, double unv[]);
    double J5(double qsq, double ml, double unv[]);
    double J6s(double qsq, double ml, double unv[]);
    double J6c(double qsq, double ml, double unv[]);
    double J7(double qsq, double ml, double unv[]);
    double J8(double qsq, double ml, double unv[]);
    double J9(double qsq, double ml, double unv[]);
    double diffWidth(double qsq, double ml, double unv[]);
    double FL(double qsq, double ml, double unv[]);
    double AFB(double qsq, double ml, double unv[]);
};

//Bd->K,ll
class BdtoKll_obserr : public BdtoKll_amperr {
public:
    double alNP(double qsq, double ml, double unv[]);
    double blNP(double qsq, double ml, double unv[]);
    double clNP(double qsq, double ml, double unv[]);
    double diffWidth(double qsq, double ml, double unv[]);
    double diffBrnch(double qsq, double ml, double unv[]);
    double diffAFB(double qsq, double ml, double unv[]);
    double diffFH(double qsq, double ml, double unv[]);
};

/*
//B->ll
class Btoll_obserr : public Btoll_amperr {
public:
    double Btoll_nf(double mBq);
    double ADeltaGammaf(double mBq, double mq, double ml1, double ml2);
    double CorrctnFctr(double mBq, double mq, double ml1, double ml2);
    double BrInst(double mBq, double mq, double ml1, double ml2);
    double BrTimeIntgratd(double mBq, double mq, double ml1, double ml2);
    double efftau(double mBq, double mq, double ml1, double ml2);
};
*/

#endif

