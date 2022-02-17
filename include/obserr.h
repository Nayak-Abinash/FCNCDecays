#ifndef obserr_h
#define obserr_h

#include "amperr.h"


//Bd->Kstr,ll:
class BdtoKstrll_obserr : public BdtoKstrll_amperr {
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

//Bs->phi,ll:
class Bstophill_obserr : public Bstophill_amperr {
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

//Bd->K,ll
class BdtoKll_obserr : public BdtoKll_amperr {
public:
    double alNP(double qsq, double ml);
    double blNP(double qsq, double ml);
    double clNP(double qsq, double ml);
    double diffWidth(double qsq, double ml);
    double diffBrnch(double qsq, double ml);
    double diffAFB(double qsq, double ml);
    double diffFH(double qsq, double ml);
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

