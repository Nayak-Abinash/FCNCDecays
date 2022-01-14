#ifndef amperr_h
#define amperr_h

#include "smpar.h"
#include "myfun.h"
#include "fferrpar.h"
#include "fferrfun.h"


//Bd->Kstr,ll
class BdtoKstrll_amperr : public BdtoKstrll_fferrfun {
public:
    double hRe(double qsq, double mq);
    double hIm(double qsq, double mq);
    double hReZ(double qsq);
    double hImZ(double qsq);
    double YRe(double qsq);
    double YIm(double qsq);
    double C9effRe(double qsq);
    double C9effIm(double qsq);
    double betal(double qsq, double ml);
    double lambda(double qsq);
    double nf(double qsq, double ml);

    double ApLRe(double qsq, double ml);
    double ApLIm(double qsq, double ml);
    double ApRRe(double qsq, double ml);
    double ApRIm(double qsq, double ml);
    double AaLRe(double qsq, double ml);
    double AaLIm(double qsq, double ml);
    double AaRRe(double qsq, double ml);
    double AaRIm(double qsq, double ml);
    double AzLRe(double qsq, double ml);
    double AzLIm(double qsq, double ml);
    double AzRRe(double qsq, double ml);
    double AzRIm(double qsq, double ml);
    double AtRe(double qsq, double ml);
    double AtIm(double qsq, double ml);
    double ASRe(double qsq, double ml);
    double ASIm(double qsq, double ml);
};

//Bs->phi,ll
class Bstophill_amperr : public Bstophill_fferrfun {
public:
    double hRe(double qsq, double mq);
    double hIm(double qsq, double mq);
    double hReZ(double qsq);
    double hImZ(double qsq);
    double YRe(double qsq);
    double YIm(double qsq);
    double C9effRe(double qsq);
    double C9effIm(double qsq);
    double betal(double qsq, double ml);
    double lambda(double qsq);
    double nf(double qsq, double ml);

    double ApLRe(double qsq, double ml);
    double ApLIm(double qsq, double ml);
    double ApRRe(double qsq, double ml);
    double ApRIm(double qsq, double ml);
    double AaLRe(double qsq, double ml);
    double AaLIm(double qsq, double ml);
    double AaRRe(double qsq, double ml);
    double AaRIm(double qsq, double ml);
    double AzLRe(double qsq, double ml);
    double AzLIm(double qsq, double ml);
    double AzRRe(double qsq, double ml);
    double AzRIm(double qsq, double ml);
    double AtRe(double qsq, double ml);
    double AtIm(double qsq, double ml);
    double ASRe(double qsq, double ml);
    double ASIm(double qsq, double ml);
};

//Bd->K,ll
class BdtoKll_amperr : public BdtoKll_fferrfun {
public:
    double hRe(double qsq, double mq);
    double hIm(double qsq, double mq);
    double hReZ(double qsq);
    double hImZ(double qsq);
    double YRe(double qsq);
    double YIm(double qsq);
    double C9effRe(double qsq);
    double C9effIm(double qsq);

    double betal(double qsq, double ml);
    double lambda(double qsq);
    double nf(double qsq, double ml);

    double xiP(double qsq);
    double TauPRe(double qsq);
    double TauPIm(double qsq);
    double FVRe(double qsq, double ml);
    double FVIm(double qsq, double ml);
    double FARe(double qsq, double ml);
    double FAIm(double qsq, double ml);
    double FSRe(double qsq, double ml);
    double FSIm(double qsq, double ml);
    double FPRe(double qsq, double ml);
    double FPIm(double qsq, double ml);
    double FTRe(double qsq, double ml);
    double FTIm(double qsq, double ml);
    double FT5Re(double qsq, double ml);
    double FT5Im(double qsq, double ml);
};

/*
//B->ll
class Btoll_amperr : public smpar {
public:
    double ampSRe(double mBq, double mq, double ml1, double ml2);
    double ampSIm(double mBq, double mq, double ml1, double ml2);
    double ampPRe(double mBq, double mq, double ml1, double ml2);
    double ampPIm(double mBq, double mq, double ml1, double ml2);
    double Btoll_lambda(double mBq, double ml1, double ml2);
    double y(double mBq);
    double tauB(double mBq);
};
*/


#endif
