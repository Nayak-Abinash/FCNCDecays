#ifndef amplitudes_mc_h
#define amplitudes_mc_h

#include "FF_functions_mc.h"


//Bd->Kstr,ll
class BdtoKstrll_amperr : public BdtoKstrll_fferrfun {
public:
    double hRe(double qsq, double mq);
    double hIm(double qsq, double mq);
    double hReZ(double qsq);
    double hImZ(double qsq);
    double YRe(double qsq, double smwc[], double npwc[], double unv[]);
    double YIm(double qsq, double smwc[], double npwc[], double unv[]);
    double C9effRe(double qsq, double smwc[], double npwc[], double unv[]);
    double C9effIm(double qsq, double smwc[], double npwc[], double unv[]);
    double betal(double qsq, double ml);
    double lambda(double qsq, double unv[]);
    double nf(double qsq, double ml, double unv[]);

    double ApLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ApLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ApRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ApRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AaLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AaLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AaRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AaRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AzLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AzLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AzRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AzRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AtRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AtIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ASRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ASIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll
class Bstophill_amperr : public Bstophill_fferrfun {
public:
    double hRe(double qsq, double mq);
    double hIm(double qsq, double mq);
    double hReZ(double qsq);
    double hImZ(double qsq);
    double YRe(double qsq, double smwc[], double npwc[], double unv[]);
    double YIm(double qsq, double smwc[], double npwc[], double unv[]);
    double C9effRe(double qsq, double smwc[], double npwc[], double unv[]);
    double C9effIm(double qsq, double smwc[], double npwc[], double unv[]);
    double betal(double qsq, double ml);
    double lambda(double qsq, double unv[]);
    double nf(double qsq, double ml, double unv[]);

    double ApLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ApLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ApRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ApRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AaLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AaLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AaRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AaRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AzLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AzLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AzRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AzRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AtRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AtIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ASRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double ASIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
class BdtoKll_amperr : public BdtoKll_fferrfun {
public:
    double hRe(double qsq, double mq);
    double hIm(double qsq, double mq);
    double hReZ(double qsq);
    double hImZ(double qsq);
    double YRe(double qsq, double smwc[], double npwc[], double unv[]);
    double YIm(double qsq, double smwc[], double npwc[], double unv[]);
    double C9effRe(double qsq, double smwc[], double npwc[], double unv[]);
    double C9effIm(double qsq, double smwc[], double npwc[], double unv[]);

    double betal(double qsq, double ml);
    double lambda(double qsq, double unv[]);
    double nf(double qsq, double ml, double unv[]);

    double xiP(double qsq, double smwc[], double npwc[], double unv[]);
    double TauPRe(double qsq, double smwc[], double npwc[], double unv[]);
    double TauPIm(double qsq, double smwc[], double npwc[], double unv[]);
    double FVRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FVIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FARe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FAIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FSRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FSIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FPRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FPIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FTRe(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FTIm(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FT5Re(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FT5Im(double qsq, double ml, double smwc[], double npwc[], double unv[]);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//B->ll
class Btoll_amperr : public smpar_mc {
public:
    double ampSRe(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
    double ampSIm(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
    double ampPRe(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
    double ampPIm(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
    double Btoll_lambda(double mBq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
    double y(double mBq, double smwc[], double npwc[], double unv[]);
    double tauB(double mBq, double smwc[], double npwc[], double unv[]);
};



#endif
