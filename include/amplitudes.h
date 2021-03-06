#ifndef amplitudes_h
#define amplitudes_h

#include "FF_functions.h"

//Bd->Kstr,ll
class BdtoKstrll_amp : public BdtoKstrll_fffun {
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
    double Renf(double qsq, double ml);
    double Imnf(double qsq, double ml);

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

    double ApbLRe(double qsq, double ml);
    double ApbLIm(double qsq, double ml);
    double ApbRRe(double qsq, double ml);
    double ApbRIm(double qsq, double ml);
    double AabLRe(double qsq, double ml);
    double AabLIm(double qsq, double ml);
    double AabRRe(double qsq, double ml);
    double AabRIm(double qsq, double ml);
    double AzbLRe(double qsq, double ml);
    double AzbLIm(double qsq, double ml);
    double AzbRRe(double qsq, double ml);
    double AzbRIm(double qsq, double ml);
    double AtbRe(double qsq, double ml);
    double AtbIm(double qsq, double ml);
    double ASbRe(double qsq, double ml);
    double ASbIm(double qsq, double ml);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll
class Bstophill_amp : public Bstophill_fffun {
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
    double Renf(double qsq, double ml);
    double Imnf(double qsq, double ml);

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

    double ApbLRe(double qsq, double ml);
    double ApbLIm(double qsq, double ml);
    double ApbRRe(double qsq, double ml);
    double ApbRIm(double qsq, double ml);
    double AabLRe(double qsq, double ml);
    double AabLIm(double qsq, double ml);
    double AabRRe(double qsq, double ml);
    double AabRIm(double qsq, double ml);
    double AzbLRe(double qsq, double ml);
    double AzbLIm(double qsq, double ml);
    double AzbRRe(double qsq, double ml);
    double AzbRIm(double qsq, double ml);
    double AtbRe(double qsq, double ml);
    double AtbIm(double qsq, double ml);
    double ASbRe(double qsq, double ml);
    double ASbIm(double qsq, double ml);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
class BdtoKll_amp : public BdtoKll_fffun {
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//B->ll
class Btoll_amp : public smpar {
public:
    double ampSRe(double mBq, double mq, double ml1, double ml2);
    double ampSIm(double mBq, double mq, double ml1, double ml2);
    double ampPRe(double mBq, double mq, double ml1, double ml2);
    double ampPIm(double mBq, double mq, double ml1, double ml2);
    double Btoll_lambda(double mBq, double ml1, double ml2);
    double y(double mBq);
    double tauB(double mBq);
};


#endif
