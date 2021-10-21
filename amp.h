#ifndef amp_h
#define amp_h

#include "smpar.h"
#include "ffpar.h"
#include "fffun.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

//Bs->phill:
class amp : public fffun {
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


#endif
