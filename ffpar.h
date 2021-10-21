#ifndef ffpar_h
#define ffpar_h

#include "smpar.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

class ffpar : public smpar {
public:
    ffpar();
    double tpd,tps,tmd,tms,tzd,tzs;
    double zd(double qsq),zs(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(),Va1(),Va2(),A0a0(),A0a1(),A0a2(),A1a0(),A1a1(),A1a2(),A12a0(),A12a1(),A12a2(),T1a0(),T1a1(),T1a2(),T2a0(),T2a1(),T2a2(),T23a0(),T23a1(),T23a2();
private:
    //Bs->phill
    double m_PSbs,m_Vbs,m_Abs;
    double  V_a0,V_a1,V_a2,A0_a0,A0_a1,A0_a2,A1_a0,A1_a1,A1_a2,A12_a0,A12_a1,A12_a2,T1_a0,T1_a1,T1_a2,T2_a0,T2_a1,T2_a2,T23_a0,T23_a1,T23_a2;
};

#endif
