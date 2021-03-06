#include "FF_functions.h"

//Bd->Kstr,ll
double BdtoKstrll_fffun::V(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( Va0()*pow(zd(qsq)-zd(0.0),0.0) + Va1()*pow(zd(qsq)-zd(0.0),1.0) + Va2()*pow(zd(qsq)-zd(0.0),2.0) );}
double BdtoKstrll_fffun::ErV(double qsq){
    double lpar_V[] = { pow(zd(qsq)-zd(0.0),0.0), pow(zd(qsq)-zd(0.0),1.0), pow(zd(qsq)-zd(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*sqrt(mtrx_tp(lpar_V,cov_V));}

double BdtoKstrll_fffun::A0(double qsq){
    return 1.0/(1.0-qsq/pow(mPSbs(),2.0))*( A0a0()*pow(zd(qsq)-zd(0.0),0.0) + A0a1()*pow(zd(qsq)-zd(0.0),1.0) + A0a2()*pow(zd(qsq)-zd(0.0),2.0) );}
double BdtoKstrll_fffun::ErA0(double qsq){
    double lpar_A0[] = { pow(zd(qsq)-zd(0.0),0.0), pow(zd(qsq)-zd(0.0),1.0), pow(zd(qsq)-zd(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mPSbs(),2.0))*sqrt(mtrx_tp(lpar_A0,cov_A0));}

double BdtoKstrll_fffun::A1(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A1a0()*pow(zd(qsq)-zd(0.0),0.0) + A1a1()*pow(zd(qsq)-zd(0.0),1.0) + A1a2()*pow(zd(qsq)-zd(0.0),2.0) );}
double BdtoKstrll_fffun::ErA1(double qsq){
    double lpar_A1[] = { pow(zd(qsq)-zd(0.0),0.0), pow(zd(qsq)-zd(0.0),1.0), pow(zd(qsq)-zd(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*sqrt(mtrx_tp(lpar_A1,cov_A1));}

double BdtoKstrll_fffun::A12(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A12a0()*pow(zd(qsq)-zd(0.0),0.0) + A12a1()*pow(zd(qsq)-zd(0.0),1.0) + A12a2()*pow(zd(qsq)-zd(0.0),2.0) );}
double BdtoKstrll_fffun::ErA12(double qsq){
    double lpar_A12[] = { pow(zd(qsq)-zd(0.0),0.0), pow(zd(qsq)-zd(0.0),1.0), pow(zd(qsq)-zd(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*sqrt(mtrx_tp(lpar_A12,cov_A12));}

double BdtoKstrll_fffun::T1(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( T1a0()*pow(zd(qsq)-zd(0.0),0.0) + T1a1()*pow(zd(qsq)-zd(0.0),1.0) + T1a2()*pow(zd(qsq)-zd(0.0),2.0) );}
double BdtoKstrll_fffun::ErT1(double qsq){
    double lpar_T1[] = { pow(zd(qsq)-zd(0.0),0.0), pow(zd(qsq)-zd(0.0),1.0), pow(zd(qsq)-zd(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*sqrt(mtrx_tp(lpar_T1,cov_T1));}

double BdtoKstrll_fffun::T2(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T2a0()*pow(zd(qsq)-zd(0.0),0.0) + T2a1()*pow(zd(qsq)-zd(0.0),1.0) + T2a2()*pow(zd(qsq)-zd(0.0),2.0) );}
double BdtoKstrll_fffun::ErT2(double qsq){
    double lpar_T2[] = { pow(zd(qsq)-zd(0.0),0.0), pow(zd(qsq)-zd(0.0),1.0), pow(zd(qsq)-zd(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*sqrt(mtrx_tp(lpar_T2,cov_T2));}

double BdtoKstrll_fffun::T23(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T23a0()*pow(zd(qsq)-zd(0.0),0.0) + T23a1()*pow(zd(qsq)-zd(0.0),1.0) + T23a2()*pow(zd(qsq)-zd(0.0),2.0) );}
double BdtoKstrll_fffun::ErT23(double qsq){
    double lpar_T23[] = { pow(zd(qsq)-zd(0.0),0.0), pow(zd(qsq)-zd(0.0),1.0), pow(zd(qsq)-zd(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*sqrt(mtrx_tp(lpar_T23,cov_T23));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll
double Bstophill_fffun::V(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( Va0()*pow(zs(qsq)-zs(0.0),0.0) + Va1()*pow(zs(qsq)-zs(0.0),1.0) + Va2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fffun::A0(double qsq){
    return 1.0/(1.0-qsq/pow(mPSbs(),2.0))*( A0a0()*pow(zs(qsq)-zs(0.0),0.0) + A0a1()*pow(zs(qsq)-zs(0.0),1.0) + A0a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fffun::A1(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A1a0()*pow(zs(qsq)-zs(0.0),0.0) + A1a1()*pow(zs(qsq)-zs(0.0),1.0) + A1a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fffun::A12(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A12a0()*pow(zs(qsq)-zs(0.0),0.0) + A12a1()*pow(zs(qsq)-zs(0.0),1.0) + A12a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fffun::T1(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( T1a0()*pow(zs(qsq)-zs(0.0),0.0) + T1a1()*pow(zs(qsq)-zs(0.0),1.0) + T1a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fffun::T2(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T2a0()*pow(zs(qsq)-zs(0.0),0.0) + T2a1()*pow(zs(qsq)-zs(0.0),1.0) + T2a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fffun::T23(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T23a0()*pow(zs(qsq)-zs(0.0),0.0) + T23a1()*pow(zs(qsq)-zs(0.0),1.0) + T23a2()*pow(zs(qsq)-zs(0.0),2.0) );}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
double BdtoKll_fffun::fz(double qsq){
    return 1.0/(1.0-qsq/pow(mSbs(),2.0))*( fza1()*pow(z(qsq)-z(0.0),1.0) + fza2()*pow(z(qsq)-z(0.0),2.0) );}
double BdtoKll_fffun::Erfz(double qsq){
    double lpar_fz[] = { pow(z(qsq)-z(0.0),1.0) , pow(z(qsq)-z(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mSbs(),2.0))*sqrt(mtrx_tp(lpar_fz,cov_fz));}

double BdtoKll_fffun::fT(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( fTa0()*pow(z(qsq)-z(0.0),0.0) + fTa1()*pow(z(qsq)-z(0.0),1.0) + fTa2()*pow(z(qsq)-z(0.0),2.0) );}
double BdtoKll_fffun::ErfT(double qsq){
    double lpar_fT[] = { pow(z(qsq)-z(0.0),0.0), pow(z(qsq)-z(0.0),1.0) , pow(z(qsq)-z(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*sqrt(mtrx_tp(lpar_fT,cov_fT));}


double BdtoKll_fffun::fp(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( fpa0()*pow(z(qsq)-z(0.0),0.0) + fpa1()*pow(z(qsq)-z(0.0),1.0) + fpa2()*pow(z(qsq)-z(0.0),2.0) );}
double BdtoKll_fffun::Erfp(double qsq){
    double lpar_fp[] = { pow(z(qsq)-z(0.0),0.0), pow(z(qsq)-z(0.0),1.0) , pow(z(qsq)-z(0.0),2.0)};
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*sqrt(mtrx_tp(lpar_fp,cov_fp));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////










