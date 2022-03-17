#include "FF_functions_mc.h"


//Bd->Kstr,ll
double BdtoKstrll_fferrfun::V(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( Va0(unv)*pow(zd(qsq, unv)-zd(0.0, unv),0.0) + Va1(unv)*pow(zd(qsq, unv)-zd(0.0, unv),1.0)
                                          + Va2(unv)*pow(zd(qsq, unv)-zd(0.0, unv),2.0) );}

double BdtoKstrll_fferrfun::A0(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mPSbs(),2.0))*( A0a0(unv)*pow(zd(qsq, unv)-zd(0.0, unv),0.0) + A0a1(unv)*pow(zd(qsq, unv)-zd(0.0, unv),1.0)
                                           + A0a2(unv)*pow(zd(qsq, unv)-zd(0.0, unv),2.0) );}

double BdtoKstrll_fferrfun::A1(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A1a0(unv)*pow(zd(qsq, unv)-zd(0.0, unv),0.0) + A1a1(unv)*pow(zd(qsq, unv)-zd(0.0, unv),1.0)
                                          + A1a2(unv)*pow(zd(qsq, unv)-zd(0.0, unv),2.0) );}

double BdtoKstrll_fferrfun::A12(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A12a0(unv)*pow(zd(qsq, unv)-zd(0.0, unv),0.0) + A12a1(unv)*pow(zd(qsq, unv)-zd(0.0, unv),1.0)
                                          + A12a2(unv)*pow(zd(qsq, unv)-zd(0.0, unv),2.0) );}

double BdtoKstrll_fferrfun::T1(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( T1a0(unv)*pow(zd(qsq, unv)-zd(0.0, unv),0.0) + T1a1(unv)*pow(zd(qsq, unv)-zd(0.0, unv),1.0)
                                          + T1a2(unv)*pow(zd(qsq, unv)-zd(0.0, unv),2.0) );}

double BdtoKstrll_fferrfun::T2(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T2a0(unv)*pow(zd(qsq, unv)-zd(0.0, unv),0.0) + T2a1(unv)*pow(zd(qsq, unv)-zd(0.0, unv),1.0)
                                          + T2a2(unv)*pow(zd(qsq, unv)-zd(0.0, unv),2.0) );}

double BdtoKstrll_fferrfun::T23(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T23a0(unv)*pow(zd(qsq, unv)-zd(0.0, unv),0.0) + T23a1(unv)*pow(zd(qsq, unv)-zd(0.0, unv),1.0)
                                          + T23a2(unv)*pow(zd(qsq, unv)-zd(0.0, unv),2.0) );}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll
double Bstophill_fferrfun::V(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( Va0(unv)*pow(zs(qsq, unv)-zs(0.0, unv),0.0) + Va1(unv)*pow(zs(qsq, unv)-zs(0.0, unv),1.0)
                                          + Va2(unv)*pow(zs(qsq, unv)-zs(0.0, unv),2.0) );}

double Bstophill_fferrfun::A0(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mPSbs(),2.0))*( A0a0(unv)*pow(zs(qsq, unv)-zs(0.0, unv),0.0) + A0a1(unv)*pow(zs(qsq, unv)-zs(0.0, unv),1.0)
                                           + A0a2(unv)*pow(zs(qsq, unv)-zs(0.0, unv),2.0) );}

double Bstophill_fferrfun::A1(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A1a0(unv)*pow(zs(qsq, unv)-zs(0.0, unv),0.0) + A1a1(unv)*pow(zs(qsq, unv)-zs(0.0, unv),1.0)
                                          + A1a2(unv)*pow(zs(qsq, unv)-zs(0.0, unv),2.0) );}

double Bstophill_fferrfun::A12(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A12a0(unv)*pow(zs(qsq, unv)-zs(0.0, unv),0.0) + A12a1(unv)*pow(zs(qsq, unv)-zs(0.0, unv),1.0)
                                          + A12a2(unv)*pow(zs(qsq, unv)-zs(0.0, unv),2.0) );}

double Bstophill_fferrfun::T1(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( T1a0(unv)*pow(zs(qsq, unv)-zs(0.0, unv),0.0) + T1a1(unv)*pow(zs(qsq, unv)-zs(0.0, unv),1.0)
                                          + T1a2(unv)*pow(zs(qsq, unv)-zs(0.0, unv),2.0) );}

double Bstophill_fferrfun::T2(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T2a0(unv)*pow(zs(qsq, unv)-zs(0.0, unv),0.0) + T2a1(unv)*pow(zs(qsq, unv)-zs(0.0, unv),1.0)
                                          + T2a2(unv)*pow(zs(qsq, unv)-zs(0.0, unv),2.0) );}

double Bstophill_fferrfun::T23(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T23a0(unv)*pow(zs(qsq, unv)-zs(0.0, unv),0.0) + T23a1(unv)*pow(zs(qsq, unv)-zs(0.0, unv),1.0)
                                          + T23a2(unv)*pow(zs(qsq, unv)-zs(0.0, unv),2.0) );}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
double BdtoKll_fferrfun::fz(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mSbs(),2.0))*( fza1(unv)*pow(z(qsq, unv)-z(0.0, unv),1.0) + fza2(unv)*pow(z(qsq, unv)-z(0.0, unv),2.0) );}

double BdtoKll_fferrfun::fT(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( fTa0(unv)*pow(z(qsq, unv)-z(0.0, unv),0.0) + fTa1(unv)*pow(z(qsq, unv)-z(0.0, unv),1.0)
                                          + fTa2(unv)*pow(z(qsq, unv)-z(0.0, unv),2.0) );}

double BdtoKll_fferrfun::fp(double qsq, double unv[]){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( fpa0(unv)*pow(z(qsq, unv)-z(0.0, unv),0.0) + fpa1(unv)*pow(z(qsq, unv)-z(0.0, unv),1.0)
                                          + fpa2(unv)*pow(z(qsq, unv)-z(0.0, unv),2.0) );}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





