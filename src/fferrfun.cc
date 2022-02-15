#include "myfun.h"
#include "smpar.h"
#include "fferrpar.h"
#include "fferrfun.h"


//Bd->Kstr,ll
double BdtoKstrll_fferrfun::V(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( Va0()*pow(zd(qsq)-zd(0.0),0.0) + Va1()*pow(zd(qsq)-zd(0.0),1.0) + Va2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fferrfun::A0(double qsq){
    return 1.0/(1.0-qsq/pow(mPSbs(),2.0))*( A0a0()*pow(zd(qsq)-zd(0.0),0.0) + A0a1()*pow(zd(qsq)-zd(0.0),1.0) + A0a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fferrfun::A1(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A1a0()*pow(zd(qsq)-zd(0.0),0.0) + A1a1()*pow(zd(qsq)-zd(0.0),1.0) + A1a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fferrfun::A12(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A12a0()*pow(zd(qsq)-zd(0.0),0.0) + A12a1()*pow(zd(qsq)-zd(0.0),1.0) + A12a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fferrfun::T1(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( T1a0()*pow(zd(qsq)-zd(0.0),0.0) + T1a1()*pow(zd(qsq)-zd(0.0),1.0) + T1a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fferrfun::T2(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T2a0()*pow(zd(qsq)-zd(0.0),0.0) + T2a1()*pow(zd(qsq)-zd(0.0),1.0) + T2a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fferrfun::T23(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T23a0()*pow(zd(qsq)-zd(0.0),0.0) + T23a1()*pow(zd(qsq)-zd(0.0),1.0) + T23a2()*pow(zd(qsq)-zd(0.0),2.0) );}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll
double Bstophill_fferrfun::V(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( Va0()*pow(zs(qsq)-zs(0.0),0.0) + Va1()*pow(zs(qsq)-zs(0.0),1.0) + Va2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fferrfun::A0(double qsq){
    return 1.0/(1.0-qsq/pow(mPSbs(),2.0))*( A0a0()*pow(zs(qsq)-zs(0.0),0.0) + A0a1()*pow(zs(qsq)-zs(0.0),1.0) + A0a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fferrfun::A1(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A1a0()*pow(zs(qsq)-zs(0.0),0.0) + A1a1()*pow(zs(qsq)-zs(0.0),1.0) + A1a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fferrfun::A12(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A12a0()*pow(zs(qsq)-zs(0.0),0.0) + A12a1()*pow(zs(qsq)-zs(0.0),1.0) + A12a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fferrfun::T1(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( T1a0()*pow(zs(qsq)-zs(0.0),0.0) + T1a1()*pow(zs(qsq)-zs(0.0),1.0) + T1a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fferrfun::T2(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T2a0()*pow(zs(qsq)-zs(0.0),0.0) + T2a1()*pow(zs(qsq)-zs(0.0),1.0) + T2a2()*pow(zs(qsq)-zs(0.0),2.0) );}

double Bstophill_fferrfun::T23(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T23a0()*pow(zs(qsq)-zs(0.0),0.0) + T23a1()*pow(zs(qsq)-zs(0.0),1.0) + T23a2()*pow(zs(qsq)-zs(0.0),2.0) );}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
double BdtoKll_fferrfun::fz(double qsq){
    return 1.0/(1.0-qsq/pow(mSbs(),2.0))*( fza1()*pow(z(qsq)-z(0.0),1.0) + fza2()*pow(z(qsq)-z(0.0),2.0) );}

double BdtoKll_fferrfun::fT(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( fTa0()*pow(z(qsq)-z(0.0),0.0) + fTa1()*pow(z(qsq)-z(0.0),1.0) + fTa2()*pow(z(qsq)-z(0.0),2.0) );}

double BdtoKll_fferrfun::fp(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( fpa0()*pow(z(qsq)-z(0.0),0.0) + fpa1()*pow(z(qsq)-z(0.0),1.0) + fpa2()*pow(z(qsq)-z(0.0),2.0) );}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





