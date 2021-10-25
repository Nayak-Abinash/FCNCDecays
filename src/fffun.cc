#include "smpar.h"
#include "ffpar.h"
#include "fffun.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

//Bd->Kstr,ll
double BdtoKstrll_fffun::V(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( Va0()*pow(zd(qsq)-zd(0.0),0.0) + Va1()*pow(zd(qsq)-zd(0.0),1.0) + Va2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fffun::A0(double qsq){
    return 1.0/(1.0-qsq/pow(mPSbs(),2.0))*( A0a0()*pow(zd(qsq)-zd(0.0),0.0) + A0a1()*pow(zd(qsq)-zd(0.0),1.0) + A0a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fffun::A1(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A1a0()*pow(zd(qsq)-zd(0.0),0.0) + A1a1()*pow(zd(qsq)-zd(0.0),1.0) + A1a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fffun::A12(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( A12a0()*pow(zd(qsq)-zd(0.0),0.0) + A12a1()*pow(zd(qsq)-zd(0.0),1.0) + A12a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fffun::T1(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs(),2.0))*( T1a0()*pow(zd(qsq)-zd(0.0),0.0) + T1a1()*pow(zd(qsq)-zd(0.0),1.0) + T1a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fffun::T2(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T2a0()*pow(zd(qsq)-zd(0.0),0.0) + T2a1()*pow(zd(qsq)-zd(0.0),1.0) + T2a2()*pow(zd(qsq)-zd(0.0),2.0) );}

double BdtoKstrll_fffun::T23(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs(),2.0))*( T23a0()*pow(zd(qsq)-zd(0.0),0.0) + T23a1()*pow(zd(qsq)-zd(0.0),1.0) + T23a2()*pow(zd(qsq)-zd(0.0),2.0) );}

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

//Bd->K,ll
double BdtoKll_fffun::lcsrfp(double qsq){
    return 1.0/(1.0-qsq/pow(lcsr_mVbs(),2.0))*( cfp0()*pow(z(qsq)-z(0.0),0.0) + cfp1()*pow(z(qsq)-z(0.0),1.0) + cfp2()*pow(z(qsq)-z(0.0),2.0) );}

double BdtoKll_fffun::lcsrfz(double qsq){
    return 1.0/(1.0-qsq/pow(lcsr_mSbs(),2.0))*( cfz0()*pow(z(qsq)-z(0.0),0.0) + cfz1()*pow(z(qsq)-z(0.0),1.0) + cfz2()*pow(z(qsq)-z(0.0),2.0) );}

double BdtoKll_fffun::lcsrft(double qsq){
    return 1.0/(1.0-qsq/pow(lcsr_mVbs(),2.0))*( cft0()*pow(z(qsq)-z(0.0),0.0) + cft1()*pow(z(qsq)-z(0.0),1.0) + cft2()*pow(z(qsq)-z(0.0),2.0) );}

double BdtoKll_fffun::latfp(double qsq){
    return 1.0/(1.0-qsq/pow(lat_mVbs(),2.0))*(b0p() + b1p()*(z(qsq)-1.0/3.0*pow(z(qsq),3.0)) + b2p()*(pow(z(qsq),2.0) + 2.0/3.0*pow(z(qsq),3.0)));}

double BdtoKll_fffun::latfz(double qsq){
    return 1.0/(1.0-qsq/pow(lat_mSbs(),2.0))*(b0z() + b1z()*z(qsq) + b2z()*pow(z(qsq),2.0));}

double BdtoKll_fffun::latft(double qsq){
    return 1.0/(1.0-qsq/pow(lat_mVbs(),2.0))*(b0t() + b1t()*(z(qsq)-1.0/3.0*pow(z(qsq),3.0)) + b2t()*(pow(z(qsq),2.0) + 2.0/3.0*pow(z(qsq),3.0)));}

double BdtoKll_fffun::fp(double qsq){
    if (qsq<4.0*pow(mc(),2.0)){
        return lcsrfp(qsq); }
    else {
        return latfp(qsq); } }

double BdtoKll_fffun::fz(double qsq){
    if (qsq<4.0*pow(mc(),2.0)){
        return lcsrfz(qsq); }
    else {
        return latfz(qsq); } }

double BdtoKll_fffun::ft(double qsq){
    if (qsq<4.0*pow(mc(),2.0)){
        return lcsrft(qsq); }
    else {
        return latft(qsq); } }










