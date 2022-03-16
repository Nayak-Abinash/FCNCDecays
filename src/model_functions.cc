#include "model_functions.h"

model_functions::model_functions(){
    rndm = new TRandom3;
    rndm->SetSeed(0);}


int model_functions::mylength(double lst[], int iter, double a, double b){
    double num=0;
    for(int i=0; i<iter ; ++i)
    {
        if(a <= lst[i] && lst[i] < b)
            num = num+1;
    }
    return num;}


double model_functions::mtrx_tp(double lpar[], double covm[][mxdm]){
    return lpar[0]*covm[0][0]*lpar[0] + lpar[0]*covm[0][1]*lpar[1] + lpar[0]*covm[0][2]*lpar[2]
            + lpar[1]*covm[1][0]*lpar[0] + lpar[1]*covm[1][1]*lpar[1] + lpar[1]*covm[1][2]*lpar[2]
            + lpar[2]*covm[2][0]*lpar[0] + lpar[2]*covm[2][1]*lpar[1] + lpar[2]*covm[2][2]*lpar[2];}


double model_functions::mnd(double mu, double sigma){
    //unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    //default_random_engine generator (seed);
    //normal_distribution<double> distribution(mu,sigma);
    //return distribution(generator);
    return rndm->Gaus(mu, sigma);}

double model_functions::mnd_default(){
    return rndm->Gaus(0.0,1.0);}

double model_functions::mnd_cov(double muv[], double cholesky[][MXdm], double untnormlv[], int i){
    double sum = 0;
    for (int j=0; j<19; j++)
        {
            sum = sum + cholesky[i][j]*untnormlv[j];
        }
    return muv[i] + sum;}

double model_functions::mnd_cov_ckm(double muv[], double cholesky[][mxdm], double untnormlv[], int i){
    double sum = 0;
    for (int j=0; j<4; j++)
        {
            sum = sum + cholesky[i][j]*untnormlv[19+j];
        }
    return muv[i] + sum;}



