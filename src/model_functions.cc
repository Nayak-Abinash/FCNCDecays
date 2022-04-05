#include "model_functions.h"

model_functions::model_functions(){
    rndm = new TRandom3;
    rndm->SetSeed(0);}

double model_functions::mean_model(double lst[], size_t n){
    double sum=0.0;
    for (int i=0;i<n;i++)
        {
            sum = sum+lst[i];
        }
    return sum/n;}

double model_functions::sd_model(double lst[], size_t n){
    double sum=0.0;
    double mn=mean_model(lst,n);
    for (int i=0;i<n;i++)
        {
            sum = sum+pow(lst[i]-mn,2);
        }
    return sqrt(sum/n);}

double model_functions::cov_model(double lst1[], double lst2[], size_t n){
    double sum=0.0;
    double mn1=mean_model(lst1,n);
    double mn2=mean_model(lst2,n);
    for (int i=0;i<n;i++)
        {
            sum = sum + (lst1[i]-mn1)*(lst2[i]-mn2);
        }
    return sum/n;}


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
    return rndm->Gaus(mu, sigma);}

double model_functions::mnd_default(){
    return rndm->Gaus(0.0,1.0);}

double model_functions::mnd_cov(double muv[], double cholesky[][MXdm], double untnormlv[], int i){
    double sum = 0;
    //int sz_muv = sizeof(muv)/sizeof(muv[0]);
    for (int j=0; j < 19; j++)
        {
            sum = sum + cholesky[i][j]*untnormlv[j];
        }
    return muv[i] + sum;}

double model_functions::mnd_cov_ckm(double muv[], double cholesky[][mxdm], double untnormlv[], int i){
    double sum = 0;
    //int sz_muv = sizeof(muv)/sizeof(muv[0]);
    for (int j=0; j < 4; j++)
        {
            sum = sum + cholesky[i][j]*untnormlv[j+19];
        }
    return muv[i] + sum;}

double model_functions::mnd_cov_wcSM(double muv[], double cholesky[][MXdm], double untnormlv[], int i){
    double sum = 0;
    //int sz_muv = sizeof(muv)/sizeof(muv[0]);
    for (int j=0; j < 13; j++)
        {
            sum = sum + cholesky[i][j]*untnormlv[j+23];
        }
    return muv[i] + sum;}

double model_functions::mnd_cov_wcNP(double muv[], double cholesky[][MXdm], double untnormlv[], int i){
    double sum = 0;
    //int sz_muv = sizeof(muv)/sizeof(muv[0]);
    for (int j=0; j < 18; j++)
        {
            sum = sum + cholesky[i][j]*untnormlv[j+36];
        }
    return muv[i] + sum;}

double model_functions::mnd_cov_smpar(double muv[], double cholesky[][MXdm], double untnormlv[], int i){
    double sum = 0;
    //int sz_muv = sizeof(muv)/sizeof(muv[0]);
    for (int j=0; j < 18; j++)
        {
            sum = sum + cholesky[i][j]*untnormlv[j+54];
        }
    return muv[i] + sum;}













