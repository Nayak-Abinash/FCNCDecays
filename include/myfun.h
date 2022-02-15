#ifndef myfun_h
#define myfun_h

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <bits/stdc++.h>
#define mxdm 5
#define MXdm 21
using namespace std;

class myfun
{
public:
    int mylength(double lst[], int iter, double a, double b);
    double mtrx_tp(double lpar[], double covm[][mxdm]);
    double mnd(double mu, double sigma);
    double mnd_cov(double mu[], double cholesky[][MXdm], int i);
};


#endif
