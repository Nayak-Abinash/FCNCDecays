#ifndef myfun_h
#define myfun_h

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <string>
#include <iostream>
using namespace std;

class myfun {
public:
    myfun();
    double mnd(double mu, double sigma);
    int mylength(double lst[], int iter, double a, double b);
    const int max_dim; //maximum allowed matrix dimension
    double Cholesky_Decomposition(double matrix[][max_dim], int n, int i);

    };


#endif
