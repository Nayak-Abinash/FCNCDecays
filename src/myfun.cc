#include "myfun.h"
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <string>
#include <iostream>
#include <bits/stdc++.h>
using namespace std;

myfun::myfun(){
    max_dim = 5;}

double myfun::mnd(double mu, double sigma){
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    normal_distribution<double> distribution (mu,sigma);
    return distribution(generator);}

int myfun::mylength(double lst[], int iter, double a, double b){
    double num=0;
    for(int i=0; i<iter ; ++i)
    {
        if(a <= lst[i] && lst[i] < b)
            num = num+1;
    }
    return num;}

double myfun::Cholesky_Decomposition(double matrix[][max_dim], int n, int i){
    double lower[n][n];
    memset(lower, 0, sizeof(lower));
    for (int i = 0; i < n; i++)
        {
        for (int j = 0; j <= i; j++)
            {
            double sum = 0;
            if (j == i)
                {
                for (int k = 0; k < j; k++)
                    sum += pow(lower[j][k], 2);
                lower[j][j] = sqrt(matrix[j][j] - sum);
                }
            else
                {
                    for (int k = 0; k < j; k++)
                    sum += (lower[i][k] * lower[j][k]);
                    lower[i][j] = (matrix[i][j] - sum) / lower[j][j];
                }
            }
        }
    return lower[i-1][i-1];}









