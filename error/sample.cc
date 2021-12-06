#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <string>
using namespace std;

double a1=0.5, b1=0.4, c1=0.3;
double a2=0.44, b2=0.33, c2=0.22;
double ma1=0.5, mb1=0.4, mc1=0.3;
double ea1=0.05, eb1=0.04, ec1=0.03;
double ma2=0.44, mb2=0.33, mc2=0.22;
double ea2=0.044, eb2=0.033, ec2=0.022;

double F1(double qsq){return a1 + b1*qsq + c1*pow(qsq,2.0);}
double F2(double qsq){return a2 + b2*qsq + c2*pow(qsq,2.0);}
double obs(double qsq){return F1(qsq)/F2(qsq);}

double mnd(double mu, double sigma){
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    normal_distribution<double> distribution (mu,sigma);
    return distribution(generator);
}

int mylength(double lst[], int iter, double a, double b){
    double num=0;
    for(int i=0; i<iter ; ++i)
    {
        if(a <= lst[i] && lst[i] < b)
            num = num+1;
    }
    return num;
}

double simF1(double qsq){return mnd(ma1,ea1) + mnd(mb1,eb1)*qsq + mnd(mc1,ec1)*pow(qsq,2.0);}
double simF2(double qsq){return mnd(ma2,ea2) + mnd(mb2,eb2)*qsq + mnd(mc2,ec2)*pow(qsq,2.0);}
double simobs(double qsq){return simF1(qsq)/simF2(qsq);}

int main(){
    double qsq;
    cout << "Enter qsq:";
    cin >> qsq;
////EquallyLikely////////////////
    int iter=100000;
    double data[iter];
    double sum=0.0, nvar=0.0;
    for(int i=0; i<iter; i++)
        {
            data[i]=simobs(qsq);
            sum=sum+data[i];
        }
    double meanv = sum/iter;
    for(int i=0; i<iter; i++)
        {
            nvar = nvar + pow(meanv-data[i],2.0);
        }
    double stdv = sqrt(nvar/iter);
////Histogram/////////////////////
    double minv = *min_element(data,data+iter);
    double maxv = *max_element(data,data+iter);
    double nbin = 50;
    double init_val = minv;
    for (int i=0; i<nbin; ++i)
        {
            int npnt = mylength(data,iter,minv+(i)*(maxv-minv)/nbin,minv+(i+1)*(maxv-minv)/nbin);
    cout << "-" << "bin-" << i << ": ";
    cout << string(npnt/(iter/50),'*') << "\t \t" << npnt <<endl;
        }

    cout << "(" << maxv << "," << minv << ")" << endl;
////OutPut///////////////////////
    cout << "The Observable value at this qsq:" << obs(qsq) << endl;
    cout << "The Mean and Standard Deviation at this qsq:" << "(" << meanv << "," << stdv <<")" << endl;
    return 0;
    }
