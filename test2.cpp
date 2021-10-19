#include<iostream>
#include<cmath>
using namespace std;

class par {
public:
    double pi,mB,mKst;
    par();
};

par :: par() {
    pi=M_PI;
    mB=5.27965,mKst=0.896;
}

int main(){
    par p;
    //p.mKst=0.5;
    cout << pow(p.mB-p.mKst,2.0) << endl;
    return 0;
}

