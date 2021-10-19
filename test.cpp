#include<iostream>
#include<cmath>
using namespace std;

class par {
public:
    double pi,mB,mKstr();
    par();
};

par :: par() {
    pi=M_PI;
    mB=5.27965;
}
double par :: mKstr() {
    return 0.896;
}

int main(){
    par p;
    cout << pow(p.mB-p.mKstr(),2.0) << endl;
    return 0;
}
