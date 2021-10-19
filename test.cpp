#include<iostream>
#include<cmath>
using namespace std;

class par{
public:
    par();
    double md();
private:
    double mdv;
};

par::par(){mdv=0.5;}

double par::md(){return mdv;}

int main(){
    par p;
    cout << pow(p.md(),2.0) << endl;
    return 0;
}
