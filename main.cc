#include "smpar.h"
#include "ffpar.h"
#include "fffun.h"
#include "amp.h"
#include "obs.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;


int main(){
    smpar p;
    obs o;
    double qsq;
    cout << "Type a value for qsq:";
    cin >> qsq;
    cout << "Observable values at this point are:"
    << o.diffWidth(qsq,p.mmu()) << "," << o.FL(qsq,p.mmu()) << "," << o.AFB(qsq,p.mmu()) << endl;
    return 0;
}

