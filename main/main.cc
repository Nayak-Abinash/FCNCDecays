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
    BdtoKstrll_obs o1;
    Bstophill_obs o2;
    Btoll_obs o3;
    BdtoKll_obs o4;
    double qsq;
    cout << "Type a value for qsq:";
    cin >> qsq;
    cout << "BdtoKstrll Observable values at this point are:"
    << o1.diffWidth(qsq,p.mmu()) << "," << o1.FL(qsq,p.mmu()) << "," << o1.AFB(qsq,p.mmu()) << endl;
    cout << "Bstophill Observable values at this point are:"
    << o2.diffWidth(qsq,p.mmu()) << "," << o2.FL(qsq,p.mmu()) << "," << o2.AFB(qsq,p.mmu()) << endl;
    cout << "Btoll Observable values are:"
    << o3.BrTimeIntgratd(p.mBs(),p.ms(),p.mmu(),p.mmu()) << "," << o3.efftau(p.mBd(),p.md(),p.me(),p.me()) << "," << o3.ADeltaGammaf(p.mBs(),p.ms(),p.me(),p.mmu()) << endl;
    cout << "BdtoKll Observable values at this point are:"
    << o4.diffBrnch(qsq,p.mmu()) << "," << o4.diffFH(qsq,p.mmu()) << "," << o4.diffAFB(qsq,p.mmu()) << endl;
    return 0;
}

