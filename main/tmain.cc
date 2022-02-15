//ProjectFiles
#include "myfun.h"
#include "smpar.h"
#include "ffpar.h"
#include "fferrpar.h"
#include "fffun.h"
#include "fferrfun.h"
#include "amp.h"
#include "amperr.h"
#include "obs.h"
#include "obserr.h"
#include "obserrDF.h"
//ROOTFiles
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TLine.h"
#include "TH1.h"

int main()
{
    //Btoll_obs o1;
    BdtoKll_obs o2; BdtoKll_obserrDF eo2;
    double qsq;
    string str;
    cout << "Type a value for qsq:";
    getline(cin,str);
    stringstream(str) >> qsq;

    //cout << o1.BrTimeIntgratd(o1.mBs(),o1.ms(),o1.mmu(),o1.mmu()) << endl;
    //cout << o1.efftau(o1.mBs(),o1.ms(),o1.mmu(),o1.mmu()) << endl;
    //cout << o1.ADeltaGammaf(o1.mBs(),o1.ms(),o1.mmu(),o1.mmu()) << endl;
    cout << o2.diffBrnch(qsq,o2.mmu()) << endl;
    cout << o2.diffAFB(qsq,o2.mmu()) << endl;
    cout << o2.diffFH(qsq,o2.mmu()) << endl;
    cout << o2.fz(qsq) << ", " << o2.Erfz(qsq) << endl;
    cout << o2.fT(qsq) << ", " << o2.ErfT(qsq) << endl;
    cout << o2.fp(qsq) << ", " << o2.Erfp(qsq) << endl;
    cout << o2.alNP(qsq,o2.mmu()) << ", " << o2.blNP(qsq,o2.mmu()) << ", " << o2.clNP(qsq,o2.mmu()) << endl;
    cout << eo2.Er_diffWidth(qsq,o2.mmu()) << endl;
    return 0;
}

