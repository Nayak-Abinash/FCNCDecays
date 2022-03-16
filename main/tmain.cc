//ProjectFiles
#include "observables.h"
#include "observables_mc.h"

int main()
{
    Btoll_obs o1;
    //cout << o1.BrTimeIntgratd(o1.mBs(),o1.ms(),o1.mmu(),o1.mmu()) << endl;
    cout << o1.ReVud() << "+ (" << o1.ImVud() << ")" << '\t' << o1.ReVus() << "+ (" << o1.ImVus() << ")" <<
            '\t' << o1.ReVub() << "+ (" << o1.ImVub() << ")" << endl;
    cout << o1.ReVcd() << "+ (" << o1.ImVcd() << ")" << '\t' << o1.ReVcs() << "+ (" << o1.ImVcs() << ")" <<
            '\t' << o1.ReVcb() << "+ (" << o1.ImVcb() << ")" << endl;
    cout << o1.ReVtd() << "+ (" << o1.ImVtd() << ")" << '\t' << o1.ReVts() << "+ (" << o1.ImVts() << ")" <<
            '\t' << o1.ReVtb() << "+ (" << o1.ImVtb() << ")" << endl;
    cout << "/////////////////////////////////////////////" << endl;
    ckm_elements_mc o2;
    double unv[] = {o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(),
                                    o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(),
                                    o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(),
                                    o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(), o1.mnd_default(),
                                    o1.mnd_default()};
    cout << o2.ReVud(unv) << "+ (" << o2.ImVud(unv) << ")" << '\t' << o2.ReVus(unv) << "+ (" << o2.ImVus(unv) << ")" <<
            '\t' << o2.ReVub(unv) << "+ (" << o2.ImVub(unv) << ")" << endl;
    cout << o2.ReVcd(unv) << "+ (" << o2.ImVcd(unv) << ")" << '\t' << o2.ReVcs(unv) << "+ (" << o2.ImVcs(unv) << ")" <<
            '\t' << o2.ReVcb(unv) << "+ (" << o2.ImVcb(unv) << ")" << endl;
    cout << o2.ReVtd(unv) << "+ (" << o2.ImVtd(unv) << ")" << '\t' << o2.ReVts(unv) << "+ (" << o2.ImVts(unv) << ")" <<
            '\t' << o2.ReVtb(unv) << "+ (" << o2.ImVtb(unv) << ")" << endl;
    cout << "/////////////////////////////////////////////" << endl;
    cout << o1.absVtbVtdStr() << '\t' << o2.absVtbVtdStr(unv) << endl;
    cout << o1.absVtbVtsStr() << '\t' << o2.absVtbVtsStr(unv) << endl;
    return 0;
}

