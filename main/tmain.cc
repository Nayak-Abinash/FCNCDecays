//ProjectFiles
#include "obs.h"
#include "obserr.h"

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
    return 0;
}

