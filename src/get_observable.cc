#include "get_observable.h"

//get_obs::get_obs(){}

void get_obs::obsval(string s, double qsq){
    //B->Vll
    BdtoKstrll_obs o1; Bstophill_obs o2;
    //B->Pll
    BdtoKll_obs o3;
    //B->ll
    Btoll_obs o4;
    //string input
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->Kll(SM) #######//
    //B0->Kmumu
    if(s=="AFB(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffAFB(qsq,o3.mmu()) << endl;}
    else if(s=="FH(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffFH(qsq,o3.mmu()) << endl;}
    else if(s=="Rmue(B0->Kll)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.mmu())/o3.diffBrnch(qsq,o3.me()) << endl;}
    else if(s=="Rtaumu(B0->Kll)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.mtau())/o3.diffBrnch(qsq,o3.mmu()) << endl;}
    else if(s=="dBR/dq2(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.mmu()) << endl;}
    //B0->Ktautau
    else if(s=="AFB(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffAFB(qsq,o3.mtau()) << endl;}
    else if(s=="FH(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffFH(qsq,o3.mtau()) << endl;}
    else if(s=="dBR/dq2(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.mtau()) << endl;}
    //B0->Kee
    else if(s=="AFB(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffAFB(qsq,o3.me()) << endl;}
    else if(s=="FH(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffFH(qsq,o3.me()) << endl;}
    else if(s=="dBR/dq2(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.me()) << endl;}


//####### B0->ll(SM) #######//
    //B0->mumu
    else if(s=="BR(B0->mumu)")
        cout << o4.BrTimeIntgratd(o4.mBd(),o4.md(),o4.mmu(),o4.mmu()) << endl;
    //B0->tautau
    else if(s=="BR(B0->tautau)")
        cout << o4.BrTimeIntgratd(o4.mBd(),o4.md(),o4.mtau(),o4.mtau()) << endl;
    //B0->ee
    else if(s=="BR(B0->ee)")
        cout << o4.BrTimeIntgratd(o4.mBd(),o4.md(),o4.me(),o4.me()) << endl;
//####### Bs->ll(SM) #######//
    //Bs->mumu
    else if(s=="ADeltaGamma(Bs->mumu)")
        cout << o4.ADeltaGammaf(o4.mBs(),o4.ms(),o4.mmu(),o4.mmu()) << endl;
    else if(s=="BR(Bs->mumu)")
        cout << o4.BrTimeIntgratd(o4.mBs(),o4.ms(),o4.mmu(),o4.mmu()) << endl;
    else if(s=="tau_mumu")
        cout << o4.efftau(o4.mBs(),o4.ms(),o4.mmu(),o4.mmu()) << endl;
    //Bs->tautau
    else if(s=="ADeltaGamma(Bs->tautau)")
        cout << o4.ADeltaGammaf(o4.mBs(),o4.ms(),o4.mtau(),o4.mtau()) << endl;
    else if(s=="BR(Bs->tautau)")
        cout << o4.BrTimeIntgratd(o4.mBs(),o4.ms(),o4.mtau(),o4.mtau()) << endl;
    else if(s=="tau_tautau")
        cout << o4.efftau(o4.mBs(),o4.ms(),o4.mtau(),o4.mtau()) << endl;
    //Bs->ee
    else if(s=="ADeltaGamma(Bs->ee)")
        cout << o4.ADeltaGammaf(o4.mBs(),o4.ms(),o4.me(),o4.me()) << endl;
    else if(s=="BR(Bs->ee)")
        cout << o4.BrTimeIntgratd(o4.mBs(),o4.ms(),o4.me(),o4.me()) << endl;
    else if(s=="tau_ee")
        cout << o4.efftau(o4.mBs(),o4.ms(),o4.me(),o4.me()) << endl;

//Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}
