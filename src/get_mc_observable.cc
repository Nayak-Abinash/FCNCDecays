#include "get_mc_observable.h"

void get_obserr::obsval(string s, double qsq, double smwc[], double npwc[]){
    //B->Vll
    BdtoKstrll_obserr eo1; Bstophill_obserr eo2;
    //B->Pll
    BdtoKll_obserr eo3;
    //B->ll
    Btoll_obserr eo4;

    int iter=100;
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->Kll #######//
    //B0->Kmumu
    if(s=="AFB(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffAFB(qsq,eo3.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="FH(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffFH(qsq,eo3.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="Rmue(B0->Kll)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.mmu(unv),smwc,npwc,unv)/eo3.diffBrnch(qsq,eo3.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="Rtaumu(B0->Kll)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.mtau(unv),smwc,npwc,unv)/eo3.diffBrnch(qsq,eo3.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="dBR/dq2(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    //B0->Ktautau
    else if(s=="AFB(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffAFB(qsq,eo3.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="FH(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffFH(qsq,eo3.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="dBR/dq2(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    //B0->Kee
    else if(s=="AFB(B0->Kee)")
        {cout << "Enter a value of qsq:";
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffAFB(qsq,eo3.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="FH(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffFH(qsq,eo3.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="dBR/dq2(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}



//####### B0->ll #######//
    //B0->mumu
    else if(s=="BR(B0->mumu)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBd(unv),eo4.md(unv),eo4.mmu(unv),eo4.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    //B0->tautau
    else if(s=="BR(B0->tautau)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBd(unv),eo4.md(unv),eo4.mtau(unv),eo4.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    //B0->ee
    else if(s=="BR(B0->ee)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBd(unv),eo4.md(unv),eo4.me(unv),eo4.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}



//####### Bs->ll #######//
    //Bs->mumu
    else if(s=="ADeltaGamma(Bs->mumu)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.ADeltaGammaf(eo4.mBs(unv),eo4.ms(unv),eo4.mmu(unv),eo4.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="BR(Bs->mumu)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBs(unv),eo4.ms(unv),eo4.mmu(unv),eo4.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="tau_mumu")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.efftau(eo4.mBs(unv),eo4.ms(unv),eo4.mmu(unv),eo4.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    //Bs->tautau
    else if(s=="ADeltaGamma(Bs->tautau)")
         {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.ADeltaGammaf(eo4.mBs(unv),eo4.ms(unv),eo4.mtau(unv),eo4.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="BR(Bs->tautau)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBs(unv),eo4.ms(unv),eo4.mtau(unv),eo4.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="tau_tautau")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.efftau(eo4.mBs(unv),eo4.ms(unv),eo4.mtau(unv),eo4.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    //Bs->ee
    else if(s=="ADeltaGamma(Bs->ee)")
         {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.ADeltaGammaf(eo4.mBs(unv),eo4.ms(unv),eo4.me(unv),eo4.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="BR(Bs->ee)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBs(unv),eo4.ms(unv),eo4.me(unv),eo4.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="tau_ee")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.efftau(eo4.mBs(unv),eo4.ms(unv),eo4.me(unv),eo4.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}



//Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}
