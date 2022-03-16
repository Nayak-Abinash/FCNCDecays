//ProjectFiles
#include "observables.h"
#include "observables_mc.h"

int main()
{
    BdtoKstrll_obs o1; BdtoKstrll_obserr eo1;
    double qsq; string str;
    cout << "Type a value for qsq:";
    getline(cin,str); stringstream(str) >> qsq;
    double cval = o1.FL(qsq,o1.mmu());
    cout << "(" << cval << "," << /*sdval <<*/ ")" << endl;
///////////Error///////////////
////68.2%ConfidenceLevel (point counting from boundary values)////////////////
    double lw_range, up_range, pe(0.3);
    if(cval > 0) lw_range = cval-pe*cval; else if(cval < 0) lw_range = cval+pe*cval; else lw_range = -1.0;
    if(cval > 0) up_range = cval+pe*cval; else if(cval < 0) up_range = cval-pe*cval; else up_range = 1.0;
    TCanvas *c1 = new TCanvas();
    TH1D *hist = new TH1D("hist", "", 100, lw_range, up_range);
//Boundary counting
    int iter(10000), sdcl(int(iter*15.9/100)); double data[iter];
    for(int i=0; i<iter; i++)
        {
            double unv[] = {eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(),
                                    eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(),
                                    eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default(),
                                    eo1.mnd_default(), eo1.mnd_default(), eo1.mnd_default()};
            data[i]= eo1.FL(qsq,o1.mmu(),unv);
            hist->Fill(data[i]);
        }
    sort(data,data+iter);
//15.9% & 84.1% points
    double llim(data[sdcl]), ulim(data[iter-sdcl-1]);
    cout << "Observable at this qsq: " << cval << "(+" << ulim- cval << ",-" << cval-llim << ")" << endl;
//plotting
    hist->GetXaxis()->SetTitle("MCDistribution");
    hist->GetYaxis()->SetTitle("Entries/bin");
    TF1 *fit = new TF1("fit","gaus", lw_range, up_range);
    //hist->SetStats(0);
    hist->Draw();
    hist->Fit("fit","R");

    TLine *ln1 = new TLine(llim,0,llim,300);
    TLine *ln2 = new TLine(ulim,0,ulim,300);
    ln1->SetLineColor(kBlue); ln1->Draw();
    ln2->SetLineColor(kBlue); ln2->Draw();
    TLegend *leg = new TLegend(0.1,0.7,0.3,0.95);
    //leg->SetHeader("Plot Legends", "C");
    leg->AddEntry(hist, "MC Simulate Data", "l");
    leg->AddEntry(fit, "Gaussian Fit", "l");
    leg->Draw();
    c1->SaveAs("obsplot.pdf");

    return 0;
}

