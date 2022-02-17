//ProjectFiles
#include "obs.h"
#include "obserr.h"
//#include "obserrDF.h"

int main()
{
    BdtoKstrll_obs o1; BdtoKstrll_obserr eo1;
    double qsq; string str;
    cout << "Type a value for qsq:";
    getline(cin,str);
    stringstream(str) >> qsq;
    double cval = o1.AFB(qsq,o1.mmu());
    cout << "(" << cval << "," << /*sdval <<*/ ")" << endl;
///////////Error///////////////
////68.2%ConfidenceLevel (point counting from boundary values)////////////////
    double lw_range;
    if(cval > 0) lw_range = cval-1.0*cval;
    else if(cval < 0) lw_range = cval+1.0*cval;
    else lw_range = -1.0;
    double up_range;
    if(cval > 0) up_range = cval+1.0*cval;
    else if(cval < 0) up_range = cval-1.0*cval;
    else up_range = 1.0;
    TCanvas *c1 = new TCanvas();
    TH1D *hist = new TH1D("hist", "", 100, lw_range, up_range);

    int iter=10000;
    int sdcl = int(iter*15.9/100);
    double data[iter];
    for(int i=0; i<iter; i++)
        {
            data[i]= eo1.AFB(qsq,o1.mmu());
            hist->Fill(data[i]);
        }
    sort(data,data+iter);
    double llim(data[sdcl]);
    double ulim(data[iter-sdcl-1]);
    cout << "Observable at this qsq: " << cval << "(+" << ulim- cval << ",-" << cval-llim << ")" << endl;

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

