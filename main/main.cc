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
    BdtoKstrll_obs o1; BdtoKstrll_obserr eo1;
    double qsq; string str;
    cout << "Type a value for qsq:";
    getline(cin,str);
    stringstream(str) >> qsq;
    //double cval = o1.AFB(qsq,o1.mmu());
    double cval = o1.V(qsq);
    double sdval = o1.ErV(qsq);
    cout << "(" << cval << "," << sdval << ")" << endl;
///////////Error///////////////
////68.2%ConfidenceLevel (point counting from boundary values)////////////////
    double lw_range(-2.0), up_range(2.0);
    //TCanvas *c1= new TCanvas("c1", "c1", 800,600);//root
    TCanvas *c1 = new TCanvas();
    TH1D *hist = new TH1D("hist", "", 100, lw_range, up_range);//root

    int iter=10000;
    int sdcl = int(iter*15.9/100);
    double data[iter];

    for(int i=0; i<iter; i++)
        {
            //data[i]= eo1.AFB(qsq,o1.mmu());
            data[i] = eo1.V(qsq);
            hist->Fill(data[i]);//root
        }
    sort(data,data+iter);
    double llim(data[sdcl]);
    double ulim(data[iter-sdcl-1]);
    cout << "Observable at this qsq: " << cval << "(+" << ulim- cval << ",-" << cval-llim << ")" << endl;

    hist->GetXaxis()->SetTitle("MCDistribution");//root
    hist->GetYaxis()->SetTitle("Entries/bin");//root
    TF1 *fit = new TF1("fit","gaus", lw_range, up_range);//root
    //hist->SetStats(0);
    hist->Draw();
    hist->Fit("gaus","R");

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

