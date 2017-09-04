#ifndef GENPARAM_H
#define GENPARAM_H

#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <iostream>
#include <string>

using namespace std;

const double ns = 1e-9;

struct globalArgs_t
{
    const char* data_folder;                /* -d option */
    const char* arg_pathToSetupFile;        /* -S option */
    TString res_folder;                     /* -o option */
    bool save_all;                          /* -a  */
    TString input;                          /* -I root  */
    double fixed_thr;                       /* -T root  */             
    int only_Vbd;                     		/* -V root  */             

} globalArgs {" "," "," ",false,"root",-1,0};

struct delayedPulse
{
	double time;
	double volt;
	string type;	// can be "DeXT" or "AP"
	bool operator<(const delayedPulse& a) const
    {
        return time < a.time;
    }
};

struct timeFitResult {
	double tau;
	double tau_error;
	double Nsig;
	double Nsig_error;
	double Nbkg;
	double Nbkg_error;
	double chi2;
	//TString name;
};

TH1 * convertGrToH(TGraph * gr)
{
    Double_t * xx  = gr->GetX();
    Double_t * yy = gr->GetY();
    double halfDT = (xx[1] - xx[0])/2.; 
    //TH1 * h = new TH1F(TString(gr->GetName())+"H",gr->GetTitle(),
    TH1 * h = new TH1F("","",
        gr->GetN(),xx[0]-halfDT,xx[gr->GetN()-1]+halfDT);
    for(int i = 1; i < gr->GetN(); i++ ) h->SetBinContent(i,yy[i]);

    return h;
}

TH1 * createHist(int npts, double * x, double * y)
{
    double halfDT = (x[1] - x[0])/2.;
    TH1 * h = new TH1D("","",npts,x[0]-halfDT,x[npts-1]+halfDT);
    for(int i(0); i < npts; i++ ) h->SetBinContent(i+1,y[i]);
    return h;
}

double computeIntegral(TGraph * g, double threshold) {
    unsigned int N = g->GetN();
    double * xx = g->GetX();
    double * yy = g->GetY();
    double integral(0);
    for(unsigned int pt(1); pt<N; ++pt) {
        if(yy[pt]>threshold) {
            //integral += 0.5 * (xx[pt-1]*yy[pt] - xx[pt]*yy[pt-1]);
            //integral += 0.5 * (xx[pt]-xx[pt-1])*(yy[pt]-yy[pt-1]);
            integral += 0.5 * (xx[pt]-xx[pt-1])*(yy[pt]+yy[pt-1]);
        }
        //integral += 0.5 * (xx[pt-1]*yy[pt] - xx[pt]*yy[pt-1]);
        //integral += 0.5 * (xx[pt]-xx[pt-1])*(yy[pt]+yy[pt-1]);
        //integral += yy[pt];
    }
    return integral;
    //return abs(integral);
    //return g->Integral();
}

TGraph * average(TGraph * gr, TGraph * gr2)
{
    static int nevts = 0;

    int npts    = gr2->GetN();
    double * x2 = gr2->GetX();
    double * y2 = gr2->GetY();

    if(!gr)
    {
        nevts = 1;
        TGraph * grout = new TGraph(npts,x2,y2);
        grout->SetName(TString(gr2->GetName())+"ForFit");
        return grout;
    }

    double avgx[npts];
    double avgy[npts];
    Double_t * time  = gr->GetX();
    Double_t * volts = gr->GetY();
    for (int row = 0; row < npts; row++) 
    {
        avgx[row] = (nevts*time[row] + x2[row])/(nevts+1);
        avgy[row] = (nevts*volts[row] + y2[row])/(nevts+1);
    }
    nevts++;

    TGraph * grout = new TGraph(npts,avgx,avgy);
    grout->SetName(gr->GetName());
    return gr;
}

TGraphErrors * average2Dhist(TH2 * h) {
	unsigned int Nx = h->GetNbinsX();
	TGraphErrors * average = new TGraphErrors(Nx);
	for(unsigned int binx(0); binx<Nx; ++binx) {
		TH1 * projy = h->ProjectionY("projy",binx+1,binx+1);
		average->SetPoint(binx, h->GetXaxis()->GetBinCenter(binx+1), projy->GetMean());
		average->SetPointError(binx, 0.5*h->GetXaxis()->GetBinWidth(binx+1), projy->GetRMS());
		delete projy;
	}
	average->SetMarkerStyle(8);
	average->SetMarkerSize(0.5);
	average->SetMarkerColor(8);
	
	return average;
}

#endif
