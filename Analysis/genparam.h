#ifndef GENPARAM_H
#define GENPARAM_H

#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
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

TH1 * createHist_and_computePulseIntegral(int npts, double * x, double * y, double& pulse_integral, double& baseline_shift)
{
    double dt = x[1] - x[0];
    double halfDT = dt/2.;
    TH1 * h = new TH1D("","",npts,x[0]-halfDT,x[npts-1]+halfDT);
    // Evaluate baseline - needed for a correct integral computation
    unsigned int bsl_npts(0);
    baseline_shift = 0;
    bool zero_reached(false);
    while(!zero_reached) {
        baseline_shift += y[bsl_npts];
        if(x[bsl_npts]>=-1*dt) zero_reached = true;
        h->SetBinContent(bsl_npts+1,y[bsl_npts]);
        bsl_npts++;
    }
    if(bsl_npts) baseline_shift = baseline_shift/bsl_npts;
    else baseline_shift = 0;
    
    pulse_integral = 0;
    for(int i(bsl_npts); i < npts; i++ ) {
        pulse_integral += 0.5 * dt *(y[i]+y[i-1] - 2*baseline_shift);
        h->SetBinContent(i+1,y[i]);
    }
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

double evaluateBaselineShift(TGraphErrors * gr) {
    // Evaluate baseline - needed for a correct long tau fit
    double baseline_shift(0);
    double x0, x1, yy;
    gr->GetPoint(0,x0,yy);
    gr->GetPoint(0,x1,yy);
    double dt = x1-x0;
    unsigned int bsl_npts(0);
    bool zero_reached(false);
    while(!zero_reached) {
        double x, y;
        gr->GetPoint(bsl_npts,x,y);
        baseline_shift += y;
        if(x>=-1*dt) zero_reached = true;
        bsl_npts++;
    }
    if(bsl_npts) baseline_shift = baseline_shift/bsl_npts;
    else baseline_shift = 0;
    return baseline_shift;
}

TGraphErrors * average(TGraphErrors * gr, TGraph * gr2)
{
    static int nevts = 0;

    const int npts    = gr2->GetN();
    double * x2 = gr2->GetX();
    double * y2 = gr2->GetY();

    if(!gr)
    {
        nevts = 1;
        double err_x[npts];
        double err_y[npts];
        err_x[0] = 0.5*(x2[1]-x2[0]); err_y[0] = 0;
        for(unsigned int row = 1; row < npts; row++) {
			err_x[row] = 0.5*(x2[row]-x2[row-1]);
			err_y[row] = 0;
		}
        TGraphErrors * grout = new TGraphErrors(npts,x2,y2, err_x,err_y);
        grout->SetName(TString(gr2->GetName())+"ForFit");
        grout->SetMarkerStyle(7);
        grout->SetMarkerColor(1);
        grout->SetLineColor(1);
        return grout;
    }

    double avgx[npts];
    double avgy[npts];
    double avgx_err[npts];
    double avgy_err[npts];
    Double_t * time  = gr->GetX();
    Double_t * volts = gr->GetY();
    Double_t * time_err  = gr->GetEX();
    Double_t * volts_err = gr->GetEY();
    for (int row = 0; row < npts; row++) 
    {
        avgx[row] = (nevts*time[row] + x2[row])/(nevts+1);
        avgy[row] = (nevts*volts[row] + y2[row])/(nevts+1);
        
        avgx_err[row] = time_err[row];
        avgy_err[row] = TMath::Sqrt(nevts/(nevts+1))*TMath::Sqrt((( (nevts*(volts_err[row]*volts_err[row] + volts[row]*volts[row])) + y2[row]*y2[row] )/(nevts+1)) - avgy[row]*avgy[row]);
    }
    nevts++;

    TGraphErrors * grout = new TGraphErrors(npts,avgx,avgy,avgx_err,avgy_err);
    grout->SetName(gr->GetName());
    return grout;
}

TGraphErrors * average2Dhist(TH2 * h, double& baseline_shift) {
	unsigned int Nx = h->GetNbinsX();
	TGraphErrors * average = new TGraphErrors(Nx);
	for(unsigned int binx(0); binx<Nx; ++binx) {
		TH1 * projy = h->ProjectionY("projy",binx+1,binx+1);
		
		//~ average->SetPoint(binx, h->GetXaxis()->GetBinCenter(binx+1), projy->GetMean());
		//~ average->SetPointError(binx, 0.5*h->GetXaxis()->GetBinWidth(binx+1), projy->GetMeanError());
		
		int binmax = projy->GetMaximumBin();
		int fit_range = 2;
		int binlow(1), binhigh(projy->GetNbinsX());
		if(binmax-fit_range > 1) binlow = binmax-fit_range;
		if(binmax+fit_range < projy->GetNbinsX()) binhigh = binmax+fit_range;
		//~ for(unsigned int biny(binlow); biny<binhigh, ++biny) {
			
		//~ }
		TF1 * gfit = new TF1("gfit","gaus",projy->GetBinLowEdge(binlow),projy->GetBinLowEdge(binhigh)+projy->GetBinWidth(binhigh));
		projy->Fit("gfit","QRF");
		double pos = h->GetXaxis()->GetBinCenter(binx+1);
		double pos_err = 0.5*h->GetXaxis()->GetBinWidth(binx+1);
		double mean = gfit->GetParameter(1);
		double error = gfit->GetParError(1);
		//if(error>0.005) { TFile* fo = new TFile("test.root","UPDATE"); projy->Write(); fo->Close();}
		average->SetPoint(binx,pos,mean);
		average->SetPointError(binx,pos_err,error);
		delete gfit;
		
		delete projy;
	}
	average->SetMarkerStyle(7);
	//~ average->SetMarkerSize(0.5);
	average->SetMarkerColor(1);
	average->SetLineColor(1);
	
	// Evaluate baseline - needed for a correct long tau fit
    baseline_shift = evaluateBaselineShift(average);
	
	return average;
}

void printValues(double values[], int n, std::string option, ostream * out = 0) {
	cout << option << " : [";
    if(out) (*out) << option << " : [";
    for(unsigned int i(0); i<n; ++i) {
        if(i<n-1) cout << values[i] << ",";
        if(i==n-1) cout << values[i] << "]" << endl;
        if(out && i<n-1) (*out) << values[i] << ",";
        if(out && i==n-1) (*out) << values[i] << "]" << endl;
    }
	return;
}

#endif
