#ifndef BREAKDOWN_VOLTAGE_H
#define BREAKDOWN_VOLTAGE_H

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooPolynomial.h"
#include "RooMsgService.h"
#include "RooFFTConvPdf.h"
#include "RooMinuit.h"

#include <TCanvas.h>
#include <iostream>
#include <TFitResultPtr.h>
#include <TPaveText.h>
#include <fstream>

#include "genparam.h"

double fitBreakdownVoltage(TGraph * Vbias_ver, TFile * file=NULL, ofstream * values=NULL)
{
    TCanvas * ca = new TCanvas("Voltage Breakdown calculation","Voltage Breakdown calculation",100,100,900,700);
    Vbias_ver->SetTitle("Voltage Breakdown calculation");
    Vbias_ver->GetYaxis()->SetTitle("Bias volatge [V]");
    //Vbias_ver->GetYaxis()->SetTitleOffset(1.2);
    Vbias_ver->GetXaxis()->SetTitle("Mean peak amplitude [V]");
    Vbias_ver->Draw("AP*");
    ca->SetGrid();
    
    TFitResultPtr fit = Vbias_ver->Fit("pol1","S");
    Double_t VBD = fit->Value(0);
    Double_t m = fit->Value(1);
    
    std::cout << "Breakdown voltage is: " << VBD << std::endl;
    std::cout << "Amplitudes from fit:" << std::endl;
    double * volt = Vbias_ver->GetY();
    cout << "PeakAmp : [";
    if(values) (*values) << "PeakAmp : [";
    for(unsigned int i(0); i<Vbias_ver->GetN(); ++i) {
		if(i<Vbias_ver->GetN()-1) cout << 1000*(volt[i] - VBD)/m << ",";
		if(i==Vbias_ver->GetN()-1) cout << 1000*(volt[i] - VBD)/m << "]" << endl;
		if(values && i<Vbias_ver->GetN()-1) (*values) << 1000*(volt[i] - VBD)/m << ",";
		if(values && i==Vbias_ver->GetN()-1) (*values) << 1000*(volt[i] - VBD)/m << "]" << endl;
	}
	cout << "VBD : " << VBD << endl;
	if(values) (*values) << "VBD : " << VBD << endl;
	cout << "VBD_error : " << fit->Error(0) << endl;
	if(values) (*values) << "VBD_error : " << fit->Error(0) << endl;

    TPaveText * pv = new TPaveText(0.2,0.65,0.5,0.75,"brNDC");
    pv->AddText(Form("V_{BD} = %2.2f#pm%2.2f V",VBD,fit->Error(0)));
    double inv_slope = 1e3/m;
    double inv_slope_err = inv_slope*(fit->Error(1)/m);
    pv->AddText(Form("Amp./#DeltaV = %2.2f#pm%2.2f mV/V",inv_slope,inv_slope_err));
    pv->SetFillColor(kWhite);
    pv->Draw();
    ca->Print(globalArgs.res_folder+"Breakdown_voltage.pdf");
    if(file) {
		file->cd("Vbd_determination");
		ca->Write();
		file->cd();
	}

    return VBD;
}

TF1 * initLongTauFunction(double pe, double baseline_shift) {
    // Init function used for long tau fit in different other functions
    
    TF1 * exp_longtau = new TF1("fit_longtau","[0]*exp(-(x-[1])/[2])+[3]",20*ns,150*ns);
    // amplitude
    exp_longtau->SetParameter(0, pe*0.2);
    exp_longtau->SetParLimits(0, 0.01*pe,1.*pe);
    // time offset
    exp_longtau->FixParameter(1, 0*ns);
    // tau
    exp_longtau->SetParameter(2, 80*ns);
    exp_longtau->SetParLimits(2, 10*ns,200*ns);
    //~ exp_longtau->FixParameter(2, 85*ns);
    // baseline shift
    exp_longtau->FixParameter(3, baseline_shift);
    
    exp_longtau->SetLineColor(kRed);
    exp_longtau->SetLineWidth(2);
    
    return exp_longtau;    
}

TF1 * fitLongTau_TGraphErrors(TGraphErrors * cleanwaves, double * amp0, double * tau, double pe, const char * vol) 
{
    // Evaluate baseline - needed for a correct long tau fit
    double baseline_shift = evaluateBaselineShift(cleanwaves);
    
    TF1 * exp_longtau = initLongTauFunction(pe, baseline_shift); // no baseline shift implemented
    TString fitname = Form("%s_TGraphErrors_%s",exp_longtau->GetName(),vol);
    exp_longtau->SetName(fitname);
        
    cleanwaves->Fit(fitname,"","",20*ns,150.*ns); // Fit boundaries for the slow component of the pulse
    
    (*amp0) = exp_longtau->GetParameter(0);
    (*tau) = exp_longtau->GetParameter(2);
    std::cout << "Long tau (from average TGraphErrors) = " << *tau << std::endl;
    
    return exp_longtau;
}

TF1 * fitLongTau_TMultiGraph(TMultiGraph * cleanwaves, double * amp0, double * tau, double pe, const char * vol, TCanvas * c, TGraph * avg) 
{
    c->cd();
	avg->SetMarkerStyle(7);
    avg->Draw("P+");
    //cleanwaves->Draw("P+");	// already drawn with drawWave

    TF1 * exp_longtau = initLongTauFunction(pe, 0); // no baseline shift implemented
    TString fitname = Form("%s_TMultiGraph_%s",exp_longtau->GetName(),vol);
    exp_longtau->SetName(fitname);
    
    cout << "Using " << cleanwaves->GetListOfGraphs()->GetSize() << " curves for long tau fit." << endl;     
    cleanwaves->Fit(fitname,"","",20*ns,150*ns); // Fit boundaries for the slow component of the pulse
    (*amp0) = exp_longtau->GetParameter(0);
    (*tau) = exp_longtau->GetParameter(2);
    std::cout << "Long tau = " << *tau << std::endl;

    c->cd();

    TPaveText * pv = new TPaveText(0.6,0.65,0.75,0.74,"brNDC");
    pv->AddText(Form("#tau_{long} = %2.1fns",1e9*(*tau)));
    pv->SetFillColor(kWhite);
    pv->Draw("SAME");
    exp_longtau->Draw("SAME");

    return exp_longtau;
}

// second version for fitting the persistence 2D histogram
TF1 * fitLongTau_TH2(TH2 * h, double * amp0, double * tau, double pe, const char * vol, TGraphErrors * &averageGraph) 
{
    double baseline_shift(0);
	averageGraph = average2Dhist(h,baseline_shift);
    averageGraph->SetName(Form("average_clean_%s",vol));
    
    TF1 * exp_longtau = initLongTauFunction(pe, baseline_shift);
    TString fitname = Form("%s_TH2_%s",exp_longtau->GetName(),vol);
    exp_longtau->SetName(fitname);
    
    averageGraph->Fit(fitname,"","",20*ns,150.*ns); // Fit boundaries for the slow component of the pulse
    (*amp0) = exp_longtau->GetParameter(0);
    (*tau) = exp_longtau->GetParameter(2);
    std::cout << "Long tau (from 2D histogram) = " << *tau << std::endl;

    return exp_longtau;
}

    // Original setting was 4ns,180ns and 30ns,100ns for QA
TF1 * fitAPTau(TGraph * APtime, double amp0, double tau, double pe, const char * vol, double * recovery, double minTime=0, double maxTime=0) 
{
    // Fit parameters and limits to calculate AP recharge
    double min_for_fit(0), max_for_fit(0);
	min_for_fit = 20*ns;
	max_for_fit = 150*ns;
    if(minTime) min_for_fit = minTime;
    if(maxTime) max_for_fit = maxTime;
    
    // function for fit is the exponential increase from the pixel recovery + long component of the pulse (exponential decrease)
    TF1 * exp = new TF1(Form("exp_%s",vol),"[0]*(1 - exp(-(x-[5])/[1])) + [2]*exp(-(x-[4])/[3])",min_for_fit,max_for_fit); // draw curve range
    //exp->SetParameter(0,pe);
    //exp->SetParLimits(0,0.2*pe,10*pe);
    exp->FixParameter(0,pe);  
    exp->SetParameter(1,50*ns);
    exp->SetParLimits(1,5*ns,120*ns);
    exp->FixParameter(2,amp0);
    //exp->SetParLimits(2,0.,1.);
    exp->FixParameter(3,tau);
    exp->FixParameter(4,0*ns);	// found to be zero so we fix it
//    exp->FixParameter(5,3*ns);	// found to be zero so we fix it QA was 0
    exp->SetParameter(5,3*ns);
    exp->SetParLimits(5,0*ns,10*ns);
    
    exp->SetLineColor(kRed);

    APtime->Fit(Form("exp_%s",vol),"","",min_for_fit,max_for_fit);
    (*recovery) = exp->GetParameter(1);
    return exp;
}

TGraph * splitAPTau(TGraph * APtime, double amp0, double tau, double pe) 
{
    unsigned int N = APtime->GetN();
    TGraph * splitG = new TGraph(N);
    splitG->SetName(TString(APtime->GetName())+"_split");
    for(unsigned int pt(0); pt<N; ++pt) {
        double x(0), y(0);
        APtime->GetPoint(pt,x,y);
        //~ y = y - amp0*TMath::Exp((-1.0*x)/tau);
        y = 1 - (y - amp0*TMath::Exp((-1.0*x)/tau))/pe;
        splitG->SetPoint(pt,x,y);
    }
    return splitG;
}

void fitTimeDist(TH1 * timeDist, TCanvas * c, timeFitResult& new_fit, const char * type, double minTime=0, double maxTime=0) {
	if(!timeDist->GetEntries()) return;
	// timeDist is expected to contain 
	// - a exponential decrease (physical distribution of time of arrivals)
	// - a constant plateau due to DCR
    
    // Pre-fit to get starting parameters
    // Total fit range:
	double minFit = timeDist->GetBinLowEdge(timeDist->GetMaximumBin());
	double maxFit = timeDist->GetBinLowEdge(timeDist->FindLastBinAbove()) + timeDist->GetBinWidth(timeDist->FindLastBinAbove());
	double expFitMax = minFit + ((maxFit-minFit)/3.0);
	// Fits first with an exponential
	TF1 * expo1 = new TF1(Form("expo1_%s",timeDist->GetName()), "expo", minFit,expFitMax);
	timeDist->Fit(Form("expo1_%s",timeDist->GetName()),"QN","",minFit,expFitMax);
	// Uses the fit results to fit the entire distribution
	double startAmp = TMath::Exp(expo1->GetParameter(0));
	double startTau = 1.0/TMath::Abs(expo1->GetParameter(1));
    double startCst = 0;
	
	c->cd();
    timeDist->Draw("E0");
    timeDist->GetXaxis()->SetTitle("Time [s]");
    timeDist->GetYaxis()->SetTitle("dN_{"+TString(type)+"}/dt");
    c->SetLogy();
	// Total fit range:
    double min_for_fit(0), max_for_fit(0);
	min_for_fit = minFit;
	max_for_fit = maxFit;
    if(minTime) min_for_fit = minTime;
    if(maxTime) max_for_fit = maxTime;
	
	TF1 * expo_and_dcr = new TF1(Form("fit_%s",timeDist->GetName()), "[0]*exp(-x/[1]) + [2]", min_for_fit,max_for_fit);
    expo_and_dcr->SetLineColor(kRed);
	expo_and_dcr->SetParameter(0,startAmp);
    expo_and_dcr->SetParLimits(0,0,timeDist->GetEntries());
	expo_and_dcr->SetParameter(1,startTau);
    expo_and_dcr->SetParLimits(1,1*ns,120*ns);
	expo_and_dcr->SetParameter(2,startCst);
    expo_and_dcr->SetParLimits(2,0,timeDist->GetEntries());
	timeDist->Fit(Form("fit_%s",timeDist->GetName()),"Q+","",min_for_fit,max_for_fit);
    
    // *********** fit results ************
	double result_tau = expo_and_dcr->GetParameter(1);
    double result_tau_error = expo_and_dcr->GetParError(1);
	double result_Nsig = (result_tau*expo_and_dcr->GetParameter(0)*( TMath::Exp(-min_for_fit/result_tau) - TMath::Exp(-max_for_fit/result_tau) ));
    double result_Nsig_error = ( (result_Nsig/result_tau) + result_tau*expo_and_dcr->GetParameter(0)*( max_for_fit*TMath::Exp(-max_for_fit/result_tau) - min_for_fit*TMath::Exp(-min_for_fit/result_tau) ) ) * max_for_fit;
    result_Nsig = result_Nsig/timeDist->GetBinWidth(1);
    result_Nsig_error = result_Nsig_error/timeDist->GetBinWidth(1);
	double result_Nbkg = ((max_for_fit-min_for_fit)/timeDist->GetBinWidth(1))*expo_and_dcr->GetParameter(2);
    double result_Nbkg_error = (expo_and_dcr->GetParError(2)/expo_and_dcr->GetParameter(2))*result_Nbkg;
	
	new_fit.tau = result_tau;	new_fit.tau_error = result_tau_error;
	new_fit.Nsig = result_Nsig;	new_fit.Nsig_error = result_Nsig_error;
	new_fit.Nbkg = result_Nbkg;	new_fit.Nbkg_error = result_Nbkg_error;
	new_fit.chi2 = expo_and_dcr->GetChisquare();
    
    cout << " *********** Fit results: ************" << endl;
	cout << "tau:\t" << result_tau << " pm " << result_tau_error << endl;
	cout << "Nsig:\t" << result_Nsig << " pm " << result_Nsig_error << endl;
	cout << "Nbkg:\t" << result_Nbkg << " pm " << result_Nbkg_error << endl;
	cout << "Tot entries:\t" << timeDist->GetEntries() << endl;
	cout << "Entries used for fit:\t" << result_Nsig+result_Nbkg << endl;
    
    TF1 * expo_comp = new TF1(Form("expo_component_%s",timeDist->GetName()), "[0]*exp(-x/[1])", min_for_fit,max_for_fit);
    expo_comp->SetParameter(0,expo_and_dcr->GetParameter(0));
    expo_comp->SetParameter(1,result_tau);
    expo_comp->SetLineColor(kRed);
    expo_comp->SetLineStyle(2);
    expo_comp->Draw("same");
    TF1 * dcr_comp = new TF1(Form("dcr_component_%s",timeDist->GetName()), "[0]", min_for_fit,max_for_fit);
    dcr_comp->SetParameter(0,expo_and_dcr->GetParameter(2));
    dcr_comp->SetLineColor(kRed);
    dcr_comp->SetLineStyle(2);
    dcr_comp->Draw("same");
    
    
    TPaveText * pv = new TPaveText(0.6,0.70,0.88,0.88,"brNDC");
    pv->AddText(Form("#tau_{%s} = %2.1f#pm%2.1f ns",type,1e9*result_tau,1e9*result_tau_error));
    pv->AddText("# DCR pulses =");
    pv->AddText(Form("%2.1f#pm%2.1f (%2.1f%s)",result_Nbkg,result_Nbkg_error,100*(result_Nbkg/(result_Nsig+result_Nbkg)),"%"));
    pv->SetFillColor(kWhite);
    pv->Draw();
    
    return;
}

void roofitTimeDist(TH1 * timeDist, TTree * tree, TCanvas * c, timeFitResult& new_fit, const char * type, double minTime=0, double maxTime=0) {
	double entries = tree->GetEntries();
	if(!timeDist->GetEntries() || !entries) return;
	// timeDist is expected to contain 
	// - a exponential decrease (physical distribution of time of arrivals)
	// - a constant plateau due to DCR
	
	c->cd();
	// Total fit range:
	double minFit = timeDist->GetBinLowEdge(timeDist->GetMaximumBin());
	double maxFit = timeDist->GetBinLowEdge(timeDist->FindLastBinAbove()) + timeDist->GetBinWidth(timeDist->FindLastBinAbove());
	double expFitMax = minFit + ((maxFit-minFit)/3.0);
	// Fits first with an exponential
	TF1 * expo1 = new TF1(Form("expo1_%s",timeDist->GetName()), "expo", minFit,expFitMax);
	timeDist->Fit(Form("expo1_%s",timeDist->GetName()),"QN","",minFit,expFitMax);
	
	// Uses the fit results to fit the entire distribution
	double startAmp = TMath::Exp(expo1->GetParameter(0));
	double startLambda = expo1->GetParameter(1);
	
	//cout << timeDist->GetName() << " " << startAmp << " " << startLambda << endl;
	// RooFit boundaries and variables
	double maxRooFit = tree->GetMaximum("time");
	double minRooFit = tree->GetMinimum("time");
    // Fit range
    double min_for_fit(minRooFit), max_for_fit(maxRooFit);
    if(minTime) min_for_fit = minTime;
    if(maxTime) max_for_fit = maxTime;
	
    //~ RooRealVar time("time","Time",minRooFit,maxRooFit,"s");
    RooRealVar time("time","Time",0,maxRooFit,"s");
    
	//Generate the dataset
	RooDataSet* data = new RooDataSet("data","data",RooArgSet(time),RooFit::Import(*tree));
    time.setRange("fit_range",min_for_fit,max_for_fit);
	
	//make the Signal model -- exponential decay PDF
	RooRealVar lambda("lambda","Lambda [Hz]",startLambda,1.5*startLambda,0.5*startLambda);
	//~ RooRealVar lambda("lambda","Lambda [Hz]",-4.2e7,1.5*-4.2e7,0.5*-4.2e7);
	RooExponential SigModel("SigModel","Exponential decay PDF",time,lambda);
	//~ RooExponential ExpModel("ExpModel","Exponential decay PDF",time,lambda);
	//~ RooRealVar mg("meangauss","meangauss",0,0,maxTime);
	//~ RooRealVar sg("sigmagauss","sigmagauss",0,0,maxTime);
	//~ RooGaussian gauss("gauss","gauss",time,mg,sg);
	//~ RooFFTConvPdf SigModel("SigModel","Exponential decay convoluted with Gaussian",time,ExpModel,gauss);
	
	//make the Background model -- constant bkg due to DCR
	RooPolynomial BkgModel("BkgModel","Constant background PDF",time,RooArgList(),1);
	
	//Yields variables
    double fraction_sig = 1.0;  // assumed fraction of signal (1 = 100% when no DCR is present)
	RooRealVar Nsig("Nsig","Number of signal event",int(entries*fraction_sig),0,entries);
	RooRealVar Nbkg("Nbkg","Number of background event",int(entries*(1-fraction_sig)),0,entries);
	
	RooAddPdf *model = new RooAddPdf("model","Sig + Bkg models", RooArgList(SigModel, BkgModel), RooArgList(Nsig,Nbkg));
	
	// Normalize log likelihood
	/*RooAbsReal * nll = model->createNLL(*data, RooFit::Extended());
	double nll_init_val = nll->getVal(nll->getVariables());
	RooFormulaVar * nll_norm = 0;
	if(nll_init_val > 0) nll_norm = new RooFormulaVar("nll_norm", ("@0-" + to_string(nll_init_val)).c_str(), RooArgList(*nll));
	else nll_norm = new RooFormulaVar("nll_norm", ("@0+" + to_string(-1*nll_init_val)).c_str(), RooArgList(*nll));
	//RooFormulaVar * nll_norm = new RooFormulaVar("nll_norm", ("@0-" + to_string(nll_init_val)).c_str(), *nll);
	cout << "Normalized log likelihood:  " << nll_init_val << " " << nll_norm->getVal(nll_norm->getVariables()) << endl;
	RooAbsReal * nll_toFit = nll_norm;
	RooMinuit m(*nll_toFit);
	//RooMinuit m(*nll_norm);
	m.setPrintLevel(-1);
	m.setWarnLevel(-1);
	m.migrad() ;
	m.hesse() ;*/
	
	//Fit the model to data
	model->fitTo(*data, RooFit::Extended(), RooFit::Range("fit_range"), RooFit::PrintLevel(-1));

	//Plot result of fit
	RooPlot* frame = time.frame(RooFit::Title(tree->GetTitle()));
    //frame->SetAxisRange(0,maxRooFit,"X");
	data->plotOn(frame);
	model->plotOn(frame);
	model->plotOn(frame,RooFit::Components(BkgModel),RooFit::LineStyle(kDashed),RooFit::LineColor(kBlue));
	model->plotOn(frame,RooFit::Components(SigModel),RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
	frame->Draw();
	
	// *********** fit results ************
	
	double result_tau = -1.0/lambda.getVal();	double result_tau_error = abs(lambda.getError()/lambda.getVal())*result_tau;
	double result_Nsig = Nsig.getVal();			double result_Nsig_error = Nsig.getError();
	double result_Nbkg = Nbkg.getVal();			double result_Nbkg_error = Nbkg.getError();
	
	new_fit.tau = result_tau;	new_fit.tau_error = result_tau_error;
	new_fit.Nsig = result_Nsig;	new_fit.Nsig_error = result_Nsig_error;
	new_fit.Nbkg = result_Nbkg;	new_fit.Nbkg_error = result_Nbkg_error;
	new_fit.chi2 = frame->chiSquare(2);
	
	//~ cout << " *********** Fit results: ************" << endl;
	//~ cout << "tau:\t" << result_tau << " pm " << result_tau_error << endl;
	//~ cout << "Nsig:\t" << result_Nsig << " pm " << result_Nsig_error << endl;
	//~ cout << "Nbkg:\t" << result_Nbkg << " pm " << result_Nbkg_error << endl;
	//~ cout << "Tot entries:\t" << entries << endl;
	
	TPaveText * pv = new TPaveText(0.6,0.70,0.88,0.88,"brNDC");
    pv->AddText(Form("#tau_{%s} = %2.1f#pm%2.1f ns",type,1e9*result_tau,1e9*result_tau_error));
    pv->AddText("# DCR pulses =");
    pv->AddText(Form("%2.1f#pm%2.1f (%2.1f%s)",result_Nbkg,result_Nbkg_error,100*(result_Nbkg/entries),"%"));
    pv->SetFillColor(kWhite);
    pv->Draw();
    
    return;
}

double fitSimple1D(TH1 * h, double& error) {
	TCanvas  * c = new TCanvas();
    // first, rough fit over large window
	double pos_maxi = h->GetBinCenter(h->GetMaximumBin());
    double maxi = h->GetMaximum();
    unsigned int bin_min(h->GetMaximumBin()), bin_max(h->GetMaximumBin());
    double fraction = 0.5;
    bool not_reach = true;
    while(bin_min>1 && not_reach) {
        --bin_min;
        double bin_content = h->GetBinContent(bin_min);
        if(bin_content<fraction*maxi) not_reach = false;
    }
    not_reach = true;
    while(bin_max<h->GetNbinsX() && not_reach) {
        ++bin_max;
        double bin_content = h->GetBinContent(bin_max);
        if(bin_content<fraction*maxi) not_reach = false;
    }
    double minf = h->GetBinCenter(bin_min);
    double maxf = h->GetBinCenter(bin_max);
    TF1 * f1 = new TF1("f1","gaus",minf,maxf);
    h->Fit("f1","RQN");
    double mean = f1->GetParameter(1);
    double sigma = f1->GetParameter(2);
    delete f1;
    // second, fit in restricted region around the peak (+/- 1 sigma)
    TF1 * f2 = new TF1("f2","gaus",mean-sigma,mean+sigma);
    f2->SetLineColor(kRed);
    h->Fit("f2","RQ+");
    mean = f2->GetParameter(1);
    error = f2->GetParError(1);
    h->SetTitle(Form("%s (<1PE>=%.2fmV.ns, <Charge>=%.2fmV.ns)",h->GetTitle(),mean*1e12,h->GetMean()*1e12));
    delete f2, c;
    return mean;
}

#endif
