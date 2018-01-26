//
//  Created by Felipe Gilberto Ortega on 25/04/16.
//  Modified by Luca Pescatore, luca.pescatore@cern.ch on 28/06/2017   
//

#ifndef Noise_analysis_largeEvent_h
#define Noise_analysis_largeEvent_h

#include "genparam.h"
#include "Thresholds.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TVirtualPad.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TF2.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TSpectrum.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TObject.h>
#include <TPave.h>
#include <TPaletteAxis.h>

using namespace std;

// Amplitude of pe calculation of raw data

TGraph * formatGr(TGraph * gr, int color, int fill_bkg, TString xtitle, TString ytitle, TString title = "", TString name = "", int marker_style=8)
{
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    gr->SetFillColor(color);
    gr->SetFillStyle(fill_bkg);
    gr->SetMarkerSize(1);
    gr->SetMarkerStyle(marker_style);
    if(title!="") gr->SetTitle(title);
    gr->GetYaxis()->SetTitle(ytitle);
    gr->GetXaxis()->SetTitle(xtitle);
    if(name!="") gr->SetName(name);

    return gr;
}

TH1 * drawWave(TGraph * waveform, int * color, TString title, TCanvas * c, float ymax, double yscale=1, double baseline_shift=0)
{
    c->cd();
    
    (*color) ++;
    int col = ((*color)%50)+50;
    //if ((*color) > 20) (*color) = 2;
    
    TH1 * waveh = convertGrToH(waveform, yscale, baseline_shift);

    waveh->SetMaximum(ymax);
    waveh->SetLineColor(col);
    waveh->SetMarkerColor(col);
    waveh->SetMarkerSize(1);
    waveh->SetTitle(title);
    waveh->GetYaxis()->SetTitle("Signal [V]");
    if(yscale!=1 && yscale!=0) waveh->GetYaxis()->SetTitle("Signal [PE]");
    waveh->GetXaxis()->SetTitle("Time [s]");
    
    if((*color)==1)
    {
        waveh->Draw("L");
        c->SetGrid();
    }
    else waveh->Draw("L SAME");

    return waveh;
}

void drawAPfit(TCanvas * c, TF1 * fit, double yscale=1) {
    c->cd();
    TF1 * exp = new TF1(TString(fit->GetName())+"_scaled","[0]*(1 - exp(-(x-[5])/[1])) + [2]*exp(-(x-[4])/[3])",0,180*ns);
    exp->SetParameter(0, fit->GetParameter(0));
    exp->SetParameter(1, fit->GetParameter(1));
    exp->SetParameter(2, fit->GetParameter(2));
    exp->SetParameter(3, fit->GetParameter(3));
    exp->SetParameter(4, fit->GetParameter(4));
    exp->SetParameter(5, fit->GetParameter(5));
    exp->SetRange(0*ns,180*ns);
    if(yscale !=1 && yscale!=0) {
        exp->SetParameter(0, fit->GetParameter(0)/yscale);
        exp->SetParameter(2, fit->GetParameter(2)/yscale);
    }
    exp->SetLineColor(kBlack);
    exp->Draw("same");
    TPaveText * pv = new TPaveText(0.55,0.78,0.85,0.85,"brNDC");
    pv->AddText(Form("#tau_{rec} = %2.1f#pm%2.1f ns",1e9*fit->GetParameter(1),1e9*fit->GetParError(1)));
    pv->SetFillColor(kWhite);
    pv->Draw();
    return;
}

void fillPersistence(TH2D * persistence, TGraph * waveform, double yscale=1) {
	double * xx = waveform->GetX();
	double * yy = waveform->GetY();
	unsigned int Npts = waveform->GetN();
    bool yscale_enabled = (yscale!=1) && (yscale!=0);
	if(yscale_enabled) for(unsigned int pt(0); pt<Npts; ++pt) persistence->Fill(xx[pt],yy[pt]/yscale);
	else for(unsigned int pt(0); pt<Npts; ++pt) persistence->Fill(xx[pt],yy[pt]);
	return;
}

TCanvas * drawPersistence(TH2D * persistence, TCanvas * c, TString title = "")
{
    c->cd();
    if(title!="") persistence->SetTitle(title);
    persistence->GetYaxis()->SetTitle("Signal [PE]");
    persistence->GetXaxis()->SetTitle("Time [s]");
    if(persistence->GetEntries()) {
	    persistence->Draw("CONT4Z");
	    c->SetLogz();
	    gPad->Update();
	    //persistence->GetListOfFunctions()->Print();
	    TPaletteAxis *palette = (TPaletteAxis*)persistence->GetListOfFunctions()->FindObject("palette");
		palette->SetX1NDC(0.91);
		palette->SetX2NDC(0.93);
		palette->SetY1NDC(0.23);
		palette->SetY2NDC(0.9);
		gPad->Modified();
		gPad->Update();
	}
    return c;
}

TF1 * fitLongTau_TH2(TH2 * h, double * amp0, double * tau, double pe, const char * vol, TGraphErrors * &averageGraph);
TF1 * fitLongTau_TGraphErrors(TGraphErrors * cleanwaves, double * amp0, double * tau, double pe, const char * vol);

TF1 * drawPersistenceWithLongTauFit(TH2D * persistence, TCanvas * c, double * amp0, double * tau, double pe, const char * vol, TString title = "", TGraphErrors * gavg = NULL)
// gavg is the average clean pulse shape measured from "very clean" events
// if gavg is not set, the long tau fit is performed directly from the 2D histogram
{
    c->cd();
    if(title!="") persistence->SetTitle(title);
    persistence->GetYaxis()->SetTitle("Signal [PE]");
    persistence->GetXaxis()->SetTitle("Time [s]");
    TF1 * exp_longtau = 0;
    TGraphErrors * average_pulse = 0;
    if(persistence->GetEntries()) {
		
        if(!gavg) {
            // fit from the 2D histogram
            exp_longtau = fitLongTau_TH2(persistence, amp0, tau, pe, vol, average_pulse);
        } else {
            // fit from the average waveform
            exp_longtau = fitLongTau_TGraphErrors(gavg, amp0, tau, pe, vol);
            average_pulse = gavg;
        }
		
	    persistence->Draw("COLZ");
	    c->SetLogz();
	    gPad->Update();
	    //persistence->GetListOfFunctions()->Print();
	    TPaletteAxis *palette = (TPaletteAxis*)persistence->GetListOfFunctions()->FindObject("palette");
		palette->SetX1NDC(0.91);
		palette->SetX2NDC(0.93);
		palette->SetY1NDC(0.25);
		palette->SetY2NDC(0.9);
		gPad->Modified();
		gPad->Update();
		
		average_pulse->Draw("P+");
		exp_longtau->Draw("SAME A*");
	    TPaveText * pv = new TPaveText(0.55,0.78,0.85,0.85,"brNDC");
	    pv->AddText(Form("#tau_{long} = %2.1f#pm%2.1f ns",1e9*(*tau),1e9*exp_longtau->GetParError(1)));
	    pv->SetFillColor(kWhite);
	    pv->Draw();
	    gPad->Update();
		
	}
    return exp_longtau;
}

Double_t Amplitude_calc(const char * vol_folder, Int_t data_size, vector<double>& minmax, string option = "root", TFile * file = NULL)
{
	cout << "Amplitude calculation for dV = " << vol_folder << endl;
    TString canvas_title = "Amplitude calculation "+TString(vol_folder);
    TH1D * volt_ampl = NULL;

    TFile * f = NULL;
    TTree * tree = NULL;
    int nsamples;
    double times[20002];
    double amps[20002];
    double amps_copy[20002];
    TGraph * waveform = NULL;
	
	double bin_size = 0.002; // bin size of amplitude histogram (PE distribution)
	//double bin_size = 0.0002; // bin size of amplitude histogram (PE distribution)
	double minval, maxval;
    double xmin(1e6), xmax(0), ymin(1e6), ymax(0);
    int nmax(0);
    if(option=="root") 
    {
        f = TFile::Open(TString(globalArgs.data_folder)+"/oscilloscope_out.root");
        tree = (TTree*)f->Get(vol_folder);
        tree->SetBranchAddress("NsampPerEv",&nsamples);
        tree->SetBranchAddress("Amps",&amps);
        tree->SetBranchAddress("Times",&times);
        //data_size = tree->GetEntries();
        // data_size defined from cfg file
        if(data_size>tree->GetEntries()) data_size = tree->GetEntries();
        // Defines amplitude max, min and bins
        
        tree->GetEntry(0);
        xmax = times[nsamples-1];
        xmin = times[0];
        nmax = nsamples;
        
        ymax = tree->GetMaximum("Amps");
        ymin = tree->GetMinimum("Amps");
        maxval = ymax;
        minval = ymin;
    } else {
		maxval = -0.05;
		minval = 0.25;
	}
    
    minmax[0] = xmin;
    minmax[1] = xmax;
    minmax[2] = ymin;
    minmax[3] = ymax;
    minmax[4] = nmax;
    
	volt_ampl = new TH1D(canvas_title, canvas_title, int(((maxval-minval)/bin_size)+0.5), minval,maxval);
    if(file) file->cd("Vbd_determination");

    //loop over every measurement on a folder

    for (int j = 0; j < data_size; j++) 
    {
        // Get the waveform
        if(option=="root") 
        {
            tree->GetEntry(j);
            waveform = new TGraph(nsamples,times,amps);
        }
        else 
        {
            TString datafilename = Form("%s%s/%i.csv",globalArgs.data_folder,vol_folder,j);
            waveform = new TGraph(datafilename,"%lg %lg","/t;,");
        }
        if (waveform->IsZombie()) continue;

        int npts         = waveform->GetN();
        Double_t * time  = waveform->GetX();
        Double_t * volts = waveform->GetY();
        Double_t  maxV   = 0.0;

        for (int pt = 0; pt < npts; pt++) 
        {
            // 1 nanosecond window for first pulse
            if (time[pt] < 0.) continue ; 
            else if (time[pt] > fall_time) break;
            if (maxV < volts[pt]) maxV = volts[pt];
        }

        volt_ampl->Fill(maxV);
        
        delete time;
        delete volts;
    }

    // Fit the first peak distribution to get pe
    
    volt_ampl->GetXaxis()->SetTitle("Pulse amplitude [V]");
    volt_ampl->GetYaxis()->SetTitle("Entries");

    double pos_maxi = volt_ampl->GetBinCenter(volt_ampl->GetMaximumBin());
    double around   = 0.5*pos_maxi;
    TF1 * f1 = new TF1("f1","gaus",pos_maxi-around,pos_maxi+around);
    f1->SetLineColor(kRed);
    // TF1 *f1 = new TF1("f1","gaus",0,1.0); // Change range for fit of MPV
    // pe bigger than 0.8 wont be detected unless changed
    volt_ampl->Fit("f1","RQ+");
    Double_t pe_volt   = f1->GetParameter(1);
    volt_ampl->SetTitle(Form("Amplitude calculation %s, <A_{pe}> = %2.3fV",vol_folder,pe_volt));
    if(option=="root") volt_ampl->Write();
    
    delete f1;
    delete volt_ampl;
    
    if(file) file->cd();

    return pe_volt;
}

// Format waveform graphs

Double_t FindTimeRange(Int_t ROWS_DATA, Double_t *time, Double_t *volts, 
    vector <Double_t> reject_time_v, vector <Double_t> time_dist_th_v, int i, int j, Double_t pe)
{ 
    TString spec = TString("Spectrum_")+Form("%d",j);
    TH1D * Spectrum = new TH1D (spec,spec,1000,0.,time[ROWS_DATA-1]);

    for (int row = 0; row < ROWS_DATA; row++) Spectrum->Fill(time[row],volts[row]);
    double tmp{0};
    double time_of_max_first_fit{0};
    double sig_max_first_fit{0};
    double sig_max_fit{0};
    
    TSpectrum * s = new TSpectrum();
    s->Search(Spectrum,1,"",0.2);

    Double_t * peaks = NULL;
    Int_t Npeaks;
    Double_t * amppeaks = NULL;

    peaks    = s->GetPositionX();
    amppeaks = s->GetPositionY();
    Npeaks   = s->GetNPeaks();

    if(Npeaks > 1 && peaks[1] > reject_time_v.at(i)*ns)
    {
        sig_max_fit=amppeaks[1];
        tmp=peaks[1];
        if(amppeaks[1] > time_dist_th_v.at(i) * pe)
        {
            time_of_max_first_fit = tmp;
            sig_max_first_fit = sig_max_fit;
        }
    }
    
    delete Spectrum;
    delete s;
    return tmp;
}

void addPointToGraph(std::map<string,TGraph *> &g, string option, double x, double y) {
    string opt = "";
    if(option=="DiXT") opt = "1:DiXT";
    else if(option=="DeXT") opt = "3:DeXT";
    else if(option=="AP") opt = "2:AP";
    else if(option=="Sec") opt = "4:Sec";
    else return;
    unsigned int N = g[opt]->GetN();
    g[opt]->SetPoint(N,x,y);
    return;
}

TCanvas * finalizeMapGraphs(std::map<string,TGraph *> &g, TString name, TString title, TString Xtitle, TString Ytitle, std::vector<TString> add_legend) {
    TCanvas * c = new TCanvas(name, name);
    double minX(1000), minY(1000), maxX(0), maxY(0);
    for ( const auto &pair : g ) {
        unsigned int N = pair.second->GetN();
        for(unsigned int pt(0); pt<N; ++pt) {
            double x,y;
            pair.second->GetPoint(pt,x,y);
            if(x<minX) minX = x;
            if(y<minY) minY = y;
            if(x>maxX) maxX = x;
            if(y>maxY) maxY = y;
        }
    }
    double distX = maxX-minX;
    double distY = maxY-minY;
    TH2D* frame = new TH2D("frame",title, 1000, minX-0.1*distX, maxX+0.1*distX, 1000, 0, maxY+0.1*distY);
    frame->GetXaxis()->SetTitle(Xtitle);
    frame->GetYaxis()->SetTitle(Ytitle);
    frame->Draw();
    
    TLegend * leg = new TLegend(0.45,0.65,0.87,0.87);
    
    vector<string> ordered_corr_noise;
    ordered_corr_noise.push_back("DiXT"); ordered_corr_noise.push_back("AP"); ordered_corr_noise.push_back("DeXT"); ordered_corr_noise.push_back("Sec");
    // Draws the graphs
    for(int i(ordered_corr_noise.size()-1); i>=0; --i) {
        for ( auto &pair : g ) {
            if(TString(pair.first).Contains(ordered_corr_noise[i])) {
                int color(1), marker(7);
                double marker_size(1);
                if(TString(pair.first).Contains("DiXT")) color = kBlue;
                if(TString(pair.first).Contains("DeXT")) color = kGreen+2;
                if(TString(pair.first).Contains("AP")) color = kOrange+7;
                if(TString(pair.first).Contains("Sec")) color = 7;
                pair.second->SetMarkerStyle(marker);
                pair.second->SetMarkerSize(marker_size);
                pair.second->SetMarkerColor(color);
                pair.second->Draw("P+");
            }
        }
    }
    // Fills the legend
    // and copy the graphs with a larger marker to be more visible in TLegend
    vector< TGraph* > g_copy(g.size(),0);
    for(unsigned int i(0); i<ordered_corr_noise.size(); ++i) {
        for ( auto &pair : g ) {
            if(TString(pair.first).Contains(ordered_corr_noise[i])) {
                TString ent("");
                if(TString(ordered_corr_noise[i]).Contains("DiXT")) ent = "Direct cross-talk";
                if(TString(ordered_corr_noise[i]).Contains("DeXT")) ent = "Delayed cross-talk";
                if(TString(ordered_corr_noise[i]).Contains("AP")) ent = "After-pulse";
                if(TString(ordered_corr_noise[i]).Contains("Sec")) ent = "Secondaries";
                TString add = "";
                
                TGraph * gcp = new TGraph();
                gcp->SetName(TString(pair.second->GetName())+"_copy");
                int color = pair.second->GetMarkerColor();
                gcp->SetMarkerStyle(8);
                gcp->SetMarkerSize(1.5);
                gcp->SetMarkerColor(color);
                g_copy[i] = gcp;
                
                if(add_legend.size()>i) add = " "+add_legend[i];
                //~ leg->AddEntry(pair.second,ent+add,"p");
                leg->AddEntry(gcp,ent+add,"p");
            }
        }
    }
    leg->Draw();
    return c;
}

double fitSimple1D(TH1 * h, double& error);
TCanvas * finalizeChargeTimeWindow(TH2D * h, TGraphErrors * r, TString name) {
    TCanvas * c = new TCanvas(name, name);
    unsigned int Nbins = h->GetNbinsX();
    for(unsigned int bin(0); bin<Nbins; ++bin) {
        TH1D * proj = h->ProjectionY("proj_time_window", bin+1, bin+1);
        
        double x = h->GetXaxis()->GetBinCenter(bin+1);
        double error(0);
        double value = fitSimple1D(proj, error);
        
        setPoint(r,  bin, x, value, 0, error);
        
        delete proj;
    }
    
    c->cd();
    h->Draw("COLZ");
    c->SetLogz();
    r->Draw("PL+3");
    
    return c;
}

vector <TString> readSetupFile(ifstream * setupFile, int * data_size, vector <Thresholds> &thresholds, double * integration_time=0, int * integration_number=0)
{
    string s;
    vector <TString> vol_folders;
    
    Thresholds thrs_code(default_thrs);
	double def_AP_minT(thrs_code.AP_minT), def_DeXT_minT(thrs_code.del_xtalk_minT), def_DiXT_maxT(thrs_code.dir_xtalk_maxT), def_AP_thrs(thrs_code.AP), def_DiXT_thrs(thrs_code.dir_xtalk), def_DeXT_thrs(thrs_code.del_xtalk);

    while (true) 
    {    
        Double_t rt, ap, delay, imme;
        
        getline((*setupFile), s);
        if ((*setupFile).eof()) break;
        
        const char* searchString = s.c_str();
        char volt [20];
        Int_t numfiles;
        double tmp_d;
        
        if (s.find("#") == 0 || s=="") continue; // Skip commented  or empty lines
        
        // Threshold defaulf values
        sscanf(searchString, "DefThrs :: AP_minT=%lf :: DeXT_minT=%lf :: DiXT_maxT=%lf :: AP_thrs=%lf :: DiXT_thrs=%lf :: DeXT_thrs=%lf", &def_AP_minT, &def_DeXT_minT, &def_DiXT_maxT, &def_AP_thrs, &def_DiXT_thrs, &def_DeXT_thrs);
        if(fall_time>def_DiXT_maxT) fall_time = def_DiXT_maxT;
        
        //Find the voltages and thresholds
        char thres_config[256] = "";
        if(sscanf(searchString, "V || %s ||%[^\t\n]", volt, thres_config)>=1)
        {
            vol_folders.push_back(volt);
            //if(thres_config == "") thresholds.push_back(Thresholds(default_thrs));
            if(thres_config == "") thresholds.push_back(defThresholds(def_AP_minT, def_DeXT_minT, def_DiXT_maxT, def_AP_thrs, def_DiXT_thrs, def_DeXT_thrs));
            else {
				// format should be: (times are given in nanoseconds and thrs in units of pe
				// AP_minT=%lf :: DeXT_minT=%lf :: DiXT_maxT=%lf :: AP_thrs=%lf :: DiXT_thrs=%lf :: DeXT_thrs=%lf
				double AP_minT, DeXT_minT, DiXT_maxT, AP_thrs, DiXT_thrs, DeXT_thrs;
				if(sscanf(thres_config, " AP_minT=%lf :: DeXT_minT=%lf :: DiXT_maxT=%lf :: AP_thrs=%lf :: DiXT_thrs=%lf :: DeXT_thrs=%lf", &AP_minT, &DeXT_minT, &DiXT_maxT, &AP_thrs, &DiXT_thrs, &DeXT_thrs)==6)
					thresholds.push_back(defThresholds(AP_minT, DeXT_minT, DiXT_maxT, AP_thrs, DiXT_thrs, DeXT_thrs));
				else thresholds.push_back(defThresholds(def_AP_minT, def_DeXT_minT, def_DiXT_maxT, def_AP_thrs, def_DiXT_thrs, def_DeXT_thrs));
			}
        }
        
        if(sscanf(searchString, "Files at each voltage || %d ||", &numfiles)==1) 
            (*data_size) = numfiles;
        
        if(sscanf(searchString, "Integration_time_step : %lf", &tmp_d)==1 && integration_time)
            (*integration_time) = tmp_d;
        if(sscanf(searchString, "Integration_time_number : %d", &numfiles)==1 && integration_number)
            (*integration_number) = numfiles;
    }
    
    for(unsigned int i(0); i<vol_folders.size(); ++i)
	    cout << "Voltage: " << vol_folders[i] << " // " << Form("AP_minT=%2.2lf :: DeXT_minT=%2.2lf :: DiXT_maxT=%2.2lf :: AP_thrs=%2.2lf :: DiXT_thrs=%2.2lf :: DeXT_thrs=%2.2lf", thresholds[i].AP_minT, thresholds[i].del_xtalk_minT, thresholds[i].dir_xtalk_maxT, thresholds[i].AP, thresholds[i].dir_xtalk, thresholds[i].del_xtalk) << endl;

    return vol_folders;
}

void setOptions(int argc, char* argv[], const char * optString = "d:S:o:T:atnVh?")
{
    int opt = getopt(argc, argv, optString);
    if(opt == -1){
        std::cerr <<  "There is no opption in the command! Type \"output -h\" for help." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    while(opt != -1){
        switch(opt){
            case 'd':
	            if(TString(optarg).EndsWith("/")) globalArgs.data_folder = optarg;
	            else globalArgs.data_folder = (string(optarg) + string("/")).c_str();
                //std::cout << "-p option path= " << globalArgs.data_folder << std::endl;
                break;
            case 'S':
                globalArgs.arg_pathToSetupFile = optarg;
                break;
            case 'o':
                if(TString(optarg).EndsWith("/")) globalArgs.res_folder = optarg;
	            else globalArgs.res_folder = (string(optarg) + string("/")).c_str();
                break;
            case 'a':
                globalArgs.save_all = true;
                break;
            case 't':
                globalArgs.save_tree = true;
                break;
            case 'n':
                globalArgs.enable_dcr = true;
                break;
            case 'T':
                globalArgs.fixed_thr = atof(optarg);
                break;
            case 'V':
                globalArgs.only_Vbd = 1;
                break;
            case 'h':
            case '?':
                std::cerr << "Usage: output -d pathToData -S pathToSetupFile -o pathToResultsFolder [-a]" << std::endl;
                std::cerr << "----------------------------------------------------------------------------------------------------" << std::endl;
                std::cerr << " '-d'+'-S'+'-o' options are necessary!" << std::endl;
                std::cerr << "-----------------------------------------------------------------------------------------------------" << std::endl;
                std::cerr << " use '-a' option afterwards to save all the plots of the analysis to further check." << std::endl;
                std::cerr << " use '-V' option afterwards to perform only Vbd measurement and quit." << std::endl;
                std::cerr << "-----------------------------------------------------------------------------------------------------" << std::endl;
                std::cerr << "Example: ./output -d /Users/Analysis_waveforms/ov_scan_pde_H2014/ -S /Users/Analysis_waveforms/config_file.txt -o /Users/Analysis_waveforms/Plots/ [-a]"<<std::endl;
                exit(EXIT_FAILURE);
                break;
            default:
                break;
        }
        opt = getopt(argc, argv, optString);
    }

    if((strncmp(globalArgs.data_folder," ",1) == 0|| strncmp(globalArgs.arg_pathToSetupFile," ",1) == 0)){
        std::cerr << "ERROR: -d or -S option is not set! Both of them has to be set correctly!"<<std::endl;
        exit(EXIT_FAILURE);
    }
    
    if(strncmp(globalArgs.res_folder," ",1) == 0){
        std::cerr << "ERROR: -o option is not set! It has to be set up correctly!"<<std::endl;
        exit(EXIT_FAILURE);
    }
}

void initOutputFile(TFile * f, vector<TString> vol_folders) {
	TDirectory *cdvbd = f->mkdir("Vbd_determination");
	for(unsigned int i(0); i<vol_folders.size(); ++i) {
		TString dirname = "pulse_shape_" + vol_folders[i];
		TDirectory *cdvbd = f->mkdir(dirname);
	}
	return;
}

#endif 
