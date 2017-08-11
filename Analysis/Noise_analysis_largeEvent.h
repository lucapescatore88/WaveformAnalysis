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

TGraph * formatGr(TGraph * gr, int color, int fill_bkg, TString xtitle, TString ytitle, TString title = "")
{
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    gr->SetFillColor(color);
    gr->SetFillStyle(fill_bkg);
    gr->SetMarkerSize(1);
    if(title!="") gr->SetTitle(title);
    gr->GetYaxis()->SetTitle(ytitle);
    gr->GetXaxis()->SetTitle(xtitle);

    return gr;
}

TH1 * drawWave(TGraph * waveform, int * color, TString title, TCanvas * c, float ymax)
{
    c->cd();

    (*color) ++;
    if ((*color) > 20) (*color) = 2;
    
    TH1 * waveh = convertGrToH(waveform);

    waveh->SetMaximum(1.8*ymax);
    waveh->SetLineColor((*color));
    waveh->SetMarkerColor((*color));
    waveh->SetMarkerSize(1);
    waveh->SetTitle(title); 
    waveh->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
    waveh->GetXaxis()->SetTitle("Time [s]");
    
    if((*color)==1)
    {
        waveh->Draw("L");
        c->SetGrid();
    }
    else waveh->Draw("L SAME");

    return waveh;
}

void fillPersistence(TH2D * persistence, TGraph * waveform) {
	double * xx = waveform->GetX();
	double * yy = waveform->GetY();
	unsigned int Npts = waveform->GetN();
	for(unsigned int pt(0); pt<Npts; ++pt) {
		persistence->Fill(xx[pt],yy[pt]);
	}
	return;
}

TCanvas * drawPersistence(TH2D * persistence, TCanvas * c, TString title = "")
{
    c->cd();
    if(title!="") persistence->SetTitle(title);
    persistence->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
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

TF1 * fitLongTau2(TH2 * h, double * amp0, double * tau, double pe, const char * vol);
TF1 * drawPersistenceWithLongTauFit(TH2D * persistence, TCanvas * c, double * amp0, double * tau, double pe, const char * vol, TString title = "")
{
    c->cd();
    if(title!="") persistence->SetTitle(title);
    persistence->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
    persistence->GetXaxis()->SetTitle("Time [s]");
    TF1 * exp_longtau = 0;
    if(persistence->GetEntries()) {
		// have to do some tricks with TPads otherwise it does not work with draw option "CONT4Z"...
		//~ TPad *pad1 = new TPad("pad1","",0,0,1,1);
		//~ TPad *pad2 = new TPad("pad2","",0,0,1,1);
		//~ pad2->SetFillStyle(0); //will be transparent
		//~ pad2->SetFillColor(0);
		//~ pad2->SetFrameFillStyle(0);
		//~ pad1->Draw();
		//~ pad1->cd();
	    //~ persistence->Draw("CONT4Z");
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
	    //~ pad1->SetLogz();
	    //~ pad1->Update();
	    //~ //persistence->GetListOfFunctions()->Print();
	    //~ TPaletteAxis *palette = (TPaletteAxis*)persistence->GetListOfFunctions()->FindObject("palette");
		//~ palette->SetX1NDC(0.91);
		//~ palette->SetX2NDC(0.93);
		//~ palette->SetY1NDC(0.25);
		//~ palette->SetY2NDC(0.9);
		//~ pad1->Modified();
		//~ pad1->Update();
		
		//~ double ymin = 0;
		//~ double ymax = 0.03;
		//~ double dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
		//~ double xmin = 0;
		//~ double xmax = 0.18e-6;
		//~ double dx = (xmax-xmin)/0.8; //10 per cent margins left and right
		//~ pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
		//~ pad2->Draw();
		//~ pad2->cd();
		exp_longtau = fitLongTau2(persistence, amp0, tau, pe, vol);
		exp_longtau->Draw("SAME A*");
	    TPaveText * pv = new TPaveText(0.6,0.65,0.75,0.74,"brNDC");
	    pv->AddText(Form("#tau_{long} = %2.1fns",1e9*(*tau)));
	    pv->SetFillColor(kWhite);
	    pv->Draw();
	    gPad->Update();
	    //~ pad2->Update();
		
	}
    return exp_longtau;
}

Double_t Amplitude_calc(const char * vol_folder, Int_t data_size, string option = "root", TFile * file = NULL)
{
    TString canvas_title = "Amplitude calculation "+TString(vol_folder);
    TH1D * volt_ampl = NULL;

    TFile * f = NULL;
    TTree * tree = NULL;
    int nsamples;
    double times[10000];
    double amps[10000];
    TGraph * waveform = NULL;
	
	double bin_size = 0.002; // bin size of amplitude histogram (PE distribution)
	double minval, maxval;
    if(option=="root") 
    {
        f = TFile::Open(TString(globalArgs.data_folder)+"/oscilloscope_out.root");
        tree = (TTree*)f->Get(vol_folder);
        tree->SetBranchAddress("NsampPerEv",&nsamples);
        tree->SetBranchAddress("Amps",&amps);
        tree->SetBranchAddress("Times",&times);
        data_size = tree->GetEntries();
        // Defines amplitude max, min and bins
        maxval = tree->GetMaximum("Amps");
        minval = tree->GetMinimum("Amps");
    } else {
		maxval = -0.05;
		minval = 0.25;
	}
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
            else if (time[pt] > 1. * ns) break;
            if (maxV < volts[pt]) maxV = volts[pt];
        }

        volt_ampl->Fill(maxV);
        delete time;
        delete volts;
    }

    // Fit the first peak distribution to get pe

    double pos_maxi = volt_ampl->GetBinCenter(volt_ampl->GetMaximumBin());
    double around   = 0.5*pos_maxi;
    TF1 * f1 = new TF1("f1","gaus",pos_maxi-around,pos_maxi+around);
    // TF1 *f1 = new TF1("f1","gaus",0,1.0); // Change range for fit of MPV
    // pe bigger than 0.8 wont be detected unless changed
    volt_ampl->Fit("f1","RQ");	// fit line in red???
    Double_t pe_volt   = f1->GetParameter(1);
    volt_ampl->SetTitle(globalArgs.res_folder+Form("Amplitude calculation %s, pe = %2.3f",vol_folder,pe_volt));
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



vector <TString> readSetupFile(ifstream * setupFile, int * data_size, vector <Thresholds> &thresholds)
{
    string s;
    vector <TString> vol_folders;

    while (true) 
    {    
        Double_t rt, ap, delay, imme;
        
        getline((*setupFile), s);
        if ((*setupFile).eof()) break;
        
        const char* searchString = s.c_str();
        char volt [20];
        Int_t numfiles;
        
        if (s.find("#") == 0 || s=="") continue; // Skip commented  or empty lines
        
        //Find the voltages
        if(sscanf(searchString, "V || %s ||", volt)==1)
        {
            vol_folders.push_back(volt);
            thresholds.push_back(Thresholds(default_thrs));
        }
        
        if(sscanf(searchString, "V/th || %s ||", volt)==1)
        {
            vol_folders.push_back(volt);
            getline((*setupFile), s);
            const char* thresholds_string = s.c_str();
            sscanf(thresholds_string, "Rej_t: %lf, AP_th: %lf, Delay_th: %lf, Imm_th: %lf", &rt,&ap,&delay,&imme);
            Thresholds defThr{rt, default_thrs.del_xtalk_minT, default_thrs.dir_xtalk_maxT, ap, imme, delay, 0.4};
            thresholds.push_back(defThr);
        }
        
        if(sscanf(searchString, "Files at each voltage || %d ||", &numfiles)==1) 
            (*data_size) = numfiles;  
    }

    return vol_folders;
}

void setOptions(int argc, char* argv[], const char * optString = "d:S:o:T:aVh?")
{
    int opt = getopt(argc, argv, optString);
    if(opt == -1){
        std::cerr <<  "There is no opption in the command! Type \"output -h\" for help." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    while(opt != -1){
        switch(opt){
            case 'd':
                globalArgs.data_folder = optarg;
                //std::cout << "-p option path= " << globalArgs.data_folder << std::endl;
                break;
            case 'S':
                globalArgs.arg_pathToSetupFile = optarg;
                break;
            case 'o':
                globalArgs.res_folder = optarg;
                break;
            case 'a':
                globalArgs.save_all = true;
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
