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
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TVirtualPad.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TSpectrum.h>

using namespace std;

//Amplitude of pe calculation of raw data

TH1 * drawWave(TGraph * waveform, int * color, TString title, TCanvas * c, float ymax)
{
    c->cd();

    (*color) ++;
    if ((*color) > 20) (*color) = 2;
    
    TH1 * waveh = convertGrToH(waveform);
    
    waveh->SetLineColor((*color));
    waveh->SetMarkerColor((*color));
    waveh->SetMarkerSize(1);
    waveh->SetTitle(title);
   
    //waveh->GetYaxis()->SetRangeUser(-0.02,0.1);
    //waveh->GetXaxis()->SetRangeUser(-40*ns,20*ns);
    waveh->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
    waveh->GetYaxis()->SetTitleOffset(1.3);
    waveh->GetXaxis()->SetTitle("Time [s]");

    if((*color)==1)
    {
        waveh->Draw("L");
        c->SetGrid();
    }
    else waveh->Draw("L SAME");

    return waveh;
}

Double_t Amplitude_calc(const char * vol_folder, Int_t data_size, string option = "root", TFile * file = NULL)
{
    TString canvas_title = "Amplitude calculation "+TString(vol_folder);
    TH1D * volt_ampl = new TH1D(canvas_title, canvas_title, 150, -0.05, 0.25);	// Adjust range of histogram, might change with amplifier used!

    cout << " -----> Amplitude calculation of pe: " << vol_folder << endl;

    TFile * f = NULL;
    TTree * tree = NULL;
    int nsamples;
    double times[10000];
    double amps[10000];
    TGraph * waveform = NULL;

    if(option=="root") 
    {
        f = TFile::Open(TString(globalArgs.data_folder)+"/oscilloscope_out.root");
        tree = (TTree*)f->Get(vol_folder);
        tree->SetBranchAddress("NsampPerEv",&nsamples);
        tree->SetBranchAddress("Amps",&amps);
        tree->SetBranchAddress("Times",&times);
        data_size = tree->GetEntries();
    }

    //loop over every measurement on a folder

    for (int j = 0; j < data_size; j++) 
    {
        // Get the waveform
        if(option=="root") 
        {
            tree->GetEntry(j);
            waveform = new TGraph(nsamples,times,amps);
            //if(j<4 && TString(vol_folder)=="54.0V") {
            //TCanvas * c1 = new TCanvas();
            //waveform->Draw("AP");
            //c1->Print(Form("%i.pdf",j)); 
            //}
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
    volt_ampl->Fit("f1","RQ");
    Double_t pe_volt   = f1->GetParameter(1);
    volt_ampl->SetTitle(globalArgs.res_folder+Form("Amplitude calculation %s, pe = %2.3f",vol_folder,pe_volt));

    delete f1;
    delete volt_ampl;

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



vector <TString> readSetupFile(ifstream * setupFile, int * data_size, vector <Thresholds_t> &thresholds)
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
            Thresholds_t defThr{4, 2, 0.5, 1.17, 0.85, 0.4};
            thresholds.push_back(defThr);
        }
        
        if(sscanf(searchString, "V/th || %s ||", volt)==1)
        {
            vol_folders.push_back(volt);
            getline((*setupFile), s);
            const char* thresholds_string = s.c_str();
            sscanf(thresholds_string, "Rej_t: %lf, AP_th: %lf, Delay_th: %lf, Imm_th: %lf", &rt,&ap,&delay,&imme);
            Thresholds_t defThr{rt, 2, ap, imme, delay, 0.4};
            thresholds.push_back(defThr);
        }
        
        if(sscanf(searchString, "Files at each voltage || %d ||", &numfiles)==1) 
            (*data_size) = numfiles;  
    }

    return vol_folders;
}

#endif 
