//
//  Noise_analysis.h
//  
//
//  Created by Felipe Gilberto Ortega on 25/04/16.
//
//

#ifndef Noise_analysis_largeEvent_h
#define Noise_analysis_largeEvent_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

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
#include <TMultiGraph.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TSpectrum.h>


using namespace std;

/* global variable declaration */

const double ns = 1e-9;

/* Auxiliar functions */

struct globalArgs_t
{
    const char* data_folder;                /* -d option */
    const char* arg_pathToSetupFile;             /* -S option */
    const char* results_folder;             /* -o option */
    int save_all;               /* -a  */
    
} globalArgs;


//Amplitude of pe calculation of raw data
Double_t Amplitude_calc(const char* vol_folder, Int_t data_size){
    
    Double_t pe_volt;
    
    Char_t canvas_title[200];
    sprintf(canvas_title,"Amplitude calculation %s",vol_folder);
    
    //TH1D* volt_ampl= new TH1D(canvas_title, canvas_title, 150, -0.05, 1.2);
    TH1D* volt_ampl= new TH1D(canvas_title, canvas_title, 150, -0.05, 0.25);	// Adjust range of histogram, might change with amplifier used!
    
    cout<<"****----->Amplitude calculation of pe: "<< vol_folder << endl;
    
    //loop over every measurement on a folder
    for (int j=0; j<data_size; j++) {
        if (j%500==0) {
            cout<<"Event: "<<j<<endl;
        }
        
        Char_t datafilename[100];
        
        //Get the waveform
        sprintf(datafilename,"%s%s/%i.csv",globalArgs.data_folder,vol_folder,j);
        TGraph* waveform = new TGraph(datafilename,"%lg %lg","/t;,");
        if (waveform->IsZombie()) continue;
        
        int ROWS_DATA = waveform->GetN();
        Double_t *time = waveform->GetX();
        Double_t *volts = waveform->GetY();
        Double_t  maxV=0.0;
        
        //loop over the time window to get height of first pulse
        // estimates the search window :
        //~ double start_time = time[0];
        //~ double end_time = time[ROWS_DATA-1];
        //~ double incr = (end_time-start_time)/(ROWS_DATA-1);
        //~ unsigned int start_row = int((-1.0*start_time)/incr)-10;
        //~ unsigned int end_row = int((-1.0*start_time + 1*ns)/incr)+10;
        for (int row = 0; row<ROWS_DATA; row++) {
        //for (int row = start_row; row<end_row; row++) {
            if (time[row] > 0 && time[row]< 1 * ns) { //1 nanosecond window for first pulse
                if (maxV<volts[row]) {
                    maxV = volts[row];
                }
            }
        }
        
        //Double_t  maxV = TMath::MaxElement(waveform->GetN(),waveform->GetY());
                
        volt_ampl->Fill(maxV);
    delete time;
    delete volts;
    waveform = NULL;


	
    time = NULL;
    volts = NULL;
 delete waveform;
    }
       

    //Fit the first peak distribution to get pe
    //TF1 *f1 = new TF1("f1","gaus",volt_ampl->GetMinimum(),volt_ampl->GetMaximum());
    double pos_maxi = volt_ampl->GetBinCenter(volt_ampl->GetMaximumBin());
    double around = 10*volt_ampl->GetBinWidth(volt_ampl->GetMaximumBin());
    //cout << "------------" << pos_maxi << "----" << around << endl;
    TF1 *f1 = new TF1("f1","gaus",pos_maxi-around,pos_maxi+around);
    //TF1 *f1 = new TF1("f1","gaus",0,1.0); //Change range for fit of MPV
                                          //pe bigger than 0.8 wont be detected unless changed
    volt_ampl->Fit("f1","R");
    pe_volt= f1->GetParameter(1);
    
    sprintf(canvas_title,"Amplitude calculation %s, pe = %2.3f",vol_folder,pe_volt);
    volt_ampl->SetTitle(canvas_title);
    
    if (globalArgs.save_all==1) {
        volt_ampl->Write();
        //To save canvas instead of histogram
        /*TCanvas* c1= new TCanvas(canvas_title,canvas_title,100,100,900,700);
        volt_ampl->Draw();
        c1->SetGrid();
        c1->Write();*/
    }


    
    return pe_volt;
    delete volt_ampl;
volt_ampl=NULL;
	delete f1;
	f1 = NULL;
}

//Amplitude of pe calculation of raw data directly from CSV file (without TGraph)
Double_t Amplitude_calc_csv(const char* vol_folder, Int_t data_size){
    
    Double_t pe_volt;
    
    Char_t canvas_title[200];
    sprintf(canvas_title,"Amplitude calculation %s",vol_folder);
    
    //TH1D* volt_ampl= new TH1D(canvas_title, canvas_title, 150, -0.05, 1.2);
    TH1D* volt_ampl= new TH1D(canvas_title, canvas_title, 150, -0.05, 0.25);	// Adjust range of histogram, might change with amplifier used!
    
    cout<<"****----->Amplitude calculation (direct from CSV) of pe: "<< vol_folder << endl;
    
    //loop over every measurement on a folder
    for (int j=0; j<data_size; j++) {
        if (j%500==0) {
            cout<<"Event: "<<j<<endl;
        }
        
        Char_t datafilename[100];
        
        //Get the waveform
        sprintf(datafilename,"%s%s/%i.csv",globalArgs.data_folder,vol_folder,j);
        
        Double_t  maxV=0.0;
		float x(0), y(0);
		if (FILE *fp = fopen(datafilename, "r")) {
		    while (fscanf(fp, "%e,%f", &x, &y) == 2) {
		        if (x > 0 && x< 1 * ns) if (maxV<y) maxV = y;
		    }
		    fclose(fp);
		} else cout << "Error opening file " << datafilename << endl;
		
		volt_ampl->Fill(maxV);
    }
       

    //Fit the first peak distribution to get pe
    //TF1 *f1 = new TF1("f1","gaus",volt_ampl->GetMinimum(),volt_ampl->GetMaximum());
    double pos_maxi = volt_ampl->GetBinCenter(volt_ampl->GetMaximumBin());
    double around = 10*volt_ampl->GetBinWidth(volt_ampl->GetMaximumBin());
    //cout << "------------" << pos_maxi << "----" << around << endl;
    TF1 *f1 = new TF1("f1","gaus",pos_maxi-around,pos_maxi+around);
    //TF1 *f1 = new TF1("f1","gaus",0,1.0); //Change range for fit of MPV
                                          //pe bigger than 0.8 wont be detected unless changed
    volt_ampl->Fit("f1","R");
    pe_volt= f1->GetParameter(1);
    
    sprintf(canvas_title,"Amplitude calculation %s, pe = %2.3f",vol_folder,pe_volt);
    volt_ampl->SetTitle(canvas_title);
    
    if (globalArgs.save_all==1) {
        volt_ampl->Write();
        //To save canvas instead of histogram
        /*TCanvas* c1= new TCanvas(canvas_title,canvas_title,100,100,900,700);
        volt_ampl->Draw();
        c1->SetGrid();
        c1->Write();*/
    }


    
    return pe_volt;
    delete volt_ampl;
volt_ampl=NULL;
	delete f1;
	f1 = NULL;
}


//Format waveform graphs
TGraph* format_graph(TGraph* waveform, char* graph_title,Double_t ymax){
    
    waveform->SetTitle(graph_title);
    waveform->GetYaxis()->SetRangeUser(-0.1,ymax);
    waveform->GetXaxis()->SetRangeUser(-10*ns,180*ns);
    waveform->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
    waveform->GetYaxis()->SetTitleOffset(1.3);
    waveform->GetXaxis()->SetTitle("Time [s]");
    
    return waveform;
 delete waveform;
}


//Int_t single_plot(const char * Voltage, const Int_t event){
//
//
//    Char_t datafilename[100];
//
//    sprintf(datafilename,"%s%s/C1H%05i.csv",data_folder.Data(),Voltage,event);
//
//    TGraph* waveform = new TGraph(datafilename,"%lg %lg","/t;,");
//
//    TCanvas* single = new TCanvas("Single signal","Single signal",100,100,900,700);
//    single->cd();
//    waveform->SetTitle(datafilename);
//    waveform->GetYaxis()->SetRangeUser(-0.1,0.5);
//    waveform->GetXaxis()->SetRangeUser(-10*ns,80*ns);
//    waveform->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
//    waveform->GetYaxis()->SetTitleOffset(1.3);
//    waveform->GetXaxis()->SetTitle("Time [s]");
//    waveform->Draw("AL");
//    single->SetGrid();
//
//    sprintf(datafilename,"C1H%05i.pdf",event);
//    single->Print(datafilename,"pdf");
//
//    return 0;
//
//
//}

Double_t FindTimeRange(Int_t ROWS_DATA, Double_t *time, Double_t *volts, vector <Double_t> reject_time_v, vector <Double_t> time_dist_th_v, int i, int j, Double_t pe){
//mari selecting peak

            TString strin; 
            strin.Form("%d",j); 
            //TString stri = "test_"+strin+".pdf";
            TString spec = "Spectrum_"+strin;
            TH1D* Spectrum = new TH1D (spec,spec,1000,0.,time[ROWS_DATA-1]);
	    double ns=0.000000001;
            for (int row = 0; row < ROWS_DATA; row++)Spectrum->Fill(time[row],volts[row]);
   
            double tmp=0;
            double time_of_max_first_fit=0;
            double sig_max_first_fit=0;
            double sig_max_fit=0;
	    
            TSpectrum *s = new TSpectrum();
            s->Search(Spectrum,1,"",0.2);
            Double_t *peaks;
            peaks=NULL;
            Int_t Npeaks;
            Double_t *amppeaks;
            amppeaks=NULL;
            peaks=s->GetPositionX();
            amppeaks=s->GetPositionY();
            Npeaks=s->GetNPeaks();
            if(Npeaks>1 && peaks[1] > reject_time_v.at(i)*ns){
            sig_max_fit=amppeaks[1];
            tmp=peaks[1];
            if(amppeaks[1]>time_dist_th_v.at(i) * pe){
            time_of_max_first_fit = tmp;
            sig_max_first_fit = sig_max_fit;}
            }
            Spectrum= NULL;
            delete Spectrum;
            s= NULL;
            delete s;

return tmp;
//mari selecting peak end

}

#endif /* Noise_analysis_h */
