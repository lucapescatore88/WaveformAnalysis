#ifndef RandomOverlap_analysis_h
#define RandomOverlap_analysis_h

#include "genparam.h"
#include "fits.h"
#include "clustering_simulation.h"

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
#include <TROOT.h>
#include <TStyle.h>
#include <TObject.h>
#include <TPave.h>
#include <TPaletteAxis.h>

using namespace std;

vector <TString> readSetupFile_ro(ifstream * setupFile, int * data_size, double * VBD, double * VBD_error, double * integration_time, int * integration_number, TString * pe_file, int * Nsim)
{
    string s;
    vector <TString> vol_folders;
    (*Nsim) = 1000;
    
    while (true) 
    {    
        
        getline((*setupFile), s);
        if ((*setupFile).eof()) break;
        
        const char* searchString = s.c_str();
        char volt [20];
        char file [256];
        Int_t numfiles;
        double tmp_d;
        
        if (s.find("#") == 0 || s=="") continue; // Skip commented  or empty lines
        
        //Find the voltages
        char pe_config[256] = "";
        if(sscanf(searchString, "V || %s ||", volt)==1)
            vol_folders.push_back(volt);
        
        if(sscanf(searchString, "Files at each voltage || %d ||", &numfiles)==1) 
            (*data_size) = numfiles;
        
        if(sscanf(searchString, "VBD : %lf", &tmp_d)==1)
            (*VBD) = tmp_d;
        if(sscanf(searchString, "VBD_error : %lf", &tmp_d)==1)
            (*VBD_error) = tmp_d;
            
        if(sscanf(searchString, "Integration_time_step : %lf", &tmp_d)==1)
            (*integration_time) = tmp_d;
        if(sscanf(searchString, "Integration_time_number : %d", &numfiles)==1)
            (*integration_number) = numfiles;
        
        if(sscanf(searchString, "pe_file : %s", file)==1)
            (*pe_file) = file;
        
        if(sscanf(searchString, "nEventSimulation : %d", &numfiles)==1)
            (*Nsim) = numfiles;
    }
    
    for(unsigned int i(0); i<vol_folders.size(); ++i)
	    cout << "Voltage: " << vol_folders[i] << endl;

    return vol_folders;
}

TGraphErrors * getPeFromFile(TString pe_file, TString vol, vector<double>& OnePECharge, vector<double>& OnePECharge_error) {
    TFile * f = new TFile(pe_file);
    //~ f->ls();
    f->cd("pulse_shape_"+vol);
    TGraphErrors * pe = (TGraphErrors*) gDirectory->Get("1PE_Charge_Time_Window_"+vol);
    if(pe->GetN() != OnePECharge.size()) {
        cout << "Number of integration steps does not match the PE charge graph size !" << endl;
        return 0;
    }
    for(unsigned int i(0); i<pe->GetN(); ++i) {
        double x;
        pe->GetPoint(i,x,OnePECharge[i]);
        OnePECharge[i] = OnePECharge[i]*1e-3*1e-9;
        OnePECharge_error[i] = pe->GetErrorY(i)*1e-3*1e-9;
    }
    f->Close();
    return pe;
}

void initChargeIntegrationTime_dist(int n, double step, vector<TH1D*>& hs, double tau_dcr, double dV, vector<double> pe_charge, vector<double> scale, TString unit="[?]") {
    hs.clear();
    int nbins = 100;
    
    for(unsigned int s(0); s<n; ++s) {
        double integration_window = (s+1)*step;
        double expected_N_pulses = integration_window/tau_dcr;
        TString name = Form("ChargeIntegratedOver%2.0lfns_dV=%2.2lfV",integration_window,dV);
        TString title = Form("Charge integrated over %2.0lf ns at #DeltaV = %2.2fV",integration_window,dV);
        TH1D * h = new TH1D(name, title, nbins, 0, 10*TMath::Sqrt(expected_N_pulses)*(pe_charge[s]/scale[s]));
        h->GetXaxis()->SetTitle("Integrated charge "+unit);
        h->GetYaxis()->SetTitle("Number of entries");
        //cout << "Created histogram with title: " << title << endl;
        hs.push_back(h);
    }
    
    return;
}

double computeIntegral_by_section(double time[], double amp[], int npts, double time_step, TH1D * h, double pe_scale=1, double baseline_shift=0, int startInd=0) {
    double tot_integral = 0.0;
    double pulse_integral = 0.0;
    double t0 = time[startInd];
    int counter(0);
    for (int row = startInd+1; row < npts; row++) {
        
        if( (time[row] - t0) > time_step*ns) {
            h->Fill(pulse_integral/pe_scale);
            tot_integral += pulse_integral;
            pulse_integral = 0.0;
            ++counter;
            t0 = time[row];
        } else {
            pulse_integral += 0.5 * (time[row]-time[row-1]) *(amp[row]+amp[row-1] - 2*baseline_shift);
        }

    }
    return (tot_integral/pe_scale)/counter;     // average charge in time window
}

void computeALLIntegral_by_section(double time[], double amp[], int npts, int n_int, double step, vector<TH1D*> hs, vector<double>& average_charge, vector<double> pe_scale, double baseline_shift=0, int startInd=0) {
    average_charge.clear();
    // integration per section
    for(unsigned int s(0); s<n_int; ++s) {
        double integration_window = (s+1)*step;
        average_charge.push_back(computeIntegral_by_section(time,amp,npts,integration_window,hs[s],pe_scale[s],baseline_shift,startInd));
    }
    return;
}

void initOutputFile_ro(TFile * f, vector<TString> vol_folders) {
	for(unsigned int i(0); i<vol_folders.size(); ++i) {
		TString dirname = "charge_analysis_" + vol_folders[i];
		TDirectory *sub = f->mkdir(dirname);
        sub->cd();
        TDirectory *subsub1 = sub->mkdir("Dir_NCR_vs_seed");
        TDirectory *subsub2 = sub->mkdir("Dir_NCR_vs_time_window");
	}
	return;
}

//~ int globcount = 0;
double evaluateBaselineShift_ro(unsigned int n, double * yy, double range) {
    TH1D * h = new TH1D("tmp_h","tmp_h",100,-1*range,range);
    h->FillN(n,yy,NULL);
    //~ if(globcount%50 == 0) {
		//~ TFile * f = new TFile("baseline_shift_eval.root","UPDATE");
	    //~ h->Write();
	    //~ f->Close();
	    //~ delete f;
	//~ }
	//~ ++globcount;
    double error, baseline;
    //~ baseline = fitSimple1D(h,error);
    // most probable value
    int mpv_bin = h->GetMaximumBin();
    double mpv = h->GetBinCenter(mpv_bin);
    double mpv_content = h->GetBinContent(mpv_bin);
    // one bin before
    bool not_found = true;
    int bin = mpv_bin;
    double prev, prev_content;
    while(not_found && bin>0) {
        --bin;
        prev = h->GetBinCenter(bin);
        prev_content = h->GetBinContent(bin);
        if(prev_content>0) not_found = false;
    }
    // one bin after
    double next, next_content;
    not_found = true;
    bin = mpv_bin;
    while(not_found && bin<=h->GetNbinsX()) {
        ++bin;
        next = h->GetBinCenter(bin);
        next_content = h->GetBinContent(bin);
        if(next_content>0) not_found = false;
    }
    //~ baseline = h->GetBinCenter(h->GetMaximumBin());
    // Ponderate with the two nearest neighbours
    baseline = (prev_content*prev + mpv_content*mpv + next_content*next) / (prev_content + mpv_content + next_content);
    delete h;
    return baseline;
}

int findIndex(double value, double * monotone_array, int size) {
    bool not_found = true;
    int ind = 0;
    while(not_found && ind<size) {
        if(monotone_array[ind] < value) ++ind;
        else not_found = false;
    }
    return ind;
}

void printIntegrals(int n_int, double step, vector<double> integrals, TString unit, ostream& out) {
    if(n_int != integrals.size()) {
        cerr << "Dimensions don't match!" << endl;
        return;
    }
    out << "Computed integrals:" << endl;
    for(unsigned int i(0); i<n_int; ++i) {
        out << "Charge for time window " << (i+1)*step << "ns " << unit << ": ";
        out << Form("%e",integrals[i]) << endl;
    }
    out << endl;
    return;
}

void printArray(double data[], int n, ostream& out) {
    out << endl << "[";
    for(unsigned int i(0); i<n-1; ++i) {
        out << data[i] << ", ";
    }
    out << data[n-1] << "]";
    out << endl;
    return;
}

// read CSV file with DCR values, measured by the oscilloscope
double getDCR_from_CSVfile(double bias, string filename) {
    ifstream file ( filename );
    
    std::string line;
    double vbias, dcr, current;
    
    bool not_found = true;
    while (std::getline(file, line) && not_found) {
        if(line[0] == '#') continue;
        if(sscanf(line.c_str(),"%lf\t%lf\t%le", &vbias, &dcr, &current) == 3)
            if(vbias == bias)
                not_found = false;
    }
    cout << Form("Found DCR = %2.0lf kHz for Vbias = %2.2lf V", dcr, bias) << endl;
    file.close();
    
    return dcr;
}

#endif 
