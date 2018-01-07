//
//  Inspired by Noise_analysis_largeEvent.C
//  Olivier Girard, olivier.girard@cern.ch on 01.2018
//

#include <TPaveStats.h>
#include "Noise_analysis_largeEvent.h"
#include "RandomOverlap_analysis.h"
#include "Thresholds.h"
#include "fits.h"
#include "lhcbStyle.h"

using namespace std;

int main(int argc, char* argv[]) 
{
    lhcbstyle();
    int fillStyle = 3003;
    //gStyle->SetOptStat(0);
    
    // Get paremeters from the command line
    setOptions(argc, argv);
    
    ifstream setupFile(globalArgs.arg_pathToSetupFile);
    if(!setupFile){
        std::cerr << "Failure: could not open file: \"" << globalArgs.arg_pathToSetupFile << "\"." << std::endl;
        std::cerr << "Please check if the path is correct or not!" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Read setup file and load input .root
    
    Int_t data_size;
    double VBD, VBD_error;
    vector<double> OnePECharge, OnePECharge_error;
    double integration_time;
    int integration_number, Nsim;
    vector <TString> vol_folders = readSetupFile_ro(&setupFile,&data_size,&VBD,&VBD_error,&integration_time,&integration_number,&OnePECharge,&OnePECharge_error,&Nsim);
    const int vol_size = vol_folders.size();
    TFile * ifile =  TFile::Open(TString(globalArgs.data_folder)+"/oscilloscope_out.root");
    
    // Define output objects

    map<string, TGraphErrors *>  results {
        {"1PE-Charge-Fit",        new TGraphErrors()} };
    
    Double_t dV;
    
    TFile * hfile = TFile::Open(TString(globalArgs.res_folder) + "randomoverlap_analysis.root","RECREATE");
    initOutputFile_ro(hfile, vol_folders);
    
    TTree * otree = 0;
    otree = new TTree("ClassifiedData","ClassifiedData");
    int npts;
    double times[500002], amps[500002];
    double Vbias, baseline_shift;
    vector<double> integral_per_section;
    
    otree->Branch("dV",&dV,"dV/D");
    otree->Branch("NsampPerEv",&npts,"NsampPerEv/I");
    otree->Branch("Amps",&amps,"Amps[NsampPerEv]/D");
    otree->Branch("Times",&times,"Times[NsampPerEv]/D");
    otree->Branch("Vbias",&Vbias,"Vbias/D");
    otree->Branch("integral_per_section",&integral_per_section,"integral_per_section[integration_number]/D");
    otree->Branch("baseline_shift",&baseline_shift,"baseline_shift/D");
    
    // Scale chosen for charge measurement
    // [V#timess], [mV#timesns] or [1PE]
    //~ TString charge_unit = "[?]";
    TString charge_unit = "[1PE]";
    //~ TString charge_unit = "[V#timess]";
    //~ TString charge_unit = "[mV#timesns]";
    
    /////////////////
    // Loop over all Voltages measured
    /////////////////
    
    cout << "\n\n-----> Random Overlap analysis *** " << endl;
    
    TGraph * waveform = NULL;
    hfile->cd();
        
    for (int i = 0; i < vol_size; i++)
    {
        const char * vol = vol_folders[i];
        cout << "\n\n----> Voltage analyzed: " << vol << endl;
    
        // Counters 
    
        unsigned int events_cnt = 0;
        
        Vbias = vol_folders[i].Atof();
        dV = Vbias - VBD;
        
        // DCR and expectations on average charge included in time window
        double dcr = getDCR_from_CSVfile(Vbias, string(globalArgs.data_folder)+"/DCR_for_RandomOverlap.txt");
        double tau_dcr = (1./dcr)*1e6;      // average time spacing in ns
        
        // Defines scale from charge unit
        double charge_scale = 1;
        if(charge_unit == "[V#timess]") charge_scale = 1;
        if(charge_unit == "[mV#timesns]") charge_scale = 1e-3*1e-9;
        if(charge_unit == "[1PE]") charge_scale = OnePECharge[i];
        
        vector<TH1D*> histograms_section_integral;
        initChargeIntegrationTime_dist(integration_number, integration_time, histograms_section_integral, tau_dcr, dV, charge_scale/OnePECharge[i], charge_unit);
        
        TGraphErrors * graphs_charge_vs_integration_time = new TGraphErrors();
        graphs_charge_vs_integration_time->SetName(Form("ChargeVsIntTime_dV=%2.2lfV",dV));
        
        // Setup input tree
        TTree * tree = NULL;
        if(globalArgs.input=="root") 
        {   
            tree = (TTree*)ifile->Get(vol);
            tree->SetBranchAddress("NsampPerEv",&npts);
            tree->SetBranchAddress("Amps",&amps);
            tree->SetBranchAddress("Times",&times);
            
            if(data_size>tree->GetEntries()) data_size = tree->GetEntries();
            
        }
        hfile->cd();
        
        // Find index for integration
        // trigger is at t=0, so we start a bit later to remove any correlation
        double timeStart = 5000*ns; // 5 us
        tree->GetEntry(0);
        int startInd = findIndex(timeStart, times, npts);
    
        // Loop over every measurement on a folder
        for (int j = 0; j < data_size; j++)
        {
            events_cnt++;
    
            if(globalArgs.input=="root") 
            {
                tree->GetEntry(j);
                waveform = new TGraph(npts,times,amps);
            }
            else  // Read from csv file
            {
                TString datafilename = Form("%s%s/%i.csv",globalArgs.data_folder,vol,j);
                waveform = new TGraph(datafilename,"%lg %lg","/t;,");
            }
            if (waveform->IsZombie()) continue;
    
            waveform->SetName(Form("%s_%i",vol,j));
    
            Double_t * time  = waveform->GetX();
            Double_t * volts = waveform->GetY();
            
            // Evaluate baseline - needed for a correct integral computation
            // baseline evaluation in time interval [mtime ; ptime]
            double mtime = -200*ns;
            double ptime = -5*ns;
            
            // Use a fit for baseline eval
            TF1 * fcste = new TF1("baseline_fit","pol0",mtime,ptime);
            waveform->Fit("baseline_fit","QNR");
            baseline_shift = fcste->GetParameter(0);
            delete fcste;
            
            // compute integral by section
            computeALLIntegral_by_section(time, volts, npts, integration_number, integration_time, histograms_section_integral, integral_per_section, baseline_shift, startInd, charge_scale);
            
            if(globalArgs.save_tree) otree->Fill();
            
            if(!(j%1000)) {
                cout << j << " / " << data_size << endl;
                printIntegrals(integration_number, integration_time, integral_per_section, charge_unit, std::cout);
            }
            delete time;
            delete volts;
        }
        
        // Integrated charge vs integration time
        for(unsigned int s(0); s<integration_number; ++s) {
            double integration_window = (s+1)*integration_time;
            double mean = histograms_section_integral[s]->GetMean();
            double mean_error = histograms_section_integral[s]->GetMeanError();
            //~ double mean_error = histograms_section_integral[s]->GetRMS();
            setPoint(graphs_charge_vs_integration_time, s, integration_window, mean, 0, mean_error);
        }
        
        // Clustering simulation
        unsigned int NSteps = 17;   // seed threshold steps (0.25PE steps, strating at 1.5PE)
        vector<TGraphErrors*> ncr_vs_seed(integration_number);
        vector<TGraphErrors*> ncr_vs_time_window(int(NSteps/4)+1);
        
        for(unsigned int s(0); s<integration_number; ++s) {
            double integration_window = (s+1)*integration_time;
            // Create cumulative distribution
            TGraph * inv_cumu = getCumulativeDistribution(histograms_section_integral[s]);
            //~ cumulative_distr.push_back(inv_cumu);
            
            // NCR vs seed plot for one time window
            TGraphErrors * g_ncr = new TGraphErrors();
            g_ncr->SetName(Form("NCR_vs_seed_%2.0lfns_%2.2lfV",integration_window,dV));
            
            // Number of clusters distribution
            vector<TH1I*> nbCluster_distr;
            for(unsigned int step_nb(0); step_nb<NSteps; ++step_nb) {
                double seed = 1.5 + step_nb*0.25;
                TString name = Form("nbClusters_%2.2lf",seed);
                TH1I * hnbC = new TH1I(name,name,100,0,100);
                nbCluster_distr.push_back(hnbC);
            }
            
            cout << "\t Generating NCR from simulation for " << Form("%2.0lfns, %2.2lfV",integration_window,dV) << endl;
            // Generate randomly a SiPM array and clusterize
            for(unsigned int n(0); n<Nsim; ++n) {
                if(n%10000==0) cout << "\t\t simulation: " << n << " / " << Nsim << endl;
                double data[128];
                generateArray(inv_cumu, data, 1);
                //~ printArray(data, 128, std::cout);
                for(unsigned int step_nb(0); step_nb<NSteps; ++step_nb) {
                    double thresholds[3];
                    thresholds[0] = 1.5 + step_nb*0.25;
                    thresholds[1] = thresholds[0] - 1.0;
                    thresholds[2] = thresholds[0] + 2.0;
                    vector<double> mean_pos, sum, size;
                    int Nclusters = Cluster_Search(data, thresholds, mean_pos, sum, size);
                    nbCluster_distr[step_nb]->Fill(Nclusters);
                }
            }
            
            // Evaluate noise cluster rate
            // Assuming 40 MHz read-out
            for(unsigned int step_nb(0); step_nb<NSteps; ++step_nb) {
                double seed = 1.5 + step_nb*0.25;
                int nClusters = 0;
                for(unsigned int bin(2); bin<=nbCluster_distr[step_nb]->GetNbinsX(); ++bin) {
                    // total number of clusters
                    nClusters += (bin-1)*nbCluster_distr[step_nb]->GetBinContent(bin);
                }
                double ncr = ((double) nClusters/(double) Nsim)*40.0;
                double ncr_error = (40.0/(double) Nsim)*TMath::Sqrt((double) nClusters);
                setPoint(g_ncr, step_nb, seed, ncr, 0, ncr_error);
                delete nbCluster_distr[step_nb];
                
                if(int(seed+0.5) == seed+0.5) {
                    if(!ncr_vs_time_window[int(step_nb/4)]) {
                        ncr_vs_time_window[int(step_nb/4)] = new TGraphErrors(integration_number);
                        ncr_vs_time_window[int(step_nb/4)]->SetName(Form("NCR_vs_time_window_%2.2lfPE_%2.2lfV",seed,dV));
                    }
                    setPoint(ncr_vs_time_window[int(step_nb/4)], s, integration_window, ncr, 0, ncr_error);
                }
            }
            ncr_vs_seed[s] = g_ncr;
        }
        TCanvas * cNCR_vs_seed = plotNCR(ncr_vs_seed, "seed");
        TCanvas * cNCR_vs_time_window = plotNCR(ncr_vs_time_window, "time_window");
        
        cout << "Saving objects" << endl;
    
        // Save/print results:
        vector<TObject*> objects_to_save;
        // Charge histogram for all integration windows
        //~ for(unsigned int hi(0); hi<histograms_section_integral.size(); ++hi) histograms_section_integral[hi]->Scale(1.0/histograms_section_integral[hi]->Integral());
        //~ for(unsigned int hi(0); hi<histograms_section_integral.size(); ++hi) objects_to_save.push_back(histograms_section_integral[hi]);
        //~ for(unsigned int hi(0); hi<histograms_section_integral.size(); ++hi) objects_to_save.push_back(histograms_section_integral[hi]->GetCumulative());
        
        formatGr(graphs_charge_vs_integration_time, kBlue, fillStyle, "#tau_{int} [ns]", "Q_{int} [1PE]", "Integrated charge vs integration time");
        objects_to_save.push_back(graphs_charge_vs_integration_time);
        
        //~ for(unsigned int c(0); c<cumulative_distr.size(); ++c) objects_to_save.push_back(cumulative_distr[c]);
        objects_to_save.push_back(cNCR_vs_seed);
        objects_to_save.push_back(cNCR_vs_time_window);
        
        TString dirname = "charge_analysis_" + TString(vol);
        for(auto obj : objects_to_save) 
        {
            hfile->cd();
            hfile->cd(dirname);
            //cout << obj->GetName() << endl;
            obj->Write();
            hfile->cd();
            delete obj;
        }
        cout << "Done" << endl;
        
    }
    
    hfile->Close();
    return 0;
}


