//
//  Created by Felipe Gilberto Ortega on 25/04/16.
//  Modified by Luca Pescatore, luca.pescatore@cern.ch on 28/06/2017
//  Modified by Olivier Girard, olivier.girard@cern.ch on 07-08/2017
//
 
#include "Noise_analysis_largeEvent.h"
#include "Thresholds.h"
#include "fits.h"

using namespace std;

int main(int argc, char* argv[]) 
{
    gROOT->ProcessLine(".x Analysis/lhcbStyle.C");
    gStyle->SetOptStat(0);
    
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
    vector <Thresholds> thrs;
    vector <TString> vol_folders = readSetupFile(&setupFile,&data_size,thrs);
    const int vol_size = vol_folders.size();
    
    TFile * ifile =  TFile::Open(TString(globalArgs.data_folder)+"/oscilloscope_out.root");
    
    // Define output objects
    
    map<string, TGraphErrors *>  results { 
        {"#DiXT",          new TGraphErrors()},
        {"#DeXT",          new TGraphErrors()},
        {"#AP",            new TGraphErrors()},
        {"#Total",         new TGraphErrors()},
        {"#DCR",           new TGraphErrors()},
        {"#SecPeaks",      new TGraphErrors()},
        {"#SecPeaksDiXT",  new TGraphErrors()},
        {"#SecPeaksDel",   new TGraphErrors()},
        {"Charge",         new TGraphErrors()},
        {"#Double",        new TGraphErrors()},
        {"#CorrectionFre", new TGraphErrors()},
        {"#CorrectionCur", new TGraphErrors()} };
    
    Char_t Category[15];
    Double_t dV;
    
    ofstream * values_for_pde = new ofstream(globalArgs.res_folder + "values_for_pde.txt",ofstream::out);
    
    TFile * hfile = TFile::Open(TString(globalArgs.res_folder) + "noiseanalysis.root","RECREATE");
    initOutputFile(hfile, vol_folders);
    
    TTree * otree = new TTree("ClassifiedData","ClassifiedData");
    int npts, noise_peaks_cnt;
    int direct_xtalk_pulse, xtalk_pulse, after_pulse, secondary_pulse;
    double times[10000], amps[10000];
    double Vbias, pe, DiXT_thr, DeXT_thr, AP_thr, pulse_integral, baseline_shift;
    
    otree->Branch("Category",Category,"Category/C");	// Categories are: DiXT, DeXT, AP, Clean
    otree->Branch("dV",&dV,"dV/D"); 
    otree->Branch("NsampPerEv",&npts,"NsampPerEv/I");
    otree->Branch("Amps",&amps,"Amps[NsampPerEv]/D");
    otree->Branch("Times",&times,"Times[NsampPerEv]/D");
    otree->Branch("NnoisePeaks",&noise_peaks_cnt,"NnoisePeaks/I");
    otree->Branch("NAP",&after_pulse,"NAP/I");
    otree->Branch("NDeXT",&xtalk_pulse,"NDeXT/I");
    otree->Branch("NSec",&secondary_pulse,"NSec/I");
    otree->Branch("Vbias",&Vbias,"Vbias/D");
    otree->Branch("pe",&pe,"pe/D");
    otree->Branch("DiXT_thr",&DiXT_thr,"DiXT_thr/D");
    otree->Branch("DeXT_thr",&DeXT_thr,"DeXT_thr/D");
    otree->Branch("AP_thr",&AP_thr,"AP_thr/D");
    otree->Branch("pulse_integral",&pulse_integral,"pulse_integral/D");
    otree->Branch("baseline_shift",&baseline_shift,"baseline_shift/D");
    
    /////////////////
    // Calculate Voltage breakdown and value of pe
    /////////////////
    
    vector <Double_t> pe_volt;
    TGraph * Vbias_vs_pe = new TGraph();
    
    std::cout << " -----> Calculation of average DCR amplitudes " << std::endl;
    for(unsigned int v(0); v<vol_size; ++v)
    {
        const char * vol = vol_folders[v];
        pe_volt.push_back(Amplitude_calc(vol, data_size, "root", hfile));
        //~ pe_volt.push_back(Amplitude_calc(vol, data_size, "root", hfile));
        Vbias_vs_pe->SetPoint(v, pe_volt.back(), TString(vol).Atof());
        std::cout << Form("Voltage: %.1fV --> Mean amplitude = %.6f", TString(vol).Atof(), pe_volt.back()) << std::endl;
    }
    
    cout << "\n\n-----> Voltage Breakdown fit" << endl;
    
    double VBD = fitBreakdownVoltage(Vbias_vs_pe, hfile, values_for_pde);
    
    if(globalArgs.only_Vbd) {hfile->Close(); values_for_pde->close(); return 0;}	// only Vbd measurement, abort other measurements
    
    
    /////////////////
    // Loop over all Voltages measured
    /////////////////
    
    cout << "\n\n-----> Noise analysis *** " << endl;
    
    TGraph * waveform = NULL;
    hfile->cd();
    
    //for (int i = vol_size-1; i < vol_size; i++) // Test on highest and more noisy voltage
    for (int i = 0; i < vol_size; i++)
    {    
        const char * vol = vol_folders[i];
        cout << "\n\n----> Voltage analyzed: " << vol << endl;
    
        map<string, int> color{{"clean",1},{"AP",1},{"DiXT",1},{"DeXT",1}};
    
        // Counters 
    
        unsigned int events_cnt = 0;
        unsigned int tot_noise_peaks_cnt = 0;
        unsigned int tot_primary_noise_peaks_cnt = 0;
        unsigned int counter_notclean = 0;
        unsigned int direct_xtalk_pulse_cnt = 0;
        unsigned int xtalk_pulse_cnt = 0;
        unsigned int after_pulse_cnt = 0;
        unsigned int secondary_pulse_cnt = 0;
        unsigned int delayed_pulse_when_primary(0), delayed_pulse_when_secondary(0);
        unsigned int nsaved(0), nsaved_DiXT(0), nsaved_DeXT(0), nsaved_AP(0);	// counters on the number of waveforms in each canvas
        unsigned int maxNwaveforms = 100; // maximum number of waveforms in the canvas
        vector<delayedPulse> delayedPulses;	// used for DeXT and AP
        double corr_noise_tot_amp(0);
        
        // Define amplitude measured at which OV
        
        Vbias = vol_folders[i].Atof();
        pe = pe_volt[i];
        dV = Vbias - VBD;
        
        // single waveforms multigraphs
        map<string,TCanvas *> canv {
            {"DiXT", new TCanvas(Form("DiXT_%s_%dwaveforms",vol,maxNwaveforms), Form("Direct cross-talk #DeltaV = %2.2f V",dV))},
            {"DeXT", new TCanvas(Form("DeXT_%s_%dwaveforms",vol,maxNwaveforms), Form("Delayed cross-talk #DeltaV = %2.2f V",dV))},
            {"AP",   new TCanvas(Form("AP_%s_%dwaveforms",vol,maxNwaveforms), Form("After-pulse #DeltaV = %2.2f V",dV))},
            {"clean",new TCanvas(Form("Clean_%s_%dwaveforms",vol,maxNwaveforms), Form("Clean #DeltaV = %2.2f V",dV))}
        };
        
        TMultiGraph * forfit    = new TMultiGraph("multigraph_for_fit","multigraph_for_fit");  
        TGraph * Expfit_AP      = new TGraph();
        TGraph * cleanforfit    = NULL;
        map<string,TH1D*> graphs {
           {"AP_arrivaltime",   new TH1D(Form("AP_arrival_times_%s",vol),"After-pulse arrival times", 140, 0, 0.2e-6)},
           {"DeXT_arrivaltime", new TH1D(Form("DeXT_arrival_times_%s",vol),"Delayed cross-talk arrival times", 140, 0, 0.2e-6)},
           {"Npeaks",           new TH1D(Form("N_peaks_%s",vol),"Number of noise peaks when not clean", 50, 0, 50)},
           {"NpeaksDiXT",       new TH1D(Form("N_peaks_DiXT_%s",vol),"Number of noise peaks when direct DiXT", 50, 0, 50)},
           {"NpeaksDel",        new TH1D(Form("N_peaks_Delayed_%s",vol),"Number of noise peaks when delayed noise", 50, 0, 50)},
           {"Charge",           new TH1D(Form("Charge_%s",vol),"Charge of all peaks", 1000, 0., ns*100*pe)}
        };
        // Canvas for time distribution fits
        TCanvas * timeDistAP   = new TCanvas(Form("canv_AP_arrival_times_%s",vol), "After-pulse arrival times");
        TCanvas * timeDistDeXT = new TCanvas(Form("canv_DeXT_arrival_times_%s",vol), "Delayed cross-talk arrival times");
        // Trees for unbinned time distribution (used with RooFit)
        TTree * timeAP   = new TTree(Form("tree_AP_arrival_times_%s",vol),"After-pulse arrival times");
        TTree * timeDeXT = new TTree(Form("tree_DeXT_arrival_times_%s",vol),"Delayed cross-talk arrival times");
        double tAP(0), tDeXT(0);
        timeAP->Branch("time",&tAP,"time/D"); 
        timeDeXT->Branch("time",&tDeXT,"time/D");
        
        // persistence waveforms
        map<string,TCanvas *> canv_persistence {
            {"DiXT", new TCanvas(Form("persistence_DiXT_%s",vol), Form("Persistence direct cross-talk #DeltaV = %2.2f V",dV))},
            {"DeXT", new TCanvas(Form("persistence_DeXT_%s",vol), Form("Persistence delayed cross-talk #DeltaV = %2.2f V",dV))},
            {"AP",   new TCanvas(Form("persistence_AP_%s",vol), Form("Persistence after-pulse #DeltaV = %2.2f V",dV))},
            {"clean",new TCanvas(Form("persistence_Clean_%s",vol), Form("Persistence Clean #DeltaV = %2.2f V",dV))}
        };
        
        double bin_size = 0.002;	// vertical bin size for persistence plot
        TH2D * persistence_clean = NULL;
        TH2D * persistence_DiXT  = NULL;
        TH2D * persistence_DeXT  = NULL;
        TH2D * persistence_AP    = NULL;
        
        double hmin, hmax, vmin, vmax;
        // Setup input tree
        TTree * tree = NULL;
        if(globalArgs.input=="root") 
        {   
            tree = (TTree*)ifile->Get(vol);
            tree->SetBranchAddress("NsampPerEv",&npts);
            tree->SetBranchAddress("Amps",&amps);
            tree->SetBranchAddress("Times",&times);
            data_size = tree->GetEntries();
            
            int Nh = tree->GetMaximum("NsampPerEv");
            hmin = tree->GetMinimum("Times");
            hmax = tree->GetMaximum("Times");
            int Nv = int(((tree->GetMaximum("Amps")-tree->GetMinimum("Amps"))/bin_size)+0.5);
            vmin = tree->GetMinimum("Amps");
            vmax = tree->GetMaximum("Amps");
            persistence_clean = new TH2D("Clean_persistence","Clean waveforms", Nh, hmin, hmax, Nv, vmin,vmax);
            persistence_DiXT  = new TH2D("DiXT_persistence","Direct cross-talk", Nh, hmin, hmax, Nv, vmin,vmax);
            persistence_DeXT  = new TH2D("DeXT_persistence","Delayed cross-talk", Nh, hmin, hmax, Nv, vmin,vmax);
            persistence_AP    = new TH2D("AP_persistence","After-pulse", Nh, hmin, hmax, Nv, vmin,vmax);
        }
        hfile->cd();
    
        Thresholds cthrs = thrs[i];
    
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
    
            /////////////////////////////////////////////////////
            // Data filtering into the different type of events
            // direct x-talk  AP   delayed x-talk
            /////////////////////////////////////////////////////
            noise_peaks_cnt = 0;
            after_pulse = 0;
            xtalk_pulse = 0;
            direct_xtalk_pulse = 0;
            secondary_pulse = 0;
            pulse_integral = 0;
            baseline_shift = 0;
            delayedPulses.clear();
            
            DiXT_thr = cthrs.dir_xtalk * pe;
            DeXT_thr = cthrs.del_xtalk * pe;
            AP_thr = cthrs.AP * pe;
            if(globalArgs.fixed_thr > 0.) 
            {
                DiXT_thr = globalArgs.fixed_thr;
                DeXT_thr = globalArgs.fixed_thr;
                AP_thr = globalArgs.fixed_thr;
            }
            double min_thrs = !(DiXT_thr<DeXT_thr)?DeXT_thr:DiXT_thr;
            min_thrs = !(min_thrs<AP_thr)?AP_thr:min_thrs;
            
            double dt = time[1] - time[0];
            // Evaluate baseline - needed for a correct integral computation
            unsigned int bsl_npts(0);
            bool zero_reached(false);
            while(!zero_reached) {
                baseline_shift += volts[bsl_npts];
                if(time[bsl_npts]>=-1*dt) zero_reached = true;
                bsl_npts++;
            }
            if(bsl_npts) baseline_shift = baseline_shift/bsl_npts;
            else baseline_shift = 0;
                
            bool done = false;
            bool veryclean = true;
            double lastPulseT = 0;	// used to apply an analysis dead time after each pulse (to avoid noise and ringing)
            //double deadTime = 5*ns;	// this is a dead time for the analysis (not of the detector)
            double deadTime = cthrs.del_xtalk_minT * ns;	// deadtime set to min delay time of DeXT (optimised to avoid ringing)
            //for (int row = 2; row < npts-2; row++) 
            for (int row = bsl_npts; row < npts-2; row++) // loop start when time[row] = 0ns
            {
                double curT = time[row];
                double curV = volts[row];
                
                // Compute pulse integral as the sum of area of trapezoids
                //if(curV-baseline_shift>0) pulse_integral += 0.5 * dt *(curV+volts[row-1] - 2*baseline_shift);
                pulse_integral += 0.5 * dt *(curV+volts[row-1] - 2*baseline_shift);
    
                // Skip if not on a maximum
                if(!( curV > volts[row-2] && curV > volts[row+2] && curV > min_thrs)) continue;
                else if (done) { done = false; continue; }
                done = true;
                
                if ( curT > 2*ns && curV > 0.4*pe) veryclean = false;
    
                // Direct x-talk: in 0-2 ns window and V > direct th.
                if( curT <= cthrs.dir_xtalk_maxT * ns && curV > DiXT_thr) 
                {
                    direct_xtalk_pulse++;
                    //noise_peaks_cnt++; // Don't count peak as it is over the DCR
                    delayedPulse dixt;
                    dixt.time = 0;
                    dixt.volt = curV;
                    dixt.type = "DiXT";
                    delayedPulses.push_back(dixt);
                    corr_noise_tot_amp += curV;
                }
                
                // we look for delayed peaks (AP, DeXT)
                // Delayed x-talk: time larger than end of DCR window and V > delayed XT th
                else if ( curT > cthrs.del_xtalk_minT * ns && curV > DeXT_thr && curT > lastPulseT+deadTime ) 
                {
                    xtalk_pulse++;
                    noise_peaks_cnt++;
                    //graphs["DeXT_arrivaltime"]->Fill(curT);
                    delayedPulse dext;
                    dext.time = curT;
                    dext.volt = curV;
                    dext.type = "DeXT";
                    delayedPulses.push_back(dext);
                    lastPulseT = curT;
                    corr_noise_tot_amp += curV;
                }
                // After-pulse: time larger 2ns and V > AP th
                else if ( curT > cthrs.AP_minT * ns && curV > AP_thr && curT > lastPulseT+deadTime ) 
                {	// if it's not a DeXT, it might be an AP
                    after_pulse++;
                    noise_peaks_cnt++;
                    //Expfit_AP->SetPoint(Expfit_AP->GetN(),curT,curV);
                    //graphs["AP_arrivaltime"]->Fill(curT);
                    delayedPulse ap;
                    ap.time = curT;
                    ap.volt = curV;
                    ap.type = "AP";
                    delayedPulses.push_back(ap);
                    lastPulseT = curT;
                    corr_noise_tot_amp += curV;
                }
                else if ( curT <= lastPulseT+deadTime ) continue; // skip pulses that might be due to ringing
                
                //cout << "Pulse at t=" << curT << "   dead time until t=" << lastPulseT+deadTime << endl;
    
            } // loop over time
            
            //------------------ Clean event ---------------------------
            if( direct_xtalk_pulse + xtalk_pulse + after_pulse == 0 ) 
            {	
                sprintf(Category,"Clean");
    
                // Get only very clean waves and make an average for long tau fit.
                if(veryclean) cleanforfit = average(cleanforfit, waveform);
    
                if(nsaved < maxNwaveforms) // Max 20 clean graphs on the plot
                {
                    nsaved++;
                    drawWave(waveform, &color["clean"], Form("Clean pulse #DeltaV = %2.2f V",dV), canv["clean"], 1.5*pe);
                    forfit->Add(waveform);
                }
                fillPersistence(persistence_clean, waveform);
            }
            
            //------------------ Event with only one primary correlated noise ------------------------
            else if( direct_xtalk_pulse + xtalk_pulse + after_pulse == 1 ) 
            {
                counter_notclean++;
                tot_primary_noise_peaks_cnt++;
                
                // They are categroized as DiXT, DeXT or AP
                if (direct_xtalk_pulse > 0)  // Primary correlated noise is DiXT
                {
                    if(direct_xtalk_pulse > 1) cout << "Attention: more then one direct DiXT found" << endl;	// secondary correlated DiXT ?
                    direct_xtalk_pulse_cnt++;
                    sprintf(Category,"DiXT");
    
                    TString graph_title = Form("Direct CrossTalk #DeltaV = %2.2f V",dV);
                    if(nsaved_DiXT < maxNwaveforms) 
                    {
                        drawWave(waveform, &color["DiXT"], graph_title, canv["DiXT"], 1.5*pe);
                        nsaved_DiXT++;
                    }
                    fillPersistence(persistence_DiXT,waveform);
                    
                    graphs["Npeaks"]->Fill(noise_peaks_cnt);
                    graphs["NpeaksDiXT"]->Fill(noise_peaks_cnt);
                }
                
                else if (xtalk_pulse > 0) // Only Delayed x-talk
                {
                    xtalk_pulse_cnt++;
                    delayed_pulse_when_primary++;
                    sprintf(Category,"DeXT");
                    
                    TString graph_title = Form("Delayed cross-talk #DeltaV = %2.2f V",dV);
                    if(nsaved_DeXT < maxNwaveforms) 
                    {
                        drawWave(waveform, &color["DeXT"], graph_title, canv["DeXT"], 1.5*pe);
                        nsaved_DeXT++;
                    }
                    fillPersistence(persistence_DeXT,waveform);
                    
                    graphs["DeXT_arrivaltime"]->Fill(delayedPulses[0].time);
                    tDeXT = delayedPulses[0].time;
                    timeDeXT->Fill();
    
                    graphs["Npeaks"]->Fill(noise_peaks_cnt);
                    graphs["NpeaksDel"]->Fill(noise_peaks_cnt-1);
                }
                
                else if (after_pulse > 0) //  Only after pulse
                {
                    after_pulse_cnt++;
                    delayed_pulse_when_primary++;
                    sprintf(Category,"AP");
    
                    TString graph_title = Form("After pulse #DeltaV = %2.2f V",dV);
                    if(nsaved_AP < maxNwaveforms) 
                    {
                        drawWave(waveform, &color["AP"], graph_title, canv["AP"], 1.5*pe);
                        nsaved_AP++;
                    }
                    fillPersistence(persistence_AP, waveform);
                    
                    graphs["AP_arrivaltime"]->Fill(delayedPulses[0].time);
                    tAP = delayedPulses[0].time;
                    timeAP->Fill();
                    Expfit_AP->SetPoint(Expfit_AP->GetN(),delayedPulses[0].time,delayedPulses[0].volt);
    
                    graphs["Npeaks"]->Fill(noise_peaks_cnt);
                    graphs["NpeaksDel"]->Fill(noise_peaks_cnt-1);
                }
            } // end of only one primary correlated noise
            
            //------------------ Event with one primary correlated noise plus secondary correlated noise(s) ------------------------
            else if( direct_xtalk_pulse + xtalk_pulse + after_pulse > 1 ) 
            {
                counter_notclean++;
                // In this case, we look for a primary correlated noise
                tot_primary_noise_peaks_cnt++;
                // The additional correlated noise are considered as secondary
                
                // Priority is given to DiXT because we are certain that it is a primary correlated noise
                if (direct_xtalk_pulse > 0)
                {
                    // Primary correlated noise is DiXT
                    if(direct_xtalk_pulse > 1) cout << "Attention: more then one direct DiXT found" << endl;	// secondary correlated DiXT ?
                    direct_xtalk_pulse_cnt++;
                    sprintf(Category,"DiXT");
    
                    TString graph_title = Form("Direct CrossTalk #DeltaV = %2.2f V",dV);
                    if(nsaved_DiXT < maxNwaveforms) 
                    {
                        drawWave(waveform, &color["DiXT"], graph_title, canv["DiXT"], 1.5*pe);
                        nsaved_DiXT++;
                    }
                    fillPersistence(persistence_DiXT,waveform);
                    
                    graphs["Npeaks"]->Fill(noise_peaks_cnt);
                    graphs["NpeaksDiXT"]->Fill(noise_peaks_cnt);
                    
                    // The other correlated noise pulses are considered as secondary
                    secondary_pulse_cnt += xtalk_pulse + after_pulse;
                    secondary_pulse = xtalk_pulse + after_pulse;
                    // Reset DeXT and AP (they are now secondary)
                    xtalk_pulse = 0;
                    after_pulse = 0;
                } 
                else 
                {
                    // looks for the first delayed pulse and consider it as primary
                    // the later coming delayed pulse(s) is secondary
                    // we don't consider DiXT on secondary pulses (approximation)
                    std::sort(delayedPulses.begin(), delayedPulses.end());
                    
                    if( delayedPulses[0].type == "DeXT" ) 
                    {
                        xtalk_pulse_cnt++;
                        delayed_pulse_when_secondary++;
                        sprintf(Category,"DeXT");
                        
                        TString graph_title = Form("Delayed cross-talk #DeltaV = %2.2f V",dV);
                        if(nsaved_DeXT < maxNwaveforms) {
                            drawWave(waveform, &color["DeXT"], graph_title, canv["DeXT"], 1.5*pe);
                            nsaved_DeXT++;
                        }
                        fillPersistence(persistence_DeXT,waveform);
                        
                        graphs["DeXT_arrivaltime"]->Fill(delayedPulses[0].time);
                        tDeXT = delayedPulses[0].time;
                        timeDeXT->Fill();
        
                        graphs["Npeaks"]->Fill(noise_peaks_cnt);
                        graphs["NpeaksDel"]->Fill(noise_peaks_cnt-1);
                        
                        // The other correlated noise pulses are considered as secondary
                        secondary_pulse_cnt += (delayedPulses.size()-1);
                        secondary_pulse = delayedPulses.size()-1;
                        // Reset AP (they are now secondary)
                        xtalk_pulse = 1;
                        after_pulse = 0;
                    }
                    
                    else if( delayedPulses[0].type == "AP" ) {
                        after_pulse_cnt++;
                        delayed_pulse_when_secondary++;
                        sprintf(Category,"AP");
        
                        TString graph_title = Form("After pulse #DeltaV = %2.2f V",dV);
                        // don't draw the waveform if there are other secondaries
                        //~ if(nsaved_AP < maxNwaveforms) {
                            //~ drawWave(waveform, &color["AP"], graph_title, canv["AP"], 1.5*pe);
                            //~ nsaved_AP++;
                        //~ }
                        fillPersistence(persistence_AP, waveform);
                        
                        graphs["AP_arrivaltime"]->Fill(delayedPulses[0].time);
                        tAP = delayedPulses[0].time;
                        timeAP->Fill();
                        // don't use it for the fit, we use only "clean" AP events without secondary
                        //~ Expfit_AP->SetPoint(Expfit_AP->GetN(),delayedPulses[0].time,delayedPulses[0].volt);
        
                        graphs["Npeaks"]->Fill(noise_peaks_cnt);
                        graphs["NpeaksDel"]->Fill(noise_peaks_cnt-1);
                        
                        // The other correlated noise pulses are considered as secondary
                        secondary_pulse_cnt += (delayedPulses.size()-1);
                        secondary_pulse = delayedPulses.size()-1;
                        // Reset DeXT (they are now secondary)
                        xtalk_pulse = 0;
                        after_pulse = 1;
                    }
                    
                    ////////////////////////////////////////////////////
                    
                }
            } // end of more than one primary correlated noise
    
            tot_noise_peaks_cnt += noise_peaks_cnt;
            graphs["Charge"]->Fill(pulse_integral);
            if(globalArgs.save_all) otree->Fill();
            //~ otree->Fill();
            
            //~ cout << Form("single event:   #NoisePeaks = %i , #DiXT = %i , #DeXT = %i , #AP = %i , #2ndOrder = %i ",
                //~ (int)noise_peaks_cnt,
                //~ (int)direct_xtalk_pulse,
                //~ (int)xtalk_pulse,
                //~ (int)after_pulse,
                //~ (int)secondary_pulse) << endl;
            //~ cout << Form("counters:       #NoisePeaks = %i , #DiXT = %i , #DeXT = %i , #AP = %i , #2ndOrder = %i ",
                //~ (int)tot_noise_peaks_cnt,
                //~ (int)direct_xtalk_pulse_cnt,
                //~ (int)xtalk_pulse_cnt,
                //~ (int)after_pulse_cnt,
                //~ (int)secondary_pulse_cnt) << endl;
            
            if(!(j%1000)) cout << j << " / " << data_size << endl;
            delete time;
            delete volts;
        }
        
        // Estimation of DCR influence on the measurement:
        // We look at the time distriubtions of AP and DeXT
        double contribution_DCR_to_AP(0), contribution_DCR_to_DeXT(0);
        TF1 * fitAPTimeDist = NULL;
        TF1 * fitDeXTTimeDist = NULL;
        timeFitResult roofitAPTimeDist;
        timeFitResult roofitDeXTTimeDist;
        // Correction for DCR influence on delayed noise pulses
        // replace 0 to 1 in order to enable DCR correction
        if(0 && graphs["AP_arrivaltime"]->GetEntries()) 
        {
            cout << "\n\n-----> AP time distribution fit ***" << endl;
            //~ fitAPTimeDist   = fitTimeDist(graphs["AP_arrivaltime"], timeDistAP, cthrs.AP_minT*ns, hmax);
            //~ contribution_DCR_to_AP = ((hmax-cthrs.AP_minT*ns)/graphs["AP_arrivaltime"]->GetBinWidth(1))*fitAPTimeDist->GetParameter(2);
            roofitAPTimeDist = roofitTimeDist(graphs["AP_arrivaltime"], timeAP, timeDistAP, cthrs.AP_minT*ns, hmax);
            contribution_DCR_to_AP = roofitAPTimeDist.Nbkg;
            cout << Form("-----> DCR contribution to AP: %2.1f pulses (%2.1f%s)",contribution_DCR_to_AP, 100*contribution_DCR_to_AP/timeAP->GetEntries(), "%") << endl;
        }
        if(0 && graphs["DeXT_arrivaltime"]->GetEntries()) 
        {
            cout << "\n\n-----> DeXT time distribution fit ***" << endl;
            //~ fitDeXTTimeDist = fitTimeDist(graphs["DeXT_arrivaltime"], timeDistDeXT, cthrs.dir_xtalk_maxT*ns, hmax);
            //~ contribution_DCR_to_DeXT = ((hmax-cthrs.dir_xtalk_maxT*ns)/graphs["DeXT_arrivaltime"]->GetBinWidth(1))*fitDeXTTimeDist->GetParameter(2);
            roofitDeXTTimeDist = roofitTimeDist(graphs["DeXT_arrivaltime"], timeDeXT, timeDistDeXT, cthrs.dir_xtalk_maxT*ns, hmax);
            contribution_DCR_to_DeXT = roofitDeXTTimeDist.Nbkg;
            cout << Form("-----> DCR contribution to DeXT: %2.1f pulses (%2.1f%s)",contribution_DCR_to_DeXT, 100*contribution_DCR_to_DeXT/timeDeXT->GetEntries(), "%") << endl;
        }
        
        double tot_DCR_contribution = contribution_DCR_to_DeXT + contribution_DCR_to_AP;
        double fraction_delayed_primary = delayed_pulse_when_primary/(delayed_pulse_when_primary + delayed_pulse_when_secondary);	// fraction of delayed primary noise vs all delayed pulses
        
        double tot_npeaks = tot_noise_peaks_cnt + events_cnt;
        //double perc_noise_peaks = (tot_primary_noise_peaks_cnt - fraction_delayed_primary*tot_DCR_contribution) / (float)events_cnt;	// number of delayed pulses (DeXT+AP) additional to the DCR on which we trigger
        
        // 1st order correlated noise ----------------------------------
        double perc_DiXT  = direct_xtalk_pulse_cnt / (float)events_cnt;		// DiXT is the number of direct pulses higher than thrs divided by the number of triggers
        // Remove DCR contribution to DeXT
        double perc_DeXT = (xtalk_pulse_cnt - contribution_DCR_to_DeXT) / (float)events_cnt;	// first order DeXT
        // Remove DCR contribution to AP
        double perc_AP  = (after_pulse_cnt - contribution_DCR_to_AP) / (float)events_cnt;	// first order AP
        // Total first order correlated noise
        double perc_noise_peaks = perc_DiXT + perc_DeXT + perc_AP;
        
        // 2nd order correlated noise ----------------------------------
        //double perc_Sec = (tot_noise_peaks_cnt - after_pulse_cnt - xtalk_pulse_cnt) / (float)tot_npeaks;
        //double perc_Sec = (secondary_pulse_cnt - (1.0-fraction_delayed_primary)*tot_DCR_contribution) / (float)events_cnt;	// DCR is removed
        double perc_Sec = secondary_pulse_cnt / (float)events_cnt;	// no subtraction of DCR because it is removed already in the 1st order correlated noise
        
        // Contribution from DCR ---------------------------------------
        double perc_DCR = (tot_DCR_contribution / (float)events_cnt);
        
        // total contribution of correlated noise in a threshold-based frequency measurement (used for PDE)
        // it must be calculated using thresholds similar to the frequency measurement: single threshold level (f.ex. thrs=0.6)
        double corr_frequency = (tot_noise_peaks_cnt - tot_DCR_contribution)/(float)(tot_npeaks - tot_DCR_contribution);
        //double corr_frequency = (tot_noise_peaks_cnt - tot_DCR_contribution)/(float)events_cnt;
        // it is the probability that a pulse chosen randomly is actually a delayed correlated pulse (the DCR is subtracted)
        
        // total contribution of correlated noise in a gain measurement (used for PDE)
        // it must be calculated using thresholds similar to the gain measurement: single threshold level (f.ex. thrs=0.6)
        double corr_current = corr_frequency + perc_DiXT;	// same as frequency measurement (delayed pulses) add contribution from DiXT
        //~ double corr_current = corr_noise_tot_amp/(corr_noise_tot_amp + events_cnt*pe);
        
        // Compute the correction for current with the charge of pulses
        // Get average charge of 1 pe clean event
        double clean_pulse_integral = fitSimple1D(graphs["Charge"]);
        //~ double corr_current = 1-(clean_pulse_integral/graphs["Charge"]->GetMean()); // This does not subtract DCR !!!
    
        cout << "\nTotal number of events: " << events_cnt << endl;
        cout << "PE: " << pe << "\n" << endl;
    
        cout << Form("Not clean events: %i [%.2f%%]",counter_notclean, (float)counter_notclean/events_cnt*100) << endl;
        cout << Form("     (DirXtalk = %i [%.2f%%], DelXtalk = %i [%.2f%%], AP = %i [%.2f%%], SecondOrder = %i [%.2f%%])",
            (int)direct_xtalk_pulse_cnt, perc_DiXT*100.,
            (int)xtalk_pulse_cnt, perc_DeXT*100,
            (int)after_pulse_cnt, perc_AP*100,
            (int)secondary_pulse_cnt, perc_Sec*100 ) << endl;
        cout << Form("Total probability of 1st order correlated noise: %.2f%%",perc_noise_peaks*100.) << endl;
        cout << Form("Total probability of 2nd order correlated noise: %.2f%%",perc_Sec*100.) << endl;
        //cout << Form("Mean peak charge: %.4f pe", graphs["Charge"]->GetMean()) << endl;
        //cout << Form("Charge correction: %.4f pe", graphs["Charge"]->Integral() / (pe*tot_npeaks) ) << endl;
        //cout << Form("Sum of amplitude of all correlated noise pulses: %.2fV", corr_noise_tot_amp) << endl;
        cout << Form("Average charge contained in clean event: %.2f mV.ns", clean_pulse_integral*1e3*1e9) << endl;
        cout << Form("Average charge contained in each event : %.2f mV.ns", graphs["Charge"]->GetMean()*1e3*1e9) << endl;
        cout << Form("Frequency correction: %.2f%% (calculated from tot_noise_peaks_cnt)", corr_frequency*100) << endl;
        cout << Form("Current correction  : %.2f%%", corr_current*100) << endl;
        
        // Final persistence plots
        double amp0, tau;
        cout << "\n\n-----> Long tau fit ***" << endl;
        TF1 * exp_longtau = drawPersistenceWithLongTauFit(persistence_clean,canv_persistence["clean"], &amp0, &tau, pe, vol,Form("Clean waveforms #DeltaV = %2.2f V",dV));
        drawPersistence(persistence_DiXT,canv_persistence["DiXT"],Form("Direct cross-talk #DeltaV = %2.2f V",dV));
        drawPersistence(persistence_DeXT,canv_persistence["DeXT"],Form("Delayed cross-talk #DeltaV = %2.2f V",dV));
        drawPersistence(persistence_AP,canv_persistence["AP"],Form("After-pulse #DeltaV = %2.2f V",dV));
        
        //~ // Fits for long tau and recovery time
        //~ cout << "\n\n-----> Long tau fit ***" << endl;
        //~ double amp0, tau;
        //~ //TF1 * exp_longtau = fitLongTau(cleanforfit, &amp0, &tau, pe, vol, canv["clean"]);
        //~ TF1 * exp_longtau = fitLongTau(forfit, &amp0, &tau, pe, vol, canv["clean"], cleanforfit);	// does not work properly
        
        cout << "\n\n-----> After pulse fit ***" << endl;
        TF1 * exp_AP = fitAPTau(Expfit_AP, amp0, tau, pe, vol, canv["AP"]);
        
        cout << "Building final results" << endl;
    
        results["#DiXT"]->SetPoint(i,dV,perc_DiXT*100.);
        double DiXT_error = TMath::Sqrt(perc_DiXT*(1.-perc_DiXT)/events_cnt)*100.;
        results["#DiXT"]->SetPointError(i,0.,DiXT_error);
        
        // Errors don't include the error on the DCR estimation yet !!!
        results["#AP"]->SetPoint(i,dV,perc_AP*100.);
        results["#AP"]->SetPointError(i,0.,TMath::Sqrt(perc_AP*(1.-perc_AP)/events_cnt)*100.);
        
        results["#DeXT"]->SetPoint(i,dV,perc_DeXT*100.);
        results["#DeXT"]->SetPointError(i,0.,TMath::Sqrt(perc_DeXT*(1.-perc_DeXT)/events_cnt)*100.);
        
        results["Charge"]->SetPoint(i,dV,graphs["Charge"]->GetMean());
        results["Charge"]->SetPointError(i,0.,graphs["Charge"]->GetMeanError());
        
        results["#CorrectionFre"]->SetPoint(i,dV,corr_frequency*100.);
        double CorrectionFre_error = 100*TMath::Sqrt(corr_frequency*(1-corr_frequency)/(float)(tot_npeaks - tot_DCR_contribution));
        results["#CorrectionFre"]->SetPointError(i,0,CorrectionFre_error);
        
        results["#CorrectionCur"]->SetPoint(i,dV,corr_current*100.);
        double CorrectionCur_error = CorrectionFre_error + DiXT_error;  // sum of error of DiXT and FreqCorr
        results["#CorrectionCur"]->SetPointError(i,0,CorrectionCur_error);
        
        results["#Total"]->SetPoint(i,dV,perc_noise_peaks*100.);
        results["#Total"]->SetPointError(i,0.,TMath::Sqrt(perc_noise_peaks*(1.-perc_noise_peaks)/(tot_noise_peaks_cnt + events_cnt))*100.);
        
        results["#DCR"]->SetPoint(i,dV,perc_DCR*100.);
        results["#DCR"]->SetPointError(i,0,TMath::Sqrt(perc_DCR*(1.-perc_DCR)/(events_cnt))*100.);
        
        results["#SecPeaks"]->SetPoint(i,dV,perc_Sec*100.);
        results["#SecPeaks"]->SetPointError(i,0.,TMath::Sqrt(perc_Sec*(1.-perc_Sec)/events_cnt)*100.);
        
        results["#SecPeaksDiXT"]->SetPoint(i,dV,graphs["NpeaksDiXT"]->GetMean());
        results["#SecPeaksDiXT"]->SetPointError(i,0.,graphs["NpeaksDiXT"]->GetMeanError());
        
        results["#SecPeaksDel"]->SetPoint(i,dV,graphs["NpeaksDel"]->GetMean());
        results["#SecPeaksDel"]->SetPointError(i,0.,graphs["NpeaksDel"]->GetMeanError());
        
        cout << "Saving objects" << endl;
    
        // Save/print results:
        vector<TObject*> objects_to_save;
        // Waveform canvas        
        for(auto const &e : canv) objects_to_save.push_back(e.second);
        for(auto const &e : graphs) objects_to_save.push_back(e.second);
        objects_to_save.push_back(timeDistAP);
        objects_to_save.push_back(timeDistDeXT);
        for(auto const &e : canv_persistence) objects_to_save.push_back(e.second);
        // fit for long tau and AP
        objects_to_save.push_back(exp_longtau);
        //objects_to_save.push_back(exp_longtau2);
        objects_to_save.push_back(exp_AP);
    
        TString dirname = "pulse_shape_" + TString(vol);
        for(auto obj : objects_to_save) 
        {
            hfile->cd();
            hfile->cd(dirname);
            obj->Write();
            hfile->cd();
            if(globalArgs.save_all)
                if(!(obj->InheritsFrom(TF1::Class()) || obj->InheritsFrom(TF2::Class())))
                    obj->SaveAs(globalArgs.res_folder+Form("%s.pdf",obj->GetName()));
        }
    
        cout << "Done" << endl;
    
        delete Expfit_AP;
        delete cleanforfit;
        for(auto const &e : graphs) delete e.second;
        for(auto const &e : canv) delete e.second;
        for(auto const &e : canv_persistence) delete e.second;
        delete persistence_clean;
        delete persistence_DiXT;
        delete persistence_DeXT;
        delete persistence_AP;
    }
    
    // Save TTree with hist of noise for each event OV and the noise classification
    hfile->cd();
    otree->Write();
    
    // Create final plot of total primary correlated noise
    TCanvas * cfinal = new TCanvas("Correlated noise","Correlated noise",100,100,900,700);
    
    Double_t tot_max_noise = TMath::MaxElement(results["#Total"]->GetN(),results["#Total"]->GetY());
    
    formatGr(results["#Total"], kBlack, 3005, "#DeltaV [V]", "Correlated noise [%]", "Correlated Noise");
    results["#Total"]->GetYaxis()->SetRangeUser(0,tot_max_noise+2);
    formatGr(results["#DiXT"], kBlue, 3005, "#DeltaV [V]", "Correlated noise [%]");
    formatGr(results["#AP"], kOrange+7, 3005, "#DeltaV [V]", "Correlated noise [%]");
    formatGr(results["#DeXT"], kGreen+2, 3005, "#DeltaV [V]", "Correlated noise [%]");
    formatGr(results["#SecPeaks"], 7, 3005, "#DeltaV [V]", "Correlated noise [%]");
    
    results["#Total"]->Draw("ALP*3");
    results["#DiXT"]->Draw("LP*3");
    results["#AP"]->Draw("LP*3");
    results["#DeXT"]->Draw("LP*3");
    
    TLegend * leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#Total"],"Total primary noise","f");
    leg->AddEntry(results["#DiXT"],"Direct cross-talk","f");
    leg->AddEntry(results["#AP"],"After-pulse","f");
    leg->AddEntry(results["#DeXT"],"Delayed cross-talk","f");
    
    leg->SetFillColor(kWhite);
    leg->Draw();
    
    cfinal->SetGrid();
    cfinal->Print(globalArgs.res_folder+"CorrelatedNoise.pdf");
    cfinal->Write();
    
    // Create final plot for primary vs secondary correlated noise and DCR contribution
    Double_t tot_max_order = TMath::MaxElement(results["#SecPeaks"]->GetN(),results["#SecPeaks"]->GetY());
    if(tot_max_order<tot_max_noise) tot_max_order = tot_max_noise;
    
    results["#SecPeaks"]->GetYaxis()->SetRangeUser(0,tot_max_order+2);
    formatGr(results["#DCR"], 2, 3005, "#DeltaV [V]", "DCR contribution");
    
    results["#SecPeaks"]->Draw("ALP*3");
    results["#Total"]->Draw("LP*3");
    results["#DCR"]->Draw("LP*3");
    
    leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#Total"],"Primary","f");
    leg->AddEntry(results["#SecPeaks"],"Secondary","f");
    leg->AddEntry(results["#DCR"],"Contribution of DCR","f");
    
    leg->SetFillColor(kWhite);
    leg->Draw();
    
    cfinal->SetName("Primary and secondary correlated noise");
    cfinal->SetTitle("Primary and secondary correlated noise");
    cfinal->Print(globalArgs.res_folder+"PrimaryVsSecondaryCorrelatedNoise.pdf");
    cfinal->Write();
    
    formatGr(results["#CorrectionFre"], kRed, 3005, "#DeltaV [V]", "Correction", "Correction");
    formatGr(results["#CorrectionCur"], kBlue,3005, "#DeltaV [V]", "Correction", "Correction");
    results["#CorrectionCur"]->Draw("ALP*3");
    results["#CorrectionFre"]->Draw("LP*3");
    leg = new TLegend(0.15,0.75,0.37,0.87);
    leg->AddEntry(results["#CorrectionCur"],"Current","f");
    leg->AddEntry(results["#CorrectionFre"],"Frequency","f");
    leg->SetFillColor(kWhite);
    leg->Draw();
    
    cfinal->SetName("Correction");
    cfinal->SetTitle("Correction");
    cfinal->Print(globalArgs.res_folder+"Correction.pdf");
    cfinal->Write();
    
    // Print corrections and DiXT in a text file for easy copy and paste into PDE config file
    // P_all_freq
    printValues(results["#CorrectionFre"]->GetY(), results["#CorrectionFre"]->GetN(), "P_all_Freq", values_for_pde);
    printValues(results["#CorrectionFre"]->GetEY(), results["#CorrectionFre"]->GetN(), "P_all_Freq_error", values_for_pde);
    // P_all_current
    printValues(results["#CorrectionCur"]->GetY(), results["#CorrectionCur"]->GetN(), "P_all_Current", values_for_pde);
    printValues(results["#CorrectionCur"]->GetEY(), results["#CorrectionCur"]->GetN(), "P_all_Current_error", values_for_pde);
    // DiXT
    printValues(results["#DiXT"]->GetY(), results["#DiXT"]->GetN(), "DiXT_Corr", values_for_pde);
    printValues(results["#DiXT"]->GetEY(), results["#DiXT"]->GetN(), "DiXT_Corr_error", values_for_pde);
    
    // ----------
    Double_t tot_max_peaks = TMath::MaxElement(results["#SecPeaksDiXT"]->GetN(),results["#SecPeaksDiXT"]->GetY());
    results["#SecPeaksDiXT"]->SetTitle("Secondary peaks after a noise peak");
    formatGr(results["#SecPeaksDiXT"], kBlue, 0, "#DeltaV [V]", "<N Secondary peaks>");
    results["#SecPeaksDiXT"]->GetYaxis()->SetRangeUser(0,tot_max_peaks*2.);
    formatGr(results["#SecPeaksDel"], kGreen+2, 0, "#DeltaV [V]", "<N Secondary peaks>");
    results["#SecPeaksDiXT"]->Draw("ALP*");
    results["#SecPeaksDel"]->Draw("LP*");
    
    leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#SecPeaksDiXT"],"After Direct Cross-Talk","lp");
    leg->AddEntry(results["#SecPeaksDel"],"After any delayed noise","lp");
    leg->SetFillColor(kWhite);
    leg->Draw();
    
    cfinal->Print(globalArgs.res_folder+"SecondaryPeaks.pdf");
    cfinal->SetName("SecondaryPeaks");
    cfinal->SetTitle("SecondaryPeaks");
    cfinal->Write();
    
    formatGr(results["Charge"], kGreen+2, 0, "#DeltaV [V]", "<Charge per counted peak>");
    results["Charge"]->Draw("AP");
    cfinal->Print(globalArgs.res_folder+"Charge.pdf");
    cfinal->SetName("Charge");
    cfinal->SetTitle("Charge");
    cfinal->Write();    
    
    // ------------
    hfile->Close();
    values_for_pde->close();
    
    return 0;
}


