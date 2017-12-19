//
//  Created by Felipe Gilberto Ortega on 25/04/16.
//  Modified by Luca Pescatore, luca.pescatore@cern.ch on 28/06/2017
//  Modified by Olivier Girard, olivier.girard@cern.ch on 07-08/2017
//

#include <TSpectrum.h>
#include <TPaveStats.h>
#include "Noise_analysis_largeEvent.h"
#include "Thresholds.h"
#include "fits.h"
#include "lhcbStyle.h"

using namespace std;

int main(int argc, char* argv[]) 
{
    lhcbstyle();

    //gROOT->ProcessLine(".x Analysis/lhcbStyle.C");

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
        {"#DiXT",                 new TGraphErrors()},
        {"#DeXT",                 new TGraphErrors()},
        {"#AP",                   new TGraphErrors()},
        {"#Total",                new TGraphErrors()},
        {"#DCR",                  new TGraphErrors()},
        {"#SecPeaks",             new TGraphErrors()},
        {"#SecPeaksDiXT",         new TGraphErrors()},
        {"#SecPeaksDel",          new TGraphErrors()},
        {"1PE-Waveform-Charge",   new TGraphErrors()},
        {"1PE-Charge-Fit",        new TGraphErrors()},
        {"Mean-Charge",           new TGraphErrors()},
        {"Mean-Charge-Corrected", new TGraphErrors()},
        {"#Double",               new TGraphErrors()},
        {"#CorrectionFre",        new TGraphErrors()},
        {"#CorrectionCur",        new TGraphErrors()} };
    
    Char_t Category[15];
    Double_t dV;
    
    ofstream * values_for_pde = new ofstream(globalArgs.res_folder + "values_for_pde.txt",ofstream::out);
    
    TFile * hfile = TFile::Open(TString(globalArgs.res_folder) + "noiseanalysis.root","RECREATE");
    initOutputFile(hfile, vol_folders);
    
    TTree * otree = 0;
    otree = new TTree("ClassifiedData","ClassifiedData");
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
    
    
    vector<double> Xmin, Xmax, Ymin, Ymax, NsamplesPerEvMax;
    Xmin.clear(); Xmax.clear(); Ymin.clear(); Ymax.clear(); NsamplesPerEvMax.clear();
    
    std::cout << " -----> Calculation of average DCR amplitudes " << std::endl;
    for(unsigned int v(0); v<vol_size; ++v)
    {
        vector<double> minmax(5,0);
        const char * vol = vol_folders[v];
        pe_volt.push_back(Amplitude_calc(vol, data_size, minmax, "root", hfile));
        //~ pe_volt.push_back(Amplitude_calc(vol, data_size, "root", hfile));
        Vbias_vs_pe->SetPoint(v, pe_volt.back(), TString(vol).Atof());
        std::cout << Form("Voltage: %.1fV --> Mean amplitude = %.6f", TString(vol).Atof(), pe_volt.back()) << std::endl;
        Xmin.push_back(minmax[0]);
        Xmax.push_back(minmax[1]);
        Ymin.push_back(minmax[2]);
        Ymax.push_back(minmax[3]);
        NsamplesPerEvMax.push_back(minmax[4]);
    }
    
    cout << "\n\n-----> Voltage Breakdown fit" << endl;
    
    double VBD = fitBreakdownVoltage(Vbias_vs_pe, hfile, values_for_pde);
    
    if(globalArgs.only_Vbd) {hfile->Close(); values_for_pde->close(); return 0;}	// only Vbd measurement, abort other measurements
    
    
    /////////////////
    // Loop over all Voltages measured
    /////////////////
    
    cout << "\n\n-----> Noise analysis *** " << endl;
    
    TGraph * waveform = NULL;
    TH1 * hist = NULL;                  // TSpectrum method
    TSpectrum * s = new TSpectrum();    // TSpectrum method
    hfile->cd();
        
    for (int i = 0; i < vol_size; i++)
    {    
        const char * vol = vol_folders[i];
        cout << "\n\n----> Voltage analyzed: " << vol << endl;
    
        map<string, int> color{{"clean",1},{"AP",1},{"DiXT",1},{"DeXT",1},{"Sec",1}};
    
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
        unsigned int nsaved(0), nsaved_DiXT(0), nsaved_DeXT(0), nsaved_AP(0), nsaved_Sec(0);	// counters on the number of waveforms in each canvas
        unsigned int maxNwaveforms = 100; // maximum number of waveforms in the canvas
        vector<delayedPulse> delayedPulses;	// used for DeXT and AP
        double corr_noise_tot_amp(0);
        
        // Define amplitude measured at which OV
        
        Vbias = vol_folders[i].Atof();
        pe = pe_volt[i];
        dV = Vbias - VBD;
        
        // single waveforms multigraphs
        map<string,TCanvas *> canv {
            {"DiXT", new TCanvas(Form("DiXT_%s_%dwaveforms",vol,maxNwaveforms), Form("Direct cross-talk #DeltaV = %2.2fV",dV))},
            {"DeXT", new TCanvas(Form("DeXT_%s_%dwaveforms",vol,maxNwaveforms), Form("Delayed cross-talk #DeltaV = %2.2fV",dV))},
            {"AP",   new TCanvas(Form("AP_%s_%dwaveforms",vol,maxNwaveforms), Form("After-pulse #DeltaV = %2.2fV",dV))},
            {"clean",new TCanvas(Form("Clean_%s_%dwaveforms",vol,maxNwaveforms), Form("Clean #DeltaV = %2.2fV",dV))},
            {"Sec",  new TCanvas(Form("Secondary_%s_%dwaveforms",vol,maxNwaveforms), Form("Waveforms with secondaries #DeltaV = %2.2fV",dV))},
            {"APtime",new TCanvas(Form("AP_%s",vol), Form("After-pulse #DeltaV = %2.2fV",dV))}
        };
        
        //TMultiGraph * forfit    = new TMultiGraph("multigraph_for_fit","multigraph_for_fit");  
        TGraph * Expfit_AP      = new TGraph();
        TGraphErrors * cleanforfit    = NULL;
        map<string,TH1D*> graphs {
           {"AP_arrivaltime",   new TH1D(Form("AP_arrival_times_%s",vol),"After-pulse arrival times", 120, 0, 0.18e-6)},
           {"DeXT_arrivaltime", new TH1D(Form("DeXT_arrival_times_%s",vol),"Delayed cross-talk arrival times", 120, 0, 0.18e-6)},
           {"Npeaks",           new TH1D(Form("N_peaks_%s",vol),"Number of noise peaks when not clean", 50, 0, 50)},
           {"NpeaksDiXT",       new TH1D(Form("N_peaks_DiXT_%s",vol),"Number of noise peaks when direct DiXT", 50, 0, 50)},
           {"NpeaksDel",        new TH1D(Form("N_peaks_Delayed_%s",vol),"Number of noise peaks when delayed noise", 50, 0, 50)},
           {"Charge",           new TH1D(Form("Charge_%s",vol),"Charge of waveforms", 1000, 0., ns*100*pe)}
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
            {"DiXT", new TCanvas(Form("persistence_DiXT_%s",vol), Form("Persistence direct cross-talk #DeltaV = %2.2fV",dV))},
            {"DeXT", new TCanvas(Form("persistence_DeXT_%s",vol), Form("Persistence delayed cross-talk #DeltaV = %2.2fV",dV))},
            {"AP",   new TCanvas(Form("persistence_AP_%s",vol), Form("Persistence after-pulse #DeltaV = %2.2fV",dV))},
            {"clean",new TCanvas(Form("persistence_Clean_%s",vol), Form("Persistence Clean #DeltaV = %2.2fV",dV))}
        };
        
        double bin_size = 0.002;	// vertical bin size for persistence plot
        TH2D * persistence_clean = NULL;
        TH2D * persistence_DiXT  = NULL;
        TH2D * persistence_DeXT  = NULL;
        TH2D * persistence_AP    = NULL;
        
        // Amplitude vs arrival time graph (DiXT, DeXT, AP, Secondaries and not classified)
        map<string,TGraph *> gAmp_vs_Time { {"1:DiXT", new TGraph()}, {"3:DeXT", new TGraph()},
			{"2:AP",   new TGraph()}, {"4:Sec",  new TGraph()} };
        
        // Setup input tree
        TTree * tree = NULL;
        if(globalArgs.input=="root") 
        {   
            tree = (TTree*)ifile->Get(vol);
            tree->SetBranchAddress("NsampPerEv",&npts);
            tree->SetBranchAddress("Amps",&amps);
            tree->SetBranchAddress("Times",&times);
            //data_size = tree->GetEntries();
            if(data_size>tree->GetEntries()) data_size = tree->GetEntries();
            // data_size set from cfg file
            
            int Nv = int(((Ymax[i]-Ymin[i])/bin_size)+0.5);
            int NvPE = int(((1.5*pe-Ymin[i])/bin_size)+0.5);
            persistence_clean = new TH2D("Clean_persistence","Clean waveforms",   NsamplesPerEvMax[i], Xmin[i], Xmax[i], NvPE, Ymin[i]/pe,1.5);
            persistence_DiXT  = new TH2D("DiXT_persistence","Direct cross-talk",  NsamplesPerEvMax[i], Xmin[i], Xmax[i], Nv,   Ymin[i]/pe,Ymax[i]/pe);
            persistence_DeXT  = new TH2D("DeXT_persistence","Delayed cross-talk", NsamplesPerEvMax[i], Xmin[i], Xmax[i], NvPE, Ymin[i]/pe,1.5);
            persistence_AP    = new TH2D("AP_persistence","After-pulse",          NsamplesPerEvMax[i], Xmin[i], Xmax[i], NvPE, Ymin[i]/pe,1.5);
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
            
            bool done = false;
            bool veryclean = true;
            double lastPulseT = -5*ns;	// used to apply an analysis dead time after each pulse (to avoid noise and ringing)
                                        // this is used to ensure a certain hysteresis (dead time also works)
                                        // set to neg value such that it works also for the DiXT region 
            double deadTime = cthrs.del_xtalk_minT * ns;	// deadtime set to min delay time of DeXT (optimised to avoid ringing)
            
            
// ---------------------------------------------------------------------
// --------------------- TSpectrum peak detection ----------------------
            //~ // Search for peaks with TSpectrum method
            //~ // Peak detection using TSpectrum
            //~ // need a TH1
            //~ hist = createHist_and_computePulseIntegral(npts,time,volts,pulse_integral,baseline_shift);
            //~ vector<double> time_peak, amp_peak;
            //~ int npeaks = 0;
            //~ 
            //~ double max_val = hist->GetMaximum();
            //~ 
            //~ if(max_val>min_thrs) {
                //~ 
                //~ //TH1 * hist = convertGrToH(waveform);
                //~ // Search( TH1 , sigma [estimated] , "" , threshold [wrt max of TH1] )
                //~ s->Search(hist,0.01,"",min_thrs/max_val);
                //~ //TFile * outtest = new TFile("tspectrum_test.root","RECREATE");
                //~ //hist->Write();
                //~ //outtest->Close();
                //~ //double x;
                //~ //cin >> x;
                //~ 
                //~ npeaks = s->GetNPeaks();
                //~ double * xpeak = s->GetPositionX();
                //~ double * ypeak = s->GetPositionY();
                //~ for(unsigned int peak(0); peak<npeaks; ++peak) {
                    //~ int bin_num = hist->FindBin(xpeak[peak]);
                    //~ double max_y(0);
                    //~ double max_x(xpeak[peak]);
                    //~ for(int b(-4); b<=4; ++b) {    // gets the position and amplitude of the peaks more precisely
                        //~ if(hist->GetBinContent(bin_num+b) > max_y) {
                            //~ max_y = hist->GetBinContent(bin_num+b);
                            //~ max_x = hist->GetBinCenter(bin_num+b);
                        //~ }
                    //~ }
                    //~ amp_peak.push_back(max_y);
                    //~ time_peak.push_back(max_x);
                //~ }
            //~ }
            //~ delete hist;
            
// ---------------------------------------------------------------------
// ------------------------ Baseline evaluation ------------------------
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
            
// ---------------------------------------------------------------------
// ------------------- Loop over the detected peaks --------------------

            for (int row = bsl_npts; row < npts-2; row++) { // loop start when time[row] = 0ns
                double curT = time[row];
                double curV = volts[row];  
                
                // Compute pulse integral as the sum of area of trapezoids
                //if(curV-baseline_shift>0) pulse_integral += 0.5 * dt *(curV+volts[row-1] - 2*baseline_shift);
                pulse_integral += 0.5 * dt *(curV+volts[row-1] - 2*baseline_shift);
    
                // Skip if not on a maximum
                // Skip if dead time condition is not fulfilled
                if(!( ((curV > volts[row-2] && curV > volts[row+2])|(curV > volts[row-3] && curV > volts[row+3])) && curV > min_thrs && curT > lastPulseT+deadTime )) continue;
                else if (done) { done = false; continue; }
                done = true;
                
            //~ for(unsigned int peak(0); peak<npeaks; ++peak) { // loop over the peaks detected by TSpectrum
                //~ double curT = time_peak[peak];
                //~ double curV = amp_peak[peak];
                
                if ( curT > 2*ns && curV > 0.4*pe) veryclean = false;
    
                // Direct x-talk: in 0-1 ns window and V > direct th.
                if( curT <= cthrs.dir_xtalk_maxT * ns && curV > DiXT_thr) 
                {
                    direct_xtalk_pulse++;
                    //noise_peaks_cnt++; // Don't count peak as it is over the DCR
                    delayedPulse dixt;
                    //~ dixt.time = 0;
                    dixt.time = curT;
                    dixt.volt = curV;
                    dixt.type = "DiXT";
                    delayedPulses.push_back(dixt);
                    corr_noise_tot_amp += curV;
                    lastPulseT = curT;
                    veryclean = false; // added
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
                    veryclean = false; // added
                }
                // After-pulse: time larger 2ns and V > AP th
                else if ( curT > cthrs.AP_minT * ns && curV > AP_thr && curT > lastPulseT+cthrs.AP_minT*ns) 
		//else if ( curT > 4 * ns && curV > DeXT_thr && curT <20*ns) 
                {	// if it's not a DeXT, it might be an AP
                    
                    if(direct_xtalk_pulse + xtalk_pulse + after_pulse >= 1) {
                        // if it is a secondary, we check if it is all the baseline which is higher than the threshold
                        // if so, we discard this peak
                        double baseline_before(0);
                        for(int t(row-10); t<row-3; ++t) if(t>0) baseline_before += volts[t];
                        baseline_before = baseline_before/7;
                        if((curV-baseline_before)/pe < 0.1) continue;
                    }
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
                    veryclean = false; // added
                }
                else if ( curT <= lastPulseT+deadTime ) continue; // skip pulses that might be due to ringing
                
                //cout << "Pulse at t=" << curT << "   dead time until t=" << lastPulseT+deadTime << endl;
    
            } // loop over time
            
            //------------------ Clean event ---------------------------
            if( direct_xtalk_pulse + xtalk_pulse + after_pulse == 0 ) 
            {
                sprintf(Category,"Clean");
    
                // Get only very clean waves and make an average for long tau fit.
                if(veryclean) cleanforfit = average(cleanforfit, waveform, pe);
                //~ cleanforfit = average(cleanforfit, waveform);
    
                if(nsaved < maxNwaveforms) // Max 20 clean graphs on the plot
                {
                    nsaved++;
                    drawWave(waveform, &color["clean"], Form("Clean pulse #DeltaV = %2.2fV",dV), canv["clean"], 2.7, pe, baseline_shift);
                    //forfit->Add(waveform);
                }
                fillPersistence(persistence_clean, waveform, pe);
            }
            
            //------------------ Event with only one primary correlated noise ------------------------
            else if( direct_xtalk_pulse + xtalk_pulse + after_pulse == 1 ) 
            {
                counter_notclean++;
                tot_primary_noise_peaks_cnt++;
                
                addPointToGraph(gAmp_vs_Time, delayedPulses[0].type, delayedPulses[0].time, delayedPulses[0].volt/pe);
                
                // They are categroized as DiXT, DeXT or AP
                if (direct_xtalk_pulse > 0)  // Primary correlated noise is DiXT
                {
                    if(direct_xtalk_pulse > 1) cout << "Attention: more then one direct DiXT found" << endl;	// secondary correlated DiXT ?
                    direct_xtalk_pulse_cnt++;
                    sprintf(Category,"DiXT");
    
                    TString graph_title = Form("Direct CrossTalk #DeltaV = %2.2fV",dV);
                    if(nsaved_DiXT < maxNwaveforms) 
                    {
                        drawWave(waveform, &color["DiXT"], graph_title, canv["DiXT"], 2.7, pe, baseline_shift);
                        nsaved_DiXT++;
                    }
                    fillPersistence(persistence_DiXT,waveform, pe);
                    
                    graphs["Npeaks"]->Fill(noise_peaks_cnt);
                    graphs["NpeaksDiXT"]->Fill(noise_peaks_cnt);
                }
                
                else if (xtalk_pulse > 0) // Only Delayed x-talk
                {
                    xtalk_pulse_cnt++;
                    delayed_pulse_when_primary++;
                    sprintf(Category,"DeXT");
                    
                    TString graph_title = Form("Delayed cross-talk #DeltaV = %2.2fV",dV);
                    if(nsaved_DeXT < maxNwaveforms) 
                    {
                        drawWave(waveform, &color["DeXT"], graph_title, canv["DeXT"], 2.7, pe, baseline_shift);
                        nsaved_DeXT++;
                    }
                    fillPersistence(persistence_DeXT,waveform, pe);
                    
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
    
                    TString graph_title = Form("After pulse #DeltaV = %2.2fV",dV);
                    if(nsaved_AP < maxNwaveforms) 
                    {
                        drawWave(waveform, &color["AP"], graph_title, canv["AP"], 2.7, pe, baseline_shift);
                        nsaved_AP++;
                    }
                    fillPersistence(persistence_AP, waveform, pe);
                    
                    graphs["AP_arrivaltime"]->Fill(delayedPulses[0].time);
                    tAP = delayedPulses[0].time;
                    timeAP->Fill();
                    Expfit_AP->SetPoint(Expfit_AP->GetN(),delayedPulses[0].time,(delayedPulses[0].volt-baseline_shift)/pe);
    
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
                
                // Sort the delayed pulses in arrival time order
                std::sort(delayedPulses.begin(), delayedPulses.end());
                // Fills the Amp vs arrival time graphs
                addPointToGraph(gAmp_vs_Time, delayedPulses[0].type, delayedPulses[0].time, delayedPulses[0].volt/pe);
                // The latest delayed noises are secondaries
                for(unsigned int pulse(1); pulse<delayedPulses.size(); ++pulse)
                    addPointToGraph(gAmp_vs_Time, "Sec", delayedPulses[pulse].time, delayedPulses[pulse].volt/pe);
                
                // Priority is given to DiXT because we are certain that it is a primary correlated noise
                if (direct_xtalk_pulse > 0)
                {
                    // Primary correlated noise is DiXT
                    if(direct_xtalk_pulse > 1) cout << "Attention: more then one direct DiXT found" << endl;	// secondary correlated DiXT ?
                    direct_xtalk_pulse_cnt++;
                    sprintf(Category,"DiXT");
    
                    TString graph_title = Form("Waveform with secondaries #DeltaV = %2.2fV",dV);
                    if(nsaved_Sec < maxNwaveforms) 
                    {
                        drawWave(waveform, &color["Sec"], graph_title, canv["Sec"], 2.7, pe, baseline_shift);
                        nsaved_Sec++;
                    }
                    fillPersistence(persistence_DiXT,waveform, pe);
                    
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
                    
                    if( delayedPulses[0].type == "DeXT" ) 
                    {
                        xtalk_pulse_cnt++;
                        delayed_pulse_when_secondary++;
                        sprintf(Category,"DeXT");
                        
                        TString graph_title = Form("Waveform with secondaries #DeltaV = %2.2fV",dV);
                        if(nsaved_Sec < maxNwaveforms) 
                        {
                            drawWave(waveform, &color["Sec"], graph_title, canv["Sec"], 2.7, pe, baseline_shift);
                            nsaved_Sec++;
                        }
                        fillPersistence(persistence_DeXT,waveform, pe);
                        
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
        
                        TString graph_title = Form("Waveform with secondaries #DeltaV = %2.2fV",dV);
                        if(nsaved_Sec < maxNwaveforms) 
                        {
                            drawWave(waveform, &color["Sec"], graph_title, canv["Sec"], 2.7, pe, baseline_shift);
                            nsaved_Sec++;
                        }
                        fillPersistence(persistence_AP, waveform, pe);
                        
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
            if(globalArgs.save_tree) otree->Fill();
            
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
        timeFitResult APTimeDist_fit_result;
        timeFitResult DeXTTimeDist_fit_result;
        // Correction for DCR influence on delayed noise pulses
        // replace 0 to 1 in order to enable DCR correction
        if(graphs["AP_arrivaltime"]->GetEntries()) 
        {
            cout << "\n\n-----> AP time distribution fit ***" << endl;
            fitTimeDist(graphs["AP_arrivaltime"], timeDistAP, APTimeDist_fit_result, "AP", 40*ns, Xmax[i]);
            //~ roofitTimeDist(graphs["AP_arrivaltime"], timeAP, timeDistAP, APTimeDist_fit_result, "AP", cthrs.AP_minT*ns, Xmax[i]);
            //~ roofitTimeDist(graphs["AP_arrivaltime"], timeAP, timeDistAP, APTimeDist_fit_result, "AP", 50*ns, Xmax[i]);   // AP mean lifetime fit starts at 50ns must be adapted maybe
            contribution_DCR_to_AP = APTimeDist_fit_result.Nbkg;
            if(!globalArgs.enable_dcr) contribution_DCR_to_AP = 0;	// Contribution of DCR sometimes wrong estimated, set to zero for QA H2017 for automated measurements!
            cout << Form("-----> DCR contribution to AP: %2.1f pulses (%2.1f%s)",contribution_DCR_to_AP, 100*contribution_DCR_to_AP/(APTimeDist_fit_result.Nsig+APTimeDist_fit_result.Nbkg), "%") << endl;
        }
        if(graphs["DeXT_arrivaltime"]->GetEntries()) 
        {
            cout << "\n\n-----> DeXT time distribution fit ***" << endl;
            fitTimeDist(graphs["DeXT_arrivaltime"], timeDistDeXT, DeXTTimeDist_fit_result, "DeXT", cthrs.dir_xtalk_maxT*ns, Xmax[i]);
            //~ roofitTimeDist(graphs["DeXT_arrivaltime"], timeDeXT, timeDistDeXT, DeXTTimeDist_fit_result, "DeXT", cthrs.del_xtalk_minT*ns, Xmax[i]);
            contribution_DCR_to_DeXT = DeXTTimeDist_fit_result.Nbkg;
            if(!globalArgs.enable_dcr) contribution_DCR_to_DeXT = 0;	// Contribution of DCR sometimes wrong estimated, set to zero for QA H2017 for automated measurements!
            cout << Form("-----> DCR contribution to DeXT: %2.1f pulses (%2.1f%s)",contribution_DCR_to_DeXT, 100*contribution_DCR_to_DeXT/(DeXTTimeDist_fit_result.Nbkg+DeXTTimeDist_fit_result.Nsig), "%") << endl;
        }
        
        double tot_DCR_contribution = contribution_DCR_to_DeXT + contribution_DCR_to_AP;
        //double fraction_delayed_primary = delayed_pulse_when_primary/(delayed_pulse_when_primary + delayed_pulse_when_secondary);	// fraction of delayed primary noise vs all delayed pulses
        
        double tot_npeaks = tot_noise_peaks_cnt + events_cnt;
        //double perc_noise_peaks = (tot_primary_noise_peaks_cnt - fraction_delayed_primary*tot_DCR_contribution) / (float)events_cnt;	// number of delayed pulses (DeXT+AP) additional to the DCR on which we trigger
        
        // -------------------------------------------------------------
        // 1st order correlated noise ----------------------------------
        // -------------------------------------------------------------
        
        double perc_DiXT  = direct_xtalk_pulse_cnt / (float)events_cnt;		// DiXT is the number of direct pulses higher than thrs divided by the number of triggers
        //~ double DiXT_error = TMath::Sqrt(perc_DiXT*(1.-perc_DiXT)/events_cnt)*100.;  // binomial error
        double DiXT_error = absoluteErrorPoisson( perc_DiXT*100., direct_xtalk_pulse_cnt, events_cnt);  // poisson error
        
        // Remove DCR contribution to DeXT
        double perc_DeXT = (xtalk_pulse_cnt - contribution_DCR_to_DeXT) / (float)events_cnt;	// first order DeXT
        //~ double DeXT_error = TMath::Sqrt(perc_DeXT*(1.-perc_DeXT)/events_cnt)*100.; // binomial error
        double DeXT_error = absoluteErrorPoisson(perc_DeXT*100., (xtalk_pulse_cnt - contribution_DCR_to_DeXT), events_cnt);  // poisson error
        
        // Remove DCR contribution to AP
        double perc_AP  = (after_pulse_cnt - contribution_DCR_to_AP) / (float)events_cnt;	// first order AP
        //~ double AP_error = TMath::Sqrt(perc_AP*(1.-perc_AP)/events_cnt)*100.;  // binomial error
        double AP_error = absoluteErrorPoisson(perc_AP*100., (after_pulse_cnt - contribution_DCR_to_AP), events_cnt);  // poisson error
        
        // Total first order correlated noise
        double perc_noise_peaks = perc_DiXT + perc_DeXT + perc_AP;
        //~ double totPrimCorr_error = TMath::Sqrt(perc_noise_peaks*(1.-perc_noise_peaks)/(tot_noise_peaks_cnt + events_cnt))*100.;  // binomial error
        double totPrimCorr_error = absoluteErrorPoisson( perc_noise_peaks*100., (direct_xtalk_pulse_cnt + xtalk_pulse_cnt - contribution_DCR_to_DeXT + after_pulse_cnt - contribution_DCR_to_AP), events_cnt);  // poisson error
        
        // -------------------------------------------------------------
        // 2nd order correlated noise ----------------------------------
        // -------------------------------------------------------------
        
        double perc_Sec = secondary_pulse_cnt / (float)events_cnt;	// no subtraction of DCR because it is removed already in the 1st order correlated noise
        //~ double SecO_error = TMath::Sqrt(perc_Sec*(1.-perc_Sec)/events_cnt)*100.;  // binomial error
        double SecO_error = absoluteErrorPoisson(perc_Sec*100.,secondary_pulse_cnt,events_cnt);  // poisson error
        
        // -------------------------------------------------------------
        // Contribution from DCR ---------------------------------------
        // -------------------------------------------------------------
        
        double perc_DCR = (tot_DCR_contribution / (float)events_cnt);
        //~ double DCR_error = TMath::Sqrt(perc_DCR*(1.-perc_DCR)/(events_cnt))*100.;  // binomial error
        double DCR_error = absoluteErrorPoisson(perc_DCR*100.,tot_DCR_contribution,events_cnt);  // poisson error
        
        // -------------------------------------------------------------
        // Corrections for PDE measurement -----------------------------
        // -------------------------------------------------------------
        
        // total contribution of correlated noise in a threshold-based frequency measurement (used for PDE)
        // it must be calculated using thresholds similar to the frequency measurement: single threshold level (f.ex. thrs=0.6)
        double corr_frequency = (tot_noise_peaks_cnt - tot_DCR_contribution)/(float)(tot_npeaks - tot_DCR_contribution);
        // it is the probability that a pulse chosen randomly is actually a delayed correlated pulse (the DCR is subtracted)
        //~ double CorrectionFre_error = 100*TMath::Sqrt(corr_frequency*(1-corr_frequency)/(float)(tot_npeaks - tot_DCR_contribution));  // binomial error
        double CorrectionFre_error = absoluteErrorPoisson(corr_frequency*100., (tot_noise_peaks_cnt - tot_DCR_contribution), (tot_npeaks - tot_DCR_contribution));  // poisson error
        
        // total contribution of correlated noise in a gain measurement (used for PDE)
        // it must be calculated using thresholds similar to the gain measurement: single threshold level (f.ex. thrs=0.6)
        double corr_current = corr_frequency + perc_DiXT;	// same as frequency measurement (delayed pulses) add contribution from DiXT
        double CorrectionCur_error = TMath::Sqrt(CorrectionFre_error*CorrectionFre_error + DiXT_error*DiXT_error);  // quadratice sum of error of DiXT and FreqCorr
        
        // -------------------------------------------------------------
        // Waveform charge from integral -------------------------------
        // -------------------------------------------------------------
        
        // Get average charge of 1 pe clean event
        double OnePEIntegral_error(0);
        double OnePEIntegral = fitSimple1D(graphs["Charge"], OnePEIntegral_error);    // from charge distribution fit
        double clean_pulse_integral_error(0);
        double clean_pulse_integral = computeIntegral(cleanforfit, evaluateBaselineShift(cleanforfit), clean_pulse_integral_error, 1./pe);    // from clean waveforms integral
        // Mean waveform charge
        double mean_waveform_charge_error = graphs["Charge"]->GetMeanError();
        double mean_waveform_charge = graphs["Charge"]->GetMean();
        // Correct the mean charge with the current correction
        double corrected_charge = mean_waveform_charge*(1-corr_current);
        double corrected_charge_error = corrected_charge*TMath::Sqrt((mean_waveform_charge_error/mean_waveform_charge)*(mean_waveform_charge_error/mean_waveform_charge) + 1e-4*(CorrectionCur_error/(1-corr_current))*(CorrectionCur_error/(1-corr_current)));
        // Compute the correction for current with the charge of pulses
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
        cout << Form("Average charge contained in clean event      : %.2f mV.ns", clean_pulse_integral*1e3*1e9) << endl;
        cout << Form("Average charge contained in each event       : %.2f mV.ns", graphs["Charge"]->GetMean()*1e3*1e9) << endl;
        cout << Form("Average charge corrected for correlated noise: %.2f mV.ns", corrected_charge*1e3*1e9) << endl;
        cout << Form("Frequency correction: %.2f%% (calculated from tot_noise_peaks_cnt)", corr_frequency*100) << endl;
        cout << Form("Current correction  : %.2f%%", corr_current*100) << endl;
        
        // Final persistence plots
        double amp0, tau;
        cout << "\n\n-----> Long tau fit ***" << endl;
        //~ TF1 * exp_longtau = drawPersistenceWithLongTauFit(persistence_clean,canv_persistence["clean"], &amp0, &tau, pe, vol,Form("Clean waveforms #DeltaV = %2.2fV, Slow component amp = %2.4f",dV,amp0));
        TF1 * exp_longtau = drawPersistenceWithLongTauFit(persistence_clean,canv_persistence["clean"], &amp0, &tau, 1, vol,Form("Clean waveforms #DeltaV = %2.2fV, Slow component amp = %2.2fPE",dV,amp0), cleanforfit);
        drawPersistence(persistence_DiXT,canv_persistence["DiXT"],Form("Direct cross-talk #DeltaV = %2.2fV",dV));
        drawPersistence(persistence_DeXT,canv_persistence["DeXT"],Form("Delayed cross-talk #DeltaV = %2.2fV",dV));
        drawPersistence(persistence_AP,canv_persistence["AP"],Form("After-pulse #DeltaV = %2.2fV",dV));
        
        // Final amp vs Time graphs
        vector<TString> add;
        add.push_back(Form("(%2.1f%%)",perc_DiXT*100)); add.push_back(Form("(%2.1f%%)",perc_AP*100)); add.push_back(Form("(%2.1f%%)",perc_DeXT*100)); add.push_back(Form("(%2.1f%%)",perc_Sec*100)); 
        TCanvas * c_amp_vs_time = finalizeMapGraphs(gAmp_vs_Time, Form("Amp_vs_Time_graph_%2.2fV",dV), "Pulse amplitude versus arrival time", "Arrival time [s]", "Pulse amplitude [PE]", add);
        double recovery; 
        //~ // Fits for long tau and recovery time
        //~ cout << "\n\n-----> Long tau fit ***" << endl;
        //~ double amp0, tau;
        //~ //TF1 * exp_longtau = fitLongTau(cleanforfit, &amp0, &tau, pe, vol, canv["clean"]);
        //~ TF1 * exp_longtau = fitLongTau(forfit, &amp0, &tau, pe, vol, canv["clean"], cleanforfit);	// does not work properly
        
        cout << "\n\n-----> After pulse fit ***" << endl;
        TF1 * exp_AP = fitAPTau(Expfit_AP, amp0, tau, 1, vol, &recovery, cthrs.AP_minT*ns, 180*ns);
        drawAPfit(canv["AP"], exp_AP, 1);
        canv["APtime"]->cd();
        formatGr(Expfit_AP, kBlue, 0, "Arrival time [s]", "Pulse amplitude [PE]", "After-pulse amplitude vs arrival time");
        Expfit_AP->SetName("APamp_vs_arrivalTime");
        Expfit_AP->Draw("AP");
        Expfit_AP->GetYaxis()->SetRangeUser(0,1);
        TPaveText * pv = new TPaveText(0.15,0.76,0.45,0.85,"brNDC");
        pv->AddText(Form("#tau_{rec} = %2.1f#pm%2.1f ns",1e9*exp_AP->GetParameter(1),1e9*exp_AP->GetParError(1)));
        pv->SetFillColor(kWhite);
        pv->Draw("SAME");
        
        cout << "Building final results" << endl;
        
        // 1st order correlated noise
        setPoint(results["#DiXT"],  i, dV, perc_DiXT*100., 0., DiXT_error);
        setPoint(results["#AP"],    i, dV, perc_AP*100.,   0., AP_error);
        setPoint(results["#DeXT"],  i, dV, perc_DeXT*100., 0., DeXT_error);
        setPoint(results["#Total"], i, dV, perc_noise_peaks*100., 0., totPrimCorr_error);
        
        // 2nd order correlated noise
        setPoint(results["#SecPeaks"], i, dV, perc_Sec*100., 0., SecO_error);
        setPoint(results["#SecPeaksDiXT"], i, dV, graphs["NpeaksDiXT"]->GetMean(), 0., graphs["NpeaksDiXT"]->GetMeanError());
        setPoint(results["#SecPeaksDel"], i, dV, graphs["NpeaksDel"]->GetMean(), 0., graphs["NpeaksDel"]->GetMeanError());
        
        // DCR contribution
        setPoint(results["#DCR"], i, dV, perc_DCR*100., 0., DCR_error);
        
        // Corrections for PDE
        setPoint(results["#CorrectionFre"], i, dV, corr_frequency*100., 0., CorrectionFre_error);
        setPoint(results["#CorrectionCur"], i, dV, corr_current*100., 0., CorrectionCur_error);
        
        // Charge
        setPoint(results["1PE-Waveform-Charge"], i, dV, clean_pulse_integral*1e3*1e9, 0., clean_pulse_integral_error*1e3*1e9);
        setPoint(results["1PE-Charge-Fit"], i, dV, OnePEIntegral*1e3*1e9, 0., OnePEIntegral_error*1e3*1e9);
        setPoint(results["Mean-Charge"], i, dV, mean_waveform_charge*1e3*1e9, 0., mean_waveform_charge_error*1e3*1e9);
        setPoint(results["Mean-Charge-Corrected"], i, dV, corrected_charge*1e3*1e9, 0., corrected_charge_error*1e3*1e9);
        
        cout << "Saving objects" << endl;
    
        // Save/print results:
        vector<TObject*> objects_to_save;
        // Thresholds
        TCanvas * c_t = listThrs(cthrs);
        objects_to_save.push_back(c_t);
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
        //objects_to_save.push_back(Expfit_AP);
        objects_to_save.push_back(c_amp_vs_time);
        //objects_to_save.push_back(splitAPTau(Expfit_AP, amp0, tau, pe));
    
        TString dirname = "pulse_shape_" + TString(vol);
        for(auto obj : objects_to_save) 
        {
            hfile->cd();
            hfile->cd(dirname);
            //cout << obj->GetName() << endl;
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
        delete tree;
        delete timeDistAP;
        delete timeDistDeXT;
        delete timeAP;
        delete timeDeXT;
        for(auto const &e : gAmp_vs_Time) delete e.second;
        delete c_t;
        //delete waveform;
    }
    
    int fillStyle = 3003;
    // Save TTree with hist of noise for each event OV and the noise classification
    hfile->cd();
    otree->Write();
    
    // Create final plot of total primary correlated noise
    TCanvas * cfinal = new TCanvas("Correlated noise","Correlated noise",100,100,900,700);
    
    Double_t tot_max_noise = TMath::MaxElement(results["#Total"]->GetN(),results["#Total"]->GetY());
    
    formatGr(results["#Total"], kBlack, fillStyle, "#DeltaV [V]", "Correlated noise [%]", "Correlated Noise");
    results["#Total"]->GetYaxis()->SetRangeUser(0,tot_max_noise+2);
    formatGr(results["#DiXT"], kBlue, fillStyle, "#DeltaV [V]", "Correlated noise [%]");
    formatGr(results["#AP"], kOrange+7, fillStyle, "#DeltaV [V]", "Correlated noise [%]");
    formatGr(results["#DeXT"], kGreen+2, fillStyle, "#DeltaV [V]", "Correlated noise [%]");
    formatGr(results["#SecPeaks"], 7, fillStyle, "#DeltaV [V]", "Correlated noise [%]");
    
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
    formatGr(results["#DCR"], 2, fillStyle, "#DeltaV [V]", "DCR contribution");
    
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
    
    formatGr(results["#CorrectionFre"], kRed, fillStyle, "#DeltaV [V]", "Correction", "Correction");
    formatGr(results["#CorrectionCur"], kBlue,fillStyle, "#DeltaV [V]", "Correction", "Correction");
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
    formatGr(results["#SecPeaksDiXT"], kBlue, fillStyle, "#DeltaV [V]", "<N Secondary peaks>");
    results["#SecPeaksDiXT"]->GetYaxis()->SetRangeUser(0,tot_max_peaks*2.);
    formatGr(results["#SecPeaksDel"], kGreen+2, fillStyle, "#DeltaV [V]", "<N Secondary peaks>");
    results["#SecPeaksDiXT"]->Draw("ALP*3");
    results["#SecPeaksDel"]->Draw("LP*3");
    
    leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#SecPeaksDiXT"],"After Direct Cross-Talk","f");
    leg->AddEntry(results["#SecPeaksDel"],"After any delayed noise","f");
    leg->SetFillColor(kWhite);
    leg->Draw();
    
    cfinal->Print(globalArgs.res_folder+"SecondaryPeaks.pdf");
    cfinal->SetName("SecondaryPeaks");
    cfinal->SetTitle("SecondaryPeaks");
    cfinal->Write();
    
    // Charge
    formatGr(results["1PE-Waveform-Charge"], kGreen+2, fillStyle, "#DeltaV [V]", "Charge [mV#timesns]");
    formatGr(results["1PE-Charge-Fit"], kBlue, fillStyle, "#DeltaV [V]", "Charge [mV#timesns]");
    formatGr(results["Mean-Charge"], kRed, fillStyle, "#DeltaV [V]", "Charge [mV#timesns]");
    formatGr(results["Mean-Charge-Corrected"], kMagenta+1, fillStyle, "#DeltaV [V]", "Charge [mV#timesns]");
    //~ results["Mean-Charge-Corrected"]->SetLineStyle(4);
    results["Mean-Charge"]->Draw("APL3");
    results["Mean-Charge"]->SetTitle("Waveform charge");
    results["1PE-Waveform-Charge"]->Draw("PL+3");
    results["1PE-Charge-Fit"]->Draw("PL+3");
    results["Mean-Charge-Corrected"]->Draw("PL+3");
    TF1 * linfit = new TF1("lin_vbd_fit","[0]*(x-[1])");
    linfit->SetParName(0,"Gain");
    linfit->SetParName(1,"V_{BD}");
    results["1PE-Charge-Fit"]->Fit("lin_vbd_fit","Q+");
    
    TPaveText * pv = new TPaveText(0.15,0.55,0.45,0.65,"brNDC");
    pv->AddText(Form("Q/#DeltaV = %2.2f#pm%2.2f mV#timesns/V",linfit->GetParameter(0),linfit->GetParError(0)));
    pv->AddText(Form("V_{BD} = %2.0f#pm%2.0f mV",linfit->GetParameter(1)*1e3,linfit->GetParError(1)*1e3));
    pv->SetFillColor(kWhite);
    pv->Draw();
    
    leg = new TLegend(0.15,0.70,0.50,0.87);
    leg->AddEntry(results["1PE-Waveform-Charge"],"1PE (clean waveforms)","f");
    leg->AddEntry(results["1PE-Charge-Fit"],"1PE (charge distr. fit)","f");
    leg->AddEntry(results["Mean-Charge"],"Mean (all waveforms)","f");
    leg->AddEntry(results["Mean-Charge-Corrected"],"Mean (corrected)","f");
    leg->SetFillColor(kWhite);
    leg->Draw();
    cfinal->Print(globalArgs.res_folder+"Charge.pdf");
    cfinal->SetName("Charge");
    cfinal->SetTitle("Charge");
    cfinal->Write();
    
    // ------------
    hfile->Close();
    values_for_pde->close();
    
    return 0;
}


