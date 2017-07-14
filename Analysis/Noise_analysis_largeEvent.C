//
//  Created by Felipe Gilberto Ortega on 25/04/16.
//  Modified by Luca Pescatore, luca.pescatore@cern.ch on 28/06/2017   
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
        {"#DeXT",         new TGraphErrors()},
        {"#AP",          new TGraphErrors()},
        {"#Total",       new TGraphErrors()},
        {"#SecPeaks",    new TGraphErrors()},
        {"#SecPeaksDiXT",  new TGraphErrors()},
        {"#SecPeaksDeXT",  new TGraphErrors()} };

    Char_t Category[15];
    Double_t dV;

    TFile * hfile = TFile::Open(TString(globalArgs.res_folder) + "noiseanalysis.root","RECREATE");
    initOutputFile(hfile, vol_folders);
    
    TTree * otree = new TTree("ClassifiedData","ClassifiedData");
    int npts, noise_peaks_cnt;
    int direct_xtalk_pulse, xtalk_pulse, after_pulse;
    double times[10000], amps[10000];
    double Vbias, pe, DiXT_thr, DeXT_thr, AP_thr;

    otree->Branch("Category",Category,"Category/C");	// Categories are: DiXT, DeXT, AP, Clean
    otree->Branch("dV",&dV,"dV/D"); 
    otree->Branch("NsampPerEv",&npts,"NsampPerEv/I");
    otree->Branch("Amps",&amps,"Amps[NsampPerEv]/D");
    otree->Branch("Times",&times,"Times[NsampPerEv]/D");
    otree->Branch("NnoisePeaks",&noise_peaks_cnt,"NnoisePeaks/I");
    otree->Branch("NAP",&after_pulse,"NAP/I");
    otree->Branch("NDeXT",&xtalk_pulse,"NDeXT/I");
    otree->Branch("Vbias",&Vbias,"Vbias/D");
    otree->Branch("pe",&pe,"pe/D");
    otree->Branch("DiXT_thr",&DiXT_thr,"DiXT_thr/D");
    otree->Branch("DeXT_thr",&DeXT_thr,"DeXT_thr/D");
    otree->Branch("AP_thr",&AP_thr,"AP_thr/D");
    //~ hfile->mkdir("Plots");
    //~ hfile->mkdir("Waves");
    //~ hfile->mkdir("Fits");

    /////////////////
    // Calculate Voltage breakdown and value of pe
    /////////////////
    
    vector <Double_t> pe_volt;
    TGraph * Vbias_vs_pe = new TGraph();

    std::cout << " -----> Calculation of average DCR amplitudes " << std::endl;
    int i = 0;
    for (auto vol : vol_folders) 
    {
        pe_volt.push_back(Amplitude_calc(vol, data_size, "root", hfile));
        Vbias_vs_pe->SetPoint(i, pe_volt.back(), vol.Atof()); i++;
        std::cout << Form("Voltage: %.1fV --> Mean amplitude = %.4f", vol.Atof(), pe_volt.back()) << std::endl;
    }

    cout << "\n\n-----> Voltage Breakdown fit" << endl;

    double VBD = fitBreakdownVoltage(Vbias_vs_pe, hfile);


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
        cout << "\n\n   --> Voltage analyzed: " << vol << endl;

        map<string, int> color{{"clean",1},{"AP",1},{"DiXT",1},{"DeXT",1}};

        // Counters 

        unsigned int events_cnt = 0;
        unsigned int tot_noise_peaks_cnt = 0;
        unsigned int tot_double_cnt = 0;
        unsigned int counter_notclean = 0;
        unsigned int direct_xtalk_pulse_cnt = 0;
        unsigned int xtalk_pulse_cnt = 0;
        unsigned int after_pulse_cnt = 0;
        unsigned int all_double_DiXT_cnt = 0;
        unsigned int nsaved(0), nsaved_DiXT(0), nsaved_DeXT(0), nsaved_AP(0);	// counters on the number of waveforms in each canvas
        unsigned int maxNwaveforms = 100; // maximum number of waveforms in the canvas
        
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
        
        // persistence waveforms
        map<string,TCanvas *> canv_persistence {
            {"DiXT", new TCanvas(Form("persistence_DiXT_%s",vol), Form("Persistence direct cross-talk #DeltaV = %2.2f V",dV))},
            {"DeXT", new TCanvas(Form("persistence_DeXT_%s",vol), Form("Persistence delayed cross-talk #DeltaV = %2.2f V",dV))},
            {"AP",   new TCanvas(Form("persistence_AP_%s",vol), Form("Persistence after-pulse #DeltaV = %2.2f V",dV))},
            {"clean",new TCanvas(Form("persistence_Clean_%s",vol), Form("Persistence Clean #DeltaV = %2.2f V",dV))}
        };
        
        TGraph * Expfit_AP      = new TGraph();
        TGraph * cleanforfit    = NULL;
        
	    TH1D * AP_arrivaltime   = new TH1D(Form("AP_arrival_time_%s",vol),"AP arrival time distribution", 140, 0, 0.2e-6);
        TH1D * DeXT_arrivaltime = new TH1D(Form("DeXT_arrival_time_%s",vol),"DeXT arrival time distribution", 140, 0, 0.2e-6);
        TH1D * Npeaks           = new TH1D(Form("Npeaks_not_clean_%s",vol),"Number of noise peaks when not clean", 50, 0, 50);
        TH1D * NpeaksDiXT       = new TH1D(Form("Npeaks_when_DiXT_%s",vol),"Number of noise peaks when DiXT", 50, 0, 50);
        TH1D * NpeaksDeXT       = new TH1D(Form("Npeaks_when_DeXT_%s",vol),"Number of noise peaks when DeXT", 50, 0, 50);
        TMultiGraph * forfit    = new TMultiGraph();        
        TH1D * Charge           = new TH1D(Form("Charge_%s",vol),"Charge of all peaks", 10, 0., 0.2);
        
        double bin_size = 0.002;	// vertical bin size for persistence plot
        TH2D * persistence_clean = NULL;
        TH2D * persistence_DiXT  = NULL;
        TH2D * persistence_DeXT  = NULL;
        TH2D * persistence_AP    = NULL;

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
            double hmin = tree->GetMinimum("Times");
            double hmax = tree->GetMaximum("Times");
            int Nv = int(((tree->GetMaximum("Amps")-tree->GetMinimum("Amps"))/bin_size)+0.5);
            double vmin = tree->GetMinimum("Amps");
            double vmax = tree->GetMaximum("Amps");
            persistence_clean = new TH2D("Clean_persistence","Clean waveforms", Nh, hmin, hmax, Nv, vmin,vmax);
            persistence_DiXT    = new TH2D("DiXT_persistence","Direct cross-talk", Nh, hmin, hmax, Nv, vmin,vmax);
            persistence_DeXT   = new TH2D("DeXT_persistence","Delayed cross-talk", Nh, hmin, hmax, Nv, vmin,vmax);
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
            all_double_DiXT_cnt = 0;
            direct_xtalk_pulse = 0;
            
            DiXT_thr = cthrs.dir_xtalk * pe;
            DeXT_thr = cthrs.del_xtalk * pe;
            AP_thr = cthrs.AP * pe;
            if(globalArgs.fixed_thr > 0.) 
            {
                DiXT_thr = globalArgs.fixed_thr;
                DeXT_thr = globalArgs.fixed_thr;
                AP_thr = globalArgs.fixed_thr;
            }

            bool done = false;
            bool veryclean = true;
            for (int row = 2; row < npts-2; row++) 
            {
                double curT = time[row];
                double curV = volts[row];

                // Skip if not on a maximum                
                if( curT < 0. || !( curV > volts[row-2] && curV > volts[row+2]) ) continue;
                else if (done) { done = false; continue; }
                done = true;

                if ( curV > AP_thr ) Charge->Fill(curV);
                if ((curT > 2*ns) && (curV > 0.4*pe)) veryclean = false;

                // Direct x-talk: in 0-2 ns window and V > direct th.
                if( curT <= cthrs.dir_xtalk_maxT * ns && curV > DiXT_thr ) 
                {
                    direct_xtalk_pulse++;
                    
                    all_double_DiXT_cnt++;	// it seems it does not count double DiXT, but simply all DiXT, no?
                    
                    //noise_peaks_cnt++; // Don't count peak as it is over the DCR?
                }

                // Delayed x-talk: time larger than end of DCR window and V > delayed XT th
                else if ( curT > cthrs.dir_xtalk_maxT * ns && curV > DeXT_thr )
                {
                    xtalk_pulse++;
                    noise_peaks_cnt++;
                    DeXT_arrivaltime->Fill(curT);
                    if(curV > DiXT_thr) all_double_DiXT_cnt++;
                }

                // After-pulse: time larger 2ns and V > AP th
                else if ( curT > cthrs.AP_minT * ns && curV > AP_thr)
                {
                    after_pulse++;
                    noise_peaks_cnt++;
                    Expfit_AP->SetPoint(Expfit_AP->GetN(),curT,curV);
                    AP_arrivaltime->Fill(curT);
                    if(curV > DiXT_thr) all_double_DiXT_cnt++;
                }

            } // loop over time
                

            if (direct_xtalk_pulse > 0)  // Check for imm x-talk and plot
            {
                if(direct_xtalk_pulse > 1) cout << "Attention: more then one direct DiXT found" << endl;
                direct_xtalk_pulse_cnt++;
                counter_notclean++;
                sprintf(Category,"DiXT");
                canv["DiXT"]->cd();

                TString graph_title = Form("Direct CrossTalk #DeltaV = %2.2f V",dV);
                if(nsaved_DiXT < maxNwaveforms) {
					drawWave(waveform, &color["DiXT"], graph_title, canv["DiXT"], 1.5*pe);
					nsaved_DiXT++;
				}
                fillPersistence(persistence_DiXT,waveform);
                
                Npeaks->Fill(noise_peaks_cnt);
                NpeaksDiXT->Fill(noise_peaks_cnt);
            }
            
            else if (xtalk_pulse > 0) // Delayed x-talk
            {
                xtalk_pulse_cnt++;
                counter_notclean++;
                sprintf(Category,"DeXT");
                canv["DeXT"]->cd();
                
                TString graph_title = Form("Delayed cross-talk #DeltaV = %2.2f V",dV);
                if(nsaved_DeXT < maxNwaveforms) {
					drawWave(waveform, &color["DeXT"], graph_title, canv["DeXT"], 1.5*pe);
					nsaved_DeXT++;
				}
                fillPersistence(persistence_DeXT,waveform);

                Npeaks->Fill(noise_peaks_cnt);
                NpeaksDeXT->Fill(noise_peaks_cnt-1);
            }
            
            else if (after_pulse > 0) //  Only after pulse
            {
                after_pulse_cnt++;
                counter_notclean++;
                sprintf(Category,"AP");

                TString graph_title = Form("After pulse #DeltaV = %2.2f V",dV);
                if(nsaved_AP < maxNwaveforms) {
					drawWave(waveform, &color["AP"], graph_title, canv["AP"], 1.5*pe);
					nsaved_AP++;
				}
                fillPersistence(persistence_AP, waveform);

                Npeaks->Fill(noise_peaks_cnt);
                NpeaksDeXT->Fill(noise_peaks_cnt-1);
            }
	           
            else    // If not noisy then it's clean
            {
                sprintf(Category,"Clean");

                // Get only very clean waves and make an average for long tau fit.
                bool veryclean = true;
                for (int row = 0; row < npts; row++) 
                    if ((time[row] > 2*ns) && (volts[row] > 0.4*pe)) 
                        veryclean = false;

                if(veryclean) cleanforfit = average(cleanforfit, waveform);	// actually this is wrong !

                if(nsaved < maxNwaveforms) // Max 20 clean graphs on the plot
                {
                    forfit->Add(waveform);
                    nsaved++;
                    drawWave(waveform, &color["clean"], Form("Clean pulse #DeltaV = %2.2f V",dV), canv["clean"], 1.5*pe);                  
                }
                fillPersistence(persistence_clean, waveform);
	        }

            tot_noise_peaks_cnt += noise_peaks_cnt;
            tot_double_cnt += all_double_DiXT_cnt;

            if(globalArgs.save_all) otree->Fill();
            
            delete time;
            delete volts;
        }
        
        double tot_npeaks = tot_noise_peaks_cnt + events_cnt;
        double perc_noise_peaks = tot_noise_peaks_cnt / (float)tot_npeaks;
        double perc_DiXT  = direct_xtalk_pulse_cnt / (float)events_cnt;
        double perc_DeXT = xtalk_pulse_cnt / (float)events_cnt;
        double perc_AP  = after_pulse_cnt / (float)events_cnt;
        double perc_double = tot_double_cnt / (float)tot_npeaks;
        double perc_Sec = (tot_noise_peaks_cnt - after_pulse_cnt - xtalk_pulse_cnt) / (float)tot_npeaks;

        cout << "\nTotal number of events: " << events_cnt << endl;
        cout << "PE: " << pe << "\n" << endl;

        cout << Form("Not clean events: %i [%.2f%%]",counter_notclean, (float)counter_notclean/events_cnt*100) << endl;
        cout << Form("     (DirXtalk = %i [%.2f%%], DelXtalk = %i [%.2f%%], AP = %i [%.2f%%])",
            (int)direct_xtalk_pulse_cnt, perc_DiXT*100.,
            (int)xtalk_pulse_cnt, perc_DeXT*100,
            (int)after_pulse_cnt, perc_AP*100 ) << endl;
        cout << Form("Total probability of secondary peaks: %.2f%%",perc_Sec*100.) << endl;
        cout << Form("Percent of noise peaks over the total (P_all): %.2f%%",perc_noise_peaks*100.) << endl;
        cout << Form("Percent of double peaks: %.2f%%",perc_double*100.) << endl;
        cout << Form("Mean peak charge: %.4f pe", Charge->GetMean()) << endl;
        cout << Form("Mean event charge: %.4f pe", Charge->Integral() / events_cnt) << endl;
        

        // Fits for long tau and recovery time

        cout << "\n\n-----> Long tau fit ***" << endl;
        double amp0, tau;
        //TF1 * exp_longtau = fitLongTau(cleanforfit, &amp0, &tau, pe, vol, canv["clean"]);
        TF1 * exp_longtau = fitLongTau(forfit, &amp0, &tau, pe, vol, canv["clean"],cleanforfit);	// does not work properly
        
        cout << "\n\n-----> After pulse fit ***" << endl;
        TF1 * exp_AP = fitAPTau(Expfit_AP, amp0, tau, pe, vol, canv["AP"]);
        
        results["#DiXT"]->SetPoint(i,dV,perc_DiXT*100.);
        results["#DiXT"]->SetPointError(i,0.,TMath::Sqrt(perc_DiXT*(1.-perc_DiXT)/events_cnt)*100.);
        results["#AP"]->SetPoint(i,dV,perc_AP*100.);
        
        results["#AP"]->SetPointError(i,0.,TMath::Sqrt(perc_AP*(1.-perc_AP)/events_cnt)*100.);
        results["#DeXT"]->SetPoint(i,dV,perc_DeXT*100.);
        results["#DeXT"]->SetPointError(i,0.,TMath::Sqrt(perc_DeXT*(1.-perc_DeXT)/events_cnt)*100.);

        results["#Total"]->SetPoint(i,dV,perc_noise_peaks*100.);
        results["#Total"]->SetPointError(i,0.,TMath::Sqrt(perc_noise_peaks*(1.-perc_noise_peaks)/(tot_noise_peaks_cnt + events_cnt))*100.);
        results["#SecPeaks"]->SetPoint(i,dV,perc_Sec*100.);
        results["#SecPeaks"]->SetPointError(i,0.,TMath::Sqrt(perc_Sec*(1.-perc_Sec)/events_cnt)*100.);
        results["#SecPeaksDiXT"]->SetPoint(i,dV,NpeaksDiXT->GetMean());
        results["#SecPeaksDiXT"]->SetPointError(i,0.,NpeaksDiXT->GetMeanError());
        results["#SecPeaksDeXT"]->SetPoint(i,dV,NpeaksDeXT->GetMean());
        results["#SecPeaksDeXT"]->SetPointError(i,0.,NpeaksDeXT->GetMeanError());
        
        // Final persistence plots
        double amp0_2, tau_2;
        TF1 * exp_longtau2 = drawPersistenceWithLongTauFit(persistence_clean,canv_persistence["clean"], &amp0_2, &tau_2, pe, vol,Form("Clean waveforms #DeltaV = %2.2f V",dV));
        drawPersistence(persistence_DiXT,canv_persistence["DiXT"],Form("Direct cross-talk #DeltaV = %2.2f V",dV));
        drawPersistence(persistence_DeXT,canv_persistence["DeXT"],Form("Delayed cross-talk #DeltaV = %2.2f V",dV));
        drawPersistence(persistence_AP,canv_persistence["AP"],Form("After-pulse #DeltaV = %2.2f V",dV));
        
        // Save/print results:
        vector<TObject*> objects_to_save;
        // Waveform canvas        
        for(auto const &e : canv) objects_to_save.push_back(e.second);
        // fit for long tau and AP
        objects_to_save.push_back(exp_longtau);
        objects_to_save.push_back(exp_longtau2);
        objects_to_save.push_back(exp_AP);
        // Time distributions and Npeak distributions
        objects_to_save.push_back(AP_arrivaltime);
 	    objects_to_save.push_back(DeXT_arrivaltime);
        objects_to_save.push_back(Npeaks);
        objects_to_save.push_back(NpeaksDiXT);
        objects_to_save.push_back(NpeaksDeXT);
        // Persistence plots
        for(auto const &e : canv_persistence) objects_to_save.push_back(e.second);
        
        for(unsigned int i(0); i<objects_to_save.size(); ++i) {
			TObject * obj = objects_to_save[i];
			hfile->cd();
			TString dirname = "pulse_shape_" + TString(vol);
			hfile->cd(dirname);
			obj->Write();
			hfile->cd();
			if(globalArgs.save_all)
				if(!(obj->InheritsFrom(TF1::Class()) || obj->InheritsFrom(TF2::Class())))
					objects_to_save[i]->SaveAs(globalArgs.res_folder+Form("%s_%s.pdf",objects_to_save[i]->GetName(),vol));
		}

	    delete AP_arrivaltime;
	    delete DeXT_arrivaltime;
        delete Npeaks;
        delete NpeaksDiXT;
        delete NpeaksDeXT;
        delete Expfit_AP;
        delete cleanforfit;
        delete exp_longtau;
        delete exp_AP;
        for(auto const &e : canv) delete e.second;
        for(auto const &e : canv_persistence) delete  e.second;
        delete persistence_clean;
		delete persistence_DiXT;
		delete persistence_DeXT;
		delete persistence_AP;
    }

    // Save TTree with hist of noise for each event OV and the noise classification
    hfile->cd();
    otree->Write();

    // Create final plot of total correlated noise
    TCanvas * cfinal = new TCanvas("Correlated Noise","Correlated Noise",100,100,900,700);
    
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
    results["#SecPeaks"]->Draw("LP*3");
    
    TLegend * leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#Total"],"Total","f");
    leg->AddEntry(results["#DiXT"],"Direct cross-talk","f");
    leg->AddEntry(results["#AP"],"After-pulse","f");
    leg->AddEntry(results["#DeXT"],"Delayed cross-talk","f");
    leg->AddEntry(results["#SecPeaks"],"Secondary noise","f");

    leg->SetFillColor(kWhite);
    leg->Draw();
    
    cfinal->SetGrid();
    cfinal->Print(globalArgs.res_folder+"CorrelatedNoise.pdf");
    cfinal->Write();

    Double_t tot_max_peaks = TMath::MaxElement(results["#SecPeaksDiXT"]->GetN(),results["#SecPeaksDiXT"]->GetY());
    results["#SecPeaksDiXT"]->SetTitle("Secondary peaks after a noise peak");
    formatGr(results["#SecPeaksDiXT"], kBlue, 0, "OverVoltage [V]", "<N Secondary peaks>");
    results["#SecPeaksDiXT"]->GetYaxis()->SetRangeUser(0,tot_max_peaks*2.);
    formatGr(results["#SecPeaksDeXT"], kGreen+2, 0, "OverVoltage [V]", "<N Secondary peaks>");
    results["#SecPeaksDiXT"]->Draw("ALP*");
    results["#SecPeaksDeXT"]->Draw("LP*");

    leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#SecPeaksDiXT"],"After Direct Cross-Talk","lp");
    leg->AddEntry(results["#SecPeaksDeXT"],"After any delayed noise","lp");
    leg->SetFillColor(kWhite);
    leg->Draw();
    
    cfinal->Print(globalArgs.res_folder+"SecondaryPeaks.pdf");
    cfinal->SetName("SecondaryPeaks");
    cfinal->Write();
    
    hfile->Close();
    return 0;
}


