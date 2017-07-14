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
        {"#CT",          new TGraphErrors()},
        {"#DCT",         new TGraphErrors()},
        {"#AP",          new TGraphErrors()},
        {"#Total",       new TGraphErrors()},
        {"#SecPeaks",    new TGraphErrors()},
        {"#SecPeaksCT",  new TGraphErrors()},
        {"#SecPeaksDT",  new TGraphErrors()},
        {"Charge",       new TGraphErrors()},
        {"#Double",      new TGraphErrors()} };

    Char_t Category[15];
    Double_t V_meas;

    TFile * hfile = TFile::Open(TString(globalArgs.res_folder) + "noiseanalysis.root","RECREATE");
    TTree * otree = new TTree("ClassifiedData","ClassifiedData");
    int npts, noise_peaks_cnt;
    int xtalk_pulse, after_pulse;
    double times[10000], amps[10000];
    double V, pe, CT_thr, DCT_thr, AP_thr;

    otree->Branch("Category",Category,"Category/C");
    otree->Branch("V_meas",&V_meas,"V_meas/D"); 
    otree->Branch("NsampPerEv",&npts,"NsampPerEv/I");
    otree->Branch("Amps",&amps,"Amps[NsampPerEv]/D");
    otree->Branch("Times",&times,"Times[NsampPerEv]/D");
    otree->Branch("NnoisePeaks",&noise_peaks_cnt,"NnoisePeaks/I");
    otree->Branch("NAfterPulses",&after_pulse,"NAfterPulses/I");
    otree->Branch("NDelayedCT",&xtalk_pulse,"NDelayedCT/I");
    otree->Branch("V",&V,"V/D");
    otree->Branch("pe",&pe,"pe/D");
    otree->Branch("CT_thr",&CT_thr,"CT_thr/D");
    otree->Branch("DCT_thr",&DCT_thr,"DCT_thr/D");
    otree->Branch("AP_thr",&AP_thr,"AP_thr/D");
    hfile->mkdir("Plots");
    hfile->mkdir("Waves");
    hfile->mkdir("Fits");

    /////////////////
    // Calculate Voltage breakdown and value of pe
    /////////////////
    
    vector <Double_t> pe_volt;
    TGraph * Vbias = new TGraph();

    std::cout << " -----> Calculation of average DCR amplitudes " << std::endl;
    int i = 0;
    for (auto vol : vol_folders)
    {
        pe_volt.push_back(Amplitude_calc(vol, data_size, "root", hfile));
        Vbias->SetPoint(i, pe_volt.back(), vol.Atof()); i++;
        std::cout << Form("Voltage: %.1fV --> Mean amplitude = %.4f", vol.Atof(), pe_volt.back()) << std::endl;
    }

    cout << "\n\n-----> Voltage Breakdown fit" << endl;

    double VBD = fitBreakdownVoltage(Vbias);


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

        map<string, int> color{{"clean",1},{"AP",1},{"CT",1},{"DCT",1}};

        // Counters 

        unsigned int events_cnt = 0;
        unsigned int tot_noise_peaks_cnt = 0;
        unsigned int tot_double_cnt = 0;
        unsigned int counter_notclean = 0;
        unsigned int direct_xtalk_pulse_cnt = 0;
        unsigned int xtalk_pulse_cnt = 0;
        unsigned int after_pulse_cnt = 0;
        unsigned int all_double_CT_cnt = 0;
        unsigned int direct_xtalk_pulse = 0;
        unsigned int nsaved = 0;
        
        // Define amplitude measured at which OV
        
        V = vol_folders[i].Atof();
        pe = pe_volt[i];
        V_meas = V - VBD;
         
        map<string,TCanvas *> canv {
            {"CT",     new TCanvas(Form("Direct CrossTalk OV = %2.2f V",V_meas))},
            {"DCT",    new TCanvas(Form("Delayed CrossTalk OV = %2.2f V",V_meas))},
            {"AP",     new TCanvas(Form("After Pulse OV = %2.2f V",V_meas))},
            {"clean",  new TCanvas(Form("Clean OV = %2.2f V",V_meas))}
        };
        
        hfile->cd("Plots");
        TMultiGraph * forfit    = new TMultiGraph();  
        TGraph * Expfit_AP      = new TGraph();
        TGraph * cleanforfit    = NULL;
        map<string,TH1D*> graphs {
	       {"AP_arrivaltime",   new TH1D(Form("AP_arrival_times_%s",vol),"AP arrival times", 140, 0, 0.2e-6)},
           {"DeXT_arrivaltime", new TH1D(Form("Del_CT_arrival_times_%s",vol),"DeCT arrival times", 140, 0, 0.2e-6)},
           {"Npeaks",           new TH1D(Form("N_peaks_%s",vol),"Number of noise peaks when not clean", 50, 0, 50)},
           {"NpeaksCT",         new TH1D(Form("N_peaks_CT_%s",vol),"Number of noise peaks when direct CT", 50, 0, 50)},
           {"NpeaksDT",         new TH1D(Form("N_peaks_Delayed_%s",vol),"Number of noise peaks when delayed noise", 50, 0, 50)},
           {"Charge",           new TH1D(Form("Charge_%s",vol),"Charge of all peaks", 10, 0., 0.2)} };

        // Setup input tree
        TTree * tree = NULL;
        if(globalArgs.input=="root") 
        {   
            tree = (TTree*)ifile->Get(vol);
            tree->SetBranchAddress("NsampPerEv",&npts);
            tree->SetBranchAddress("Amps",&amps);
            tree->SetBranchAddress("Times",&times);
            data_size = tree->GetEntries();
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
            all_double_CT_cnt = 0;
            direct_xtalk_pulse = 0;
            
            CT_thr = cthrs.dir_xtalk * pe;
            DCT_thr = cthrs.del_xtalk * pe;
            AP_thr = cthrs.AP * pe;
            if(globalArgs.fixed_thr > 0.) 
            {
                CT_thr = globalArgs.fixed_thr;
                DCT_thr = globalArgs.fixed_thr;
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

                if ( curV > AP_thr ) graphs["Charge"]->Fill(curV);
                if ( curT > 2*ns && curV > 0.4*pe) veryclean = false;

                // Direct x-talk: in 0-2 ns window and V > direct th.
                if( curT <= cthrs.dir_xtalk_maxT * ns && curV > CT_thr ) 
                {
                    direct_xtalk_pulse++;
                    all_double_CT_cnt++;
                    //noise_peaks_cnt++; // Don't count peak as it is over the DCR?
                }

                // Delayed x-talk: time larger than end of DCR window and V > delayed CT th
                else if ( curT > cthrs.dir_xtalk_maxT * ns && curV > DCT_thr )
                {
                    xtalk_pulse++;
                    noise_peaks_cnt++;
                    graphs["DeXT_arrivaltime"]->Fill(curT);
                    if(curV > CT_thr) all_double_CT_cnt++;
                }

                // After-pulse: time larger 2ns and V > AP th
                else if ( curT > cthrs.AP_minT * ns && curV > AP_thr)
                {
                    after_pulse++;
                    noise_peaks_cnt++;
                    Expfit_AP->SetPoint(Expfit_AP->GetN(),curT,curV);
                    graphs["AP_arrivaltime"]->Fill(curT);
                    if(curV > CT_thr) all_double_CT_cnt++;
                }

            } // loop over time
                

            if (direct_xtalk_pulse > 0)  // Check for imm x-talk and plot
            {
                if(direct_xtalk_pulse > 1) cout << "Attention: more then one direct CT found" << endl;
                direct_xtalk_pulse_cnt++;
                counter_notclean++;
                sprintf(Category,"Direct_Crosstalk");
                canv["CT"]->cd();

                TString graph_title = Form("Direct CrossTalk OV = %2.2f V",V_meas);
                drawWave(waveform, &color["CT"], graph_title, canv["CT"], 1.5*pe);
                
                graphs["Npeaks"]->Fill(noise_peaks_cnt);
                graphs["NpeaksCT"]->Fill(noise_peaks_cnt);
            }
            
            else if (xtalk_pulse > 0) // Delayed x-talk
            {
                xtalk_pulse_cnt++;
                counter_notclean++;
                sprintf(Category,"Delayed_Crosstalk");
                canv["DCT"]->cd();
                
                TString graph_title = Form("Delayed cross-talk OV = %2.2f V",V_meas);
                drawWave(waveform, &color["DCT"], graph_title, canv["DCT"], 1.5*pe);

                graphs["Npeaks"]->Fill(noise_peaks_cnt);
                graphs["NpeaksDT"]->Fill(noise_peaks_cnt-1);
            }
            
            else if (after_pulse > 0) //  Only after pulse
            {
                after_pulse_cnt++;
                counter_notclean++;
                sprintf(Category,"After_Pulse");

                TString graph_title = Form("After pulse OV = %2.2f V",V_meas);
                drawWave(waveform, &color["AP"], graph_title, canv["AP"], 1.5*pe);

                graphs["Npeaks"]->Fill(noise_peaks_cnt);
                graphs["NpeaksDT"]->Fill(noise_peaks_cnt-1);
            }
	           
            else    // If not noisy then it's clean
            {
                sprintf(Category,"Clean");

                // Get only very clean waves and make an average for long tau fit.
                if(veryclean) cleanforfit = average(cleanforfit, waveform);

                if(nsaved < 5) // Max 20 clean graphs on the plot
                {
                    forfit->Add(waveform);
                    nsaved++;
                    drawWave(waveform, &color["clean"], Form("Clean pulse OV = %2.2f V",V_meas), canv["clean"], 1.5*pe);                  
                }
	        }

            tot_noise_peaks_cnt += noise_peaks_cnt;
            tot_double_cnt += all_double_CT_cnt;

            if(globalArgs.save_all) otree->Fill();
            
            delete time;
            delete volts;
        }
       
        double tot_npeaks = tot_noise_peaks_cnt + events_cnt;
        double perc_noise_peaks = tot_noise_peaks_cnt / (float)tot_npeaks;
        double perc_CT  = direct_xtalk_pulse_cnt / (float)events_cnt;
        double perc_DCT = xtalk_pulse_cnt / (float)events_cnt;
        double perc_AP  = after_pulse_cnt / (float)events_cnt;
        double perc_double = tot_double_cnt / (float)tot_npeaks;
        double perc_Sec = (tot_noise_peaks_cnt - after_pulse_cnt - xtalk_pulse_cnt) / (float)tot_npeaks;

        cout << "\nTotal number of events: " << events_cnt << endl;
        cout << "PE: " << pe << "\n" << endl;

        cout << Form("Not clean events: %i [%.2f%%]",counter_notclean, (float)counter_notclean/events_cnt*100) << endl;
        cout << Form("     (DirXtalk = %i [%.2f%%], DelXtalk = %i [%.2f%%], AP = %i [%.2f%%])",
            (int)direct_xtalk_pulse_cnt, perc_CT*100.,
            (int)xtalk_pulse_cnt, perc_DCT*100,
            (int)after_pulse_cnt, perc_AP*100 ) << endl;
        cout << Form("Total probability of secondary peaks: %.2f%%",perc_Sec*100.) << endl;
        cout << Form("Percent of noise peaks over the total (P_all): %.2f%%",perc_noise_peaks*100.) << endl;
        cout << Form("Percent of double peaks: %.2f%%",perc_double*100.) << endl;
        cout << Form("Mean peak charge: %.4f pe", graphs["Charge"]->GetMean()) << endl;
        cout << Form("Mean event charge: %.4f pe", graphs["Charge"]->Integral() / events_cnt) << endl;
        

        // Fits for long tau and recovery time

        hfile->cd("Fits");

        cout << "\n\n-----> Long tau fit ***" << endl;
        double amp0, tau;
        //fitLongTau(cleanforfit, &amp0, &tau, pe, vol, canv["clean"]);
        fitLongTau(forfit, &amp0, &tau, pe, vol, canv["clean"],cleanforfit);

        cout << "\n\n-----> After pulse fit ***" << endl;
        fitAPTau(Expfit_AP, amp0, tau, pe, vol, canv["AP"])->Write();


        // Final result: Correlated noise + save/print everything

        results["#CT"]->SetPoint(i,V_meas,perc_CT*100.);
        results["#CT"]->SetPointError(i,0.,TMath::Sqrt(perc_CT*(1.-perc_CT)/events_cnt)*100.);
        results["#AP"]->SetPoint(i,V_meas,perc_AP*100.);
        results["#AP"]->SetPointError(i,0.,TMath::Sqrt(perc_AP*(1.-perc_AP)/events_cnt)*100.);
        results["#DCT"]->SetPoint(i,V_meas,perc_DCT*100.);
        results["#DCT"]->SetPointError(i,0.,TMath::Sqrt(perc_DCT*(1.-perc_DCT)/events_cnt)*100.);
        results["#Double"]->SetPoint(i,V_meas,perc_double*100.);
        results["#Double"]->SetPointError(i,0.,TMath::Sqrt(perc_double*(1.-perc_double)/events_cnt)*100.);
        results["Charge"]->SetPoint(i,V_meas,graphs["Charge"]->GetMean());
        results["Charge"]->SetPointError(i,0.,graphs["Charge"]->GetMeanError());

        results["#Total"]->SetPoint(i,V_meas,perc_noise_peaks*100.);
        results["#Total"]->SetPointError(i,0.,TMath::Sqrt(perc_noise_peaks*(1.-perc_noise_peaks)/(tot_noise_peaks_cnt + events_cnt))*100.);
        results["#SecPeaks"]->SetPoint(i,V_meas,perc_Sec*100.);
        results["#SecPeaks"]->SetPointError(i,0.,TMath::Sqrt(perc_Sec*(1.-perc_Sec)/events_cnt)*100.);
        results["#SecPeaksCT"]->SetPoint(i,V_meas,graphs["NpeaksCT"]->GetMean());
        results["#SecPeaksCT"]->SetPointError(i,0.,graphs["NpeaksCT"]->GetMeanError());
        results["#SecPeaksDT"]->SetPoint(i,V_meas,graphs["NpeaksDT"]->GetMean());
        results["#SecPeaksDT"]->SetPointError(i,0.,graphs["NpeaksDT"]->GetMeanError());
        
        // Save/print reults:
        hfile->cd("Waves");
        for(auto const &e : canv) e.second->Write();
	    hfile->cd("Plots");
        for(auto const &e : graphs) e.second->Write();
        
        canv["CT"]->Print(globalArgs.res_folder+Form("Immcrosstalk_%s.pdf",vol));
        canv["DCT"]->Print(globalArgs.res_folder+Form("Delcrosstalk_%s.pdf",vol));
        
        if(globalArgs.save_all) 
        {
            TCanvas * tmpc = new TCanvas();
            for(auto const &e : graphs) 
            {
                e.second->Draw();
                tmpc->Print(globalArgs.res_folder+Form((e.first+"_%s.pdf").c_str(),vol));
            }
            /*
            DeXT_arrivaltime->Draw();
            tmpc->Print(globalArgs.res_folder+Form("DeXT_arrivetime_%s.pdf",vol));
            Npeaks->Draw("HIST");
            tmpc->Print(globalArgs.res_folder+Form("Npeaks_%s.pdf",vol));
            NpeaksCT->Draw("HIST");
            tmpc->Print(globalArgs.res_folder+Form("Npeaks_whenCT_%s.pdf",vol));
            NpeaksDT->Draw("HIST");
            tmpc->Print(globalArgs.res_folder+Form("Npeaks_whenDelayedNoise_%s.pdf",vol));
            */
            delete tmpc;
        }

        delete Expfit_AP;
        delete cleanforfit;
        for(auto const &e : graphs) delete e.second;
        for(auto const &e : canv) delete e.second;
    }

    // Save TTree with hist of noise for each event OV and the noise classification
    hfile->cd();
    otree->Write();

    // Create final plot of total correlated noise
    TCanvas * cfinal = new TCanvas("Correlated Noise","Correlated Noise",100,100,900,700);
    
    Double_t tot_max_noise = TMath::MaxElement(results["#Total"]->GetN(),results["#Total"]->GetY());
    
    formatGr(results["#Total"], kBlack, 3005, "OverVoltage [V]", "Noise [%]", "Correlated Noise");
    results["#Total"]->GetYaxis()->SetRangeUser(0,tot_max_noise+2);
    formatGr(results["#CT"], kBlue, 3005, "OverVoltage [V]", "Noise [%]");
    formatGr(results["#AP"], kOrange+7, 3005, "OverVoltage [V]", "Noise [%]");
    formatGr(results["#DCT"], kGreen+2, 3005, "OverVoltage [V]", "Noise [%]");
    formatGr(results["#SecPeaks"], 7, 3005, "OverVoltage [V]", "Noise [%]");
    formatGr(results["#Double"], 8, 3005, "OverVoltage [V]", "Noise [%]");

    results["#Total"]->Draw("ALP*3");
    results["#CT"]->Draw("LP*3");
    results["#AP"]->Draw("LP*3");
    results["#DCT"]->Draw("LP*3");
    results["#SecPeaks"]->Draw("LP*3");
    results["#Double"]->Draw("LP*3");
    
    TLegend * leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#Total"],"Total","f");
    leg->AddEntry(results["#CT"],"Direct Cross-Talk","f");
    leg->AddEntry(results["#AP"],"After Pulse","f");
    leg->AddEntry(results["#DCT"],"Delayed Cross-Talk","f");
    leg->AddEntry(results["#SecPeaks"],"Secondary noise","f");
    leg->AddEntry(results["#Double"],"Double hight peaks","f");

    leg->SetFillColor(kWhite);
    leg->Draw();
    
    cfinal->SetGrid();
    cfinal->Print(globalArgs.res_folder+"CorrelatedNoise.pdf");
    cfinal->Write();

    Double_t tot_max_peaks = TMath::MaxElement(results["#SecPeaksCT"]->GetN(),results["#SecPeaksCT"]->GetY());
    results["#SecPeaksCT"]->SetTitle("Secondary peaks after a noise peak");
    formatGr(results["#SecPeaksCT"], kBlue, 0, "OverVoltage [V]", "<N Secondary peaks>");
    results["#SecPeaksCT"]->GetYaxis()->SetRangeUser(0,tot_max_peaks*2.);
    formatGr(results["#SecPeaksDT"], kGreen+2, 0, "OverVoltage [V]", "<N Secondary peaks>");
    results["#SecPeaksCT"]->Draw("ALP*");
    results["#SecPeaksDT"]->Draw("LP*");

    leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#SecPeaksCT"],"After Direct Cross-Talk","lp");
    leg->AddEntry(results["#SecPeaksDT"],"After any delayed noise","lp");
    leg->SetFillColor(kWhite);
    leg->Draw();
    
    cfinal->Print(globalArgs.res_folder+"SecondaryPeaks.pdf");
    cfinal->SetName("SecondaryPeaks");
    cfinal->Write();

    results["Charge"]->Draw("AP");
    cfinal->Print(globalArgs.res_folder+"Charge.pdf");
    cfinal->SetName("Charge");
    cfinal->Write();    
  
    return 0;
}


