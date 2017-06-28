
//
//  Noise_analysis.C
//  
//
//  Created by Felipe Gilberto Ortega on 25/04/16.
//
//
 
#include "Noise_analysis_largeEvent.h"
#include "Thresholds.h"
#include "fits.h"
#include "genparam.h"

static const char * optString = "d:S:o:ah?";

using namespace std;

int main(int argc, char* argv[]) 
{    
    // Get paremeters from the command line

    int opt = getopt(argc, argv, optString);
    if(opt == -1){
        std::cerr <<  "There is no opption in the command! Type \"output -h\" for help." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    while(opt != -1){
        switch(opt){
            case 'd':
                globalArgs.data_folder = optarg;
                std::cout << "-p option path= " << globalArgs.data_folder << std::endl;
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
            
    ifstream setupFile(globalArgs.arg_pathToSetupFile);
    if(!setupFile){
        std::cerr << "Failure: could not open file: \"" << globalArgs.arg_pathToSetupFile << "\"." << std::endl;
        std::cerr << "Please check if the path is correct or not!" << std::endl;
        exit(EXIT_FAILURE);
    }
    

    ////////////////
    //Define thresholds
    ////////////////

    // Read setup file and load input .root

    Int_t data_size;
    vector <Thresholds_t> thrs;
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
        {"#SecPeaksDCT", new TGraphErrors()} };

    Char_t Category[15];
    Double_t V_meas;

    TFile * hfile = TFile::Open(TString(globalArgs.res_folder) + "noiseanalysis.root","RECREATE");
    TTree * otree = new TTree("T","Noise Analysis");
    otree->Branch("Category",Category,"Category/C");
    otree->Branch("V_meas",&V_meas,"V_meas/D"); 


    /////////////////
    // Calculate Voltage breakdown and value of pe
    /////////////////
    
    vector <Double_t> pe_volt;
    TGraph * Vbias_ver = new TGraph();

    int i = 0;
    for (auto vol : vol_folders) 
    {
        pe_volt.push_back(Amplitude_calc(vol, data_size, "root"));
        Vbias_ver->SetPoint(i, pe_volt.back(), vol.Atof()); i++;
        std::cout << Form("(Voltge,Mean amplitude) = (%f,%f)", vol.Atof(), pe_volt.back()) << "\n" << std::endl;
    }

    cout << "\n\n-----> Voltage Breakdown fit" << endl;

    double VBD = fitBreakdownVoltage(Vbias_ver);


    /////////////////
    // Loop over all Voltages measured
    /////////////////

    cout << "\n\n-----> Noise analysis *** " << endl;

    TGraph * waveform = NULL;
    hfile->cd();

    //for (int i = vol_size-1; i < vol_size; i++) // Test on hihest and more noisy voltage
    for (int i = 0; i < vol_size; i++)
    {    
        const char * vol = vol_folders[i];
        cout << "\n\n-----> Voltage analyzed: " << vol << endl;

        map<string, int> color{{"clean",1},{"AP",1},{"CT",1},{"DCT",1}};

        // Counters 

        unsigned int events_cnt = 0;
        unsigned int counter_notclean = 0;
        unsigned int direct_xtalk_pulse_cnt = 0;
        unsigned int xtalk_pulse_cnt = 0;
        unsigned int after_pulse_cnt = 0;
        unsigned int nsaved = 0;
        
        // Define amplitude measured at which OV
        
        Double_t pe = pe_volt[i];
        V_meas = vol_folders[i].Atof() - VBD;
         
        map<string,TCanvas *> canv {
            {"CT",    new TCanvas(Form("Direct CrossTalk OV = %2.2f V",V_meas),Form("Direct CrossTalk OV = %2.2f V",V_meas),100,100,900,700)},
            {"DCT",   new TCanvas(Form("Delayed CrossTalk OV = %2.2f V",V_meas),Form("Delayed CrossTalk OV = %2.2f V",V_meas),100,100,900,700)},
            {"AP",    new TCanvas(Form("After Pulse OV = %2.2f V",V_meas),Form("After Pulse OV = %2.2f V",V_meas),100,100,900,700)},
            {"clean", new TCanvas(Form("Clean OV = %2.2f V",V_meas),Form("Clean OV = %2.2f V",V_meas),100,100,900,700)}
        };
        
        TGraph * Expfit_AP      = new TGraph();
        TGraph * cleanforfit    = NULL;
	    TH1D * AP_arrivaltime   = new TH1D("Histo_AP","AP arrival time", 140, 0, 0.2e-6);
        TH1D * DeXT_arrivaltime = new TH1D("Histo_DeXT","DeXT arrival time", 140, 0, 0.2e-6);
        TH1D * Npeaks           = new TH1D("Histo_Npeaks","Number of noise peaks when not clean", 50, 0, 50);
        TH1D * NpeaksCT         = new TH1D("Histo_NpeaksCT","Number of noise peaks when direct CT", 50, 0, 50);
        TH1D * NpeaksDCT        = new TH1D("Histo_NpeaksDCT","Number of noise peaks when delayed CT", 50, 0, 50);

        // Setup input tree
        TTree * tree = NULL;
        int npts;
        double times[10000];
        double amps[10000];
        if(globalArgs.input=="root") 
        {   
            tree = (TTree*)ifile->Get(vol);
            tree->SetBranchAddress("NsampPerEv",&npts);
            tree->SetBranchAddress("Amps",&amps);
            tree->SetBranchAddress("Times",&times);
            data_size = tree->GetEntries();
        }

        // Loop over every measurement on a folder
        for (int j = 0; j < data_size; j++)
        {
            events_cnt++;

            if(globalArgs.input=="root") 
            {
                tree->GetEntry(j);
                waveform = new TGraph(npts,times,amps);
            }
            //else  // Read from csv file
            //{
            //   TString datafilename = Form("%s%s/%i.csv",globalArgs.data_folder,vol,j);
            //    waveform = new TGraph(datafilename,"%lg %lg","/t;,");
            //}
            //if (waveform->IsZombie()) continue;

            waveform->SetName(Form("%s_%i",vol,j));
    
            Double_t * time  = waveform->GetX();
            Double_t * volts = waveform->GetY();

            /////////////////////////////////////////////////////
            // Data filtering into the different type of events
            // direct x-talk  AP   delayed x-talk
            /////////////////////////////////////////////////////
            unsigned int after_pulse = 0;
            unsigned int xtalk_pulse = 0;
            unsigned int direct_xtalk_pulse = 0;
            unsigned int time_of_pulse = 0;
            unsigned int sig_max = 0;
            unsigned int noise_peaks_cnt = 0;
            double sig_max_first = -1;
            double time_of_max_first = -1;
            double time_of_max_DeXT = 0;
            
            for (int row = 0; row < npts; row++) 
            {
                /////////////////////////////////////////////////////
                // direct x-talk
                if ( time[row] > 0 && volts[row] > thrs[i].direct_xtalk * pe ) // time larger 0ns
                    direct_xtalk_pulse++;  

                /////////////////////////////////////////////////////
                // after-pulse threshold
                if ( time[row]>thrs[i].reject_time*ns && volts[row] > thrs[i].after_pulse * pe ) // time larger 2ns and ap_th
                    after_pulse++;

                /////////////////////////////////////////////////////
                // delayed x-talk
                if ( time[row] > thrs[i].delxtalk_reject_time*ns && volts[row] > thrs[i].xtalk * pe ) // time larger 4ns and larger xtalk_th
                { 
                    xtalk_pulse++;
                    time_of_max_DeXT = time[row];
                }

                /////////////////////////////////////////////////////////////////////
                // Detect peaks in data after 4ns, count the number of maxima and
                // measure the time of arrival of first maxima, used later for AP exp fit
                /////////////////////////////////////////////////////////////////////

		        if (time[row] > thrs[i].reject_time*ns    // time larger 4ns
                        && volts[row] >= volts[row-1] && volts[row] >= volts[row-2] 
                        && volts[row] >= volts[row+1] && volts[row] >= volts[row+2] 
                        && volts[row] > thrs[i].after_pulse * pe) 
                {
                    noise_peaks_cnt++;

                    if (sig_max_first < 0)
                    {
				        sig_max_first     = volts[row];
				        time_of_max_first = time[row];
                    }
                }

            } // loop over time
                

            if (direct_xtalk_pulse > 0)  // Check for imm x-talk and plot
            {
                direct_xtalk_pulse_cnt++;
                counter_notclean++;
                sprintf(Category,"ImmCrosstalk");
                canv["CT"]->cd();

                TString graph_title = Form("Direct CrossTalk OV = %2.2f V",V_meas);
                drawWave(waveform, &color["CT"], graph_title, canv["CT"], 1.5*pe);
                
                Npeaks->Fill(noise_peaks_cnt);
                NpeaksCT->Fill(noise_peaks_cnt);
            }
            
            else if (xtalk_pulse > 0) // Only delayed x-talk
            {
                xtalk_pulse_cnt++;
                counter_notclean++;
                sprintf(Category,"DelCrosstalk");
                canv["DCT"]->cd();
                
                TString graph_title = Form("Delayed cross-talk OV = %2.2f V",V_meas);
                drawWave(waveform, &color["DCT"], graph_title, canv["DCT"], 1.5*pe);

                DeXT_arrivaltime->Fill(time_of_max_DeXT);
                Npeaks->Fill(noise_peaks_cnt);
                NpeaksDCT->Fill(noise_peaks_cnt);
            }
            
            else if (after_pulse > 0) //  Only after pulse
            {
                after_pulse_cnt++;
                counter_notclean++;
                sprintf(Category,"AfterPulse");

                TString graph_title = Form("After pulse OV = %2.2f V",V_meas);
                drawWave(waveform, &color["AP"], graph_title, canv["AP"], 1.5*pe);

                Expfit_AP->SetPoint(after_pulse_cnt-1,time_of_max_first,sig_max_first);
                AP_arrivaltime->Fill(time_of_max_first); 
                Npeaks->Fill(noise_peaks_cnt);
            }
	           
            else    // If not noisy then it's clean
            {
                sprintf(Category,"Clean");

                // Get only very clean waves and make an average for long tau fit.
                bool veryclean = true;
                for (int row = 0; row < npts; row++) 
                    if ((time[row] > 2) && (volts[row] > 0.5*pe)) 
                        veryclean = false;

                if(veryclean) cleanforfit = average(cleanforfit, waveform);

                if(nsaved < 20) // Max 20 clean graphs on the plot
                {
                    nsaved++;
                    drawWave(waveform, &color["clean"], Form("Clean pulse OV = %2.2f V",V_meas), canv["clean"], 1.5*pe);                  
		        }
	        }

            otree->Fill();
            
            delete time;
            delete volts;
        }
       
        cout << "----- Total number of events: " << events_cnt << endl;
        cout << Form("----- Not clean events: %i [%.2f%%]",counter_notclean, (float)(counter_notclean/events_cnt*100)) << endl;
        cout << Form("     (DirXtalk = %i [%.2f%%], DelXtalk = %i [%.2f%%], AP = %i [%.2f%%])",
            (int)direct_xtalk_pulse_cnt, (float)(direct_xtalk_pulse_cnt/events_cnt*100),
            (int)xtalk_pulse_cnt, (float)(xtalk_pulse_cnt/events_cnt*100),
            (int)after_pulse_cnt, (float)(after_pulse_cnt/events_cnt*100) ) << endl;

        cout << "-----> Long tau fit ***" << endl;
        double amp0, tau;
        fitLongTau(cleanforfit, &amp0, &tau, pe, vol, canv["clean"]);
   
        cout << "\n\n-----> After pulse fit ***" << endl;
        fitAPTau(Expfit_AP, amp0, tau, pe, vol, canv["AP"])->Write();

        // Final result: Correlated noise
        double perc_CT  = (float)direct_xtalk_pulse_cnt/events_cnt;
        double perc_DCT = (float)xtalk_pulse_cnt/events_cnt;
        double perc_AP  = (float)after_pulse_cnt/events_cnt;
        double perc_tot = perc_CT + perc_DCT + perc_AP;

        results["#CT"]->SetPoint(i,V_meas,perc_CT*100.);
        results["#CT"]->SetPointError(i,0.,TMath::Sqrt(perc_CT*(1.-perc_CT)/events_cnt)*100.);
        results["#AP"]->SetPoint(i,V_meas,perc_AP*100.);
        results["#AP"]->SetPointError(i,0.,TMath::Sqrt(perc_AP*(1.-perc_AP)/events_cnt)*100.);
        results["#DCT"]->SetPoint(i,V_meas,perc_DCT*100.);
        results["#DCT"]->SetPointError(i,0.,TMath::Sqrt(perc_DCT*(1.-perc_DCT)/events_cnt)*100.);
        results["#Total"]->SetPoint(i,V_meas,perc_tot*100.);
        results["#Total"]->SetPointError(i,0.,TMath::Sqrt(perc_tot*(1.-perc_tot)/events_cnt)*100.);
        results["#SecPeaks"]->SetPoint(i,V_meas,Npeaks->GetMean());
        results["#SecPeaks"]->SetPointError(i,0.,Npeaks->GetMeanError());
        results["#SecPeaksCT"]->SetPoint(i,V_meas,NpeaksCT->GetMean());
        results["#SecPeaksCT"]->SetPointError(i,0.,NpeaksCT->GetMeanError());
        results["#SecPeaksDCT"]->SetPoint(i,V_meas,NpeaksDCT->GetMean());
        results["#SecPeaksDCT"]->SetPointError(i,0.,NpeaksDCT->GetMeanError());
        
        // Save/print reults:
        
        for(auto const &e : canv) e.second->Write();
	    AP_arrivaltime->Write();
 	    DeXT_arrivaltime->Write();
        Npeaks->Write();
        NpeaksCT->Write();
        NpeaksDCT->Write();
        
        canv["CT"]->Print(globalArgs.res_folder+Form("Immcrosstalk_%s.pdf",vol));
        canv["DCT"]->Print(globalArgs.res_folder+Form("Delcrosstalk_%s.pdf",vol));
        
        if(globalArgs.save_all) 
        {
            TCanvas * tmpc = new TCanvas();
            AP_arrivaltime->Draw();
            tmpc->Print(globalArgs.res_folder+Form("AP_arrivetime_%s.pdf",vol));
            DeXT_arrivaltime->Draw();
            tmpc->Print(globalArgs.res_folder+Form("DeXT_arrivetime_%s.pdf",vol));
            Npeaks->Draw("HIST");
            tmpc->Print(globalArgs.res_folder+Form("Npeaks_%s.pdf",vol));
            NpeaksCT->Draw("HIST");
            tmpc->Print(globalArgs.res_folder+Form("Npeaks_whenCT_%s.pdf",vol));
            NpeaksDCT->Draw("HIST");
            tmpc->Print(globalArgs.res_folder+Form("Npeaks_whenDCT_%s.pdf",vol));
            delete tmpc;
        }

	    delete AP_arrivaltime;
	    delete DeXT_arrivaltime;
        delete Npeaks;
        delete NpeaksCT;
        delete NpeaksDCT;
        delete Expfit_AP;
        delete cleanforfit;
        for(auto const &e : canv) delete e.second;
    }

    // Save TTree with hist of noise for each event OV and the noise classification
    otree->Write();
    
    // Create final plot of total correlated noise
    TCanvas * cfinal = new TCanvas("Correlated Noise","Correlated Noise",100,100,900,700);
    
    Double_t tot_max_noise = TMath::MaxElement(results["#Total"]->GetN(),results["#Total"]->GetY());
    
    results["#Total"]->SetTitle("Correlated Noise");
    results["#Total"]->SetMarkerColor(kBlack);
    results["#Total"]->SetLineColor(kBlack);
    results["#Total"]->GetYaxis()->SetRangeUser(0,tot_max_noise+2);
    results["#Total"]->GetYaxis()->SetTitle("Noise [%]");
    results["#Total"]->GetXaxis()->SetTitle("OverVoltage [V]");
    results["#Total"]->SetFillColor(kBlack);
    results["#Total"]->SetFillStyle(3005);
    results["#Total"]->Draw("ALP*3");
    
    results["#CT"]->SetLineColor(kBlue);
    results["#AP"]->SetLineColor(kOrange+7);
    results["#DCT"]->SetLineColor(kGreen+2);
    results["#CT"]->SetFillColor(kBlue);
    results["#AP"]->SetFillColor(kOrange+7);
    results["#DCT"]->SetFillColor(kGreen+2);
    results["#CT"]->SetFillStyle(3005);
    results["#AP"]->SetFillStyle(3005);
    results["#DCT"]->SetFillStyle(3005);
    results["#CT"]->SetMarkerColor(kBlue);
    results["#AP"]->SetMarkerColor(kOrange+7);
    results["#DCT"]->SetMarkerColor(kGreen+2);
    results["#CT"]->Draw("LP*3");
    results["#AP"]->Draw("LP*3");
    results["#DCT"]->Draw("LP*3");
    
    TLegend * leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#Total"],"Total","f");
    leg->AddEntry(results["#CT"],"Direct Cross-Talk","f");
    leg->AddEntry(results["#AP"],"After Pulse","f");
    leg->AddEntry(results["#DCT"],"Delayed Cross-Talk","f");
    leg->Draw();
    
    cfinal->SetGrid();
    cfinal->Print(globalArgs.res_folder+"CorrelatedNoise.pdf");
    cfinal->Write();

    Double_t tot_max_peaks = TMath::MaxElement(results["#SecPeaks"]->GetN(),results["#SecPeaks"]->GetY());
    results["#SecPeaks"]->SetTitle("Secondary peaks after a noise peak");
    results["#SecPeaks"]->SetMarkerColor(kBlack);
    results["#SecPeaks"]->SetLineColor(kBlack);
    results["#SecPeaks"]->GetYaxis()->SetRangeUser(0,tot_max_peaks*2.);
    results["#SecPeaks"]->GetYaxis()->SetTitle("<N Secondary peaks>");
    results["#SecPeaks"]->GetXaxis()->SetTitle("OverVoltage [V]");
    results["#SecPeaks"]->Draw("ALP*");

    results["#SecPeaksCT"]->SetLineColor(kBlue);
    results["#SecPeaksDCT"]->SetLineColor(kGreen+2);
    results["#SecPeaksCT"]->SetMarkerColor(kBlue);
    results["#SecPeaksDCT"]->SetMarkerColor(kGreen+2);
    results["#SecPeaksCT"]->Draw("LP*");
    results["#SecPeaksDCT"]->Draw("LP*");

    leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(results["#SecPeaks"],"After any noise","lp");
    leg->AddEntry(results["#SecPeaksCT"],"After Direct Cross-Talk","lp");
    leg->AddEntry(results["#SecPeaksDCT"],"After Delayed Cross-Talk","lp");
    leg->Draw();
    
    cfinal->SetGrid();
    cfinal->Print(globalArgs.res_folder+"SecondaryPeaks.pdf");
  
    return 0;
}


