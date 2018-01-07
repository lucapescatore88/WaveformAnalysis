#ifndef clustering_simulation_h
#define clustering_simulation_h

#include "genparam.h"
#include "fits.h"

#include <vector>
#include <map>

#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
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
#include <TRandom3.h>

TGraph * getCumulativeDistribution(TH1D * distr) {
    distr->Scale(1.0/distr->Integral());
    TH1 * cumu = distr->GetCumulative();
    unsigned int Nx = cumu->GetNbinsX();
	TGraph * inv_cumu = new TGraph(Nx);
	for(unsigned int binx(0); binx<Nx; ++binx) {
        double bin_center = cumu->GetBinCenter(binx+1);
        double bin_content = cumu->GetBinContent(binx+1);
        inv_cumu->SetPoint(binx, bin_content, bin_center);
    }
    return inv_cumu;
}

void generateArray(TGraph * cumu, double (&event_array)[128], double yscale=1) {
    
    TRandom3* rnd_gen = new TRandom3(0);
    double data[128];
    rnd_gen->RndmArray(128, data);
    for (int channel=0; channel<128; channel++){
        event_array[channel] = cumu->Eval(data[channel]) / yscale;
    }
    
    return;
}

// Search Cluster Algorithm (3 Thresholds: Seed, Neighb, Sum)
int Cluster_Search(double event_array[128], double threshold[3], vector<double>& mean_position_vec, vector<double>& sum_vec, vector<double>& clusterSize_vec){
	
    mean_position_vec.clear();
    sum_vec.clear();
    clusterSize_vec.clear();
	
    bool    EndOfCluster = false;
    int     clusterSize = 0;
    int     clusterSizeMAX = 15;
    int     nbCluster = 0;
    int     first_ch_in_cluster = 0;
    double  sum = 0;
    
    double   SeedCutPE = threshold[0];
    double   NeighborCutPE = threshold[1];
    double   SumCutPE = threshold[2];

    double mean_position = 0.0;

    double mean_position_num = 0.0;
    double mean_position_den = 0.0;
    
    double event_array_copy[128];
    
    for (int channel=0; channel<128; channel++){
        event_array_copy[channel] = event_array[channel];
    }
        

    // cluster search, start to find channels passing the seed cut
    for (int channel=0; channel<128; channel++) {

        if (event_array_copy[channel] >= SeedCutPE){                                             // Search for seed threshold
            ++clusterSize;
        }
        else if ((clusterSize >= 1 )){                                                  // no more seeds in cluster
            EndOfCluster=true;
            first_ch_in_cluster = channel - clusterSize;
        }

        if ((clusterSize == 1) && (channel == 128 - 1)){                     // last channel reached (Check for single channel Cluster)
            if(event_array_copy[channel] >= SeedCutPE) {
                first_ch_in_cluster = channel;
            }
            EndOfCluster=true;
        }

        if ((clusterSize > 1 ) && (channel == 128 - 1)){                     // Looking for Cluster with the last channel included
            EndOfCluster=true;
            if(event_array_copy[channel] >= SeedCutPE) {
                first_ch_in_cluster = channel - clusterSize + 1;
            }
        }

        if ((clusterSize == clusterSizeMAX - 2) && (!EndOfCluster))  {                  // maximum size reached
            ++channel;                                                                  // go to next channel to search for neighbour threshold
            first_ch_in_cluster = channel - clusterSize;                                // NOTE: Cluster with the size greater than clusterSizeMAX will be breaked
            EndOfCluster=true;
        }

        // cluster finished so we need to add the neighbors
        if (EndOfCluster){

            if(first_ch_in_cluster > 0){                                                // check left neighbour
                if (event_array_copy[first_ch_in_cluster - 1] >= NeighborCutPE){
                    ++clusterSize;
                    first_ch_in_cluster = first_ch_in_cluster - 1;
                }
            }
            if(channel == 128 - 1){  
                if(channel == first_ch_in_cluster+clusterSize+1) {
                    if (event_array_copy[channel] >= NeighborCutPE){
                        ++clusterSize;
                    }
                }
            } else if(channel < 128 - 1){                                           // check right neighbour
                if (event_array_copy[channel] >= NeighborCutPE){
                    ++clusterSize;
                }
            }
            
            // calculate and check the sum threshold
            for (int cluster_ch = 0; cluster_ch < clusterSize; cluster_ch++){
                sum += event_array_copy[first_ch_in_cluster + cluster_ch];
            }

            if(sum >= SumCutPE){

                nbCluster++;
                clusterSize_vec.push_back(clusterSize);
                sum_vec.push_back(sum);
                
                // Weighted Mean calculation
                for(int itt=0; itt < clusterSize; itt++) {
                    mean_position_num += event_array_copy[first_ch_in_cluster+itt]*(double)itt;
                    mean_position_den += event_array_copy[first_ch_in_cluster+itt];
                    event_array_copy[first_ch_in_cluster+itt] = 0;   // Clean values already evaluated to avoid duplicate channel in two clusters
                }
                
                mean_position = mean_position_num / mean_position_den;
                mean_position_vec.push_back(first_ch_in_cluster + mean_position);
            }
            
            EndOfCluster=false;
            clusterSize=0;
            mean_position = 0.0;
            mean_position_num = 0.0;
            mean_position_den = 0.0;
            sum = 0.0;
            first_ch_in_cluster = 0;
            
        }
    }
    
    return nbCluster;
}

TCanvas * plotNCR(vector<TGraphErrors*> gs, TString option="seed") {
    if(!gs.size()) return NULL;
    double dV(0), para(0);
    TString title, xtitle, header;
    string pattern(""), leg_pat;
    int Ncol;
    if(option == "seed") { pattern = "NCR_vs_"+string(option)+"_%lfns_%lfV"; leg_pat="%2.0lf"; title="NCR vs seed threshold"; xtitle="Seed threshold [PE]"; header="Integration time [ns]"; Ncol=3; }
    if(option == "time_window") { pattern = "NCR_vs_"+string(option)+"_%lfPE_%lfV"; leg_pat="%2.1lf"; title="NCR vs integration time"; xtitle="#tau_{int} [ns]"; header="Seed threshold [PE]"; Ncol=2; }
    
    double tmp1, tmp2;
    if(sscanf(gs[0]->GetName(), pattern.c_str(), &tmp1, &tmp2) == 2) {
        para = tmp1;
        dV = tmp2;
    }
    
    double xmin(1e6), xmax(0), ymin(1e6), ymax(0);
    for(unsigned int i(0); i<gs.size(); ++i) {
        double x = TMath::MinElement(gs[i]->GetN(),gs[i]->GetX());
        double X = TMath::MaxElement(gs[i]->GetN(),gs[i]->GetX());
        double y = TMath::MinElement(gs[i]->GetN(),gs[i]->GetY());
        double Y = TMath::MaxElement(gs[i]->GetN(),gs[i]->GetY());
        if(x<xmin) xmin = x;
        if(X>xmax) xmax = X;
        if(y<ymin) ymin = y;
        if(Y>ymax) ymax = Y;
    }
    if(ymin==0) ymin = 0.01;
    TCanvas * cNCR = new TCanvas(Form("NCR_vs_"+option+"_%2.2lfV",dV), Form("NCR at #DeltaV=%2.2lfV",dV), 900, 700);
    TLegend * leg;
    if(option == "seed") leg = new TLegend(0.50,0.60,0.89,0.89);
    if(option == "time_window") leg = new TLegend(0.16,0.60,0.4,0.89);
    leg->SetNColumns(Ncol);
    int colour = 1;
    TH2D * frame = new TH2D("frame","frame",200,xmin-0.5,xmax+0.5, 200,0.5*ymin,1.5*ymax);
    frame->GetYaxis()->SetTitle("NCR [MHz]");
    frame->SetTitle(title);
    frame->GetXaxis()->SetTitle(xtitle);
    leg->SetHeader(header);
    cNCR->cd();
    frame->Draw();
    for(unsigned int i(0); i<gs.size(); ++i) {
        if(sscanf(gs[i]->GetName(), pattern.c_str(), &tmp1, &tmp2) == 2) {
            para = tmp1;
            dV = tmp2;
        }
        if(gs.size()>=10) colour+=3;
        else colour+=8;
        int col = (colour%50)+50;
        formatGr(gs[i], col, 3003, xtitle, "NCR [MHz]");
        leg->AddEntry(gs[i], Form(leg_pat.c_str(),para), "f");
        cNCR->cd();
        gs[i]->Draw("PL+3");
    }
    leg->Draw();
    cNCR->SetLogy();
    cNCR->SetGrid();
    gStyle->SetGridColor(17);
    return cNCR;
}












#endif
