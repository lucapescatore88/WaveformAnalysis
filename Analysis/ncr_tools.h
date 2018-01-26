#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <map>
#include <dirent.h>

#include <TROOT.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TString.h>
#include <TFile.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>

using namespace std;

void lhcbstyle();

vector<TString> files;
TString pulse_shape_analysis_file = "";
map<double, double> correlated_noise;	// dV, tot_corr_noise
map<TString, TString> correlated_noise_str;	// dV, tot_corr_noise
map< TString, TGraphErrors* > gNCRvsDCR_corr_noise;	// dV, NCR_vs_DCR_plot

void loadFiles() {
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir ("./")) != NULL) {
	  /* print all the files and directories within directory */
	  while ((ent = readdir (dir)) != NULL) {
		  TString ele = TString(ent->d_name);
		  if(ele.Contains("_pulse_shape")) pulse_shape_analysis_file = ele + "/noiseanalysis.root";
		  if(ele.Contains("_ro_")) files.push_back(ele+"/");	// no sorting
	  }
	  closedir (dir);
	} else {
	  /* could not open directory */
	  perror ("");
	  return;
	}
	cout << "Found pulse analysis file: " << pulse_shape_analysis_file << endl;
	cout << "Found ro files: " << files.size() << endl;
}

map< double, int > colours {
	{1.0, 1},
	{1.5, kCyan-3},
	{2.0, kTeal+2},
	{2.5, kBlue},
	{3.0, kViolet},
	{3.5, kRed},
	{4.0, kOrange+7},
	{4.5, kOrange-8},
	{5.0, kOrange+3},
	{5.5, 1}
};

map< double, int > style {
	{1.0, 20},
	{1.5, 23},
	{2.0, 20},
	{2.5, 21},
	{3.0, 22},
	{3.5, 34},
	{4.0, 29},
	{4.5, 33},
	{5.0, 39},
	{5.5, 20},
};

map<TString, unsigned int> getCorrelatedNoise(map<double, double>& m, map<TString, TString>& m_str, map< TString, TGraphErrors* >& mg, TString filename, TString plotName="tot_prim_corr_noise") {
	m.clear();
	m_str.clear();
	mg.clear();
	map<TString, unsigned int> counter;
	counter.clear();
	
	TFile* f = new TFile(filename);
	TCanvas * c = (TCanvas*) f->Get("Correlated noise");
	TGraphErrors * g = (TGraphErrors*) c->FindObject(plotName);
	if(!g) {
		cout << "Problem getting graph " << plotName << endl;
		return counter;
	}
	
	for(int pt(0); pt<g->GetN(); ++pt) {
		double x,y;
		g->GetPoint(pt,x,y);
		m.insert ( std::pair<double,double>(x,y) );
		m_str.insert ( std::pair<TString,TString>(Form("%.2lf",x),Form("%.1lf",y)) );
		
		TGraphErrors * gNCR_vs_DCR = new TGraphErrors();
		TString name = Form("dV_%.2lf_corr_%.1lf",x,y);
		gNCR_vs_DCR->SetName(name);
		double rounddV = (floor((x*2)+0.5)/2);
		gNCR_vs_DCR->SetLineColor(colours[rounddV]);
		gNCR_vs_DCR->SetMarkerColor(colours[rounddV]);
		gNCR_vs_DCR->SetMarkerStyle(style[rounddV]);
		gNCR_vs_DCR->SetMarkerSize(1.5);
		mg.insert ( std::pair<TString,TGraphErrors*>(Form("%.2lf",x),gNCR_vs_DCR) );
		
		counter.insert ( std::pair<TString,unsigned int>(Form("%.2lf",x),0) );
	}
	f->Close();
	return counter;
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
    //cout << Form("Found DCR = %2.0lf kHz for Vbias = %2.2lf V", dcr, bias) << endl;
    file.close();
    
    return dcr;
}

int NCRvsDCR(double seedThrs = 2.5, double time_window = 100) {
	
	lhcbstyle();
	
	loadFiles();
	
	// Get correlated noise from pulse shape file
	map<TString, unsigned int> counter = getCorrelatedNoise(correlated_noise, correlated_noise_str, gNCRvsDCR_corr_noise, pulse_shape_analysis_file, "tot_prim_corr_noise");
	//~ cout << correlated_noise.size() << " " << correlated_noise_str.size() << " " << gNCRvsDCR_corr_noise.size() << endl;
	//~ for(auto const& g: gNCRvsDCR_corr_noise) cout << g.second->GetName() << endl;
	
	// Init NCR vs DCR plot
	TLegend * leg = new TLegend(0.55,0.18,1.0,0.58);
	leg->SetHeader("#DeltaV  /  total correlated noise");
	leg->SetFillColor(0);
	
	double minX(1e9), minY(1e9), maxX(0), maxY(0);
	
	for(unsigned int j(0); j<files.size(); ++j) {
		
		TString dcr_filename = files[j]+"DCR_for_RandomOverlap.txt";
		
		// Get NCR vs seed plot
		TString filename = files[j]+"randomoverlap_analysis.root";
		cout << filename << endl;
		TFile * f = new TFile(filename);
		// Loop over the directories (Vbias)
		TList * l = f->GetListOfKeys();
		TIter next(l);
		TObject* object = 0;
		while ((object = next())) {	// loop over Vbias
			//cout << "HERE: " << object->GetName() << endl;
			TString dirname = TString(object->GetName())+"/Dir_NCR_vs_seed";
			double Vbias;
			if(sscanf(object->GetName(), "charge_analysis_%lfV", &Vbias)!=1) {
				cout << "problem scanning the file" << endl;
				return 0;
			}
			f->cd(dirname);
			TString graphname = "";
			//gDirectory->GetListOfKeys()->Print();
			TIter next(gDirectory->GetListOfKeys());
			TObject* obj = 0;
			double dV(0), time(0);
			while ((obj = next())) {
				//cout << "HERE2: " << obj->GetName() << endl;
				if(sscanf(obj->GetName(), "NCR_vs_seed_%lfns_%lfV", &time, &dV)==2) {
					if(time == time_window) graphname = TString(obj->GetName());
				}
			}
			delete obj;
			if(graphname == "") {
				cout << "problem getting the ncr_vs_seed plot" << endl;
				return 0;
			}
			TGraphErrors * ncr_vs_seed = (TGraphErrors*) gDirectory->Get(graphname);
			
			// Get NCR from plot
			double ncr_value(0.0), ncr_error(0.0);
			const int Npoints = ncr_vs_seed->GetN();
			double ax[Npoints],ay[Npoints];
			double ax_err[Npoints],ay_err[Npoints];
			for(int pt=0; pt<Npoints; pt++) {
				ncr_vs_seed->GetPoint(pt,ax[pt],ay[pt]);
				ax_err[pt] = ncr_vs_seed->GetErrorX(pt);
				ay_err[pt] = ncr_vs_seed->GetErrorY(pt);
				if(ax[pt] == seedThrs) {
					ncr_value = ay[pt];
					ncr_error = ay_err[pt];
				}
				//~ cout << ax[pt] << " " << ay[pt] << " " << ax_err[pt] << " " << ay_err[pt] << endl;
			}
			
			// Get DCR from txt file
			double dcr_value = getDCR_from_CSVfile(Vbias, dcr_filename.Data());
			
			cout << "Vbias = " << Vbias << "V, dV = " << dV << "V, DCR = " << Form("%.2lf",dcr_value/1000.0) << "MHz, NCR = " << Form("%.2lf +/- %.2lf",ncr_value,ncr_error) << "MHz" << endl;
			
			
			
			// Builds NCR vs DCR plot
			TString dVstr = Form("%.2lf",dV);
			if(ncr_value!=0 && dcr_value!=0) {
				gNCRvsDCR_corr_noise[dVstr]->SetPoint(counter[dVstr], dcr_value/1000.0, ncr_value);
				gNCRvsDCR_corr_noise[dVstr]->SetPointError(counter[dVstr], 0, ncr_error);
				++counter[dVstr];
				if(dcr_value/1000.0 < minX) minX = dcr_value/1000.0;
				if(ncr_value < minY) minY = ncr_value;
				if(dcr_value/1000.0 > maxX) maxX = dcr_value/1000.0;
				if(ncr_value > maxY) maxY = ncr_value;
			}
			
			
			f->cd();
		}
		delete object;
		f->Close();
		delete f;
	}
	
	TCanvas * c1 = new TCanvas("NCR_vs_DCR","NCR vs DCR",100,100,900,700);
	double zoom = 0.5;
	TH2D * frame = new TH2D("frame","frame",1000,zoom*minX,(1./zoom)*maxX,1000,zoom*minY,(1./zoom)*maxY);
	TString tit1 = Form("NCR, seed = %.1lf PE, #tau_{int} = %.0lf ns", seedThrs, time_window);
	frame->SetTitle(tit1);
	frame->Draw();
	frame->GetXaxis()->SetTitle("DCR per channel [MHz]");
	frame->GetYaxis()->SetTitle("NCR per array [MHz]");
	
	gStyle->SetGridColor(17);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
	for(auto const& g: gNCRvsDCR_corr_noise) {
		TString ent = g.first + "V  /  " + correlated_noise_str[g.first] + "%";
		g.second->Draw("P+");
		leg->AddEntry(g.second,ent,"P");
		
		// Fit
		//~ TString fitname = TString(g.second->GetName()) + "_fit";
		//~ TF1 * fit = new TF1(fitname,"([0]^2)*x**([1]^2)",minX,maxX);
		//~ double x1, x2, y1, y2;
		//~ g.second->GetPoint(0,x1,y1);
		//~ g.second->GetPoint(g.second->GetN()-1,x2,y2);
		//~ double A, B;
		//~ B = TMath::Log(y1/y2)/TMath::Log(x1/x2);
		//~ A = y2/(TMath::Power(x2,B));
		//~ fit->SetParameter(0, TMath::Sqrt(A));
		//~ fit->SetParameter(1, TMath::Sqrt(B));
		//~ fit->SetLineColor(g.second->GetMarkerColor());
		//~ g.second->Fit(fitname,"QR+");
	}
	
	gPad->SetGrid(1);
    gPad->SetLogy();
    gPad->SetLogx();
    leg->Draw();
    
    return 0;
}




// all users - please change the name of this file to lhcbStyle.C
// Commits to lhcbdocs svn of .C files are not allowed
void lhcbstyle(){
 
  // define names for colours
  Int_t black  = 1;
  Int_t red    = 2;
  Int_t green  = 3;
  Int_t blue   = 4;
  Int_t yellow = 5;
  Int_t magenta= 6;
  Int_t cyan   = 7;
  Int_t purple = 9;
 
 
////////////////////////////////////////////////////////////////////
// PURPOSE:
//
// This macro defines a standard style for (black-and-white)
// "publication quality" LHCb ROOT plots.
//
// USAGE:
//
// Include the lines
//   gROOT->ProcessLine(".L lhcbstyle.C");
//   lhcbStyle();
// at the beginning of your root macro.
//
// Example usage is given in myPlot.C
//
// COMMENTS:
//
// Font:
//
// The font is chosen to be 132, this is Times New Roman (like the text of
//  your document) with precision 2.
//
// "Landscape histograms":
//
// The style here is designed for more or less square plots.
// For longer histograms, or canvas with many pads, adjustements are needed.
// For instance, for a canvas with 1x5 histograms:
//  TCanvas* c1 = new TCanvas("c1", "L0 muons", 600, 800);
//  c1->Divide(1,5);
//  Adaptions like the following will be needed:
//  gStyle->SetTickLength(0.05,"x");
//  gStyle->SetTickLength(0.01,"y");
//  gStyle->SetLabelSize(0.15,"x");
//  gStyle->SetLabelSize(0.1,"y");
//  gStyle->SetStatW(0.15);
//  gStyle->SetStatH(0.5);
//
// Authors: Thomas Schietinger, Andrew Powell, Chris Parkes, Niels Tuning
// Maintained by Editorial board member (currently Niels)
///////////////////////////////////////////////////////////////////
 
  // Use times new roman, precision 2
  Int_t lhcbFont        = 132;  // Old LHCb style: 62;
  // Line thickness
  Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
  // Text size
  Double_t lhcbTSize    = 0.06;
 
  // use plain black on white colors
  gROOT->SetStyle("Plain");
  TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");
 
  //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X
 
  lhcbStyle->SetFillColor(1);
  lhcbStyle->SetFillStyle(1001);   // solid
  lhcbStyle->SetFrameFillColor(0);
  lhcbStyle->SetFrameBorderMode(0);
  lhcbStyle->SetPadBorderMode(0);
  lhcbStyle->SetPadColor(0);
  lhcbStyle->SetCanvasBorderMode(0);
  lhcbStyle->SetCanvasColor(0);
  lhcbStyle->SetStatColor(0);
  lhcbStyle->SetLegendBorderSize(0);
  lhcbStyle->SetLegendFont(132);
 
  // If you want the usual gradient palette (blue -> red)
  lhcbStyle->SetPalette(1);
  // If you want colors that correspond to gray scale in black and white:
  //int colors[8] = {0,5,7,3,6,2,4,1};
  //lhcbStyle->SetPalette(8,colors);
 
  // set the paper & margin sizes
  lhcbStyle->SetPaperSize(20,26);
  lhcbStyle->SetPadTopMargin(0.1);
  lhcbStyle->SetPadRightMargin(0.1); // adjusted for TH2 colz
  lhcbStyle->SetPadBottomMargin(0.16);
  lhcbStyle->SetPadLeftMargin(0.14);
 
  // use large fonts
  lhcbStyle->SetTextFont(lhcbFont);
  lhcbStyle->SetTextSize(lhcbTSize);
  lhcbStyle->SetLabelFont(lhcbFont,"x");
  lhcbStyle->SetLabelFont(lhcbFont,"y");
  lhcbStyle->SetLabelFont(lhcbFont,"z");
  lhcbStyle->SetLabelSize(lhcbTSize,"x");
  lhcbStyle->SetLabelSize(lhcbTSize,"y");
  lhcbStyle->SetLabelSize(lhcbTSize,"z");
  lhcbStyle->SetTitleFont(lhcbFont);
  lhcbStyle->SetTitleFont(lhcbFont,"x");
  lhcbStyle->SetTitleFont(lhcbFont,"y");
  lhcbStyle->SetTitleFont(lhcbFont,"z");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");
 
  // use medium bold lines and thick markers
  lhcbStyle->SetLineWidth(lhcbWidth);
  lhcbStyle->SetFrameLineWidth(lhcbWidth);
  lhcbStyle->SetHistLineWidth(lhcbWidth);
  lhcbStyle->SetFuncWidth(lhcbWidth);
  lhcbStyle->SetGridWidth(lhcbWidth);
  lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  lhcbStyle->SetMarkerStyle(20);
  lhcbStyle->SetMarkerSize(1.0);
 
  // label offsets
  lhcbStyle->SetLabelOffset(0.012,"X");
  lhcbStyle->SetLabelOffset(0.010,"Y");
 
  // by default, do not display histogram decorations:
  lhcbStyle->SetOptStat(0);  
  //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
  lhcbStyle->SetOptTitle(1);
  lhcbStyle->SetOptFit(0);
  //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  lhcbStyle->SetTitleOffset(0.95,"X");
  lhcbStyle->SetTitleOffset(0.95,"Y");
  lhcbStyle->SetTitleOffset(1.2,"Z");
  lhcbStyle->SetTitleFillColor(0);
  lhcbStyle->SetTitleStyle(0);
  lhcbStyle->SetTitleBorderSize(0);
  lhcbStyle->SetTitleFont(lhcbFont,"title");
  lhcbStyle->SetTitleX(0.0);
  lhcbStyle->SetTitleY(1.0);
  lhcbStyle->SetTitleW(1.0);
  lhcbStyle->SetTitleH(0.05);
 
  // look of the statistics box:
  lhcbStyle->SetStatBorderSize(0);
  lhcbStyle->SetStatFont(lhcbFont);
  lhcbStyle->SetStatFontSize(0.05);
  lhcbStyle->SetStatX(0.9);
  lhcbStyle->SetStatY(0.9);
  lhcbStyle->SetStatW(0.25);
  lhcbStyle->SetStatH(0.15);
 
  // put tick marks on top and RHS of plots
  lhcbStyle->SetPadTickX(1);
  lhcbStyle->SetPadTickY(1);
 
  // histogram divisions: only 5 in x to avoid label overlaps
  lhcbStyle->SetNdivisions(505,"x");
  lhcbStyle->SetNdivisions(510,"y");
 
  gROOT->SetStyle("lhcbStyle");
  gROOT->ForceStyle();
 
  // add LHCb label
  TPaveText *lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                           0.87 - gStyle->GetPadTopMargin(),
                           gStyle->GetPadLeftMargin() + 0.20,
                           0.95 - gStyle->GetPadTopMargin(),
                           "BRNDC");
  lhcbName->AddText("LHCb");
  lhcbName->SetFillColor(0);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);
 
  TText *lhcbLabel = new TText();
  lhcbLabel->SetTextFont(lhcbFont);
  lhcbLabel->SetTextColor(1);
  lhcbLabel->SetTextSize(lhcbTSize);
  lhcbLabel->SetTextAlign(12);
 
  TLatex *lhcbLatex = new TLatex();
  lhcbLatex->SetTextFont(lhcbFont);
  lhcbLatex->SetTextColor(1);
  lhcbLatex->SetTextSize(lhcbTSize);
  lhcbLatex->SetTextAlign(12);
 
 // cout << "-------------------------" << endl;  
 // cout << "Set LHCb Style - Feb 2012" << endl;
 // cout << "-------------------------" << endl;  
 
}
