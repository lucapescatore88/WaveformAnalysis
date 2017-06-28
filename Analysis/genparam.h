#ifndef GENPARAM_H
#define GENPARAM_H

#include "TString.h"
#include "TGraph.h"
#include "TH1.h"
#include <iostream>

using namespace std;

const double ns = 1e-8;

struct globalArgs_t
{
    const char* data_folder;                /* -d option */
    const char* arg_pathToSetupFile;        /* -S option */
    TString res_folder;                     /* -o option */
    bool save_all;                           /* -a  */
    TString input;                          /* -I root  */             

} globalArgs {" "," "," ",false,"root"};

TH1 * convertGrToH(TGraph * gr)
{
    Double_t * xx  = gr->GetX();
    Double_t * yy = gr->GetY();
    double halfDT = (xx[1] - xx[0])/2.; 
    TH1 * h = new TH1F(TString(gr->GetName())+"H",gr->GetTitle(),
        gr->GetN(),xx[0]-halfDT,xx[gr->GetN()-1]+halfDT);
    for(int i = 1; i < gr->GetN(); i++ ) h->SetBinContent(i,yy[i]);

    return h;
}

TGraph * average(TGraph * gr, TGraph * gr2)
{
    static int nevts = 0;

    int npts    = gr2->GetN();
    double * x2 = gr2->GetX();
    double * y2 = gr2->GetY();

    if(!gr)
    {
        nevts = 1;
        TGraph * grout = new TGraph(npts,x2,y2);
        grout->SetName(TString(gr2->GetName())+"ForFit");
        return grout;
    }

    double avgx[npts];
    double avgy[npts];
    Double_t * time  = gr->GetX();
    Double_t * volts = gr->GetY();
    for (int row = 0; row < npts; row++) 
    {
        avgx[row] = (nevts*time[row] + x2[row])/(nevts+1);
        avgy[row] = (nevts*volts[row] + y2[row])/(nevts+1);
    }
    nevts++;

    TGraph * grout = new TGraph(npts,avgx,avgy);
    grout->SetName(gr->GetName());
    return gr;
}



#endif