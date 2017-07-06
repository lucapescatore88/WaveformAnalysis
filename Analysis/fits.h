#ifndef BREAKDOWN_VOLTAGE_H
#define BREAKDOWN_VOLTAGE_H

#include <TCanvas.h>
#include <iostream>
#include <TFitResultPtr.h>
#include <TPaveText.h>

#include "genparam.h"

double fitBreakdownVoltage(TGraph * Vbias_ver)
{
    TCanvas * ca = new TCanvas("Voltage Breakdown calculation","Voltage Breakdown calculation",100,100,900,700);
    Vbias_ver->SetTitle("Voltage Breakdown calculation");
    Vbias_ver->GetYaxis()->SetTitle("Bias Volatge [V]");
    //Vbias_ver->GetYaxis()->SetTitleOffset(1.2);
    Vbias_ver->GetXaxis()->SetTitle("Mean peak amplitude [V]");
    Vbias_ver->Draw("AP*");
    ca->SetGrid();
    
    TFitResultPtr fit = Vbias_ver->Fit("pol1","S");
    Double_t VBD = fit->Value(0);
    
    std::cout << "Breakdown voltage is: " << VBD << std::endl;

    TPaveText * pv = new TPaveText(0.2,0.65,0.35,0.74,"brNDC");
    pv->AddText(Form("V_{BD} = %2.2f",VBD));
    pv->SetFillColor(kWhite);
    pv->Draw();
    ca->Print(globalArgs.res_folder+"Breakdown_voltage.pdf");

    return VBD;
}

TF1 * fitLongTau(TGraph * cleanwaves, double * amp0, double * tau, double pe, const char * vol, TCanvas * c) 
{
    TH1 * waveh = convertGrToH(cleanwaves);

    TCanvas * ctmp = new TCanvas(Form("LongTauFit_%s",vol));

    // Fit parameters and limits to calculate slow component of the pulse
    TF1 * exp_longtau = new TF1(Form("exptau_%s",vol),"[0]*exp(-x/[1])",0.,180.*ns);
    exp_longtau->SetParameter(0,pe*0.2);
    exp_longtau->SetParLimits(0,0.01*pe,1.*pe);
    exp_longtau->SetParameter(1,  80*ns);
    exp_longtau->SetParLimits(1,10*ns,200*ns); 
        
    waveh->Fit(Form("exptau_%s",vol),"","",4*ns,180.*ns); // Fit boundaries for the slow component of the pulse
    (*amp0) = exp_longtau->GetParameter(0);
    (*tau) = exp_longtau->GetParameter(1);
    std::cout << "Long tau = " << *tau << std::endl;
    ctmp->Write();
    delete ctmp;

    c->cd();

    TPaveText * pv = new TPaveText(0.2,0.65,0.35,0.74,"brNDC");
    pv->AddText(Form("#tau_{long} = %2.2e ",*tau));
    pv->SetFillColor(kWhite);
    pv->Draw();
    exp_longtau->Draw("SAME");

    c->Print(globalArgs.res_folder+Form("Clean_%s.pdf",vol));

    return exp_longtau;
}

TF1 * fitAPTau(TGraph * APtime, double amp0, double tau, double pe, const char * vol, TCanvas * c) 
{
    // Fit parameters and limits to calculate AP recharge

    TCanvas * ctmp = new TCanvas(Form("APTauFit_%s",vol));
    APtime->Draw("AP*");
    TF1 * exp = new TF1(Form("exp_%s",vol),"(([0])/(exp(-4E-9/[1])))*(exp(-4E-9/[1])-exp(-x/[1]))+[2]*exp(-x/[3])",4*ns,180 * ns);
    exp->SetParameter(0,pe);
    //exp->SetParLimits(0,0.2*pe,10*pe);    
    exp->FixParameter(0,pe);  
    exp->SetParameter(1,30*ns); // 30
    exp->SetParLimits(1,10*ns,500*ns); // 4 500
    exp->FixParameter(2,amp0);
    exp->FixParameter(3,tau);
    APtime->Fit(Form("exp_%s",vol));
    exp->Draw("same");
    //ctmp->Print(Form("APFit_%s.pdf",vol));
    ctmp->Write();
    delete ctmp;
    
    c->cd();
    
    TPaveText * pv = new TPaveText(0.2,0.65,0.35,0.74,"brNDC");
    pv->AddText(Form("Recovery time = %2.2e ",exp->GetParameter(1)));
    pv->SetFillColor(kWhite);
    pv->Draw();
    exp->Draw("SAME");

    c->Print(globalArgs.res_folder+Form("Afterpulse_%s.pdf",vol));

    return exp;
}



#endif