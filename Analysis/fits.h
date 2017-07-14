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
    TF1 * exp_longtau = new TF1(Form("exptau_%s",vol),"[0]*exp(-(x-[1])/[2])",20*ns,120*ns);
    exp_longtau->SetParameter(0, pe*0.02);
    exp_longtau->SetParLimits(0, 0.005*pe,0.2*pe);
    //exp_longtau->SetParameter(1, 2*ns);
    //exp_longtau->SetParLimits(1, 0*ns,5*ns);
    exp_longtau->FixParameter(1, 0*ns);
    exp_longtau->SetParameter(2, 50*ns);
    exp_longtau->SetParLimits(2, 30*ns,70*ns);
    //exp_longtau->SetLineColor(2);
        
    waveh->Fit(Form("exptau_%s",vol),"","",5*ns,140*ns); // Fit boundaries for the slow component of the pulse
    (*amp0) = exp_longtau->GetParameter(0);
    (*tau) = exp_longtau->GetParameter(2);
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

TF1 * fitLongTau(TMultiGraph * cleanwaves, double * amp0, double * tau, double pe, const char * vol, TCanvas * c, TGraph * avg) 
{
    TCanvas * ctmp = new TCanvas(Form("LongTauFit_%s",vol));

    avg->Draw("AP");
    //cleanwaves->Draw("P SAME");

    // Fit parameters and limits to calculate slow component of the pulse
    TF1 * exp_longtau = new TF1(Form("exptau_%s",vol),"[0]*exp(-(x-[1])/[2])",15*ns,120*ns);
    exp_longtau->SetParameter(0, pe*0.02);
    exp_longtau->SetParLimits(0, 0.005*pe,0.2*pe);
    //exp_longtau->SetParameter(1, 2*ns);
    //exp_longtau->SetParLimits(1, 0*ns,5*ns);
    exp_longtau->FixParameter(1, 0*ns);
    exp_longtau->SetParameter(2, 50*ns);
    exp_longtau->SetParLimits(2, 30*ns,70*ns);
    //exp_longtau->SetLineColor(2);
        
    cleanwaves->Fit(Form("exptau_%s",vol),"R","",50*ns,120*ns); // Fit boundaries for the slow component of the pulse
    (*amp0) = exp_longtau->GetParameter(0);
    (*tau) = exp_longtau->GetParameter(2);
    std::cout << "Long tau = " << *tau << std::endl;
    ctmp->Write();
    delete ctmp;

    c->cd();

    TPaveText * pv = new TPaveText(0.4,0.65,0.55,0.74,"brNDC");
    pv->AddText(Form("#tau_{long} = %2.2e #pm %2.2e",*tau, exp_longtau->GetParError(2)));
    pv->SetFillColor(kWhite);
    pv->Draw("SAME");
    exp_longtau->Draw("SAME");

    c->Print(globalArgs.res_folder+Form("Clean_%s.pdf",vol));

    return exp_longtau;
}

TF1 * fitAPTau(TGraph * APtime, double amp0, double tau, double pe, const char * vol, TCanvas * c) 
{
    // Fit parameters and limits to calculate AP recharge

    TCanvas * ctmp = new TCanvas(Form("APTauFit_%s",vol));
    APtime->Draw("AP*");
    TF1 * exp = new TF1(Form("exp_%s",vol),"[0]*(1 - exp(-(x-[5])/[1])) + [2]*exp(-(x-[4])/[3])",4*ns,180 * ns);
    //exp->SetParameter(0,pe);
    //exp->SetParLimits(0,0.2*pe,10*pe);    
    exp->FixParameter(0,pe);  
    exp->SetParameter(1,50*ns);
    exp->SetParLimits(1,5*ns,100*ns);
    exp->FixParameter(2,amp0);
    //exp->SetParLimits(2,0.,1.);
    exp->FixParameter(3,tau);
    exp->FixParameter(4,0*ns);
    exp->FixParameter(5,0*ns);
    //exp->SetParameter(5,8*ns);
    //exp->SetParLimits(5,0*ns,15*ns);

    APtime->Fit(Form("exp_%s",vol),"","",30*ns,100*ns);
    exp->Draw("same");
    ctmp->Write();
    delete ctmp;
    
    c->cd();
    
    TPaveText * pv = new TPaveText(0.4,0.65,0.55,0.74,"brNDC");
    pv->AddText(Form("Recovery time = %2.2e #pm %2.2e",exp->GetParameter(1), exp->GetParError(1)));
    pv->SetFillColor(kWhite);
    pv->Draw();
    exp->Draw("SAME");

    c->Print(globalArgs.res_folder+Form("Afterpulse_%s.pdf",vol));

    return exp;
}



#endif