#ifndef THRESHOLDS_H
#define THRESHOLDS_H

#include "TCanvas.h"
#include "TPaveText.h"

struct Thresholds
{
    // In  nanoseconds
    const double AP_minT;
    const double del_xtalk_minT;
    const double dir_xtalk_maxT;
    
    // In  % of pe
    const double AP;
    const double dir_xtalk;
    const double del_xtalk;

    //const double time_dist; // this is not used, I remove it. (Olivier)

//} default_thrs {4., 2., 2., 0.4, 1.17, 0.8, 0.4};	// H2017 Luca
//} default_thrs {25., 2., 2., 0.4, 1.17, 0.8};	// H2017 Olivier
//} default_thrs {25., 2., 2., 0.5, 1.17, 0.8};	// H2017 Olivier
} default_thrs {25., 2., 2., 0.6, 1.17, 0.6};	// H2017 Olivier - PDE single threshold
//} default_thrs {4., 2., 2., 0.6, 1.2, 0.6};
//} default_thrs {10., 5., 5., 0.5, 1.2, 0.8};	//FBK
//} default_thrs {10., 5., 5., 0.6, 1.2, 0.6};

Thresholds defThresholds(double AP_minT, double DeXT_minT, double DiXT_maxT, double AP_thrs, double DiXT_thrs, double DeXT_thrs) {
	Thresholds thrs{ AP_minT, DeXT_minT, DiXT_maxT, AP_thrs, DiXT_thrs, DeXT_thrs };
	return thrs;
}

TCanvas * listThrs(Thresholds thrs) {
    TCanvas * c = new TCanvas("Thresholds","Thresholds",300,300);
    TPaveText * pv = new TPaveText(0.1,0.1,0.9,0.9,"brNDC");
    pv->SetTextAlign(12);
    pv->AddText("DiXT:");                                                       //(pv->GetListOfLines()->Last())->SetTextColor(2);
    pv->AddText(Form("    - Time: [ %.1f, %.1f ] ns",0.,thrs.dir_xtalk_maxT));  //(pv->GetListOfLines()->Last())->SetTextColor(2);
    pv->AddText(Form("    - Ampl: > %.2f PE",thrs.dir_xtalk));                  //(pv->GetListOfLines()->Last())->SetTextColor(2);
    pv->AddText("DeXT:");                                                       //(pv->GetListOfLines()->Last())->SetTextColor(880);
    pv->AddText(Form("    - Time: > %.1f ns",thrs.del_xtalk_minT));             //(pv->GetListOfLines()->Last())->SetTextColor(880);
    pv->AddText(Form("    - Ampl: > %.2f PE",thrs.del_xtalk));                  //(pv->GetListOfLines()->Last())->SetTextColor(880);
    pv->AddText("AP:");                                                         //(pv->GetListOfLines()->Last())->SetTextColor(807);
    pv->AddText(Form("    - Time: > %.1f ns",thrs.AP_minT));                    //(pv->GetListOfLines()->Last())->SetTextColor(807);
    pv->AddText(Form("    - Ampl: [ %.2f, %.2f [ PE",thrs.AP,thrs.del_xtalk));  //(pv->GetListOfLines()->Last())->SetTextColor(807);
    pv->SetFillColor(kWhite);
    pv->Draw();
    return c;
}

#endif
