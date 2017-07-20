import ROOT as R
import sys, os
from array import array
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('datafile',default=None)
parser.add_argument('--volt',type=float,default=None)
opts = parser.parse_args()

f = R.TFile.Open(opts.datafile)
t = f.Get('ClassifiedData')
print "Tree entries:",t.GetEntries()

c = R.TCanvas()
for ev in t :

    if opts.volt is not None and abs(ev.Vbias - opts.volt) > 0.1: continue
    
    T, A = [],[]
    for i in range(int(ev.NsampPerEv)):
        T.append(ev.Times[i])
        A.append(ev.Amps[i])
    gr = R.TGraph(ev.NsampPerEv,array('d',T),array('d',A))
    
    gr.SetMaximum(1.5*ev.pe)
    gr.SetLineColor(1)
    gr.Draw("APL")
    
    # draws the thresholds
    pe_line = R.TGraph(2,array('d',[T[0],T[-1]]),array('d',[ev.pe,ev.pe]));
    pe_line.SetLineStyle(9);
    pe_line.SetLineWidth(2);
    pe_line.SetLineColor(4);
    pe_line.Draw("l+");
    
    DiXT_line = R.TGraph(2,array('d',[-10e-9,10e-9]),array('d',[ev.DiXT_thr,ev.DiXT_thr]));
    DiXT_line.SetLineStyle(9);
    DiXT_line.SetLineWidth(2);
    DiXT_line.SetLineColor(2);
    DiXT_line.Draw("l+");
    
    DeXT_line = R.TGraph(2,array('d',[T[0],T[-1]]),array('d',[ev.DeXT_thr,ev.DeXT_thr]));
    DeXT_line.SetLineStyle(9);
    DeXT_line.SetLineWidth(2);
    DeXT_line.SetLineColor(880);
    DeXT_line.Draw("l+");
    
    AP_line = R.TGraph(2,array('d',[T[0],T[-1]]),array('d',[ev.AP_thr,ev.AP_thr]));
    AP_line.SetLineStyle(9);
    AP_line.SetLineWidth(2);
    AP_line.SetLineColor(807);
    AP_line.Draw("l+");
    
    pv = R.TPaveText(0.84,0.74,0.99,0.99,"brNDC");
    pv.AddText(ev.Category.replace("_"," "));
    pv.AddText("#DeltaV = {:2.2f}V".format(ev.dV));
    pv.AddText("N_{{noise peaks}} = {:.0f}".format(ev.NnoisePeaks));
    pv.AddText("N_{{AP}} = {:.0f}".format(ev.NAP));
    pv.AddText("N_{{DeXT}} = {:.0f}".format(ev.NDeXT));
    pv.AddText("N_{{Sec.peaks}} = {:.0f}".format(ev.NSec));
    pv.AddText("PE = {:.3f}".format(ev.pe));	(pv.GetListOfLines().Last()).SetTextColor(4);
    pv.SetFillColor(R.kWhite);
    pv.Draw();
    
    pv2 = R.TPaveText(0.84,0.59,0.99,0.72,"brNDC");
    pv2.AddText("Thresholds used");
    pv2.AddText("Direct CT = {:.03f}".format(ev.DiXT_thr));	(pv2.GetListOfLines().Last()).SetTextColor(2);
    pv2.AddText("Delayed CT = {:.3f}".format(ev.DeXT_thr));	(pv2.GetListOfLines().Last()).SetTextColor(880);
    pv2.AddText("After pulse = {:.03f}".format(ev.AP_thr));	(pv2.GetListOfLines().Last()).SetTextColor(807);
    pv2.SetFillColor(R.kWhite);
    pv2.Draw();

    c.Print('wave.pdf')

    out = raw_input("Press Enter to go to the next event or q to quit    ")
    if out == 'q' : break

print "Events finished"
