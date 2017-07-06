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

c = R.TCanvas()
for ev in t :

    if opts.volt is not None and abs(ev.V - opts.volt) > 0.1: continue

    T, A = [],[]
    for i in range(int(ev.NsampPerEv)):
        T.append(ev.Times[i])
        A.append(ev.Amps[i])
    gr = R.TGraph(ev.NsampPerEv,array('d',T),array('d',A))
    
    gr.SetMaximum(1.5*ev.pe)
    gr.SetLineColor(1)
    gr.Draw("APL")

    pv = R.TPaveText(0.6,0.6,0.75,0.85,"brNDC");
    pv.AddText(ev.Category.replace("_"," "));
    pv.AddText("V_{{meas}} = {:2.2f}".format(ev.V_meas));
    pv.AddText("N_{{noise peaks}} = {:.0f}".format(ev.NnoisePeaks));
    pv.AddText("N_{{after pulses}} = {:.0f}".format(ev.NAfterPulses));
    pv.AddText("N_{{delayed CT}} = {:.0f}".format(ev.NDelayedCT));
    pv.AddText("PE = {:.3f}".format(ev.pe));
    pv.SetFillColor(R.kWhite);
    pv.Draw();
    
    pv2 = R.TPaveText(0.6,0.45,0.75,0.58,"brNDC");
    pv2.AddText("Thrsholds used");
    pv2.AddText("Direct CT = {:.03f}".format(ev.CT_thr));
    pv2.AddText("Delayed CT = {:.3f}".format(ev.DCT_thr));
    pv2.AddText("After pulse = {:.03f}".format(ev.AP_thr));
    pv2.SetFillColor(R.kWhite);
    pv2.Draw();

    c.Print('wave.pdf')

    out = raw_input("Press Enter to go to the next event or q to quit    ")
    if out == 'q' : break

print "Events finished"
