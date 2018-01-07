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

#~ R.gROOT.ProcessLine(".x lhcbstyle.h+")
import LHCbStyle

c = R.TCanvas()
for ev in t :

    if opts.volt is not None and abs(ev.Vbias - opts.volt) > 0.1: continue
    
    T, A = [],[]
    max_val = 0
    for i in range(int(ev.NsampPerEv)):
        T.append(ev.Times[i])
        A.append(ev.Amps[i])
        if max_val < ev.Amps[i]:
            max_val = ev.Amps[i]
    gr = R.TGraph(ev.NsampPerEv,array('d',T),array('d',A))
    gr.SetTitle("Event display")
    gr.GetXaxis().SetTitle("Time [s]")
    gr.GetYaxis().SetTitle("Signal [PE]")
    
    gr.SetLineColor(1)
    gr.SetLineWidth(1)
    gr.Draw("AL")
    gr.GetXaxis().SetRangeUser(ev.Times[0],ev.Times[ev.NsampPerEv-1])
    c.SetGrid(1)
    
    # draw baseline
    baseline = R.TGraph(2,array('d',[T[0],T[-1]]),array('d',[ev.baseline_shift,ev.baseline_shift]));
    baseline.SetLineStyle(3);
    baseline.SetLineWidth(2);
    baseline.SetLineColor(40);
    baseline.Draw("l+");

    c.Print('wave_ro.pdf')

    out = raw_input("Press Enter to go to the next event or q to quit    ")
    if out == 'q' : break

print "Events finished"
