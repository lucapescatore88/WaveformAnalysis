import ROOT as R
import sys, os
from array import array
import argparse

# program meant to view raw data from the oscilloscope
# saved as binary and converted with "readbinary.py"

parser = argparse.ArgumentParser()
parser.add_argument('datafile',default=None)
parser.add_argument('--volt',type=str,default=None)
opts = parser.parse_args()

f = R.TFile.Open(opts.datafile)
print "Selected voltage:", opts.volt
t = f.Get(opts.volt)
print "Tree entries:",t.GetEntries()

import LHCbStyle

c = R.TCanvas()
for ev in t :
    
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
    gr.GetYaxis().SetTitle("Signal [V]")
    
    gr.SetLineColor(1)
    gr.SetLineWidth(1)
    gr.Draw("AL")
    gr.GetXaxis().SetRangeUser(ev.Times[0],ev.Times[ev.NsampPerEv-1])
    c.SetGrid(1)

    c.Print('raw_wave.pdf')

    out = raw_input("Press Enter to go to the next event or q to quit    ")
    if out == 'q' : break

print "Events finished"
