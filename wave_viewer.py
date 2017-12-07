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
    halfDT = (ev.Times[1] - ev.Times[0])/2.
    hist = R.TH1F("waveform","waveform", ev.NsampPerEv,ev.Times[0]-halfDT,ev.Times[ev.NsampPerEv-1]+halfDT)   
    max_val = 0
    min_thrs = min(ev.DiXT_thr,min(ev.DeXT_thr,ev.AP_thr))
    artificial_offset = 0.00
    for i in range(int(ev.NsampPerEv)):
        T.append(ev.Times[i])
        A.append(ev.Amps[i])
        hist.SetBinContent(i,ev.Amps[i]+artificial_offset)
        if max_val < ev.Amps[i]:
            max_val = ev.Amps[i]
    gr = R.TGraph(ev.NsampPerEv,array('d',T),array('d',A))
    gr.GetXaxis().SetTitle("Time [s]")
    gr.GetYaxis().SetTitle("Signal [V]")
    
    gr.SetMaximum(1.5*ev.pe+artificial_offset)
    gr.SetLineColor(1)
    gr.SetLineWidth(1)
    gr.Draw("APL")
    c.SetGrid(1)
    
    hist.Draw("same L")
    hist.GetYaxis().SetRangeUser(-0.2*ev.pe,1.5*ev.pe+artificial_offset)
    hist.SetStats(0)
    s = R.TSpectrum()
    sigma = 0.01
    s.Search(hist,sigma,"",min_thrs/max_val)
    
    npeaks = s.GetNPeaks()
    xpeaks_float = s.GetPositionX()
    ypeaks_float = s.GetPositionY()
    for peak in range(npeaks):
        bin_num = hist.FindBin(xpeaks_float[peak])
        max_y = 0
        for b in range(-4,5):
            if max_y < hist.GetBinContent(bin_num+b):
                max_y = hist.GetBinContent(bin_num+b)
        print "peak", peak, " max_y =", max_y
    
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
    
    # draw baseline
    baseline = R.TGraph(2,array('d',[T[0],T[-1]]),array('d',[ev.baseline_shift,ev.baseline_shift]));
    baseline.SetLineStyle(3);
    baseline.SetLineWidth(2);
    baseline.SetLineColor(40);
    baseline.Draw("l+");
    
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
    
    pv3 = R.TPaveText(0.84,0.45,0.99,0.57,"brNDC");
    pv3.AddText("Charge =");
    pv3.AddText("{:.02f} mV#timesns".format(1e3*1e9*ev.pulse_integral));
    pv3.SetFillColor(R.kWhite);
    pv3.Draw();

    c.Print('wave.pdf')

    out = raw_input("Press Enter to go to the next event or q to quit    ")
    if out == 'q' : break

print "Events finished"
