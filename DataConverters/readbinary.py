import sys, struct
import readTrc
import os, glob, argparse
from ROOT import TFile, TTree, vector
from array import array

parser = argparse.ArgumentParser()
parser.add_argument('folder',default=None)
parser.add_argument('--tocsv',action="store_true",default=None)
opts = parser.parse_args()

voltages = glob.glob(opts.folder+'/*')               # Input
out = TFile("oscilloscope_out.root","recreate")      # Output

for volt in voltages :                               # Loop on folders
    
    print "Analysing", volt 
    
    files = glob.glob(volt+'/*.trc')                 
    datX, datY, m = readTrc.readTrc( files[0] )      # Decode first file to get metadata
    #print m

    nsamples = m['WAVE_ARRAY_COUNT']
    nsubsamples = m['NOM_SUBARRAY_COUNT']
    samplesPerEvent =  m['WAVE_ARRAY_COUNT'] / m['NOM_SUBARRAY_COUNT']

    print "N files: ", len(files)
    print "N samples per event", samplesPerEvent
    print "Total N events", nsubsamples*len(files)

    # Define tree structure    

    nsamps = array('i',[samplesPerEvent])
    times  = array('d',[0.0]*samplesPerEvent)
    ampls  = array('d',[0.0]*samplesPerEvent)

    tname = volt[len(volt)-(volt[::-1]).find('/'):]
    t = TTree(tname,tname) 

    t.Branch('NsampPerEv',nsamps,'NsampPerEv/I') 
    t.Branch('Times',times,'Times[NsampPerEv]/D')
    t.Branch('Amps',ampls,'Amps[NsampPerEv]/D'.format(nsamps=samplesPerEvent))

    # Fill tree

    delta = 0.
    shift = -1e9
    csvf = None

    for fi,ff in enumerate(files) :

        datX, datY, m = readTrc.readTrc( ff )
        shift = datX[0]
            
        for i in range(len(datX)) :
            if i%samplesPerEvent == 0 : 
                if opts.tocsv :
                    if csvf is not None : csvf.close()
                    csvf = open(volt+"/{0}.csv".format(fi*nsubsamples+i/samplesPerEvent),"w")
                if i > 0 : t.Fill()
                
                delta = datX[i] - shift
                
            if datX[i] - delta < -20e-9 : print "WRONG"
            times[int(i%samplesPerEvent)] = datX[i] - delta
            ampls[int(i%samplesPerEvent)] = datY[i]

            if opts.tocsv : csvf.write("{:.10f},{:.10f}\n".format(datX[i],datY[i]))
        
        t.Fill()
        if csvf is not None : csvf.close()

    t.Write()

out.Close()
    
    





