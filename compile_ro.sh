g++ -o random_overlap.o Analysis/RandomOverlap_analysis.C -c `root-config --cflags` -lSpectrum -IAnalysis -lRooFit -lRooFitCore -lMinuit
g++ -o random_overlap random_overlap.o `root-config --libs` -lSpectrum -lRooFit -lRooFitCore -lMinuit
