g++ -o output_large.o Analysis/Noise_analysis_largeEvent.C -c `root-config --cflags` -lSpectrum -IAnalysis
g++ -o output_large output_large.o `root-config --libs` -lSpectrum 
