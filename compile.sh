g++ -o output_large.o Noise_analysis_largeEvent.C -c `root-config --cflags` -lSpectrum 
g++ -o output_large output_large.o `root-config --libs` -lSpectrum 
