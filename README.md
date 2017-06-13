# Analysis_waveforms
Original programmed by Felipe Ortega 2016

==============================

The code reads out waveforms taken with an oscilloscope
and characterizes them in order to quantify the noise
 of a SiPM.

To run here is an example:
```
./output -d /Users/Analysis_waveforms/ov_example/ -S /Users/Analysis_waveforms/config_file.txt -o /Users/Analysis_waveforms/Results/ -a
```

Compilation:
sh compile.sh

Convention:
Make a data capture on the oscilloscope (eg in 0.5V steps)
Make a directory eg Serial20_ch22
Create directories in this with name the bias voltage values "56.5V", an empty template 
is located on the osc
Make a zip file for transport the data 
Place the data directory in eg H2015 
Copy a cfg.txt file to the directory, eg in Serial20_ch22
Exit the results of the analysis in Serial20_ch22

Run the code with:
cd Serial20_ch22
root@lphe1pc68 Serial20_ch22 ../../output -d ./ -S cfg.txt -o ./ -a 

There is an already compiled version of the program in this repository,
but if desired the file "sh compile.sh" together with the cpp files
can be used to compile it. The program ROOT is necessary: -> root.cern.ch

There is a folder with waveforms examples to be analyzed.
The folder Results has the output of the analyzed waveforms
using the configure file also in found in this repository.
