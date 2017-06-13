# Analysis_waveforms
Original programmed by Felipe Ortega 2016
======

The code reads out waveforms taken with an oscilloscope
and characterizes them in order to quantify the noise
 of a SiPM.

## To run here is an example:
1. Extract the folder inside Analysis_waveforms 
2. Convert the files:
modify the file ```Reformat_all.sh``` in such a waz to point to the folder of interest: ```Analysis_waveforms/H2017/H2017_XX_chYY``` where XX is the detector number and YY is the channel number
3. Run ```source Reformat_all.sh```
4. Create a cfg.txt file inside ```/production2/Analysis_waveforms/H2017/H2017_XX_chYY``` on the basis of ```/production2/Analysis_waveforms/cfg_example.txt```
5. Run the analysis:

```
./output -d /production2/Analysis_waveforms/H2017/H2017_XX_chYY/ -S /production2/Analysis_waveforms/H2017/H2017_XX_chYY/cfg.txt -o /production2/Analysis_waveforms/H2017/H2017_XX_chYY/ -a
```
 
## Compilation:
1. source compile.sh (on a 64 bit platform, otherwise compile_32.sh)

## Convention:
1. Make a data capture on the oscilloscope (eg in 0.5V steps)
2. Make a directory eg H2017_XX_chYY
3. Create directories in this with name the bias voltage values "56.5V", an empty template is located on the osc
4. Make a zip file for transport the data 
5. Place the data directory in eg H2017 
6. Run the analysis

The program ROOT is necessary: -> root.cern.ch
