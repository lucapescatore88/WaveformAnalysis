#Import of needed librarys the pywin32 library
import time
import serial
import sys
import numpy as np
import re
import win32com.client
from KeithleyControl import KeithleyControl

##############################################################################
#Creation of the folders following the configfile

scope=win32com.client.Dispatch("LeCroy.ActiveDSOCtrl.1")  #creates instance of the ActiveDSO control
scope.MakeConnection("IP:128.178.89.125") #Connects to the oscilloscope.  Substitute your IP address

#Loads the config file
from utils import get_checked_cfg
cfgfile="C:\Users\lphe\Desktop\Oscilloscope_Control\Record_Frequency\cfg.yml"
cfg = get_checked_cfg(cfgfile)

#Configuration of the txt file with the measurement results
text_file = open("C:\Users\lphe\Desktop\Oscilloscope_Control\Record_Frequency\PDE.txt","w")
Header = "Bias Voltage[V]      PDE Frequency[*10kHz] \n"
text_file.write(Header)

#Recall internal setup
scope.WriteString("VBS app.SaveRecall.Setup.RecallInternal4",1)

for j in range(0,cfg.NumberOfBias):

	#Load the parameters from the config file for each bias point
	VerticalScale= cfg.VerticalScale[j]*0.001
	VerticalOffset= cfg.VerticalOffset[j]*0.001
	EdgeDetectionThreshold = float(cfg.EdgeThr)
	AmplitudeTreshold=cfg.PeakAmp[j]*EdgeDetectionThreshold*0.001
	
	#Send the amplitude scale and offset to the oscilloscope
	scope.WriteString("""VBS 'app.Acquisition.C1.VerScale = "%f"'"""%VerticalScale,1)
	scope.WriteString("""VBS 'app.Acquisition.C1.VerOffset = "%f"'"""%VerticalOffset,1)

	#Set the Edge at level parameters like an absolut threshold level (Amplitude * Edge detection threshold)
	scope.WriteString("""VBS 'app.Measure.P2.Operator.LevelType="Absolute"'""",1)
	scope.WriteString("""VBS app.Measure.P2.Operator.AbsLevel= "%f"'"""%AmplitudeTreshold,1)

	#Set the Keithley to the right voltage
	KeithleyControl(cfg.Voltage[j],cfg.Port)

	#Wait time to have enough statistics on the mean value of number of peak detected
	time.sleep(1)#30

	#Access the mean value and reads it
	scope.WriteString("""VBS? 'return=app.Measure.P2.Statistics("histo").Result.Mean' """,1) #Queries the P1 parameter
	value = scope.ReadString(80)#reads a maximum of 80 bytes

	#Writte the text file 
	text_file.write(cfg.Voltage[j])
	text_file.write("                  ")
	text_file.write(value.encode('utf8'))
	text_file.write("\n")
	text_file.flush()
	
text_file.close()
	

#######################################################################################


scope.Disconnect() #Disconnects from the oscilloscope


