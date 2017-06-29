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

from utils import get_checked_cfg
cfgfile="C:\Users\lphe\Desktop\Oscilloscope_Control\config_file.yml"
cfg = get_checked_cfg(cfgfile)

scope.WriteString("VBS app.SaveRecall.Setup.RecallInternal3",1)
for j in range(0,cfg.NumberOfBias):
	FileName= "D:\\HardCopy\\H\\"+cfg.SiPMName+"\\"+ str(cfg.Voltage[j])+"V"
	scope.WriteString("""VBS 'app.SaveRecall.Utilities.Directory = "%s" '""" %FileName,1)
	scope.WriteString("VBS app.SaveRecall.Utilities.CreateDir",1)
	#Loading default OscilloscopeParameter and set threshold and amplitude
	VerticalScale  = cfg.VerticalScale[j]*0.001
	VerticalOffset = cfg.VerticalOffset[j]*0.001
	scope.WriteString("""VBS 'app.Acquisition.C1.VerScale = "%f"'"""%VerticalScale,1)
	scope.WriteString("""VBS 'app.Acquisition.C1.VerOffset = "%f"'"""%VerticalOffset,1)
	KeithleyControl(cfg.Voltage[j],cfg.Port)
	TriggerLevel=cfg.TriggerLevel[j]*0.001
	scope.WriteString("""VBS 'app.Acquisition.Trigger.C1.Level = "%f"'"""%TriggerLevel,1)
	time.sleep(1)
	scope.WriteString("""VBS 'app.SaveRecall.Waveform.SaveTo = "File"'""",1)
	scope.WriteString("""VBS 'app.SaveRecall.Waveform.WaveformDir = "%s" '"""%FileName ,1)
	for l in range (0,20):
		#scope.WriteString("""VBS 'app.Acquisition.TriggerMode = "Single"'""",1)
		#scope.WriteString("VBS? 'return=app.Measure.P4.Out.Result.Value' ",1) #Queries the P1 parameter
		scope.WriteString("VBS app.SaveRecall.Waveform.DoSave",1)
		#value = scope.ReadString(80)#reads a maximum of 80 bytes
		#print (value) #Print value to Interactive Window
		time.sleep(0.3)

#######################################################################################

scope.Disconnect() #Disconnects from the oscilloscope
