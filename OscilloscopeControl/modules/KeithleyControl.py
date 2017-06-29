import time
import serial
import sys
import argparse

def KeithleyControl(voltage, port='COM3'):

	ser = serial.Serial(port, 9600, timeout=0)
	ser.isOpen()
	print 'Serial port opened ...'

	# Reset and clear register of KY2400 
	ser.write("*RST" + '\r')
	ser.write("*CLS" + '\r')
	ser.write("*IDN?" + '\r')

	out = '' 
	time.sleep(4)

	while ser.inWaiting() > 0:
		out += ser.read(1)
	if out != '': print out

	if voltage < 0.01 :
		print "OUTPUT off"
		ser.write("OUTP OFF" + '\r')
	else :
		ser.write("FORM:ELEM VOLT,CURR" + '\r')
		ser.write("SOUR:VOLT:RANG:AUTO 1" + '\r')

		vcommand = "SOUR:VOLT:LEV:IMM:AMPL " + str(voltage)
		ser.write(vcommand + '\r')

		ser.write("OUTP ON" + '\r')
		ser.write("INIT" + '\r')

		ser.write("READ?" + '\r')

		out = '' 
		time.sleep(1)
		while ser.inWaiting() > 0:
			out += ser.read(1)
		if out != '': print out

		time.sleep(3)
		ser.write("SYST:LOC" + '\r')

	ser.close()
	print 'Serial port closed ...'


if __name__ == '__main__' :

    parser = argparse.ArgumentParser()
    parser.add_argument('voltage', help="Volts to set or 0 to turn off the output")
    opts = parser.parse_args()

    KeithleyControl(opts.voltage)

