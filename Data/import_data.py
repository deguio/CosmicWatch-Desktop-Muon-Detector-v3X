# This script requires one library:
# pyserial
# to install, type: >> pip install pyserial

from __future__ import print_function
import serial 
import time
import glob
import sys
import os
import os.path
import signal
from datetime import datetime
from multiprocessing import Process
import numpy as np
import math
import random
import platform

# NEW: check for -p flag (print to screen)
PRINT_TO_SCREEN = ('-p' in sys.argv)

print("\n========================================================")
print("  CosmicWatch: Desktop Muon Detector – Data Acquisition")
print("  Initialize your detectors, then select serial port(s).")
print('  Operating System: ',platform.system())
if PRINT_TO_SCREEN:
    print("  [Option] Printing data to screen enabled (-p)")
print("========================================================")


def signal_handler(signal, frame):
    print('You pressed Ctrl+C!')
    ComPort.close()     
    file.close() 
    sys.exit(0)

def serial_ports():
    if sys.platform.startswith('win'):
        ports = ['COM%s' % (i + 1) for i in range(256)]
    elif sys.platform.startswith('linux') or sys.platform.startswith('cygwin'):
        # this excludes your current terminal "/dev/tty"
        ports = glob.glob('/dev/tty[A-Za-z]*')
    elif sys.platform.startswith('darwin'):
        ports = glob.glob('/dev/tty.*')
    else:
        raise EnvironmentError('Unsupported platform')
        sys.exit(0)
    result = []
    for port in ports:
        try: 
            s = serial.Serial(port)
            s.close()
            result.append(port)
        except (OSError, serial.SerialException):
            pass
    return result

t1 = time.time()
port_list = serial_ports()
if (time.time()-t1)>2:
    print('Listing ports is taking unusually long...')

print('\nWhich ports do you want to read from?')
for i in range(len(port_list)):
    print('  ['+str(i+1)+'] ' + str(port_list[i]))




# Account for Python 2 and Python 3 syntax
if sys.version_info[:3] > (3,0):
    ArduinoPort = input("Select port: ")
    ArduinoPort = ArduinoPort.split(',')

elif sys.version_info[:3] > (2,5,2):
    ArduinoPort = raw_input("Select port(s): ")
    ArduinoPort = ArduinoPort.split(',')

nDetectors = len(ArduinoPort)


port_name_list = []
for i in range(len(ArduinoPort)):
	port_name_list.append(str(port_list[int(ArduinoPort[i])-1]))

# Ask for file name:
cwd = os.getcwd()
print('')
default_fname = cwd+"/CW_data.txt"
if sys.version_info[:3] > (3,0):
    fname = input("Enter file name (press Enter for default: "+default_fname+"):")
elif sys.version_info[:3] > (2,5,2):
    fname = raw_input("Enter file name (press Enter for default: "+default_fname+"):")
if fname == '':
    fname = default_fname
# If the input is just a file name (no path separators), prepend cwd
elif '/' not in fname and '\\' not in fname:
    fname = os.path.join(cwd, fname)
print(' -- Saving data to: '+fname)

# Assign a unique software name to each selected USB port. This name replaces
# the identical firmware default (for example "AxLab") in the output file.
detector_name_list = []
used_detector_names = set()
print('\nAssign a name to each detector (names must be unique):')
for i, port in enumerate(port_name_list):
    default_name = 'Detector_%d' % (i + 1)
    while True:
        prompt = "  Name for %s (press Enter for %s): " % (port, default_name)
        if sys.version_info[:3] > (3,0):
            detector_name = input(prompt)
        else:
            detector_name = raw_input(prompt)
        detector_name = detector_name.strip() or default_name
        detector_name = detector_name.replace('\t', '_').replace('\r', '').replace('\n', '')
        if detector_name in used_detector_names:
            print("  Name '%s' is already in use. Choose a different name." % detector_name)
            continue
        detector_name_list.append(detector_name)
        used_detector_names.add(detector_name)
        break

print('  Detector mapping:')
for port, detector_name in zip(port_name_list, detector_name_list):
    print('    %s -> %s' % (port, detector_name))


print()
for i in range(nDetectors):
    time.sleep(0.1)
    port = port_name_list[i]
    baudrate = 115200
    globals()['Det%s' % str(i)] = serial.Serial(port,baudrate)
    time.sleep(0.1)
file = open(fname, "w")

# Get list of names, using 5 seconds of data.
'''
print('')
print('Acquiring detector names')
det_names = []
t1 = time.time()
while (time.time()-t1) < 5:
    for i in range(nDetectors):
        if globals()['Det%s' % str(i)].inWaiting():
            data = globals()['Det%s' % str(i)].readline().decode().replace('\r\n','')    # Wait and read data 
            data = data.split("\t")
            det_names.append(data[-1])
            
#print("\nHere is a list of the detectors I see:")
det_names = list(set(det_names))
print(det_names)
for i in range(len(det_names)):
    print("  "+str(i+1)+') '+det_names[i])
'''
# Start recording data to file.
print("Taking data ...")
if platform.system() == "Windows":
    print("ctrl+break to termiante process")
else:
    print("Press ctl+c to terminate process")

file.write("###########################################################################################################################################################\n")

file.write("#                                                          CosmicWatch: The Desktop Muon Detector v3X\n")
file.write("#                                                                   Questions? saxani@udel.edu\n")
file.write("# Event  Timestamp[s]  Coincident[bool]  ADC[12b]  SiPM[mV]  Deadtime[s]  Name  Time  Date\n")
file.write("###########################################################################################################################################################\n")

while True:
    for i in range(nDetectors):
        if globals()['Det%s' % str(i)].inWaiting():
            
            raw_data = globals()['Det%s' % str(i)].readline().decode(errors='replace').rstrip('\r\n\t')
            serial_data = raw_data.split("\t")

            # A valid event needs the six core fields. Temperature, pressure,
            # accelerometer, gyroscope and the firmware-provided name are ignored.
            try:
                float(serial_data[0]); float(serial_data[1])
                float(serial_data[3]); float(serial_data[4]); float(serial_data[5])
            except (ValueError, IndexError):
                continue
            
            ti = str(datetime.now()).split(" ")
            comp_time = ti[-1]
            comp_date = ti[0].split('-')
            date_string = comp_date[2] + '/' + comp_date[1] + '/' + comp_date[0]

            # The selected USB port determines the detector identity in the merged file.
            data = serial_data[:6] + [detector_name_list[i], comp_time, date_string]
            file.write('\t'.join(data) + '\n')
            
            # NEW: optionally also print to screen
            if PRINT_TO_SCREEN:
                print('\t'.join(data))


            event_number = int(data[0])
            if event_number % 1 ==0:
                file.flush() 


#for i in range(nDetectors):
globals()['Det%s' % str(0)].close()     
file.close()  


