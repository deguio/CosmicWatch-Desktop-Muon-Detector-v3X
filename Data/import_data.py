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
import threading
import queue

# NEW: check for -p flag (print to screen)
PRINT_TO_SCREEN = ('-p' in sys.argv)

print("\n========================================================")
print("  CosmicWatch: Desktop Muon Detector – Data Acquisition")
print("  Initialize your detectors, then select serial port(s).")
print('  Operating System: ',platform.system())
if PRINT_TO_SCREEN:
    print("  [Option] Printing data to screen enabled (-p)")
print("========================================================")


def current_time_ns():
    """Return wall-clock time with the best resolution available."""
    if hasattr(time, 'time_ns'):
        return time.time_ns()
    return int(time.time() * 1_000_000_000)


def format_pc_timestamp(timestamp_ns):
    """Convert a nanosecond Unix timestamp to the existing Time/Date columns."""
    seconds, nanoseconds = divmod(timestamp_ns, 1_000_000_000)
    timestamp = datetime.fromtimestamp(seconds)
    # Keep all nine fractional digits reported by time_ns(); plot.py accepts
    # arbitrary decimal precision in the Time column.
    comp_time = timestamp.strftime('%H:%M:%S') + '.%09d' % nanoseconds
    comp_date = timestamp.strftime('%d/%m/%Y')
    return comp_time, comp_date


def build_output_row(raw_bytes, detector_name, timestamp_ns):
    """Validate one serial event and add its software name and PC timestamp."""
    raw_data = raw_bytes.decode(errors='replace').rstrip('\r\n\t')
    serial_data = raw_data.split('\t')
    try:
        float(serial_data[0]); float(serial_data[1])
        float(serial_data[3]); float(serial_data[4]); float(serial_data[5])
    except (ValueError, IndexError):
        return None

    comp_time, comp_date = format_pc_timestamp(timestamp_ns)
    return serial_data[:6] + [detector_name, comp_time, comp_date]


def read_detector(detector_index, connection, event_queue, stop_event):
    """Continuously read one USB device and timestamp each complete line immediately."""
    while not stop_event.is_set():
        try:
            raw_bytes = connection.readline()
            received_ns = current_time_ns()
        except (OSError, serial.SerialException) as error:
            if not stop_event.is_set():
                event_queue.put(('error', detector_index, error))
            return
        if raw_bytes:
            event_queue.put(('event', detector_index, received_ns, raw_bytes))

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
detectors = []
for i in range(nDetectors):
    time.sleep(0.1)
    port = port_name_list[i]
    baudrate = 115200
    # A short timeout lets the reader threads stop promptly on Ctrl+C while each
    # thread otherwise blocks efficiently waiting for its own USB device.
    detectors.append(serial.Serial(port, baudrate, timeout=0.1))
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
file.write("# PC timestamp assigned immediately after receipt of each complete serial line\n")
for port, detector_name in zip(port_name_list, detector_name_list):
    file.write("# Detector alias: %s = %s\n" % (detector_name, port))
file.write("# Event  Timestamp[s]  Coincident[bool]  ADC[12b]  SiPM[mV]  Deadtime[s]  Name  Time  Date\n")
file.write("###########################################################################################################################################################\n")
file.flush()

# One dedicated reader per USB device. The main thread only serializes already
# timestamped rows to disk, so disk I/O never delays timestamp assignment.
event_queue = queue.Queue()
stop_event = threading.Event()
reader_threads = []
for detector_index, connection in enumerate(detectors):
    reader = threading.Thread(
        target=read_detector,
        args=(detector_index, connection, event_queue, stop_event),
        name='CosmicWatchReader-%d' % (detector_index + 1),
        daemon=True,
    )
    reader.start()
    reader_threads.append(reader)

last_flush = time.monotonic()

def write_queued_event(item):
    if item[0] == 'error':
        detector_index, error = item[1], item[2]
        print(
            'Serial connection lost for %s: %s'
            % (detector_name_list[detector_index], error),
            file=sys.stderr,
        )
        return

    detector_index, received_ns, raw_bytes = item[1], item[2], item[3]
    data = build_output_row(
        raw_bytes, detector_name_list[detector_index], received_ns
    )
    if data is None:
        return
    line = '\t'.join(data)
    file.write(line + '\n')
    if PRINT_TO_SCREEN:
        print(line)


try:
    while True:
        try:
            queued_item = event_queue.get(timeout=0.25)
            write_queued_event(queued_item)
        except queue.Empty:
            pass

        # Flushing periodically instead of after every event greatly reduces disk
        # overhead while limiting unwritten buffered data to about half a second.
        if time.monotonic() - last_flush >= 0.5:
            file.flush()
            last_flush = time.monotonic()
except KeyboardInterrupt:
    print('\nStopping acquisition ...')
finally:
    stop_event.set()
    for connection in detectors:
        connection.close()
    for reader in reader_threads:
        reader.join(timeout=1.0)

    # Preserve events that were timestamped just before Ctrl+C.
    while True:
        try:
            write_queued_event(event_queue.get_nowait())
        except queue.Empty:
            break

    file.flush()
    file.close()
    print('Data file closed: ' + fname)
