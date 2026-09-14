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


def build_output_row(raw_bytes, timestamp_ns):
    """Validate an event, preserve its sensors and add PC metadata.

    USB rows end with the detector name configured on the microSD card.
    """
    raw_data = raw_bytes.decode(errors='replace').rstrip('\r\n\t')
    serial_data = raw_data.split('\t')
    try:
        float(serial_data[0]); float(serial_data[1])
        float(serial_data[3]); float(serial_data[4]); float(serial_data[5])
    except (ValueError, IndexError):
        return None

    # Current firmware emits 7/9/11 USB fields: six core values,
    # zero/two/four sensor values, and its configured detector name.
    if len(serial_data) not in (7, 9, 11):
        return None
    sensor_values = serial_data[6:-1]
    firmware_name = serial_data[-1].strip()
    if not firmware_name:
        return None

    if len(sensor_values) == 0:
        sensor_names = ()
    elif len(sensor_values) == 2:
        if all(':' in value for value in sensor_values):
            sensor_names = ('accel', 'gyro')
        else:
            sensor_names = ('temperature', 'pressure')
    elif len(sensor_values) == 4:
        sensor_names = ('temperature', 'pressure', 'accel', 'gyro')
    else:
        return None

    try:
        if 'temperature' in sensor_names:
            float(sensor_values[sensor_names.index('temperature')])
            float(sensor_values[sensor_names.index('pressure')])
        for vector_name in ('accel', 'gyro'):
            if vector_name in sensor_names:
                components = sensor_values[sensor_names.index(vector_name)].split(':')
                if len(components) != 3:
                    return None
                for component in components:
                    float(component)
    except ValueError:
        return None

    comp_time, comp_date = format_pc_timestamp(timestamp_ns)
    output_row = (
        serial_data[:6] + sensor_values
        + [firmware_name, comp_time, comp_date]
    )
    return output_row, sensor_names, firmware_name


def output_header(sensor_names):
    """Build a column header matching the sensor schema received over USB."""
    labels = [
        'Event', 'Timestamp[s]', 'Coincident[bool]', 'ADC[12b]',
        'SiPM[mV]', 'Deadtime[s]',
    ]
    sensor_labels = {
        'temperature': 'Temp[C]',
        'pressure': 'Pressure[Pa]',
        'accel': 'Accel(X:Y:Z)[g]',
        'gyro': 'Gyro(X:Y:Z)[deg/sec]',
    }
    labels.extend(sensor_labels[name] for name in sensor_names)
    labels.extend(['Name', 'Time', 'Date'])
    return '# ' + '  '.join(labels)


def read_detector(detector_index, connection, event_queue, stop_event):
    """Read one USB device and timestamp the arrival of each line's first byte."""
    while not stop_event.is_set():
        try:
            first_byte = connection.read(1)
            if not first_byte:
                continue
            received_ns = current_time_ns()
            raw_bytes = first_byte + connection.read_until(b'\n')
        except (OSError, serial.SerialException) as error:
            if not stop_event.is_set():
                event_queue.put(('error', detector_index, error))
            return
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
# Device names and sensor columns are announced together after every selected
# USB port has delivered its first valid event.
print("Waiting for data from the selected devices ...")

file.write("###########################################################################################################################################################\n")

file.write("#                                                          CosmicWatch: The Desktop Muon Detector v3X\n")
file.write("#                                                                   Questions? saxani@udel.edu\n")
file.write("# PC timestamp assigned immediately after receipt of each serial line's first byte\n")
file.flush()

# One dedicated reader per USB device timestamps the first byte, then completes
# the line. The main thread only serializes already timestamped rows to disk.
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
record_separator = "#" * 155 + "\n"
record_schema = None
reported_schema_mismatches = set()
firmware_name_by_detector = {}
detector_index_by_firmware_name = {}
reported_name_errors = set()
acquisition_announced = False

def write_queued_event(item):
    global record_schema, acquisition_announced
    if item[0] == 'error':
        detector_index, error = item[1], item[2]
        print(
            'Serial connection lost for %s: %s'
            % (port_name_list[detector_index], error),
            file=sys.stderr,
        )
        return

    detector_index, received_ns, raw_bytes = item[1], item[2], item[3]
    parsed_event = build_output_row(raw_bytes, received_ns)
    if parsed_event is None:
        return
    data, sensor_names, firmware_name = parsed_event

    known_name = firmware_name_by_detector.get(detector_index)
    if known_name is None:
        existing_detector = detector_index_by_firmware_name.get(firmware_name)
        if existing_detector is not None and existing_detector != detector_index:
            error_key = (detector_index, firmware_name)
            if error_key not in reported_name_errors:
                print(
                    "Skipping %s: firmware name '%s' is already used by %s. "
                    "Assign unique names in each microSD config file."
                    % (
                        port_name_list[detector_index], firmware_name,
                        port_name_list[existing_detector],
                    ),
                    file=sys.stderr,
                )
                reported_name_errors.add(error_key)
            return
        firmware_name_by_detector[detector_index] = firmware_name
        detector_index_by_firmware_name[firmware_name] = detector_index
    elif firmware_name != known_name:
        error_key = (detector_index, firmware_name)
        if error_key not in reported_name_errors:
            print(
                "Skipping %s: firmware name changed from '%s' to '%s'"
                % (port_name_list[detector_index], known_name, firmware_name),
                file=sys.stderr,
            )
            reported_name_errors.add(error_key)
        return

    if record_schema is None:
        record_schema = sensor_names
        file.write(output_header(record_schema) + '\n')
        file.write(record_separator)
        file.flush()
    elif sensor_names != record_schema:
        mismatch = (firmware_name, sensor_names)
        if mismatch not in reported_schema_mismatches:
            print(
                'Skipping events from %s: sensor columns %s do not match %s'
                % (
                    firmware_name,
                    ', '.join(sensor_names) or 'none',
                    ', '.join(record_schema) or 'none',
                ),
                file=sys.stderr,
            )
            reported_schema_mismatches.add(mismatch)
        return

    if (not acquisition_announced
            and len(firmware_name_by_detector) == nDetectors):
        print('\nDetected USB devices:')
        for selected_index, port in enumerate(port_name_list):
            print(
                '  %s -> %s'
                % (port, firmware_name_by_detector[selected_index])
            )
        print(
            'USB sensor columns: '
            + (', '.join(record_schema) if record_schema else 'none')
        )
        print('Taking data ...')
        if platform.system() == 'Windows':
            print('Press Ctrl+Break to terminate process')
        else:
            print('Press Ctrl+C to terminate process')
        acquisition_announced = True

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
