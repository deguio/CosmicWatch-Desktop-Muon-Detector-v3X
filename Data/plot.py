#***********************************************************************************
# Master import
#***********************************************************************************

import sys, os, time, warnings, argparse
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np


warnings.filterwarnings('ignore')


font = {
    "family": "serif",
    "serif": "Computer Modern Roman",
    "weight": 200,
    "size": 15
}
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

# Define your own color palette
mycolors = ['#c70039','#ff5733','#ff8d1a','#ffc300','#eddd53','#add45c','#57c785',
               '#00baad','#2a7b9b','#3d3d6b','#511849','#900c3f','#900c3f']

def time_coincidence_analysis(timestamps_s, detector_names, window_ms=10.0):
    """Return the coincident-event mask and number of cross-detector pairs."""
    timestamps_s = np.asarray(timestamps_s, dtype=float)
    detector_names = np.asarray(detector_names, dtype=str)
    if len(timestamps_s) != len(detector_names):
        raise ValueError('Timestamps and detector names must have the same length')
    if window_ms <= 0:
        raise ValueError('The coincidence window must be greater than zero')
    if len(set(detector_names)) < 2:
        raise ValueError(
            'Time-based coincidences require at least two different detector names'
        )

    window_s = window_ms / 1000.0
    coincident = np.zeros(len(timestamps_s), dtype=bool)
    pair_count = 0
    order = np.argsort(timestamps_s, kind='mergesort')
    sorted_times = timestamps_s[order]
    sorted_names = detector_names[order]

    left = 0
    for right in range(len(order)):
        while sorted_times[right] - sorted_times[left] > window_s:
            left += 1
        for candidate in range(left, right):
            if sorted_names[candidate] != sorted_names[right]:
                coincident[order[candidate]] = True
                coincident[order[right]] = True
                pair_count += 1

    return coincident, pair_count


def coincidences_from_computer_time(timestamps_s, detector_names, window_ms=10.0):
    """Mark events from different detectors separated by at most ``window_ms``."""
    return time_coincidence_analysis(timestamps_s, detector_names, window_ms)[0]

class CWClass():
    def __init__(self, fname, bin_size=60, coincidence_source='device',
                 coincidence_window_ms=10.0):
        self.name = fname.split('/')[-1]
        self.bin_size = bin_size
        self.coincidence_source = coincidence_source
        self.coincidence_window_ms = coincidence_window_ms
        self.coincidence_pair_count = None
        if coincidence_source not in ('device', 'time'):
            raise ValueError("coincidence_source must be 'device' or 'time'")
        
        fileHandle = open(fname,"r" )
        lineList = fileHandle.readlines()
        fileHandle.close()
        # Sensor columns are optional. Determine the layout from valid data rows,
        # ignoring comments, incomplete final rows and trailing tab characters.
        supported_columns = (6, 9, 10, 13)
        column_counts = []
        for line in lineList[-200:]:
            stripped = line.rstrip('\t\r\n')
            if stripped and not stripped.lstrip().startswith('#'):
                n_columns = len(stripped.split('\t'))
                if n_columns in supported_columns:
                    column_counts.append(n_columns)
        if not column_counts:
            raise ValueError('No valid CosmicWatch event rows found in file')

        number_of_columns = max(set(column_counts), key=column_counts.count)
        print('Number of columns in file: ', number_of_columns)

        self.file_from_computer = False
        self.file_from_sdcard   = False
        self.has_MPU6050 = False
        self.has_BMP280 = False
        
        if number_of_columns in (9, 13):
            self.file_from_computer = True
            print('  -> File from Computer')
            metadata_columns = (6, 7, 8) if number_of_columns == 9 else (10, 11, 12)
            selected_columns = (0, 1, 2, 3, 4, 5, *metadata_columns)
            data = np.genfromtxt(
                fname, dtype=str, delimiter='\t', comments='#',
                usecols=selected_columns, invalid_raise=False
            )
            data = np.atleast_2d(data)
            event_number = data[:,0].astype(float) #first column of data
            PICO_timestamp_s = data[:,1].astype(float)
            coincident = np.array([
                value.strip().lower() not in ('0', 'false', '') for value in data[:,2]
            ])
            adc = data[:,3].astype(int)
            sipm = data[:,4].astype(float)
            deadtime = data[:,5].astype(float)
            detName = data[:,6]
            comp_time = data[:,7]
            comp_date = data[:,8]
        
        elif number_of_columns in (6, 10):
            print('  -> File from MicroSD Card')
            self.file_from_sdcard = True
            data = np.genfromtxt(
                fname, dtype=str, delimiter='\t', comments='#',
                usecols=(0, 1, 2, 3, 4, 5), invalid_raise=False
            )
            data = np.atleast_2d(data)
            event_number = data[:,0].astype(float)#first column of data
            PICO_timestamp_s = data[:,1].astype(float)
            coincident = np.array([
                value.strip().lower() not in ('0', 'false', '') for value in data[:,2]
            ])
            adc = data[:,3].astype(int)
            sipm = data[:,4].astype(float)
            deadtime = data[:,5].astype(float)

        deadtime = deadtime - min(deadtime)

        # Defaults for files recorded without environmental or motion sensors.
        missing_sensor_data = np.full(len(event_number), np.nan, dtype=float)
        temperature = missing_sensor_data.copy()
        pressure = missing_sensor_data.copy()
        accel_x = missing_sensor_data.copy()
        accel_y = missing_sensor_data.copy()
        accel_z = missing_sensor_data.copy()
        gyro_x = missing_sensor_data.copy()
        gyro_y = missing_sensor_data.copy()
        gyro_z = missing_sensor_data.copy()

        # Read sensor values only from the two legacy full layouts.
        if number_of_columns in (10, 13):
            sensor_data = np.genfromtxt(
                fname, dtype=str, delimiter='\t', comments='#',
                usecols=(6, 7, 8, 9), invalid_raise=False
            )
            sensor_data = np.atleast_2d(sensor_data)
            temperature = sensor_data[:,0].astype(float)
            pressure = sensor_data[:,1].astype(float)

            def split_xyz(values):
                return np.asarray([value.split(':') for value in values], dtype=float).T

            accel_x, accel_y, accel_z = split_xyz(sensor_data[:,2])
            gyro_x, gyro_y, gyro_z = split_xyz(sensor_data[:,3])
            self.has_MPU6050 = True
            self.has_BMP280 = True

        # Convert the computer time to an absolute time (MJD).
        if self.file_from_computer:
            time_stamp = []
            for i in range(len(comp_date)):
                
                day  = int(comp_date[i].split('/')[0])
                month = int(comp_date[i].split('/')[1])
                year   = int(comp_date[i].split('/')[2])
                hour  = int(comp_time[i].split(':')[0])
                mins  = int(comp_time[i].split(':')[1])
                sec   = int(np.floor(float(comp_time[i].split(':')[2])))
                try:  
                    decimal = float('0.'+str(comp_time[i].split('.')[-1]))
                except:
                    decimal = 0.0
                time_stamp.append(float(time.mktime((year, month, day, hour, mins, sec, 0, 0, 0))) + decimal) 


            self.time_stamp_s     = np.asarray(time_stamp) -  min(np.asarray(time_stamp))       # The absolute time of an event in seconds
            self.time_stamp_ms    = self.time_stamp_s*1000  # The absolute time of an event in miliseconds   
            self.total_time_s     = max(self.time_stamp_s) -  min(self.time_stamp_s)     # The absolute time of an event in seconds
            self.detector_name    = detName                                
            self.n_detector       = len(set(detName))

            if coincidence_source == 'time':
                coincident, self.coincidence_pair_count = time_coincidence_analysis(
                    self.time_stamp_s, detName, coincidence_window_ms
                )
                print(
                    '  -> Coincidences calculated from computer Date/Time '
                    'with a %.3f ms window' % coincidence_window_ms
                )

        elif coincidence_source == 'time':
            raise ValueError(
                "Time-based coincidences require a computer file with Name, Time and Date columns"
            )

        if coincidence_source == 'device':
            print('  -> Coincidences read from the device flag column')

        # Convert cumulative deadtime separately for every detector. Interleaving
        # cumulative counters from multiple USB devices would otherwise create
        # artificial negative deadtimes.
        event_deadtime_s = np.zeros(len(deadtime), dtype=float)
        detector_total_deadtimes = []
        if self.file_from_computer:
            for detector_name in set(detName):
                detector_indices = np.where(detName == detector_name)[0]
                detector_deadtime = deadtime[detector_indices]
                detector_deltas = np.diff(np.append(detector_deadtime[0], detector_deadtime))
                event_deadtime_s[detector_indices] = np.maximum(detector_deltas, 0.0)
                detector_total_deadtimes.append(
                    max(detector_deadtime) - min(detector_deadtime)
                )
        else:
            event_deadtime_s = np.maximum(
                np.diff(np.append(deadtime[0], deadtime)), 0.0
            )
            detector_total_deadtimes.append(max(deadtime) - min(deadtime))

        # The RP Pico absolute time isn't great. Over the course of a few hours, it will be off by several seconds. 
        # The computer will give you accurate time down to about 1ms. Reading from the serial port has ~ms scale uncertainty.
        # The RP Pico can give you a precise measurement (down to 1us), but the absolute time will drift. Expect it to be off by roughly 1min per day.
        #self.PICO_time_ms      = PICO_time_ms
        self.PICO_timestamp_s       = PICO_timestamp_s
        
        self.PICO_total_time_s = max(self.PICO_timestamp_s) - min(self.PICO_timestamp_s)
        self.PICO_total_time_ms= self.PICO_total_time_s * 1000.

        self.event_number     = np.asarray(event_number)  # an arrray of the event numbers
        self.total_counts     = len(event_number)
        self.select_coincident        = coincident         # an arrray of the measured event ADC value

        self.adc              = adc         # an arrray of the measured event ADC value
        self.sipm             = sipm        # an arrray of the measured event SiPM value
        
        self.temperature      = temperature         # an arrray of the measured event ADC value
        self.pressure        = pressure         # an arrray of the measured event ADC value

        self.accel_x        = accel_x         # an arrray event acceleration x data
        self.accel_y        = accel_y        # an arrray event acceleration x data
        self.accel_z        = accel_z        # an arrray event acceleration x data

        self.gyro_x        = gyro_x         # an arrray event acceleration x data
        self.gyro_y        = gyro_y        # an arrray event acceleration x data
        self.gyro_z        = gyro_z        # an arrray event acceleration x data    

        self.event_deadtime_s   = event_deadtime_s    # an array of the measured event deadtime in seconds
        #print(self.event_deadtime_s)
        self.event_deadtime_ms  = self.event_deadtime_s*1000            # an array of the measured event deadtime in miliseconds
        # For a merged stream, use the mean detector deadtime against the common
        # wall-clock duration. This keeps the aggregate rate denominator physical.
        self.total_deadtime_s   = float(np.mean(detector_total_deadtimes))
        self.total_deadtime_ms  = self.total_deadtime_s*1000. # The total deadtime in seconds
                
         
        # The time between events is well described by the PICO timestamp. 
        # The 'diff' command takes the difference between each element in the array.
        self.PICO_event_livetime_s = np.diff(np.append([0],self.PICO_timestamp_s)) - self.event_deadtime_s
        
        def round(x, err):
            """Round x and err based on the first significant digit of err."""
            if err == 0:
                return x, err  # Avoid division by zero
            # Find order of magnitude of error
            order_of_magnitude = int(np.floor(np.log10(err)))
            # Find the first significant digit of err
            first_digit = int(err / (10 ** order_of_magnitude))
            # Round both values to the first significant digit of err
            rounded_x = np.round(x, -order_of_magnitude+1)
            rounded_err = np.round(err, -order_of_magnitude+1)#first_digit * (10 ** (order_of_magnitude)) 
            return rounded_x, rounded_err

        '''
        if self.file_from_computer:
            
            self.live_time_s        = (self.total_time_s - self.total_deadtime_s)
            self.weights          = np.ones(len(event_number)) / self.live_time_s
            self.count_rate       = self.total_counts/self.live_time_s 
            self.count_rate_err   = np.sqrt(self.total_counts)/self.live_time_s 

            n = 4
            print("    -- Total Count Rate: ", np.round(self.total_counts/self.live_time_s,n),"+/-",
                    np.round(np.sqrt(self.total_counts)/self.live_time_s,n),"Hz")

            self.count_rate, self.count_rate_err = round(
                    self.total_counts/self.live_time_s, 
                    np.sqrt(self.total_counts)/self.live_time_s)

            bins = range(0,int(max(self.time_stamp_s)), self.bin_size)
            counts, binEdges       = np.histogram(self.time_stamp_s, bins = bins)
            bin_livetime, binEdges = np.histogram(self.time_stamp_s, bins = bins, weights = self.PICO_event_livetime_s)
        
            # Bin the amount of deadtime
            sum_deadtime, _ = np.histogram(self.time_stamp_s, bins=bins, weights=self.pressure)
            count_pressure, _ = np.histogram(self.time_stamp_s, bins=bins)
            self.binned_pressure = sum_pressure / np.maximum(count_pressure, 1)  # Avoid division by zero

            # Bin the pressure by taking the average pressure in each bin
            sum_pressure, _ = np.histogram(self.time_stamp_s, bins=bins, weights=self.pressure)
            count_pressure, _ = np.histogram(self.time_stamp_s, bins=bins)
            self.binned_pressure = sum_pressure / np.maximum(count_pressure, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_temperature, _ = np.histogram(self.time_stamp_s, bins=bins, weights=self.temperature)
            count_temperature, _ = np.histogram(self.time_stamp_s, bins=bins)
            self.binned_temperature = sum_temperature / np.maximum(count_temperature, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_accel_x, _ = np.histogram(self.time_stamp_s, bins=bins, weights=self.accel_x)
            count_accel_x, _ = np.histogram(self.time_stamp_s, bins=bins)
            self.binned_accel_x = sum_accel_x / np.maximum(count_accel_x, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_accel_y, _ = np.histogram(self.time_stamp_s, bins=bins, weights=self.accel_y)
            count_accel_y, _ = np.histogram(self.time_stamp_s, bins=bins)
            self.binned_accel_y = sum_accel_y / np.maximum(count_accel_y, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_accel_z, _ = np.histogram(self.time_stamp_s, bins=bins, weights=self.accel_z)
            count_accel_z, _ = np.histogram(self.time_stamp_s, bins=bins)
            self.binned_accel_z = sum_accel_z / np.maximum(count_accel_z, 1)  # Avoid division by zero
        
        if self.file_from_computer:
            self.live_time_s        = (self.total_time_s - self.total_deadtime_s)
            self.live_time_ms        = self.live_time_s/1000.
            self.weights          = np.ones(len(event_number)) / self.live_time_s

            n = 4
            print("    -- Total Count Rate: ", np.round(self.total_counts/self.live_time_s,n),"+/-",
                    np.round(np.sqrt(self.total_counts)/self.live_time_s,n),"Hz")

            self.count_rate, self.count_rate_err = round(
                    self.total_counts/self.live_time_s, 
                    np.sqrt(self.total_counts)/self.live_time_s)
            
            

            bins = range(int(min(self.PICO_timestamp_s)),int(max(self.PICO_timestamp_s)),self.bin_size)
            counts, binEdges = np.histogram(self.PICO_timestamp_s, bins = bins)
            bin_livetime, binEdges = np.histogram(self.PICO_timestamp_s, bins = bins, weights = self.PICO_event_livetime_s)

            self.bin_size          = bin_size
            self.binned_counts     = counts
            self.binned_counts_err = np.sqrt(counts)
            self.binned_count_rate = counts/bin_livetime
            self.binned_count_rate_err = np.sqrt(counts)/bin_livetime

            counts_coincident, binEdges      = np.histogram(self.PICO_timestamp_s[self.select_coincident], bins = bins)
            bin_deadtime, binEdges      = np.histogram(self.PICO_timestamp_s, bins = bins, weights = self.event_deadtime_s)

            self.total_coincident = len(self.PICO_timestamp_s[self.select_coincident])
            
            print("    -- Count Rate Coincident (coincident): ",np.round(self.total_coincident/self.live_time_s,n),"+/-" ,
                        np.round(np.sqrt(self.total_coincident)/self.live_time_s,n),"Hz")

            self.count_rate_coincident, self.count_rate_err_coincident = round(
                    self.total_coincident/self.live_time_s, 
                    np.sqrt(self.total_coincident)/self.live_time_s)
            
            
			# Bin the amount of deadtime
            self.binned_deadtime_percentage = bin_deadtime/bin_size * 100
            self.binned_counts_coincident     = counts_coincident
            self.binned_counts_err_coincident = np.sqrt(counts_coincident)
            self.binned_count_rate_coincident = counts_coincident/(bin_size-bin_deadtime)
            self.binned_count_rate_err_coincident = np.sqrt(counts_coincident)/(bin_size-bin_deadtime)

            counts_non_coincident, binEdges      = np.histogram(self.PICO_timestamp_s[~self.select_coincident], bins = bins)
            bin_deadtime, binEdges      = np.histogram(self.PICO_timestamp_s, bins = bins, weights = self.event_deadtime_s)
            self.total_non_coincident = len(self.PICO_timestamp_s[~self.select_coincident])
            self.binned_counts_non_coincident     = counts_non_coincident
            self.binned_counts_err_non_coincident = np.sqrt(counts_non_coincident)
            self.binned_count_rate_non_coincident = counts_non_coincident/(bin_size-bin_deadtime)
            self.binned_count_rate_err_non_coincident = np.sqrt(counts_non_coincident)/(bin_size-bin_deadtime)

            print("    -- Count Rate Non-Coincident: ",np.round(self.total_non_coincident/self.live_time_s,n),"+/-",
                        np.round(np.sqrt(self.total_non_coincident)/self.live_time_s,n),"Hz")

            self.count_rate_non_coincident, self.count_rate_err_non_coincident = round(
                    self.total_non_coincident/self.live_time_s, 
                    np.sqrt(self.total_non_coincident)/self.live_time_s)

            sum_pressure, _ = np.histogram(self.PICO_timestamp_s, bins=bins, weights=self.pressure)
            count_pressure, _ = np.histogram(self.PICO_timestamp_s, bins=bins)
            self.binned_pressure = sum_pressure / np.maximum(count_pressure, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_temperature, _ = np.histogram(self.PICO_timestamp_s, bins=bins, weights=self.temperature)
            count_temperature, _ = np.histogram(self.PICO_timestamp_s, bins=bins)
            self.binned_temperature = sum_temperature / np.maximum(count_temperature, 1)  # Avoid division by zero

			
            if self.has_MPU6050:
                # Bin the temperature by taking the average temperature in each bin
                sum_accel_x, _ = np.histogram(self.PICO_timestamp_s, bins=bins, weights=self.accel_x)
                count_accel_x, _ = np.histogram(self.PICO_timestamp_s, bins=bins)
                self.binned_accel_x = sum_accel_x / np.maximum(count_accel_x, 1)  # Avoid division by zero

                # Bin the temperature by taking the average temperature in each bin
                sum_accel_y, _ = np.histogram(self.PICO_timestamp_s, bins=bins, weights=self.accel_y)
                count_accel_y, _ = np.histogram(self.PICO_timestamp_s, bins=bins)
                self.binned_accel_y = sum_accel_y / np.maximum(count_accel_y, 1)  # Avoid division by zero

                # Bin the temperature by taking the average temperature in each bin
                sum_accel_z, _ = np.histogram(self.PICO_timestamp_s, bins=bins, weights=self.accel_z)
                count_accel_z, _ = np.histogram(self.PICO_timestamp_s, bins=bins)
                self.binned_accel_z = sum_accel_z / np.maximum(count_accel_z, 1)  # Avoid division by zero

                # Bin the temperature by taking the average temperature in each bin
                sum_gyro_x, _ = np.histogram(self.PICO_timestamp_s, bins=bins, weights=self.gyro_x)
                count_gyro_x, _ = np.histogram(self.PICO_timestamp_s, bins=bins)
                self.binned_gyro_x = sum_gyro_x / np.maximum(count_gyro_x, 1)  # Avoid division by zero

                # Bin the temperature by taking the average temperature in each bin
                sum_gyro_y, _ = np.histogram(self.PICO_timestamp_s, bins=bins, weights=self.gyro_y)
                count_gyro_y, _ = np.histogram(self.PICO_timestamp_s, bins=bins)
                self.binned_gyro_y = sum_gyro_y / np.maximum(count_gyro_y, 1)  # Avoid division by zero

                # Bin the temperature by taking the average temperature in each bin
                sum_gyro_z, _ = np.histogram(self.PICO_timestamp_s, bins=bins, weights=self.gyro_z)
                count_gyro_z, _ = np.histogram(self.PICO_timestamp_s, bins=bins)
                self.binned_gyro_z = sum_gyro_z / np.maximum(count_gyro_z, 1)  # Avoid division by zero

            
            # Coincident binned data
        '''

        if self.file_from_computer:
            self.live_time_s        = (self.total_time_s - self.total_deadtime_s)
            analysis_timestamp_s = self.time_stamp_s
        elif self.file_from_sdcard:
            self.live_time_s        = (self.PICO_total_time_s - self.total_deadtime_s)
            analysis_timestamp_s = self.PICO_timestamp_s - min(self.PICO_timestamp_s)
        self.live_time_s = max(self.live_time_s, np.finfo(float).eps)
        self.live_time_ms = self.live_time_s * 1000.0
        self.analysis_timestamp_s = np.asarray(analysis_timestamp_s)
        self.weights          = np.ones(len(event_number)) / self.live_time_s

        self.detector_count_rates = {}
        self.coincidence_pair_rate = None
        self.accidental_pair_rate = None
        self.corrected_pair_rate = None
        if self.file_from_computer:
            detector_names = sorted(set(detName))
            self.detector_count_rates = {
                name: np.count_nonzero(detName == name) / self.live_time_s
                for name in detector_names
            }
            if coincidence_source == 'time':
                self.coincidence_pair_rate = self.coincidence_pair_count / self.live_time_s
                window_s = coincidence_window_ms / 1000.0
                self.accidental_pair_rate = sum(
                    2.0 * self.detector_count_rates[detector_names[i]]
                    * self.detector_count_rates[detector_names[j]] * window_s
                    for i in range(len(detector_names))
                    for j in range(i + 1, len(detector_names))
                )
                self.corrected_pair_rate = max(
                    self.coincidence_pair_rate - self.accidental_pair_rate, 0.0
                )

        n = 4
        print("    -- Total Count Rate: ", np.round(self.total_counts/self.live_time_s,n),"+/-",
                np.round(np.sqrt(self.total_counts)/self.live_time_s,n),"Hz")

        self.count_rate, self.count_rate_err = round(
                self.total_counts/self.live_time_s, 
                np.sqrt(self.total_counts)/self.live_time_s)
        
        

        max_analysis_time = float(max(self.analysis_timestamp_s))
        if max_analysis_time <= 0:
            binEdges = np.asarray([0.0, float(self.bin_size)])
        else:
            binEdges = np.arange(0.0, max_analysis_time, float(self.bin_size))
            if len(binEdges) == 0 or binEdges[0] != 0:
                binEdges = np.insert(binEdges, 0, 0.0)
            binEdges = np.append(binEdges, max_analysis_time)

        counts, binEdges = np.histogram(self.analysis_timestamp_s, bins=binEdges)
        bin_deadtime_total, _ = np.histogram(
            self.analysis_timestamp_s, bins=binEdges, weights=self.event_deadtime_s
        )
        detector_count = self.n_detector if self.file_from_computer else 1
        effective_bin_deadtime = bin_deadtime_total / max(detector_count, 1)
        bin_widths = np.diff(binEdges)
        bin_livetime = np.maximum(
            bin_widths - effective_bin_deadtime, np.finfo(float).eps
        )

        self.bin_size          = bin_size
        self.binned_counts     = counts
        self.binned_counts_err = np.sqrt(counts)
        self.binned_count_rate = counts/bin_livetime
        self.binned_count_rate_err = np.sqrt(counts)/bin_livetime

        counts_coincident, _ = np.histogram(
            self.analysis_timestamp_s[self.select_coincident], bins=binEdges
        )

        self.total_coincident = int(np.count_nonzero(self.select_coincident))
        
        print("    -- Count Rate Coincident (coincident): ",np.round(self.total_coincident/self.live_time_s,n),"+/-" ,
                    np.round(np.sqrt(self.total_coincident)/self.live_time_s,n),"Hz")

        self.count_rate_coincident, self.count_rate_err_coincident = round(
                self.total_coincident/self.live_time_s, 
                np.sqrt(self.total_coincident)/self.live_time_s)
        
        
        # Bin the amount of deadtime
        self.binned_deadtime_percentage = effective_bin_deadtime/bin_widths * 100
        self.binned_counts_coincident     = counts_coincident
        self.binned_counts_err_coincident = np.sqrt(counts_coincident)
        self.binned_count_rate_coincident = counts_coincident/bin_livetime
        self.binned_count_rate_err_coincident = np.sqrt(counts_coincident)/bin_livetime

        counts_non_coincident, _ = np.histogram(
            self.analysis_timestamp_s[~self.select_coincident], bins=binEdges
        )
        self.total_non_coincident = int(np.count_nonzero(~self.select_coincident))
        self.binned_counts_non_coincident     = counts_non_coincident
        self.binned_counts_err_non_coincident = np.sqrt(counts_non_coincident)
        self.binned_count_rate_non_coincident = counts_non_coincident/bin_livetime
        self.binned_count_rate_err_non_coincident = np.sqrt(counts_non_coincident)/bin_livetime

        print("    -- Count Rate Non-Coincident: ",np.round(self.total_non_coincident/self.live_time_s,n),"+/-",
                    np.round(np.sqrt(self.total_non_coincident)/self.live_time_s,n),"Hz")

        self.count_rate_non_coincident, self.count_rate_err_non_coincident = round(
                self.total_non_coincident/self.live_time_s, 
                np.sqrt(self.total_non_coincident)/self.live_time_s)

        if self.has_BMP280:
            sum_pressure, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges, weights=self.pressure)
            count_pressure, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges)
            self.binned_pressure = sum_pressure / np.maximum(count_pressure, 1)

            # Bin the temperature by taking the average temperature in each bin
            sum_temperature, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges, weights=self.temperature)
            count_temperature, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges)
            self.binned_temperature = sum_temperature / np.maximum(count_temperature, 1)

        
        if self.has_MPU6050:
            # Bin the temperature by taking the average temperature in each bin
            sum_accel_x, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges, weights=self.accel_x)
            count_accel_x, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges)
            self.binned_accel_x = sum_accel_x / np.maximum(count_accel_x, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_accel_y, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges, weights=self.accel_y)
            count_accel_y, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges)
            self.binned_accel_y = sum_accel_y / np.maximum(count_accel_y, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_accel_z, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges, weights=self.accel_z)
            count_accel_z, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges)
            self.binned_accel_z = sum_accel_z / np.maximum(count_accel_z, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_gyro_x, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges, weights=self.gyro_x)
            count_gyro_x, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges)
            self.binned_gyro_x = sum_gyro_x / np.maximum(count_gyro_x, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_gyro_y, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges, weights=self.gyro_y)
            count_gyro_y, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges)
            self.binned_gyro_y = sum_gyro_y / np.maximum(count_gyro_y, 1)  # Avoid division by zero

            # Bin the temperature by taking the average temperature in each bin
            sum_gyro_z, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges, weights=self.gyro_z)
            count_gyro_z, _ = np.histogram(self.analysis_timestamp_s, bins=binEdges)
            self.binned_gyro_z = sum_gyro_z / np.maximum(count_gyro_z, 1)  # Avoid division by zero

            
            # Coincident binned data

        #else:
        #    print('Error')
        
        bincenters = 0.5*(binEdges[1:]+ binEdges[:-1])
        self.binned_time_s     = bincenters
        self.binned_time_m     = bincenters/60.
        self.weights           = np.ones(len(event_number)) / self.live_time_s  

 
def plusSTD(n,array):
    xh = np.add(n,np.sqrt(np.abs(array)))
    return xh

def subSTD(n,array):
    xl = np.subtract(n,np.sqrt(np.abs(array)))
    return xl


def fill_between_steps(x, y1, y2=0, h_align='mid', ax=None, lw=2, **kwargs):
    # If no Axes object given, grab the current one:
    if ax is None:
        ax = plt.gca()
    
    # First, duplicate the x values
    xx = np.ravel(np.column_stack((x, x)))[1:]
    
    # Now: calculate the average x bin width
    xstep = np.ravel(np.column_stack((x[1:] - x[:-1], x[1:] - x[:-1])))
    xstep = np.concatenate(([xstep[0]], xstep, [xstep[-1]]))
    
    # Now: add one step at the end of the row
    xx = np.append(xx, xx.max() + xstep[-1])

    # Adjust step alignment
    if h_align == 'mid':
        xx -= xstep / 2.
    elif h_align == 'right':
        xx -= xstep

    # Also, duplicate each y coordinate in both arrays
    y1 = np.ravel(np.column_stack((y1, y1)))
    if isinstance(y2, np.ndarray):
        y2 = np.ravel(np.column_stack((y2, y2)))

    # Plotting
    ax.fill_between(xx, y1, y2=y2, lw=lw, **kwargs)
    return ax


class NPlot():
    def __init__(self, 
                 data,
                 weights,
                 colors,
                 labels,
                 xmin,xmax,ymin,ymax,
                 figsize = [8,6],fontsize = 15,nbins = 101, alpha = 0.85,fit_gaussian=False,
                 xscale = 'log',yscale = 'log',xlabel = '',loc = 1,pdf_name='',lw=2, title=''):

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True, 
                                       gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.1})  # `hspace=0` removes space

        # --- Automatically determine axis limits if not provided ---
        if xmin is None:
            xmin = min(np.nanmin(d) for d in data if len(d) > 0)*0.9
        if xmax is None:
            xmax = max(np.nanmax(d) for d in data if len(d) > 0)*1.1

        if ymin is None or ymax is None:
            all_y = []
            for d, w in zip(data, weights):
                counts, edges = np.histogram(d[~np.isnan(d)], bins=nbins, weights=w[~np.isnan(w)])
                all_y.append(counts)
            all_y = np.concatenate(all_y)
            print(np.nanmax(all_y))
            if yscale == 'log':
                ymin = np.nanmin(all_y[all_y > 0]) * 0.9 if ymin is None else ymin
                ymax = np.nanmax(all_y) * 2.4 if ymax is None else ymax
            elif yscale == 'linear':
                ymin = np.nanmin(all_y[all_y > 0]) * 0.9 if ymin is None else ymin
                ymax = np.nanmax(all_y) * 0.2 if ymax is None else ymax
            #print(ymax)

        # --- Choose bin spacing automatically ---
        if xscale == 'log':
            bins = np.logspace(np.log10(xmin), np.log10(xmax), nbins)
        else:
            bins = np.linspace(xmin, xmax, nbins)
    
        ax1.set_axisbelow(True)
        ax1.grid(which='both', linestyle='--', alpha=0.5, zorder=0)
        ax1.set_title(title, fontsize=fontsize + 1)

        hist_data = []
        std = []
        bin_centers = []

        # Define the Gaussian function
        def gaussian(x, a, mu, sigma):
            return a * np.exp(-(x - mu)**2 / (2 * sigma**2))

        for i in range(len(data)):
            valid_data = data[i][~np.isnan(data[i])]
            valid_weights = weights[i][~np.isnan(weights[i])]

            counts, bin_edges = np.histogram(valid_data, bins=bins, weights=valid_weights)
            bin_center = 0.5 * (bin_edges[1:] + bin_edges[:-1])

            sum_weights_sqrd, _ = np.histogram(valid_data, bins=bins, weights=np.power(valid_weights, 2))

            hist_data.append(counts)
            upper_value = plusSTD(counts,sum_weights_sqrd)
            lower_value = subSTD(counts,sum_weights_sqrd)
            std.append([upper_value,lower_value])
            bin_centers.append(bin_center)
            fill_between_steps(bin_center, upper_value,lower_value,  color = colors[i],alpha = alpha,lw=lw,ax=ax1)
            ax1.plot([1e14,1e14], label = labels[i],color = colors[i],alpha = alpha,linewidth = 2)
            
            

        ax1.set_yscale(yscale)
        ax1.set_xscale(xscale)
        ax1.legend(fontsize=fontsize - 2, loc=loc, fancybox=True, frameon=True)
        ax1.set_ylabel(r'Rate/bin [s$^{-1}$]', size=fontsize)
        ax1.set_xlim(xmin, xmax)
        ax1.set_ylim(ymin, ymax)

        
        ax1.tick_params(axis='both', which='major', labelsize=fontsize-3)
        ax1.tick_params(axis='both', which='minor', labelsize=fontsize-3) 

        # --- Ratio Plot ---
        reference_hist = hist_data[0]
        for i in range(1, len(hist_data)):
            #ratio = np.divide(hist_data[i], reference_hist, out=np.zeros_like(hist_data[i]), where=reference_hist != 0)
            upper_value = np.divide(std[i][0], reference_hist, out=np.zeros_like(std[i][0]), where=reference_hist != 0)
            lower_value = np.divide(std[i][1], reference_hist, out=np.zeros_like(std[i][1]), where=reference_hist != 0)
            fill_between_steps(bin_centers[0], upper_value, lower_value, ax=ax2, color=colors[i], alpha=0.85, lw=lw)
            #ax2.plot(bin_centers[0], ratio, marker='.', linestyle='-', color=colors[i], alpha=alpha, label=f'{labels[i]} / {labels[0]}')

        ax2.set_yscale('linear')

        ax2.axhline(1.0, color='black', linestyle='--', linewidth=1)  # Reference line at 1
        ax2.set_ylabel("Ratio", size=fontsize)
        ax2.set_xlabel(xlabel, labelpad=10, size=fontsize)
        ax2.set_ylim(0., 1.)  # Adjust as needed
        ax2.grid(which='both', linestyle='--', alpha=0.5)

        ax2.tick_params(axis='both', which='major', labelsize=fontsize - 3)
        ax2.tick_params(axis='both', which='major', labelsize=fontsize - 3)

        plt.tight_layout()
        
        if pdf_name != '':
            print('Saving Figure to: '+os.getcwd() +  '/'+pdf_name)
            plt.savefig(pdf_name, format='pdf',transparent =True)
        plt.show()

class ratePlot():
    def __init__(self,
                 time,
                 count_rates,
                 count_rates_err,
                 colors,
                 labels,
                 xmin,xmax,ymin,ymax,fmt,
                 figsize = [8,8],fontsize = 16, alpha = 0.9,
                 xscale = 'linear',yscale = 'linear',
                 xlabel = '',ylabel = '',
                 loc = 2, pdf_name='',title = '', legend_extra=None):
        
        f = plt.figure(figsize=(figsize[0], figsize[1])) 
        ax1 = f.add_subplot(111)

        ax1.set_axisbelow(True)
        ax1.grid(which='both', linestyle='--', alpha=0.5, zorder=0)
        if len(fmt)!=len(time):
            fmt = ['ko']*len(time)

        if len(alpha)!=len(time):
            alpha = [1]*len(time)

        for i in range(len(count_rates)):
            plt.errorbar(time[i], 
                           count_rates[i],
                           xerr=0, yerr=count_rates_err[i],alpha = alpha[i],
                           fmt= fmt[i],label = labels[i], linewidth = 2, ecolor = colors[i], markersize = 1.5)

        plt.yscale(yscale)
        plt.xscale(xscale)
        plt.ylabel(ylabel,size=fontsize)
        plt.xlabel(xlabel,size=fontsize)
        plt.axis([xmin, xmax, ymin,ymax])
        
        ax1.tick_params(axis='both', which='major', labelsize=fontsize-3)
        ax1.tick_params(axis='both', which='minor', labelsize=fontsize-3) 
        ax1.xaxis.labelpad = 0 

        for extra_label in legend_extra or []:
            ax1.plot([], [], linestyle='none', marker='', label=extra_label)

        plt.legend(fontsize=fontsize-3,loc = loc,  fancybox = True,frameon=True)
        
        plt.title(title,fontsize=fontsize+1)
        plt.tight_layout()
        if pdf_name != '':
            print('Saving Figure to: '+os.getcwd() +  '/'+pdf_name)
            plt.savefig(pdf_name, format='pdf',transparent =True)
        plt.show()
        
def main():
    parser = argparse.ArgumentParser(description="Process CosmicWatch data.")
    parser.add_argument('-i', '--input', required=True, help="Input file name or full path")
    parser.add_argument('-b', '--bin_width', required=False, help="The width of bins for rate vs time plot in seconds", type=int,default=60)
    parser.add_argument(
        '-c', '--coincidence-source', choices=('device', 'time'), default='device',
        help="Use the device flag or calculate coincidences from computer Time/Date (default: device)"
    )
    parser.add_argument(
        '-w', '--coincidence-window-ms', type=float, default=10.0,
        help="Time coincidence window in milliseconds when --coincidence-source=time (default: 10)"
    )

    args = parser.parse_args()

    infile_name = args.input.split('/')[-1].split('.')[0]
    print("Plotting infile name: ", infile_name + '.txt')
    if os.path.isfile(os.path.join(os.getcwd(), args.input)):
        file_path = os.path.join(os.getcwd(), args.input)
    elif os.path.isfile(args.input):
        file_path = args.input  # Full path provided by user
    else:
        print(f"Error: File '{args.input}' not found in the current directory.")
        print(f"    -- Example Usage: >> python plot.py -i ExampleData/AxLab_000.txt")
        sys.exit(1)

    file_location = os.path.dirname(file_path)
    pdf_file_location = os.path.join(file_location, 'Figures')

    # Create the directory if it doesn't exist
    os.makedirs(pdf_file_location, exist_ok=True)
    
    # Load the data file, set the binsize for the rate as a function of time plot.
    f1 = CWClass(
        file_path,
        bin_size=args.bin_width,
        coincidence_source=args.coincidence_source,
        coincidence_window_ms=args.coincidence_window_ms,
    )

    
    # Plot the ADC values from the coincident and non-coincident events
    c = NPlot(
        data=[ f1.adc,f1.adc[~f1.select_coincident],f1.adc[f1.select_coincident]],
        weights=[f1.weights,f1.weights[~f1.select_coincident],f1.weights[f1.select_coincident]],
        colors=[mycolors[7], mycolors[3],mycolors[1]],
        labels=[r'All Events:  ' + str(f1.count_rate) + '+/-' + str(f1.count_rate_err) +' Hz',
                r'Non-Coincident:  ' + str(f1.count_rate_non_coincident) + '+/-' + str(f1.count_rate_err_non_coincident) +' Hz',
                r'Coincident: ' + str(f1.count_rate_coincident) + '+/-' + str(f1.count_rate_err_coincident) +' Hz'],
        xmin=None, xmax=None, ymin=None, ymax=None,nbins=101,yscale='log',xscale='log',
        xlabel='Meausred ADC peak value [0-4095]',
        pdf_name=pdf_file_location+'/'+infile_name+'_ADC.pdf',title = '')
    
    # Plot the Calculated SiPM Peak voltages coincident and non-coincident events
    c = NPlot(
        data=[f1.sipm, f1.sipm[~f1.select_coincident],f1.sipm[f1.select_coincident]],
        weights=[f1.weights, f1.weights[~f1.select_coincident],f1.weights[f1.select_coincident]],
        colors=[mycolors[7], mycolors[3],mycolors[1]],
        labels=[r'All Events:  ' + str(f1.count_rate) + '+/-' + str(f1.count_rate_err) +' Hz',
                r'Non-Coincident:  ' + str(f1.count_rate_non_coincident) + '+/-' + str(f1.count_rate_err_non_coincident) +' Hz',
                r'Coincident: ' + str(f1.count_rate_coincident) + '+/-' + str(f1.count_rate_err_coincident) +' Hz'],
        xmin=None, xmax=None, ymin=None, ymax=None,xscale='log',nbins = 51,
        xlabel='SiPM Peak Voltage [mV]',fit_gaussian=True,
        pdf_name=pdf_file_location+'/'+infile_name+'_SiPM_peak_voltage.pdf',title = '',)
    
    
    rate_legend_extra = []
    if f1.coincidence_source == 'time':
        rate_legend_extra.extend(
            '%s rate: %.4f Hz' % (name, rate)
            for name, rate in f1.detector_count_rates.items()
        )
        rate_legend_extra.extend([
            'Pair rate: %.4f Hz' % f1.coincidence_pair_rate,
            r'$R_{acc}$ ($\pm$%.3g ms): %.4f Hz'
            % (f1.coincidence_window_ms, f1.accidental_pair_rate),
            'Corrected pair rate: %.4f Hz' % f1.corrected_pair_rate,
        ])

    # Plot rate as a function of time
    c = ratePlot(time = [f1.binned_time_m,f1.binned_time_m,f1.binned_time_m],
        count_rates = [f1.binned_count_rate,f1.binned_count_rate_non_coincident,f1.binned_count_rate_coincident],
        count_rates_err = [f1.binned_count_rate_err,f1.binned_count_rate_err_non_coincident,f1.binned_count_rate_err_coincident], 
        colors=[mycolors[7], mycolors[3], mycolors[1]],
        labels=[r'All Events: ' + str(f1.count_rate) + '+/-' + str(f1.count_rate_err) +' Hz', 
                r'Non-Coincident:  ' + str(f1.count_rate_non_coincident) + '+/-' + str(f1.count_rate_err_non_coincident) +' Hz',
                r'Coincident:  ' + str(f1.count_rate_coincident) + '+/-' + str(f1.count_rate_err_coincident) +' Hz'],
        xmin = min(f1.binned_time_m), xmax = max(f1.binned_time_m),ymin = 0,ymax = 1.35*max(f1.binned_count_rate),
        figsize = [7,5],fmt = ['ko'],
        fontsize = 16,alpha = [1],
        xscale = 'linear',yscale = 'linear',xlabel = 'Time [min]',ylabel = r'Rate [s$^{-1}$]',
        loc = 1, pdf_name=pdf_file_location+'/'+infile_name+'_rate.pdf',title = '',
        legend_extra=rate_legend_extra)

    
    c = ratePlot(time = [f1.binned_time_m,],
        count_rates = [f1.binned_deadtime_percentage],
        count_rates_err = [np.zeros(len(f1.binned_time_m))], # Uncertainty on pressure is 100 Pa
        colors =[mycolors[6]],
        xmin = min(f1.binned_time_m),xmax = max(f1.binned_time_m),ymin = 0.04,ymax =4,
        figsize = [7,5],labels=['Deadtime Percentage'],
        fontsize = 16,alpha = [1],fmt = ['ko'],
        xscale = 'linear',yscale = 'log',xlabel = 'Time [min]',ylabel = r'Deadtime Percentage [%]',
        loc = 2,pdf_name=pdf_file_location+'/'+infile_name+'_deadtime.pdf',title = '')


    if f1.has_BMP280:
        c = ratePlot(time = [f1.binned_time_m,],
            count_rates = [f1.binned_pressure],
            count_rates_err = [np.ones(len(f1.binned_pressure)) * 100],
            colors =[mycolors[6]],
            xmin = min(f1.binned_time_m),xmax = max(f1.binned_time_m),ymin = min(f1.binned_pressure)-500,ymax =max(f1.binned_pressure)+500,
            figsize = [7,5],labels=['Pressure Data'],
            fontsize = 16,alpha = [1],fmt = ['ko'],
            xscale = 'linear',yscale = 'linear',xlabel = 'Time [min]',ylabel = r'Pressure [Pa]',
            loc = 3,pdf_name=pdf_file_location+'/'+infile_name+'_pressure.pdf',title = '')

        c = ratePlot(time = [f1.binned_time_m,],
            count_rates = [f1.binned_temperature],
            count_rates_err = [np.ones(len(f1.binned_temperature))*0.1],
            colors =[mycolors[5]],
            xmin = min(f1.binned_time_m),xmax = max(f1.binned_time_m),ymin = min(f1.binned_temperature)-0.4,ymax = max(f1.binned_temperature)+0.4,
            figsize = [7,5],fmt = ['ko'],
            fontsize = 16,alpha = [1],labels=['Temperature Data'],
            xscale = 'linear',yscale = 'linear',xlabel = 'Time [min]',ylabel = r'Temperature [$^{\circ}$C]',
            loc = 3,pdf_name=pdf_file_location+'/'+infile_name+'_temperature.pdf',title = '')

    if f1.has_MPU6050:
        c = ratePlot(time = [f1.binned_time_m,f1.binned_time_m,f1.binned_time_m,],
            count_rates = [f1.binned_accel_z,f1.binned_accel_y,f1.binned_accel_x],
            count_rates_err = [np.ones(len(f1.binned_accel_z))*0.001,np.ones(len(f1.binned_accel_z))*0.001,np.ones(len(f1.binned_accel_z))*0.001], # Uncertainty on pressure is 0.1C
            colors=[mycolors[7], mycolors[3], mycolors[1]],
            xmin = min(f1.binned_time_m),xmax = max(f1.binned_time_m),ymin = -1.3,ymax = 1.3,#ymin = min(f1.binned_accel_z)-0.001,ymax = max(f1.binned_accel_z)+0.001,
            figsize = [7,5], fmt=['go-','ro-','bo-'],
            fontsize = 16,alpha = [0.3,0.3,0.3],labels=['Acceleration Z Data','Acceleration Y Data','Acceleration X Data'],
            xscale = 'linear',yscale = 'linear',xlabel = 'Time [min]',ylabel = "Averaged Instantaneous \n Linear Acceleration [g]",
            loc = 4,pdf_name=pdf_file_location+'/'+infile_name+'_accel.pdf',title = '')

        c = ratePlot(time = [f1.binned_time_m,f1.binned_time_m,f1.binned_time_m,],
            count_rates = [f1.binned_gyro_z,f1.binned_gyro_y,f1.binned_gyro_x],
            count_rates_err = [np.ones(len(f1.binned_gyro_z))*0.1,np.ones(len(f1.binned_gyro_y))*0.1,np.ones(len(f1.binned_gyro_x))*0.1], # Uncertainty on pressure is 0.1C
            colors=[mycolors[7], mycolors[3], mycolors[1]],
            xmin = min(f1.binned_time_m),xmax = max(f1.binned_time_m),ymin = -100,ymax = 100,
            figsize = [7,5],fmt=['go-','ro-','bo-'],
            fontsize = 16,alpha = [0.3,0.3,0.3],labels=['Gyro Z Data','Gyro Y Data','Gyro X Data'],
            xscale = 'linear',yscale = 'linear',xlabel = 'Time [min]',ylabel = "Averaged Instantaneous \n Angular Velocity [deg/s]",
            loc = 4,pdf_name=pdf_file_location+'/'+infile_name+'_gyro.pdf',title = '')

    t = f1.binned_time_m  # shared x-axis (Time [min])

    
    '''fig, axes = plt.subplots(
        nrows=5, ncols=1,
        sharex=True,
        figsize=(5, 8),
        gridspec_kw={'hspace': 0.1}
    )
    '''

    plot_rows = 1 + (2 if f1.has_BMP280 else 0) + (2 if f1.has_MPU6050 else 0)
    fig, axes = plt.subplots(
        nrows=plot_rows, ncols=1,
        sharex=True,
        figsize=(5, max(3, 1.6 * plot_rows)),
        constrained_layout=True,   # <— replaces tight_layout
        gridspec_kw={'hspace': 0.1}
    )
    axes = np.atleast_1d(axes)

    # 1) Total rate
    axes[0].plot(t, f1.binned_count_rate, color=mycolors[7], label='All Events')
    axes[0].plot(t, f1.binned_count_rate_non_coincident, color=mycolors[3], label='Non-Coincident')
    axes[0].plot(t, f1.binned_count_rate_coincident, color=mycolors[1], label='Coincident')
    axes[0].set_ylabel('Rate [Hz]')
    axes[0].legend(loc='upper right', fontsize=6)
    axes[0].grid(True, which='both', linestyle='--', alpha=0.5)


    axis_index = 1
    if f1.has_BMP280:
        axes[axis_index].plot(t, f1.binned_pressure/1000., 'o-', color=mycolors[6])
        axes[axis_index].set_ylabel('Pressure [kPa]')
        axis_index += 1

        axes[axis_index].plot(t, f1.binned_temperature, 'o-', color=mycolors[5])
        axes[axis_index].set_ylabel('Temperature [°C]')
        axis_index += 1

    # 5) Acceleration (if present)
    if f1.has_MPU6050:
        accel_err = np.full_like(f1.binned_accel_x, 0.001)  # Example 1 mg uncertainty
        axes[axis_index].errorbar(t, f1.binned_accel_x, yerr=accel_err, fmt='o-', color=mycolors[7], alpha=0.7, markersize=2, label='Ax')
        axes[axis_index].errorbar(t, f1.binned_accel_y, yerr=accel_err, fmt='o-', color=mycolors[3], alpha=0.7, markersize=2, label='Ay')
        axes[axis_index].errorbar(t, f1.binned_accel_z, yerr=accel_err, fmt='o-', color=mycolors[1], alpha=0.7, markersize=2, label='Az')
        axes[axis_index].set_ylabel('Accel [g]')
        axes[axis_index].legend(loc='upper right', fontsize=8)
        axis_index += 1

    # 6) Angular velocity (if present)
    if f1.has_MPU6050:
        gyro_err = np.full_like(f1.binned_gyro_x, 0.1)  # Example 0.1°/s uncertainty
        axes[axis_index].errorbar(t, f1.binned_gyro_x, yerr=gyro_err, fmt='o-', color=mycolors[7], alpha=0.7, markersize=2, label='ωx')
        axes[axis_index].errorbar(t, f1.binned_gyro_y, yerr=gyro_err, fmt='o-', color=mycolors[3], alpha=0.7, markersize=2, label='ωy')
        axes[axis_index].errorbar(t, f1.binned_gyro_z, yerr=gyro_err, fmt='o-', color=mycolors[1], alpha=0.7, markersize=2, label='ωz')
        axes[axis_index].set_ylabel('Gyro [°/s]')
        axes[axis_index].legend(loc='upper right', fontsize=8)

    for axis in axes:
        axis.grid(True, which='both', linestyle='--', alpha=0.5)

    # common x-label
    axes[-1].set_xlabel('Time [min]')

    #plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
