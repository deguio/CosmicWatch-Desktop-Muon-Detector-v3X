#***********************************************************************************
# Master import
#***********************************************************************************

import sys, os, time, warnings, argparse
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, least_squares
from scipy.signal import fftconvolve
from scipy.stats import landau
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

# Compact legends leave more of the measured distributions visible.
LEGEND_FONTSIZE = 9
DENSE_LEGEND_FONTSIZE = 8
LEGEND_STYLE = {
    'fancybox': True,
    'frameon': True,
    'borderpad': 0.35,
    'labelspacing': 0.3,
    'handlelength': 1.8,
    'handletextpad': 0.5,
}

# Mode of scipy.stats.landau in its standardized parametrization.  Including
# this offset makes the fitted ``mpv`` parameter the actual distribution peak.
LANDAU_STANDARD_MODE = -0.4293145383

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


def unique_time_coincidence_pairs(timestamps_s, detector_names, window_ms=10.0):
    """Build closest one-to-one matches independently for each detector pair."""
    timestamps_s = np.asarray(timestamps_s, dtype=float)
    detector_names = np.asarray(detector_names, dtype=str)
    if len(timestamps_s) != len(detector_names):
        raise ValueError('Timestamps and detector names must have the same length')
    if window_ms <= 0:
        raise ValueError('The coincidence window must be greater than zero')

    window_s = window_ms / 1000.0
    names = sorted(set(detector_names))
    unique_pairs = []
    for name_index, name_a in enumerate(names):
        indices_a = np.where(detector_names == name_a)[0]
        indices_a = indices_a[np.argsort(timestamps_s[indices_a], kind='mergesort')]
        for name_b in names[name_index + 1:]:
            indices_b = np.where(detector_names == name_b)[0]
            indices_b = indices_b[np.argsort(timestamps_s[indices_b], kind='mergesort')]
            times_b = timestamps_s[indices_b]
            candidates = []
            for index_a in indices_a:
                time_a = timestamps_s[index_a]
                first = np.searchsorted(times_b, time_a - window_s, side='left')
                last = np.searchsorted(times_b, time_a + window_s, side='right')
                for position_b in range(first, last):
                    index_b = int(indices_b[position_b])
                    candidates.append((
                        abs(timestamps_s[index_b] - time_a),
                        int(index_a), index_b,
                    ))

            # Events cannot be reused within the same detector pair. They may
            # still participate in a different detector pair (a valid triple or
            # higher-order coincidence).
            used_a = set()
            used_b = set()
            for delta_t, index_a, index_b in sorted(candidates):
                if index_a in used_a or index_b in used_b:
                    continue
                used_a.add(index_a)
                used_b.add(index_b)
                unique_pairs.append((index_a, index_b, delta_t))
    return unique_pairs


def strict_all_detector_coincidence_groups(timestamps_s, detector_names,
                                           window_ms=10.0):
    """Return disjoint groups containing one event from every detector.

    A group is accepted only when its full span, ``max(time) - min(time)``, is
    no larger than ``window_ms``. Events are consumed chronologically and can
    belong to at most one group.
    """
    timestamps_s = np.asarray(timestamps_s, dtype=float)
    detector_names = np.asarray(detector_names, dtype=str)
    if len(timestamps_s) != len(detector_names):
        raise ValueError('Timestamps and detector names must have the same length')
    if window_ms <= 0:
        raise ValueError('The coincidence window must be greater than zero')

    names = sorted(set(detector_names))
    if len(names) < 2:
        raise ValueError(
            'Strict time coincidences require at least two different detector names'
        )

    window_s = window_ms / 1000.0
    indices_by_name = {
        name: np.flatnonzero(detector_names == name)[
            np.argsort(timestamps_s[detector_names == name], kind='mergesort')
        ]
        for name in names
    }
    positions = {name: 0 for name in names}
    groups = []

    while all(positions[name] < len(indices_by_name[name]) for name in names):
        group_indices = tuple(
            int(indices_by_name[name][positions[name]]) for name in names
        )
        group_times = np.asarray(
            [timestamps_s[index] for index in group_indices], dtype=float
        )
        earliest_time = float(np.min(group_times))
        latest_time = float(np.max(group_times))

        if latest_time - earliest_time <= window_s:
            groups.append((group_indices, latest_time - earliest_time))
            for name in names:
                positions[name] += 1
        else:
            # Every future event is at least as late as the current head of its
            # detector. Therefore an earliest head outside this span can never
            # participate in a valid all-detector group and may be discarded.
            for name, index in zip(names, group_indices):
                if timestamps_s[index] == earliest_time:
                    positions[name] += 1

    return groups


def cross_detector_time_differences(timestamps_s, detector_names, max_abs_ms=20.0):
    """Return signed time differences for every unordered detector pair.

    For a pair ``A, B`` the sign convention is ``t_B - t_A``. All event
    combinations inside ``+/- max_abs_ms`` are retained, so the prompt peak
    and the approximately flat accidental background can both be inspected.
    """
    timestamps_s = np.asarray(timestamps_s, dtype=float)
    detector_names = np.asarray(detector_names, dtype=str)
    if len(timestamps_s) != len(detector_names):
        raise ValueError('Timestamps and detector names must have the same length')
    if max_abs_ms <= 0:
        raise ValueError('The delta-t histogram range must be greater than zero')

    max_abs_s = max_abs_ms / 1000.0
    names = sorted(set(detector_names))
    differences = {}
    for i, name_a in enumerate(names):
        times_a = np.sort(timestamps_s[detector_names == name_a])
        for name_b in names[i + 1:]:
            times_b = np.sort(timestamps_s[detector_names == name_b])
            pair_differences = []
            for time_a in times_a:
                first = np.searchsorted(times_b, time_a - max_abs_s, side='left')
                last = np.searchsorted(times_b, time_a + max_abs_s, side='right')
                if last > first:
                    pair_differences.extend((times_b[first:last] - time_a) * 1000.0)
            differences[(name_a, name_b)] = np.asarray(pair_differences, dtype=float)
    return differences


def plot_delta_t_histogram(timestamps_s, detector_names, window_ms, max_abs_ms,
                           pdf_name, nbins=201):
    """Plot cross-detector delta-t distributions and the coincidence window."""
    differences = cross_detector_time_differences(
        timestamps_s, detector_names, max_abs_ms=max_abs_ms
    )
    if not differences:
        print('Skipping delta-t histogram: fewer than two detectors are present')
        return

    fig, ax = plt.subplots(figsize=(8, 5.5))
    bins = np.linspace(-max_abs_ms, max_abs_ms, nbins)
    plotted = False
    for index, ((name_a, name_b), values) in enumerate(differences.items()):
        if len(values) == 0:
            continue
        ax.hist(
            values, bins=bins, histtype='step', linewidth=1.6,
            color=mycolors[index % len(mycolors)], label=f'{name_b} - {name_a}'
        )
        plotted = True

    if not plotted:
        plt.close(fig)
        print('Skipping delta-t histogram: no event pairs fall inside the requested range')
        return

    ax.axvline(-window_ms, color='black', linestyle='--', linewidth=1.2)
    ax.axvline(window_ms, color='black', linestyle='--', linewidth=1.2,
               label=r'Coincidence window: $\pm$%.3g ms' % window_ms)
    ax.set_xlabel(r'$\Delta t = t_B-t_A$ [ms]')
    ax.set_ylabel('Event pairs / bin')
    ax.set_xlim(-max_abs_ms, max_abs_ms)
    ax.grid(which='both', linestyle='--', alpha=0.5)
    ax.legend(fontsize=DENSE_LEGEND_FONTSIZE, **LEGEND_STYLE)
    fig.tight_layout()
    print('Saving Figure to: ' + os.getcwd() + '/' + pdf_name)
    fig.savefig(pdf_name, format='pdf', transparent=True)
    plt.show()

class CWClass():
    def __init__(self, fname, bin_size=60, coincidence_source='device',
                 coincidence_window_ms=10.0,
                 time_coincidence_mode='pairwise'):
        self.name = fname.split('/')[-1]
        self.bin_size = bin_size
        self.coincidence_source = coincidence_source
        self.coincidence_window_ms = coincidence_window_ms
        self.time_coincidence_mode = time_coincidence_mode
        self.coincidence_pair_count = None
        self.unique_coincidence_pairs = []
        self.coincidence_groups = []
        self.coincidence_group_count = None
        if coincidence_source not in ('device', 'time'):
            raise ValueError("coincidence_source must be 'device' or 'time'")
        if time_coincidence_mode not in ('pairwise', 'all'):
            raise ValueError(
                "time_coincidence_mode must be 'pairwise' or 'all'"
            )
        if coincidence_source != 'time' and time_coincidence_mode != 'pairwise':
            raise ValueError(
                "time_coincidence_mode='all' requires coincidence_source='time'"
            )
        
        # Take one immutable snapshot of the input.  The acquisition process may
        # keep appending to the original file, but every column parsed below must
        # come from this same set of lines.
        with open(fname, "r") as fileHandle:
            lineList = fileHandle.readlines()
        # Sensor columns are optional. Determine the layout from valid data rows,
        # ignoring comments, incomplete final rows and trailing tab characters.
        # The firmware may independently enable BMP280 (Temp/Press) and MPU6050
        # (Accel/Gyro), so valid files are not limited to the old 6/10 layouts.
        column_counts = []
        sample_fields = None
        for line in lineList[-200:]:
            stripped = line.rstrip('\t\r\n')
            if stripped and not stripped.lstrip().startswith('#'):
                n_columns = len(stripped.split('\t'))
                if 6 <= n_columns <= 13:
                    column_counts.append(n_columns)
        if not column_counts:
            raise ValueError('No valid CosmicWatch event rows found in file')

        number_of_columns = max(set(column_counts), key=column_counts.count)
        print('Number of columns in file: ', number_of_columns)

        for line in lineList:
            stripped = line.rstrip('\t\r\n')
            if stripped and not stripped.lstrip().startswith('#'):
                fields = stripped.split('\t')
                if len(fields) == number_of_columns:
                    sample_fields = fields
                    break

        event_headers = [
            line for line in lineList
            if line.lstrip().startswith('# Event') and 'Timestamp' in line
        ]
        event_header = event_headers[-1] if event_headers else ''

        def looks_like_computer_metadata(fields):
            if fields is None or len(fields) < 9:
                return False
            try:
                time_parts = fields[-2].split(':')
                date_parts = fields[-1].split('/')
                return len(time_parts) == 3 and len(date_parts) == 3
            except (AttributeError, ValueError):
                return False

        file_from_computer = looks_like_computer_metadata(sample_fields)
        data_column_count = number_of_columns - (3 if file_from_computer else 0)
        sensor_count = data_column_count - 6
        if sensor_count < 0:
            raise ValueError('CosmicWatch event rows contain fewer than 6 columns')

        # Recover the sensor order from the firmware-generated header. This
        # distinguishes, for example, an 8-column Temp/Press file from an
        # 8-column Accel/Gyro file.
        sensor_markers = {
            'temperature': ('Temp[', 'Temperature['),
            'pressure': ('Press[', 'Pressure['),
            'accel': ('Accel(', 'Acceleration('),
            'gyro': ('Gyro(', 'Gyroscope('),
        }
        sensor_positions = []
        for sensor_name, markers in sensor_markers.items():
            positions = [event_header.find(marker) for marker in markers]
            positions = [position for position in positions if position >= 0]
            if positions:
                sensor_positions.append((min(positions), sensor_name))
        sensor_names = [name for _, name in sorted(sensor_positions)]

        # Header-less legacy files remain supported through the shape/content
        # of their optional columns.
        if len(sensor_names) != sensor_count:
            sensor_values = sample_fields[6:data_column_count]
            if sensor_count == 0:
                sensor_names = []
            elif sensor_count == 4:
                sensor_names = ['temperature', 'pressure', 'accel', 'gyro']
            elif sensor_count == 2 and all(':' in value for value in sensor_values):
                sensor_names = ['accel', 'gyro']
            elif sensor_count == 2:
                sensor_names = ['temperature', 'pressure']
            else:
                raise ValueError(
                    'Unable to identify the %d optional sensor column(s)'
                    % sensor_count
                )

        if sensor_names:
            print('  -> Sensor columns: ' + ', '.join(sensor_names))
        else:
            print('  -> No sensor columns')

        def valid_snapshot_row(fields):
            """Reject incomplete/malformed rows captured while the file grows."""
            if len(fields) != number_of_columns:
                return False
            try:
                float(fields[0])
                float(fields[1])
                float(fields[3])
                float(fields[4])
                float(fields[5])

                for column_offset, sensor_name in enumerate(sensor_names):
                    value = fields[6 + column_offset]
                    if sensor_name in ('temperature', 'pressure'):
                        float(value)
                    else:
                        components = value.split(':')
                        if len(components) != 3:
                            return False
                        for component in components:
                            float(component)

                if file_from_computer:
                    time_parts = fields[-2].split(':')
                    date_parts = fields[-1].split('/')
                    if len(time_parts) != 3 or len(date_parts) != 3:
                        return False
                    int(time_parts[0])
                    int(time_parts[1])
                    float(time_parts[2])
                    int(date_parts[0])
                    int(date_parts[1])
                    int(date_parts[2])
            except (TypeError, ValueError, IndexError):
                return False
            return True

        snapshot_rows = []
        for line in lineList:
            stripped = line.rstrip('\t\r\n')
            if not stripped or stripped.lstrip().startswith('#'):
                continue
            fields = stripped.split('\t')
            if valid_snapshot_row(fields):
                snapshot_rows.append(fields)

        if not snapshot_rows:
            raise ValueError('No complete CosmicWatch event rows found in file snapshot')
        data = np.asarray(snapshot_rows, dtype=str)

        self.file_from_computer = False
        self.file_from_sdcard   = False
        self.has_MPU6050 = False
        self.has_BMP280 = False
        
        if file_from_computer:
            self.file_from_computer = True
            print('  -> File from Computer')
            event_number = data[:,0].astype(float) #first column of data
            PICO_timestamp_s = data[:,1].astype(float)
            coincident = np.array([
                value.strip().lower() not in ('0', 'false', '') for value in data[:,2]
            ])
            adc = data[:,3].astype(int)
            sipm = data[:,4].astype(float)
            deadtime = data[:,5].astype(float)
            detName = data[:,data_column_count]
            comp_time = data[:,data_column_count + 1]
            comp_date = data[:,data_column_count + 2]
        
        else:
            print('  -> File from MicroSD Card')
            self.file_from_sdcard = True
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

        if sensor_count:
            sensor_data = data[:, 6:data_column_count]
            sensor_columns = {
                name: sensor_data[:, index]
                for index, name in enumerate(sensor_names)
            }

            def split_xyz(values):
                components = [value.split(':') for value in values]
                if any(len(component) != 3 for component in components):
                    raise ValueError('Invalid three-axis sensor value in data file')
                return np.asarray(components, dtype=float).T

            if 'temperature' in sensor_columns:
                temperature = sensor_columns['temperature'].astype(float)
            if 'pressure' in sensor_columns:
                pressure = sensor_columns['pressure'].astype(float)
            if 'accel' in sensor_columns:
                accel_x, accel_y, accel_z = split_xyz(sensor_columns['accel'])
            if 'gyro' in sensor_columns:
                gyro_x, gyro_y, gyro_z = split_xyz(sensor_columns['gyro'])

            self.has_BMP280 = (
                'temperature' in sensor_columns and 'pressure' in sensor_columns
            )
            self.has_MPU6050 = (
                'accel' in sensor_columns and 'gyro' in sensor_columns
            )

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
                if time_coincidence_mode == 'pairwise':
                    coincident, self.coincidence_pair_count = time_coincidence_analysis(
                        self.time_stamp_s, detName, coincidence_window_ms
                    )
                    self.unique_coincidence_pairs = unique_time_coincidence_pairs(
                        self.time_stamp_s, detName, coincidence_window_ms
                    )
                    print(
                        '  -> Pairwise coincidences calculated from computer '
                        'Date/Time with a %.3f ms window' % coincidence_window_ms
                    )
                    print(
                        '  -> %d unique one-to-one detector-pair matches '
                        '(%d total pair associations)'
                        % (len(self.unique_coincidence_pairs), self.coincidence_pair_count)
                    )
                else:
                    self.coincidence_groups = strict_all_detector_coincidence_groups(
                        self.time_stamp_s, detName, coincidence_window_ms
                    )
                    self.coincidence_group_count = len(self.coincidence_groups)
                    coincident = np.zeros(len(self.time_stamp_s), dtype=bool)
                    for group_indices, _ in self.coincidence_groups:
                        coincident[list(group_indices)] = True

                    # Keep the pair representation available to the existing
                    # per-detector spectrum code. All pairs here belong to a
                    # rigorously accepted N-detector group.
                    for group_indices, _ in self.coincidence_groups:
                        for first_position in range(len(group_indices)):
                            for second_position in range(
                                    first_position + 1, len(group_indices)):
                                first = group_indices[first_position]
                                second = group_indices[second_position]
                                delta_t = abs(
                                    self.time_stamp_s[second]
                                    - self.time_stamp_s[first]
                                )
                                self.unique_coincidence_pairs.append(
                                    (first, second, delta_t)
                                )
                    print(
                        '  -> Strict all-detector coincidences calculated from '
                        'computer Date/Time with a %.3f ms maximum span'
                        % coincidence_window_ms
                    )
                    print(
                        '  -> %d disjoint %d-fold groups; no event reused'
                        % (self.coincidence_group_count, self.n_detector)
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
        self.coincident_sipm_by_detector = {}
        if self.file_from_computer and self.unique_coincidence_pairs:
            paired_indices = np.asarray([
                index
                for first, second, _ in self.unique_coincidence_pairs
                for index in (first, second)
            ], dtype=int)
            for detector_name in sorted(set(detName)):
                detector_indices = np.unique(
                    paired_indices[detName[paired_indices] == detector_name]
                )
                self.coincident_sipm_by_detector[detector_name] = self.sipm[detector_indices]
        
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
        self.coincidence_group_rate = None
        self.coincidence_group_rate_err = None
        self.accidental_group_rate = None
        self.corrected_group_rate = None
        if self.file_from_computer:
            detector_names = sorted(set(detName))
            self.detector_count_rates = {
                name: np.count_nonzero(detName == name) / self.live_time_s
                for name in detector_names
            }
            if coincidence_source == 'time':
                window_s = coincidence_window_ms / 1000.0
                if time_coincidence_mode == 'pairwise':
                    self.coincidence_pair_rate = (
                        self.coincidence_pair_count / self.live_time_s
                    )
                    self.accidental_pair_rate = sum(
                        2.0 * self.detector_count_rates[detector_names[i]]
                        * self.detector_count_rates[detector_names[j]] * window_s
                        for i in range(len(detector_names))
                        for j in range(i + 1, len(detector_names))
                    )
                    self.corrected_pair_rate = max(
                        self.coincidence_pair_rate - self.accidental_pair_rate,
                        0.0,
                    )
                else:
                    multiplicity = len(detector_names)
                    self.coincidence_group_rate = (
                        self.coincidence_group_count / self.live_time_s
                    )
                    self.coincidence_group_rate_err = (
                        np.sqrt(self.coincidence_group_count) / self.live_time_s
                    )
                    # For independent, low-occupancy Poisson streams, the
                    # relative-time volume satisfying max(t)-min(t) <= W is
                    # N*W**(N-1). This reduces to 2*R1*R2*W for N=2.
                    self.accidental_group_rate = (
                        multiplicity
                        * window_s ** (multiplicity - 1)
                        * np.prod(list(self.detector_count_rates.values()))
                    )
                    self.corrected_group_rate = max(
                        self.coincidence_group_rate
                        - self.accidental_group_rate,
                        0.0,
                    )

        n = 4
        print("    -- Total Count Rate: ", np.round(self.total_counts/self.live_time_s,n),"+/-",
                np.round(np.sqrt(self.total_counts)/self.live_time_s,n),"Hz")

        self.count_rate, self.count_rate_err = round(
                self.total_counts/self.live_time_s, 
                np.sqrt(self.total_counts)/self.live_time_s)

        if (coincidence_source == 'time'
                and time_coincidence_mode == 'all'):
            print(
                '    -- %d-fold group rate: %.6g Hz'
                % (self.n_detector, self.coincidence_group_rate)
            )
            print(
                '    -- Accidental %d-fold group rate: %.6g Hz'
                % (self.n_detector, self.accidental_group_rate)
            )
            print(
                '    -- Corrected %d-fold group rate: %.6g Hz'
                % (self.n_detector, self.corrected_group_rate)
            )
        
        

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
        self.binned_coincidence_group_rate = None
        self.binned_coincidence_group_rate_err = None
        if (coincidence_source == 'time'
                and time_coincidence_mode == 'all'):
            group_timestamps = np.asarray([
                np.mean(self.analysis_timestamp_s[list(group_indices)])
                for group_indices, _ in self.coincidence_groups
            ], dtype=float)
            group_counts, _ = np.histogram(group_timestamps, bins=binEdges)
            self.binned_coincidence_group_rate = group_counts / bin_livetime
            self.binned_coincidence_group_rate_err = (
                np.sqrt(group_counts) / bin_livetime
            )

        self.total_coincident = int(np.count_nonzero(self.select_coincident))
        
        coincident_description = (
            '%d-fold event rows' % self.n_detector
            if coincidence_source == 'time' and time_coincidence_mode == 'all'
            else 'coincident'
        )
        print("    -- Count Rate Coincident (%s): " % coincident_description,np.round(self.total_coincident/self.live_time_s,n),"+/-" ,
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

        # Environmental and motion data are binned independently for every
        # detector.  The combined curve is then the arithmetic mean of those
        # detector means, so a higher event rate cannot give one device more
        # weight than another.
        if self.file_from_computer:
            sensor_detector_masks = {
                name: detName == name for name in sorted(set(detName))
            }
        else:
            sensor_detector_masks = {
                os.path.splitext(self.name)[0]: np.ones(len(event_number), dtype=bool)
            }
        self.sensor_detector_names = list(sensor_detector_masks)

        def binned_sensor_by_detector(sensor_values):
            sensor_values = np.asarray(sensor_values, dtype=float)
            result = {}
            for detector_name, detector_mask in sensor_detector_masks.items():
                valid = detector_mask & np.isfinite(sensor_values)
                sums, _ = np.histogram(
                    self.analysis_timestamp_s[valid], bins=binEdges,
                    weights=sensor_values[valid],
                )
                samples, _ = np.histogram(
                    self.analysis_timestamp_s[valid], bins=binEdges,
                )
                means = np.full(len(samples), np.nan, dtype=float)
                np.divide(sums, samples, out=means, where=samples > 0)
                result[detector_name] = means
            return result

        def equal_detector_mean(series_by_detector):
            stacked = np.asarray(list(series_by_detector.values()), dtype=float)
            contributors = np.sum(np.isfinite(stacked), axis=0)
            mean = np.full(stacked.shape[1], np.nan, dtype=float)
            np.divide(
                np.nansum(stacked, axis=0), contributors,
                out=mean, where=contributors > 0,
            )
            return mean

        self.binned_pressure_by_detector = {}
        self.binned_temperature_by_detector = {}
        self.binned_accel_x_by_detector = {}
        self.binned_accel_y_by_detector = {}
        self.binned_accel_z_by_detector = {}
        self.binned_gyro_x_by_detector = {}
        self.binned_gyro_y_by_detector = {}
        self.binned_gyro_z_by_detector = {}

        if self.has_BMP280:
            self.binned_pressure_by_detector = binned_sensor_by_detector(self.pressure)
            self.binned_temperature_by_detector = binned_sensor_by_detector(self.temperature)
            self.binned_pressure = equal_detector_mean(
                self.binned_pressure_by_detector
            )
            self.binned_temperature = equal_detector_mean(
                self.binned_temperature_by_detector
            )

        if self.has_MPU6050:
            self.binned_accel_x_by_detector = binned_sensor_by_detector(self.accel_x)
            self.binned_accel_y_by_detector = binned_sensor_by_detector(self.accel_y)
            self.binned_accel_z_by_detector = binned_sensor_by_detector(self.accel_z)
            self.binned_gyro_x_by_detector = binned_sensor_by_detector(self.gyro_x)
            self.binned_gyro_y_by_detector = binned_sensor_by_detector(self.gyro_y)
            self.binned_gyro_z_by_detector = binned_sensor_by_detector(self.gyro_z)

            self.binned_accel_x = equal_detector_mean(self.binned_accel_x_by_detector)
            self.binned_accel_y = equal_detector_mean(self.binned_accel_y_by_detector)
            self.binned_accel_z = equal_detector_mean(self.binned_accel_z_by_detector)
            self.binned_gyro_x = equal_detector_mean(self.binned_gyro_x_by_detector)
            self.binned_gyro_y = equal_detector_mean(self.binned_gyro_y_by_detector)
            self.binned_gyro_z = equal_detector_mean(self.binned_gyro_z_by_detector)

            
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
                 figsize = [8,6],fontsize = 15,nbins = 101, alpha = 0.85,
                 fit_gaussian=False, fit_landau=False, landau_data_index=0,
                 landau_fit_range=(15.0, 80.0), landau_initial_mpv=40.0,
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
        hist_uncertainties = []
        hist_edges = []
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
            hist_uncertainties.append(np.sqrt(sum_weights_sqrd))
            hist_edges.append(bin_edges)
            upper_value = plusSTD(counts,sum_weights_sqrd)
            lower_value = subSTD(counts,sum_weights_sqrd)
            std.append([upper_value,lower_value])
            bin_centers.append(bin_center)
            fill_between_steps(bin_center, upper_value,lower_value,  color = colors[i],alpha = alpha,lw=lw,ax=ax1)
            ax1.plot([1e14,1e14], label = labels[i],color = colors[i],alpha = alpha,linewidth = 2)

        if fit_landau:
            requested_fit_indices = np.atleast_1d(landau_data_index).astype(int)
            for requested_index in requested_fit_indices:
                fit_index = min(max(requested_index, 0), len(hist_data) - 1)
                centers = bin_centers[fit_index]
                values = hist_data[fit_index]
                uncertainties = hist_uncertainties[fit_index]
                edges = hist_edges[fit_index]
                fit_min, fit_max = landau_fit_range
                fit_mask = (
                    np.isfinite(centers) & np.isfinite(values) & (values >= 0)
                    & (centers >= fit_min) & (centers <= fit_max)
                )

                try:
                    if np.count_nonzero(fit_mask) < 5:
                        raise RuntimeError('not enough populated bins in the fit interval')
                    fit_x = centers[fit_mask]
                    fit_y = values[fit_mask]
                    fit_low_edges = edges[:-1][fit_mask]
                    fit_high_edges = edges[1:][fit_mask]
                    fit_sigma = uncertainties[fit_mask]
                    positive_sigma = fit_sigma[fit_sigma > 0]
                    if len(positive_sigma) > 0:
                        fit_sigma = np.where(
                            fit_sigma > 0, fit_sigma, np.min(positive_sigma)
                        )
                    else:
                        fit_sigma = None

                    # Integrate the true Landau PDF over each logarithmic
                    # histogram bin. ``mpv`` is the actual peak position.
                    def landau_bin_model(_x, area, mpv, width):
                        upper = (
                            (fit_high_edges - mpv) / width + LANDAU_STANDARD_MODE
                        )
                        lower = (
                            (fit_low_edges - mpv) / width + LANDAU_STANDARD_MODE
                        )
                        return area * (landau.cdf(upper) - landau.cdf(lower))

                    initial_width = 10.0
                    initial_amplitude = max(np.sum(fit_y), np.finfo(float).eps)
                    parameters, covariance = curve_fit(
                        landau_bin_model, fit_x, fit_y,
                        p0=(initial_amplitude, landau_initial_mpv, initial_width),
                        sigma=fit_sigma,
                        absolute_sigma=fit_sigma is not None,
                        bounds=(
                            (0.0, fit_min, 0.05),
                            (np.inf, fit_max, fit_max - fit_min),
                        ),
                        maxfev=20000,
                    )
                    _, fitted_mpv, fitted_width = parameters
                    mpv_error = float(np.sqrt(max(covariance[1, 1], 0.0)))
                    fitted_values = landau_bin_model(fit_x, *parameters)
                    if fit_sigma is not None:
                        chi_squared = np.sum(((fit_y - fitted_values) / fit_sigma) ** 2)
                        ndf = max(len(fit_y) - len(parameters), 1)
                        reduced_chi_squared = chi_squared / ndf
                    else:
                        reduced_chi_squared = np.nan
                    fit_name = labels[fit_index].split(':', 1)[0]
                    fit_name = fit_name.replace('Coincident ', '')
                    ax1.plot(
                        fit_x, fitted_values,
                        color=colors[fit_index], linestyle='--', linewidth=1.8,
                        label=(r'Landau %s: MPV=%.2f $\pm$ %.2f mV, '
                               r'w=%.2f mV'
                               % (fit_name, fitted_mpv, mpv_error, fitted_width)),
                    )
                    print(
                        '    -- Landau %s fit: MPV = %.3f +/- %.3f mV, '
                        'width = %.3f mV, chi2/ndf = %.2f'
                        % (fit_name, fitted_mpv, mpv_error, fitted_width,
                           reduced_chi_squared)
                    )
                except (RuntimeError, ValueError, FloatingPointError) as error:
                    print('Warning: SiPM Landau fit failed for %s: %s'
                          % (labels[fit_index], error))

        ax1.set_yscale(yscale)
        ax1.set_xscale(xscale)
        legend_fontsize = (
            DENSE_LEGEND_FONTSIZE if len(labels) > 3 else LEGEND_FONTSIZE
        )
        ax1.legend(fontsize=legend_fontsize, loc=loc, **LEGEND_STYLE)
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


def plot_coincident_sipm_langauss(values, live_time_s, detector_name, pdf_name,
                                  fit_range=(25.0, 200.0), nbins=51,
                                  color='#c70039', singles_values=None,
                                  accidental_rate_hz=0.0, x_range=None):
    """Plot one detector's unique coincident pulses and fit Landau (x) Gaussian."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if len(values) < 10:
        print('Skipping Langau fit for %s: not enough coincident events'
              % detector_name)
        return None

    if x_range is None:
        xmin = max(float(np.min(values)) * 0.9, 0.1)
        xmax = float(np.max(values)) * 1.1
    else:
        xmin, xmax = map(float, x_range)
        if not np.isfinite(xmin) or not np.isfinite(xmax) or xmin <= 0 or xmax <= xmin:
            raise ValueError('The shared SiPM x range must satisfy 0 < xmin < xmax')
    bins = np.logspace(np.log10(xmin), np.log10(xmax), nbins)
    counts, edges = np.histogram(values, bins=bins)
    rates = counts / live_time_s
    errors = np.sqrt(counts) / live_time_s
    centers = 0.5 * (edges[:-1] + edges[1:])

    # The accidental pulse-height shape is the singles spectrum. Normalize it
    # to the analytic accidental rate involving this detector and keep it fixed
    # while fitting the Langau signal component.
    accidental_rates = np.zeros_like(rates)
    if singles_values is not None and accidental_rate_hz > 0:
        singles_values = np.asarray(singles_values, dtype=float)
        singles_values = singles_values[
            np.isfinite(singles_values) & (singles_values > 0)
        ]
        if len(singles_values) > 0:
            singles_counts, _ = np.histogram(singles_values, bins=edges)
            singles_rate_hz = len(singles_values) / live_time_s
            accidental_rates = (
                singles_counts / live_time_s
                * accidental_rate_hz / singles_rate_hz
            )
    fit_min, fit_max = fit_range
    fit_mask = (centers >= fit_min) & (centers <= fit_max)
    fit_x = centers[fit_mask]
    fit_y = rates[fit_mask]
    fit_counts = counts[fit_mask]
    fit_accidental_rates = accidental_rates[fit_mask]
    fit_low_edges = edges[:-1][fit_mask]
    fit_high_edges = edges[1:][fit_mask]
    fit_sigma = errors[fit_mask]
    fit_sigma = np.where(fit_sigma > 0, fit_sigma, 1.0 / live_time_s)

    # A fixed fine grid makes the numerical convolution stable during fitting.
    grid_step = 0.5  # mV
    grid_min = fit_min - 1200.0
    grid_max = fit_max + 2000.0
    grid = np.arange(grid_min, grid_max + grid_step, grid_step)
    if len(grid) % 2 == 0:
        grid = np.append(grid, grid[-1] + grid_step)
    kernel_x = (np.arange(len(grid)) - len(grid) // 2) * grid_step

    def convolved_pdf_and_cdf(mpv, landau_width, gaussian_sigma):
        standardized = (
            (grid - mpv) / landau_width + LANDAU_STANDARD_MODE
        )
        landau_pdf = landau.pdf(standardized) / landau_width
        gaussian_pdf = np.exp(-0.5 * (kernel_x / gaussian_sigma) ** 2)
        gaussian_pdf /= np.sqrt(2.0 * np.pi) * gaussian_sigma
        convolved_pdf = fftconvolve(
            landau_pdf, gaussian_pdf, mode='same'
        ) * grid_step
        convolved_pdf = np.maximum(convolved_pdf, 0.0)
        convolved_cdf = np.zeros_like(convolved_pdf)
        convolved_cdf[1:] = np.cumsum(
            0.5 * (convolved_pdf[1:] + convolved_pdf[:-1]) * grid_step
        )
        return convolved_pdf, convolved_cdf

    def langauss_bin_model(_x, area, mpv, landau_width, gaussian_sigma):
        _, convolved_cdf = convolved_pdf_and_cdf(
            mpv, landau_width, gaussian_sigma
        )
        upper = np.interp(fit_high_edges, grid, convolved_cdf)
        lower = np.interp(fit_low_edges, grid, convolved_cdf)
        return area * (upper - lower) + fit_accidental_rates

    fit_result = None
    try:
        if np.count_nonzero(fit_mask) < 6:
            raise RuntimeError('not enough bins in the fit interval')
        initial_parameters = np.asarray([
            max(np.sum(fit_y - fit_accidental_rates), np.finfo(float).eps),
            40.0, 7.0, 5.0,
        ])

        def poisson_deviance_residuals(parameters):
            expected_counts = np.maximum(
                langauss_bin_model(fit_x, *parameters) * live_time_s,
                np.finfo(float).eps,
            )
            deviance_terms = expected_counts - fit_counts
            positive = fit_counts > 0
            deviance_terms[positive] += fit_counts[positive] * np.log(
                fit_counts[positive] / expected_counts[positive]
            )
            signs = np.sign(fit_counts - expected_counts)
            return signs * np.sqrt(np.maximum(2.0 * deviance_terms, 0.0))

        optimization = least_squares(
            poisson_deviance_residuals,
            x0=initial_parameters,
            bounds=(
                (0.0, fit_min, 0.1, 0.1),
                (np.inf, fit_max, 200.0, 200.0),
            ),
            max_nfev=30000,
        )
        if not optimization.success:
            raise RuntimeError(optimization.message)
        parameters = optimization.x
        covariance = np.linalg.pinv(optimization.jac.T @ optimization.jac)
        parameter_errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        fitted_values = langauss_bin_model(fit_x, *parameters)
        poisson_deviance = np.sum(
            poisson_deviance_residuals(parameters) ** 2
        )
        ndf = max(len(fit_y) - len(parameters), 1)
        reduced_chi_squared = poisson_deviance / ndf
        area, mpv, landau_width, gaussian_sigma = parameters
        convolved_pdf, _ = convolved_pdf_and_cdf(
            mpv, landau_width, gaussian_sigma
        )
        langauss_peak = float(grid[np.argmax(convolved_pdf)])
        fit_result = {
            'mpv': mpv,
            'mpv_error': parameter_errors[1],
            'landau_width': landau_width,
            'landau_width_error': parameter_errors[2],
            'gaussian_sigma': gaussian_sigma,
            'gaussian_sigma_error': parameter_errors[3],
            'peak': langauss_peak,
            'chi2_ndf': reduced_chi_squared,
            'accidental_rate': accidental_rate_hz,
        }
        print(
            '    -- Langau %s: Landau MPV = %.3f +/- %.3f mV, '
            'Landau width = %.3f +/- %.3f mV, Gaussian sigma = '
            '%.3f +/- %.3f mV, convolved peak = %.3f mV, deviance/ndf = %.2f'
            % (detector_name, mpv, parameter_errors[1], landau_width,
               parameter_errors[2], gaussian_sigma, parameter_errors[3],
               langauss_peak, reduced_chi_squared)
        )
    except (RuntimeError, ValueError, FloatingPointError) as error:
        print('Warning: Langau fit failed for %s: %s'
              % (detector_name, error))

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 6), sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.08}
    )
    lower_rates = np.maximum(rates - errors, 0.0)
    upper_rates = rates + errors
    fill_between_steps(
        centers, upper_rates, lower_rates, color=color, alpha=0.35, ax=ax1
    )
    ax1.step(
        edges[:-1], rates, where='post', color=color, linewidth=1.8,
        label='Coincident %s: %.4f Hz'
        % (detector_name, len(values) / live_time_s)
    )

    if np.any(accidental_rates > 0):
        ax1.step(
            edges[:-1], accidental_rates, where='post', color='gray',
            linestyle=':', linewidth=1.6,
            label=r'Accidental template: %.4f Hz' % accidental_rate_hz
        )

    if fit_result is not None:
        if (
            fit_result['gaussian_sigma_error']
            > 2.0 * fit_result['gaussian_sigma']
            or fit_result['gaussian_sigma'] <= 0.2
        ):
            gaussian_label = r'$\sigma_G$ unconstrained'
            print('Warning: Gaussian resolution is not constrained for %s'
                  % detector_name)
        else:
            gaussian_label = (
                r'$\sigma_G$=%.2f $\pm$ %.2f mV'
                % (fit_result['gaussian_sigma'],
                   fit_result['gaussian_sigma_error'])
            )
        ax1.plot(
            fit_x, fitted_values, color='black', linestyle='--', linewidth=2,
            label=(r'Langau + accidental: MPV$_L$=%.2f $\pm$ %.2f mV, '
                   '\n' r'%s, deviance/ndf=%.1f'
                   % (fit_result['mpv'], fit_result['mpv_error'],
                      gaussian_label, fit_result['chi2_ndf']))
        )
        expected_counts = np.maximum(
            fitted_values * live_time_s, np.finfo(float).eps
        )
        pulls = (fit_counts - expected_counts) / np.sqrt(expected_counts)
        ax2.axhline(0.0, color='black', linewidth=1)
        ax2.plot(fit_x, pulls, marker='o', linestyle='none', color=color,
                 markersize=3)

    positive_rates = rates[rates > 0]
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xmin, xmax)
    if len(positive_rates) > 0:
        ax1.set_ylim(np.min(positive_rates) * 0.5, np.max(positive_rates) * 3.0)
    ax1.set_ylabel(r'Rate/bin [s$^{-1}$]')
    ax1.set_title('Coincident SiPM spectrum - %s' % detector_name)
    ax1.grid(which='both', linestyle='--', alpha=0.5)
    ax1.legend(fontsize=DENSE_LEGEND_FONTSIZE, **LEGEND_STYLE)

    ax2.set_xscale('log')
    ax2.set_ylabel('Pull')
    ax2.set_xlabel('SiPM Peak Voltage [mV]')
    ax2.set_ylim(-6.0, 6.0)
    ax2.grid(which='both', linestyle='--', alpha=0.5)
    fig.tight_layout()
    print('Saving Figure to: ' + os.getcwd() + '/' + pdf_name)
    fig.savefig(pdf_name, format='pdf', transparent=True)
    plt.show()
    return fit_result

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

        legend_entry_count = len(labels) + len(legend_extra or [])
        legend_fontsize = (
            DENSE_LEGEND_FONTSIZE
            if legend_entry_count > 4 else LEGEND_FONTSIZE
        )
        # A dense rate legend otherwise hides the highest-rate series when the
        # caller requests the traditional upper-right position.
        dense_upper_right = legend_entry_count > 4 and loc in (1, 'upper right')
        legend_loc = 'center right' if dense_upper_right else loc
        legend_position = {'bbox_to_anchor': (1.0, 0.38)} if dense_upper_right else {}
        plt.legend(
            fontsize=legend_fontsize,
            loc=legend_loc,
            **legend_position,
            **LEGEND_STYLE,
        )
        
        plt.title(title,fontsize=fontsize+1)
        plt.tight_layout()
        if pdf_name != '':
            print('Saving Figure to: '+os.getcwd() +  '/'+pdf_name)
            plt.savefig(pdf_name, format='pdf',transparent =True)
        plt.show()


def draw_scalar_detector_curves(ax, time_minutes, series_by_detector,
                                equal_weight_mean, quantity_label,
                                value_scale=1.0):
    """Draw one scalar sensor curve per detector plus their equal-weight mean."""
    for detector_index, (detector_name, values) in enumerate(
            series_by_detector.items()):
        ax.plot(
            time_minutes, np.asarray(values) * value_scale,
            marker='o', markersize=2, linewidth=1.1, alpha=0.7,
            color=mycolors[detector_index % len(mycolors)],
            label='%s %s' % (quantity_label, detector_name),
        )
    if len(series_by_detector) > 1:
        ax.plot(
            time_minutes, np.asarray(equal_weight_mean) * value_scale,
            color='black', linewidth=3.0, zorder=10,
            label='Mean (equal detector weight)',
        )


def draw_vector_detector_curves(ax, time_minutes, component_series,
                                component_means, quantity_label):
    """Draw X/Y/Z for each detector and thick black component means."""
    component_styles = {'X': '-', 'Y': '--', 'Z': ':'}
    detector_names = list(next(iter(component_series.values())).keys())
    for detector_index, detector_name in enumerate(detector_names):
        detector_color = mycolors[detector_index % len(mycolors)]
        for component_name, series_by_detector in component_series.items():
            ax.plot(
                time_minutes, series_by_detector[detector_name],
                color=detector_color,
                linestyle=component_styles[component_name],
                linewidth=1.1, alpha=0.65,
                label='%s %s %s'
                      % (quantity_label, component_name, detector_name),
            )
    if len(detector_names) > 1:
        for component_name, mean_values in component_means.items():
            ax.plot(
                time_minutes, mean_values, color='black',
                linestyle=component_styles[component_name],
                linewidth=3.0, zorder=10,
                label='Mean %s (equal detector weight)' % component_name,
            )


def save_sensor_figure(draw_function, ylabel, pdf_name, ylim=None):
    """Apply common styling and save an environmental/motion sensor figure."""
    fig, ax = plt.subplots(figsize=(8, 5.5))
    draw_function(ax)
    ax.set_xlabel('Time [min]')
    ax.set_ylabel(ylabel)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.grid(which='both', linestyle='--', alpha=0.5)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles, labels,
        fontsize=DENSE_LEGEND_FONTSIZE,
        ncol=2 if len(labels) > 6 else 1,
        **LEGEND_STYLE,
    )
    fig.tight_layout()
    print('Saving Figure to: ' + os.getcwd() + '/' + pdf_name)
    fig.savefig(pdf_name, format='pdf', transparent=True)
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
    parser.add_argument(
        '--time-coincidence-mode', choices=('pairwise', 'all'),
        default='pairwise',
        help=(
            "With --coincidence-source=time, find independent detector pairs "
            "or require one event from every detector inside the same window "
            "(default: pairwise)"
        ),
    )
    parser.add_argument(
        '--delta-t-range-ms', type=float, default=20.0,
        help="Half-range of the cross-detector delta-t histogram in ms (default: 20)"
    )

    args = parser.parse_args()
    if (args.time_coincidence_mode == 'all'
            and args.coincidence_source != 'time'):
        parser.error(
            '--time-coincidence-mode=all requires --coincidence-source=time'
        )

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
        time_coincidence_mode=args.time_coincidence_mode,
    )
    strict_all_mode = (
        f1.coincidence_source == 'time'
        and f1.time_coincidence_mode == 'all'
    )
    coincident_event_label = (
        '%d-fold event rows' % f1.n_detector
        if strict_all_mode else 'Coincident'
    )

    
    # Plot the ADC values from the coincident and non-coincident events
    c = NPlot(
        data=[ f1.adc,f1.adc[~f1.select_coincident],f1.adc[f1.select_coincident]],
        weights=[f1.weights,f1.weights[~f1.select_coincident],f1.weights[f1.select_coincident]],
        colors=[mycolors[7], mycolors[3],mycolors[1]],
        labels=[r'All Events:  ' + str(f1.count_rate) + '+/-' + str(f1.count_rate_err) +' Hz',
                r'Non-Coincident:  ' + str(f1.count_rate_non_coincident) + '+/-' + str(f1.count_rate_err_non_coincident) +' Hz',
                coincident_event_label + ': ' + str(f1.count_rate_coincident)
                + '+/-' + str(f1.count_rate_err_coincident) +' Hz'],
        xmin=None, xmax=None, ymin=None, ymax=None,nbins=101,yscale='log',xscale='log',
        xlabel='Meausred ADC peak value [0-4095]',
        pdf_name=pdf_file_location+'/'+infile_name+'_ADC.pdf',title = '')
    
    # For time coincidences, do not merge the two detector measurements into a
    # single doubled spectrum. Use unique one-to-one pairs and retain one
    # conditional pulse-height spectrum per detector.
    sipm_data = [f1.sipm, f1.sipm[~f1.select_coincident]]
    sipm_weights = [f1.weights, f1.weights[~f1.select_coincident]]
    sipm_colors = [mycolors[7], mycolors[3]]
    sipm_labels = [
        r'All Events:  ' + str(f1.count_rate) + '+/-' + str(f1.count_rate_err) +' Hz',
        r'Non-Coincident:  ' + str(f1.count_rate_non_coincident) + '+/-'
        + str(f1.count_rate_err_non_coincident) +' Hz',
    ]
    if f1.coincidence_source == 'time' and f1.coincident_sipm_by_detector:
        detector_plot_colors = [mycolors[1], mycolors[10], mycolors[5], mycolors[8]]
        for detector_index, (detector_name, values) in enumerate(
                f1.coincident_sipm_by_detector.items()):
            sipm_data.append(values)
            sipm_weights.append(np.ones(len(values)) / f1.live_time_s)
            sipm_colors.append(
                detector_plot_colors[detector_index % len(detector_plot_colors)]
            )
            sipm_labels.append(
                'Coincident %s: %.4f Hz'
                % (detector_name, len(values) / f1.live_time_s)
            )
        landau_fit_indices = list(range(2, len(sipm_data)))
    else:
        sipm_data.append(f1.sipm[f1.select_coincident])
        sipm_weights.append(f1.weights[f1.select_coincident])
        sipm_colors.append(mycolors[1])
        sipm_labels.append(
            coincident_event_label + ': ' + str(f1.count_rate_coincident) + '+/-'
            + str(f1.count_rate_err_coincident) +' Hz'
        )
        landau_fit_indices = [2]

    # Plot the calculated SiPM peak voltages.
    c = NPlot(
        data=sipm_data,
        weights=sipm_weights,
        colors=sipm_colors,
        labels=sipm_labels,
        xmin=None, xmax=None, ymin=None, ymax=None,xscale='log',nbins = 51,
        xlabel='SiPM Peak Voltage [mV]', fit_landau=True,
        landau_data_index=landau_fit_indices, landau_fit_range=(25.0, 200.0),
        landau_initial_mpv=40.0,
        pdf_name=pdf_file_location+'/'+infile_name+'_SiPM_peak_voltage.pdf',title = '',)

    # Produce one dedicated, deduplicated Langau plot for each detector in the
    # time-coincidence sample.
    if f1.coincidence_source == 'time' and f1.coincident_sipm_by_detector:
        # Define one common logarithmic range from the union of all coincident
        # detector samples.  Passing it to every plot also gives them identical
        # histogram bin edges, not only matching displayed axis limits.
        shared_coincident_values = np.concatenate([
            np.asarray(detector_values, dtype=float)
            for detector_values in f1.coincident_sipm_by_detector.values()
            if len(detector_values) > 0
        ])
        shared_coincident_values = shared_coincident_values[
            np.isfinite(shared_coincident_values) & (shared_coincident_values > 0)
        ]
        shared_sipm_x_range = None
        if len(shared_coincident_values) > 0:
            shared_xmin = max(float(np.min(shared_coincident_values)) * 0.9, 0.1)
            shared_xmax = float(np.max(shared_coincident_values)) * 1.1
            if shared_xmax > shared_xmin:
                shared_sipm_x_range = (shared_xmin, shared_xmax)

        for detector_index, (detector_name, values) in enumerate(
                f1.coincident_sipm_by_detector.items()):
            safe_detector_name = ''.join(
                character if character.isalnum() or character in '-_'
                else '_'
                for character in detector_name
            )
            detector_rate = f1.detector_count_rates[detector_name]
            if strict_all_mode:
                # Every accidental N-fold group contributes exactly one event
                # to each detector's conditional spectrum.
                detector_accidental_rate = f1.accidental_group_rate
            else:
                detector_accidental_rate = sum(
                    2.0 * detector_rate * other_rate
                    * f1.coincidence_window_ms / 1000.0
                    for other_name, other_rate in f1.detector_count_rates.items()
                    if other_name != detector_name
                )
            plot_coincident_sipm_langauss(
                values=values,
                live_time_s=f1.live_time_s,
                detector_name=detector_name,
                fit_range=(25.0, 200.0),
                singles_values=f1.sipm[f1.detector_name == detector_name],
                accidental_rate_hz=detector_accidental_rate,
                x_range=shared_sipm_x_range,
                color=detector_plot_colors[
                    detector_index % len(detector_plot_colors)
                ],
                pdf_name=(
                    pdf_file_location+'/'+infile_name
                    +'_SiPM_coincident_'+safe_detector_name+'_Langau.pdf'
                ),
            )

    if f1.file_from_computer and f1.n_detector >= 2:
        plot_delta_t_histogram(
            f1.time_stamp_s, f1.detector_name,
            window_ms=args.coincidence_window_ms,
            max_abs_ms=args.delta_t_range_ms,
            pdf_name=pdf_file_location+'/'+infile_name+'_coincidence_delta_t.pdf',
        )
    
    
    rate_legend_extra = []
    if f1.coincidence_source == 'time':
        rate_legend_extra.extend(
            '%s rate: %.4f Hz' % (name, rate)
            for name, rate in f1.detector_count_rates.items()
        )
        if strict_all_mode:
            rate_legend_extra.extend([
                r'%d-fold $R_{acc}$ (span $\leq$ %.3g ms): %.4g Hz'
                % (f1.n_detector, f1.coincidence_window_ms,
                   f1.accidental_group_rate),
                'Corrected %d-fold rate: %.4f Hz'
                % (f1.n_detector, f1.corrected_group_rate),
            ])
        else:
            rate_legend_extra.extend([
                'Pair rate: %.4f Hz' % f1.coincidence_pair_rate,
                r'$R_{acc}$ ($\pm$%.3g ms): %.4f Hz'
                % (f1.coincidence_window_ms, f1.accidental_pair_rate),
                'Corrected pair rate: %.4f Hz' % f1.corrected_pair_rate,
            ])

    if strict_all_mode:
        plotted_coincidence_rate = f1.binned_coincidence_group_rate
        plotted_coincidence_rate_err = f1.binned_coincidence_group_rate_err
        plotted_coincidence_label = (
            '%d-fold groups: %.4f+/-%.4f Hz'
            % (f1.n_detector, f1.coincidence_group_rate,
               f1.coincidence_group_rate_err)
        )
    else:
        plotted_coincidence_rate = f1.binned_count_rate_coincident
        plotted_coincidence_rate_err = f1.binned_count_rate_err_coincident
        plotted_coincidence_label = (
            'Coincident:  ' + str(f1.count_rate_coincident) + '+/-'
            + str(f1.count_rate_err_coincident) + ' Hz'
        )

    # Plot rate as a function of time
    c = ratePlot(time = [f1.binned_time_m,f1.binned_time_m,f1.binned_time_m],
        count_rates = [f1.binned_count_rate,f1.binned_count_rate_non_coincident,plotted_coincidence_rate],
        count_rates_err = [f1.binned_count_rate_err,f1.binned_count_rate_err_non_coincident,plotted_coincidence_rate_err],
        colors=[mycolors[7], mycolors[3], mycolors[1]],
        labels=[r'All Events: ' + str(f1.count_rate) + '+/-' + str(f1.count_rate_err) +' Hz', 
                r'Non-Coincident:  ' + str(f1.count_rate_non_coincident) + '+/-' + str(f1.count_rate_err_non_coincident) +' Hz',
                plotted_coincidence_label],
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
        save_sensor_figure(
            lambda ax: draw_scalar_detector_curves(
                ax, f1.binned_time_m,
                f1.binned_pressure_by_detector, f1.binned_pressure,
                'Pressure',
            ),
            ylabel='Pressure [Pa]',
            pdf_name=pdf_file_location+'/'+infile_name+'_pressure.pdf',
        )

        save_sensor_figure(
            lambda ax: draw_scalar_detector_curves(
                ax, f1.binned_time_m,
                f1.binned_temperature_by_detector, f1.binned_temperature,
                'Temperature',
            ),
            ylabel=r'Temperature [$^{\circ}$C]',
            pdf_name=pdf_file_location+'/'+infile_name+'_temperature.pdf',
        )

    if f1.has_MPU6050:
        save_sensor_figure(
            lambda ax: draw_vector_detector_curves(
                ax, f1.binned_time_m,
                {
                    'X': f1.binned_accel_x_by_detector,
                    'Y': f1.binned_accel_y_by_detector,
                    'Z': f1.binned_accel_z_by_detector,
                },
                {
                    'X': f1.binned_accel_x,
                    'Y': f1.binned_accel_y,
                    'Z': f1.binned_accel_z,
                },
                'Acceleration',
            ),
            ylabel='Linear acceleration [g]',
            pdf_name=pdf_file_location+'/'+infile_name+'_accel.pdf',
            ylim=(-1.3, 1.3),
        )

        save_sensor_figure(
            lambda ax: draw_vector_detector_curves(
                ax, f1.binned_time_m,
                {
                    'X': f1.binned_gyro_x_by_detector,
                    'Y': f1.binned_gyro_y_by_detector,
                    'Z': f1.binned_gyro_z_by_detector,
                },
                {
                    'X': f1.binned_gyro_x,
                    'Y': f1.binned_gyro_y,
                    'Z': f1.binned_gyro_z,
                },
                'Gyro',
            ),
            ylabel='Angular velocity [deg/s]',
            pdf_name=pdf_file_location+'/'+infile_name+'_gyro.pdf',
            ylim=(-100, 100),
        )

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
    axes[0].plot(
        t, plotted_coincidence_rate, color=mycolors[1],
        label=('%d-fold groups' % f1.n_detector
               if strict_all_mode else 'Coincident')
    )
    axes[0].set_ylabel('Rate [Hz]')
    axes[0].legend(loc='upper right', fontsize=6)
    axes[0].grid(True, which='both', linestyle='--', alpha=0.5)


    axis_index = 1
    if f1.has_BMP280:
        draw_scalar_detector_curves(
            axes[axis_index], t,
            f1.binned_pressure_by_detector, f1.binned_pressure,
            'Pressure', value_scale=1.0/1000.0,
        )
        axes[axis_index].set_ylabel('Pressure [kPa]')
        axes[axis_index].legend(fontsize=6, **LEGEND_STYLE)
        axis_index += 1

        draw_scalar_detector_curves(
            axes[axis_index], t,
            f1.binned_temperature_by_detector, f1.binned_temperature,
            'Temperature',
        )
        axes[axis_index].set_ylabel('Temperature [°C]')
        axes[axis_index].legend(fontsize=6, **LEGEND_STYLE)
        axis_index += 1

    # 5) Acceleration (if present)
    if f1.has_MPU6050:
        draw_vector_detector_curves(
            axes[axis_index], t,
            {
                'X': f1.binned_accel_x_by_detector,
                'Y': f1.binned_accel_y_by_detector,
                'Z': f1.binned_accel_z_by_detector,
            },
            {
                'X': f1.binned_accel_x,
                'Y': f1.binned_accel_y,
                'Z': f1.binned_accel_z,
            },
            'Acceleration',
        )
        axes[axis_index].set_ylabel('Accel [g]')
        axes[axis_index].set_ylim(-1.3, 1.3)
        axes[axis_index].legend(
            loc='upper right', fontsize=6, ncol=2, **LEGEND_STYLE
        )
        axis_index += 1

    # 6) Angular velocity (if present)
    if f1.has_MPU6050:
        draw_vector_detector_curves(
            axes[axis_index], t,
            {
                'X': f1.binned_gyro_x_by_detector,
                'Y': f1.binned_gyro_y_by_detector,
                'Z': f1.binned_gyro_z_by_detector,
            },
            {
                'X': f1.binned_gyro_x,
                'Y': f1.binned_gyro_y,
                'Z': f1.binned_gyro_z,
            },
            'Gyro',
        )
        axes[axis_index].set_ylabel('Gyro [°/s]')
        axes[axis_index].set_ylim(-100, 100)
        axes[axis_index].legend(
            loc='upper right', fontsize=6, ncol=2, **LEGEND_STYLE
        )

    for axis in axes:
        axis.grid(True, which='both', linestyle='--', alpha=0.5)

    # common x-label
    axes[-1].set_xlabel('Time [min]')

    #plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
