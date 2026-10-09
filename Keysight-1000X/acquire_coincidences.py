#!/usr/bin/env python3
"""Acquire two-channel muon coincidences from a Keysight InfiniiVision 1000X (e.g. DSOX1202A).

Every trigger is read out in binary over VISA and written in the instrument's own CSV export
layout, so analyze_pulses.py reads the files unchanged. events.csv logs PC time, live time and
peak amplitudes of every trigger, rejected ones included; session_*.json keeps the settings.
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
import time
import warnings
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyvisa


RUN_BIT = 8  # :OPERegister:CONDition? bit 3, set while an acquisition is running
MAX_SAMPLE_RATE = 2e9  # DSOX1202A, with both channels on (measured)
LOG_NAME = 'events.csv'
LOG_FIELDS = ['index', 'file', 'accepted', 'pc_time', 'utc_time', 'live_s', 'readout_s',
              'amplitude1_mV', 'amplitude2_mV', 'adc_limit1', 'adc_limit2', 'points', 'dt_ns']
SETTINGS = [f':CHANnel{n}:{item}?' for n in (1, 2)
            for item in ('DISPlay', 'SCALe', 'OFFSet', 'COUPling', 'BWLimit', 'PROBe')] + [
    ':TIMebase:SCALe?', ':TIMebase:POSition?', ':TIMebase:REFerence?',
    ':ACQuire:TYPE?', ':ACQuire:SRATe?', ':ACQuire:POINts?',
    ':TRIGger:MODE?', ':TRIGger:SWEep?', ':TRIGger:EDGE:SOURce?', ':TRIGger:EDGE:SLOPe?',
    ':TRIGger:EDGE:LEVel? CHANnel1', ':TRIGger:EDGE:LEVel? CHANnel2',
    ':TRIGger:PATTern?', ':TRIGger:PATTern:QUALifier?', ':TRIGger:HOLDoff?',
    ':WAVeform:POINts:MODE?', ':WAVeform:POINts?']
# Communication failures worth reopening the session for during a long run.
RECOVERABLE = (pyvisa.errors.Error, OSError, ValueError, subprocess.SubprocessError)


class ConfigurationError(RuntimeError):
    """The instrument rejected a setting: reconnecting cannot fix it."""


class ReadoutError(ValueError):
    """Inconsistent instrument state or data; the session is reopened."""


def per_channel(parser, name, values):
    if values is None:
        return None
    if len(values) not in (1, 2):
        parser.error(f'{name} takes one value (both channels) or two (C1 C2).')
    return list(values) * (3 - len(values))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help='Run directory; an existing run is resumed without overwriting events')
    connection = parser.add_mutually_exclusive_group()
    connection.add_argument('--usb', action='store_true', help='Use the first Keysight oscilloscope found on USB')
    connection.add_argument('--resource', help='Explicit VISA resource, e.g. USB0::10893::903::<serial>::0::INSTR or TCPIP::192.168.1.20::INSTR')
    connection.add_argument('--ssh-bastion', help='SSH destination (user@host or ~/.ssh/config alias) to tunnel the SCPI socket through')
    parser.add_argument('--host', help='Oscilloscope IP/hostname: VXI-11 directly, or as seen from the bastion')
    parser.add_argument('--instrument-port', type=int, default=5025, help='SCPI socket port for --ssh-bastion (default: 5025)')
    parser.add_argument('--local-port', type=int, default=15026, help='Loopback port of the SSH tunnel (default: 15026)')
    parser.add_argument('--visa-library', default='@py', help="PyVISA backend: '@py' (pyvisa-py, default) or '@ivi' for vendor VISA")
    parser.add_argument('--timeout-ms', type=int, default=10000, help='Timeout for each VISA operation (default: 10000)')

    # Fixed for the 3 m run (positive pulses, onset near the trigger, ~190 ns decay).
    # 25 mV/div with the baseline half a division above the screen bottom: the ADC
    # stays valid ~5.1 div above center, so up to ~215 mV; at 20 mV/div one event
    # in 26 saturated at 173 mV. 100 ns/div centered at +300 ns: -200..+800 ns,
    # ~190 ns of baseline and the tail below 5% of the peak, at 2 GSa/s.
    setup = parser.add_argument_group('instrument settings (defaults: fixed values for the 3 m run)')
    setup.add_argument('--scale-mv', type=float, nargs='+', default=[25.], metavar='MV',
                       help='Vertical scale per division, C1 [C2] (default: 25)')
    setup.add_argument('--offset-mv', type=float, nargs='+', default=[87.5], metavar='MV',
                       help='Voltage at screen center, C1 [C2] (default: 87.5, i.e. -12.5..+187.5 mV on screen)')
    setup.add_argument('--timebase-ns', type=float, default=100., help='Horizontal scale per division (default: 100)')
    setup.add_argument('--delay-ns', type=float, default=300.,
                       help='Time of screen center after the trigger, reference CENTer (default: 300)')
    setup.add_argument('--trigger', choices=('keep', 'pattern', 'edge'), default='keep',
                       help='pattern: hardware coincidence, both channels beyond their level; '
                            'edge: single channel, coincidence only in software; keep: front panel (default)')
    setup.add_argument('--trigger-level-mv', type=float, nargs='+', metavar='MV', help='Trigger level, C1 [C2]')
    setup.add_argument('--edge-source', type=int, choices=(1, 2), default=1, help='Channel for --trigger edge (default: 1)')
    setup.add_argument('--polarity', choices=('positive', 'negative'), default='positive',
                       help='Pulse polarity, for trigger and amplitudes (default: positive)')
    setup.add_argument('--points-mode', choices=('raw', 'normal', 'maximum'), default='raw',
                       help='raw: full acquisition record at the real sample rate (default)')
    setup.add_argument('--points', type=int,
                       help='Requested points (default: all available); fewer than the raw record decimates it')

    run = parser.add_argument_group('selection and run control')
    run.add_argument('--window-ns', type=float, nargs=2, default=[-250., 800.], metavar=('START', 'STOP'),
                     help='Time window saved around the trigger (default: -250 800)')
    run.add_argument('--full-record', action='store_true', help='Save the whole transferred record instead of --window-ns')
    run.add_argument('--baseline-end-ns', type=float, default=-20.,
                     help='Baseline for the logged amplitudes ends here (default: -20); '
                          'with fewer than 20 samples before it, the first 15%% of the record is used')
    run.add_argument('--min-amplitude-mv', type=float, nargs='+', metavar='MV',
                     help='Software coincidence: peak above baseline required in C1 [C2]; default saves every trigger')
    run.add_argument('--save-rejected', action='store_true', help='Also save triggers failing the cut, in rejected/')
    run.add_argument('--max-events', type=int, help='Stop after this many saved coincidences')
    run.add_argument('--duration-h', type=float, help='Stop after this many hours')
    run.add_argument('--poll-s', type=float, default=.25, help='Polling interval while armed (default: 0.25)')
    run.add_argument('--status-s', type=float, default=300.,
                     help='Status line and session file update, also while waiting (default: 300)')
    run.add_argument('--max-reconnects', type=int,
                     help='Consecutive failed reconnections before giving up (default: unlimited)')
    run.add_argument('--check', action='store_true', help='Apply the requested settings, print the readback and exit')
    args = parser.parse_args(argv)

    if not (args.usb or args.resource or args.host):
        parser.error('Give --usb, --resource or --host (with or without --ssh-bastion).')
    if args.ssh_bastion and (args.ssh_bastion.startswith('-') or any(c.isspace() for c in args.ssh_bastion)):
        parser.error('--ssh-bastion must be user@host or an SSH alias, not a shell command.')
    if not all(1 <= port <= 65535 for port in (args.local_port, args.instrument_port)):
        parser.error('TCP ports must be between 1 and 65535.')
    for name in ('timeout_ms', 'points', 'poll_s', 'status_s', 'timebase_ns', 'max_events', 'duration_h'):
        value = getattr(args, name)
        if value is not None and not value > 0:
            parser.error(f'--{name.replace("_", "-")} must be positive.')
    args.scale_mv = per_channel(parser, '--scale-mv', args.scale_mv)
    args.offset_mv = per_channel(parser, '--offset-mv', args.offset_mv)
    args.trigger_level_mv = per_channel(parser, '--trigger-level-mv', args.trigger_level_mv)
    args.min_amplitude_mv = per_channel(parser, '--min-amplitude-mv', args.min_amplitude_mv)
    if args.scale_mv and min(args.scale_mv) <= 0:
        parser.error('--scale-mv must be positive.')
    if args.trigger != 'keep' and args.trigger_level_mv is None:
        parser.error(f'--trigger {args.trigger} requires --trigger-level-mv.')
    if args.window_ns[0] >= args.window_ns[1]:
        parser.error('--window-ns START must be before STOP.')
    return args


@contextmanager
def instrument_resource(args):
    """Yield the VISA resource; own and close only the SSH tunnel opened here."""
    if args.usb:
        yield find_usb(args.visa_library)
        return
    if not args.ssh_bastion:
        yield args.resource or f'TCPIP::{args.host}::INSTR'
        return
    # OpenSSH handles keys, agent, MFA and ~/.ssh/config; the private control
    # socket closes this tunnel without touching other SSH sessions.
    with tempfile.TemporaryDirectory(prefix='cwssh-', dir='/tmp') as directory:
        control = str(Path(directory) / 'control')
        forward = f'127.0.0.1:{args.local_port}:{args.host}:{args.instrument_port}'
        command = ['ssh', '-M', '-S', control, '-fNT', '-o', 'ControlPersist=no',
                   '-o', 'ExitOnForwardFailure=yes', '-o', 'ConnectTimeout=15',
                   '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=3',
                   '-L', forward, '--', args.ssh_bastion]
        print(f'Opening SSH tunnel via {args.ssh_bastion}: {forward}', flush=True)
        try:
            subprocess.run(command, check=True)
            yield f'TCPIP::127.0.0.1::{args.local_port}::SOCKET'
        finally:
            subprocess.run(['ssh', '-S', control, '-O', 'exit', '--', args.ssh_bastion],
                           check=False, capture_output=True, timeout=10)


def find_usb(visa_library):
    manager = pyvisa.ResourceManager(visa_library)
    try:
        # pyvisa-py also scans the network and warns that psutil/zeroconf would
        # widen that scan: irrelevant for USB.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='TCPIP', category=UserWarning)
            resources = manager.list_resources('USB?*INSTR')
        # Keysight USB vendor ID 0x2A8D (10893); 0x0957 (2391) on older Agilent firmware.
        found = [r for r in resources if r.split('::')[1].upper() in ('10893', '0X2A8D', '2391', '0X0957')]
    finally:
        manager.close()
    if not found:
        raise ConfigurationError('No Keysight oscilloscope found by PyVISA on USB: ' + usb_diagnosis())
    return found[0]


def usb_diagnosis():
    """Why pyvisa-py lists no Keysight device: it skips USB without pyusb/libusb and
    silently drops devices it cannot open (Linux without a udev rule)."""
    try:
        import usb.core
    except ImportError:
        return 'pyusb is missing in this Python environment: pip install pyusb libusb-package'
    vendors = (0x2A8D, 0x0957)
    try:
        try:
            import libusb_package  # the same libusb lookup order as pyvisa-py
            devices = list(libusb_package.find(find_all=True, custom_match=lambda d: d.idVendor in vendors))
        except ImportError:
            devices = list(usb.core.find(find_all=True, custom_match=lambda d: d.idVendor in vendors))
    except usb.core.NoBackendError:
        return 'libusb not found: pip install libusb-package'
    if not devices:
        return ('no Keysight device on the USB bus: is the cable in the rear USB device port '
                '(square type-B connector) and the oscilloscope on?')
    device = devices[0]
    node = f'/dev/bus/usb/{device.bus:03d}/{device.address:03d}'
    if os.path.exists(node) and not os.access(node, os.R_OK | os.W_OK):
        return (f'the device is on the bus ({device.idVendor:04x}:{device.idProduct:04x}) but {node} '
                'is not writable by this user: add the udev rule from README_acquisition.md and replug the cable')
    try:
        device.serial_number
    except (usb.core.USBError, ValueError) as error:
        return (f'the device is on the bus but cannot be opened ({error}); on Linux the usbtmc kernel '
                'driver may hold it (see README_acquisition.md), or another program is using it')
    return 'the device is on the bus and readable, but pyvisa-py does not list it: run python -m pyvisa info'


@contextmanager
def open_scope(args):
    with instrument_resource(args) as resource:
        print(f'Connecting to {resource} ...', flush=True)
        manager = pyvisa.ResourceManager(args.visa_library)
        scope = None
        try:
            # USB: start from a clean USBTMC state whatever a previous session left.
            usbtmc_clear(resource)
            scope = manager.open_resource(resource)
            scope.timeout = args.timeout_ms
            scope.read_termination = scope.write_termination = '\n'
            if not resource.upper().startswith('USB'):
                drain(scope)  # reading with nothing pending is not harmless on USBTMC
            error_queue(scope)  # e.g. -310 left by an interrupted transfer
            yield scope
        finally:
            if scope is not None:
                scope.close()
            manager.close()


def drain(scope):
    """Discard a reply left unread by an interrupted session (USB has no device clear)."""
    timeout, scope.timeout = scope.timeout, 200
    try:
        while scope.read_raw():
            pass
    except pyvisa.errors.VisaIOError:
        pass
    finally:
        scope.timeout = timeout


def as_int(reply):
    return int(float(reply.strip()))


def error_queue(scope, limit=50):
    errors = []
    for _ in range(limit):
        reply = scope.query(':SYSTem:ERRor?').strip()
        if as_int(reply.split(',', 1)[0]) == 0:
            break
        errors.append(reply)
    return errors


def configure(scope, args):
    """Apply the requested settings plus those the readout loop depends on."""
    error_queue(scope)  # discard errors left by front-panel use or earlier sessions
    positive = args.polarity == 'positive'
    # Auto sweep would fire without any pulse; averaging would merge different muons.
    commands = [':TRIGger:SWEep NORMal', ':ACQuire:TYPE NORMal',
                ':CHANnel1:DISPlay 1', ':CHANnel2:DISPlay 1']
    for n, scale in enumerate(args.scale_mv or [], 1):
        commands.append(f':CHANnel{n}:SCALe {scale/1e3:.6g}')
    for n, offset in enumerate(args.offset_mv or [], 1):
        commands.append(f':CHANnel{n}:OFFSet {offset/1e3:.6g}')
    if args.timebase_ns is not None:
        commands += [':TIMebase:MODE MAIN', f':TIMebase:SCALe {args.timebase_ns/1e9:.6g}']
    if args.delay_ns is not None:
        commands += [':TIMebase:REFerence CENTer', f':TIMebase:POSition {args.delay_ns/1e9:.6g}']
    levels = args.trigger_level_mv
    if args.trigger == 'edge':
        source = f'CHANnel{args.edge_source}'
        commands += [':TRIGger:MODE EDGE', f':TRIGger:EDGE:SOURce {source}', ':TRIGger:EDGE:COUPling DC',
                     f':TRIGger:EDGE:SLOPe {"POSitive" if positive else "NEGative"}',
                     f':TRIGger:EDGE:LEVel {levels[args.edge_source-1]/1e3:.6g},{source}']
    elif args.trigger == 'pattern':
        # Fires when both channels are beyond their own level at the same time:
        # the ~170 ns wide pulses overlap despite the ~10 ns flight time.
        commands += [':TRIGger:MODE PATTern', ':TRIGger:PATTern:FORMat ASCii',
                     pattern_command(scope, '1' if positive else '0'), ':TRIGger:PATTern:QUALifier ENTered']
        commands += [f':TRIGger:EDGE:LEVel {level/1e3:.6g},CHANnel{n}' for n, level in enumerate(levels, 1)]
    commands += [':WAVeform:FORMat BYTE', ':WAVeform:UNSigned 1',
                 f':WAVeform:POINts:MODE {args.points_mode.upper()}']
    if args.points:  # by default the whole record; a larger number is an error, not a clamp
        commands.append(f':WAVeform:POINts {args.points}')
    send(scope, commands)


def send(scope, commands):
    for command in commands:
        scope.write(command)
        # One check per command, so the rejected one can be named.
        errors = error_queue(scope)
        if errors:
            hint = ' Set the trigger from the front panel and use --trigger keep.' if ':TRIG' in command else ''
            raise ConfigurationError(f'Instrument rejected {command!r}: {"; ".join(errors)}.{hint}')


def pattern_command(scope, state):
    """Both channels at `state`, every other pattern input (e.g. EXT) don't care."""
    scope.write(':TRIGger:PATTern:FORMat ASCii')
    current = scope.query(':TRIGger:PATTern?').split(',')[0].strip().strip('"')
    if len(current) < 2:
        raise ConfigurationError(f'Unexpected pattern readback: {current!r}')
    # The DSOX1202A reports "C1 C2 EXT", e.g. "11X".
    return f':TRIGger:PATTern "{state * 2}{"X" * (len(current) - 2)}"'


def readback(scope):
    """Query the settings for the session record; unsupported queries become None."""
    values, timeout = {}, scope.timeout
    scope.timeout = min(timeout, 3000)
    try:
        for query in SETTINGS:
            try:
                values[query] = scope.query(query).strip()
            except pyvisa.errors.VisaIOError:
                values[query] = None
    finally:
        scope.timeout = timeout
    error_queue(scope)  # errors from unsupported queries are not a reason to stop
    return values


def setting_warnings(settings):
    notes = []
    for n in (1, 2):
        probe = settings.get(f':CHANnel{n}:PROBe?')
        if probe is not None and abs(float(probe) - 1) > 1e-6:
            notes.append(f'C{n} probe factor is {float(probe):g}: volts include it; use 1 for a direct cable.')
        if (settings.get(f':CHANnel{n}:COUPling?') or '').upper().startswith('AC'):
            notes.append(f'C{n} is AC coupled.')
        if (settings.get(f':CHANnel{n}:BWLimit?') or '0').strip() in ('1', 'ON'):
            notes.append(f'C{n} has the 20 MHz bandwidth limit on: rise times are slowed.')
    if (settings.get(':TRIGger:SWEep?') or 'NORM').upper()[:4] != 'NORM':
        notes.append('Trigger sweep is not NORMal.')
    rate = settings.get(':ACQuire:SRATe?')
    if rate is not None and float(rate) < MAX_SAMPLE_RATE:
        notes.append(f'Sample rate {float(rate)/1e9:g} GSa/s, below the {MAX_SAMPLE_RATE/1e9:g} GSa/s maximum: '
                     'shorten --timebase-ns.')
    return notes


def wait_for_trigger(scope, poll_s, deadline=None, arm_timeout_s=5., on_idle=None):
    """Arm one acquisition and wait for it; return the armed (live) time, None if none.

    on_idle(armed_s) is called at every poll, e.g. for a heartbeat during hours-long waits.
    """
    scope.query(':AER?')  # arm and trigger event registers clear on read
    scope.query(':TER?')
    scope.write(':SINGle')
    start = time.monotonic()
    # The run bit may still show the previous state: wait for the new arm first.
    while not as_int(scope.query(':AER?')):
        if time.monotonic() - start > arm_timeout_s:
            raise ReadoutError('The oscilloscope did not arm after :SINGle')
        time.sleep(.01)
    armed = time.monotonic()
    while as_int(scope.query(':OPERegister:CONDition?')) & RUN_BIT:
        if deadline is not None and time.monotonic() >= deadline:
            scope.write(':STOP')
            return None
        if on_idle is not None:
            on_idle(time.monotonic() - armed)
        time.sleep(poll_s)
    live = time.monotonic() - armed
    # Stopped from the front panel, for example, rather than triggered.
    return live if as_int(scope.query(':TER?')) else None


def read_channel(scope, channel):
    scope.write(f':WAVeform:SOURce CHANnel{channel}')
    preamble = [float(value) for value in scope.query(':WAVeform:PREamble?').split(',')]
    if len(preamble) != 10 or preamble[0] != 0:
        raise ReadoutError(f'Unexpected preamble for C{channel}: {preamble}')
    codes = np.asarray(scope.query_binary_values(':WAVeform:DATA?', datatype='B', container=np.array,
                                                 expect_termination=True), dtype=float)
    if len(codes) < 30:
        raise ReadoutError(f'Only {len(codes)} points from C{channel}')
    xinc, xorigin, xref, yinc, yorigin, yref = preamble[4:]
    t = (np.arange(len(codes)) - xref) * xinc + xorigin
    return t * 1e9, ((codes - yref) * yinc + yorigin) * 1e3, codes


def read_event(scope, window_ns=None):
    """Both channels of the stopped acquisition, cropped around the trigger (ns, mV)."""
    (t, y1, c1), (t2, y2, c2) = read_channel(scope, 1), read_channel(scope, 2)
    dt = t[1] - t[0]
    if len(t) != len(t2) or abs(t2[0] - t[0]) > 1e-3 * dt or abs((t2[1] - t2[0]) - dt) > 1e-6 * dt:
        raise ReadoutError('C1 and C2 records have different time axes')
    keep = np.ones(len(t), bool) if window_ns is None else (t >= window_ns[0]) & (t <= window_ns[1])
    if keep.sum() < 30:
        raise ReadoutError(f'Fewer than 30 samples in the window: record spans {t[0]:.1f} to {t[-1]:.1f} ns')
    # Byte codes 0/255 are holes or clipped samples beyond the ADC range.
    limits = [int(np.sum((c[keep] <= 0) | (c[keep] >= 255))) for c in (c1, c2)]
    return t[keep], np.column_stack([y1[keep], y2[keep]]), limits


def peak_amplitudes(t, y, polarity, baseline_end=None):
    """Peak beyond the pre-trigger baseline; analyze_pulses needs the same --baseline-end-ns."""
    # A short pre-trigger window puts the first 15% of the record on the rising edge.
    mask = t <= baseline_end if baseline_end is not None else np.zeros(len(t), bool)
    if mask.sum() < 20:
        mask = t <= t[0] + .15 * np.ptp(t)
    base = y[mask].mean(axis=0)
    sign = 1 if polarity == 'positive' else -1
    return np.max(sign * (y - base), axis=0)


def write_keysight_csv(path, t_ns, y_mV):
    """Same layout as a front-panel CSV export, published only once complete."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f'.{path.name}.part')
    with partial.open('w', newline='') as stream:
        stream.write('x-axis,1,2\nsecond,Volt,Volt\n')
        np.savetxt(stream, np.column_stack([t_ns * 1e-9, y_mV * 1e-3]),
                   fmt=['%+.9E', '%+.6E', '%+.6E'], delimiter=',')
    os.replace(partial, path)


class EventLog:
    """Append-only events.csv; numbering resumes after the last logged event."""

    def __init__(self, directory):
        self.directory = directory
        self.path = directory / LOG_NAME
        self.next_index = 0
        if self.path.exists() and self.path.stat().st_size:
            with self.path.open(newline='') as stream:
                reader = csv.DictReader(stream)
                if reader.fieldnames != LOG_FIELDS:
                    raise ValueError(f'{self.path} has unexpected columns: {reader.fieldnames}')
                self.next_index = max((int(row['index']) for row in reader), default=-1) + 1
            new = False
        else:
            new = True
        self.stream = self.path.open('a', newline='')
        self.writer = csv.DictWriter(self.stream, LOG_FIELDS)
        if new:
            self.writer.writeheader()
            self.stream.flush()

    def claim_name(self):
        """Never overwrite an event file, even if the log was edited or lost."""
        while any((self.directory / sub / f'scope_{self.next_index:06d}.csv').exists()
                  for sub in ('.', 'rejected')):
            self.next_index += 1
        index, self.next_index = self.next_index, self.next_index + 1
        return index, f'scope_{index:06d}.csv'

    def append(self, row):
        self.writer.writerow(row)
        self.stream.flush()

    def close(self):
        self.stream.close()


def usbtmc_clear(resource):
    """USBTMC device clear: INITIATE_CLEAR, CHECK_CLEAR_STATUS, then clear the bulk halts.

    pyvisa-py has no device clear on USB, and a transfer left half done (timeout,
    Ctrl-C) leaves the DSOX1202A refusing every later command with write timeouts;
    a USB port reset does not recover it, this does (verified 9 Oct 2026).
    The acquisition in progress on the oscilloscope is not touched.
    """
    parts = (resource or '').split('::')
    if len(parts) < 4 or not parts[0].upper().startswith('USB'):
        return
    try:
        import usb.core
        import usb.util
        vendor, product, serial = int(parts[1], 0), int(parts[2], 0), parts[3]
        try:
            import libusb_package
            devices = list(libusb_package.find(find_all=True, idVendor=vendor, idProduct=product))
        except ImportError:
            devices = list(usb.core.find(find_all=True, idVendor=vendor, idProduct=product))
        for device in devices:
            try:
                if device.serial_number != serial:
                    continue
                interface = next(i for i in device.get_active_configuration()
                                 if i.bInterfaceClass == 0xFE and i.bInterfaceSubClass == 3)
                number = interface.bInterfaceNumber
                usb.util.claim_interface(device, number)
                # Class requests, recipient interface: 5 INITIATE_CLEAR, 6 CHECK_CLEAR_STATUS (2 = pending).
                device.ctrl_transfer(0xA1, 5, 0, number, 1, timeout=5000)
                for _ in range(100):
                    if device.ctrl_transfer(0xA1, 6, 0, number, 2, timeout=5000)[0] != 2:
                        break
                    time.sleep(.05)
                for endpoint in interface:
                    if usb.util.endpoint_type(endpoint.bmAttributes) == usb.util.ENDPOINT_TYPE_BULK:
                        device.clear_halt(endpoint.bEndpointAddress)
                usb.util.release_interface(device, number)
            finally:
                usb.util.dispose_resources(device)
    except Exception as error:  # best effort: opening the session is attempted anyway
        print(f'USBTMC clear not possible: {error!r}', flush=True)


def pending_trigger(scope, unread=False):
    """True if the stopped oscilloscope holds a triggered acquisition not yet saved:
    one that arrived while the connection was down, or whose readout failed (unread)."""
    if as_int(scope.query(':OPERegister:CONDition?')) & RUN_BIT:
        return False
    triggered = bool(as_int(scope.query(':TER?')))  # read anyway: it clears on read
    return unread or triggered


def prepare(scope, args, session):
    idn = scope.query('*IDN?').strip()
    print(f'Connected: {idn}', flush=True)
    configure(scope, args)
    settings = readback(scope)
    if 'settings' not in session:
        session.update(idn=idn, settings=settings, warnings=setting_warnings(settings))
        for warning in session['warnings']:
            print(f'WARNING: {warning}', flush=True)
    elif settings != session['settings']:
        session.setdefault('settings_after_reconnect', []).append(
            {'time': datetime.now().astimezone().isoformat(timespec='seconds'), 'settings': settings})
        print('WARNING: settings differ from the start of the session (see session JSON).', flush=True)


def save_setup(scope, path):
    """Binary front-panel setup, restorable with :SYSTem:SETup; optional."""
    try:
        path.write_bytes(scope.query_binary_values(':SYSTem:SETup?', datatype='B', container=bytes))
    except RECOVERABLE as error:
        print(f'Setup not saved ({error})', flush=True)
        error_queue(scope)


def acquire(args, connect=open_scope):
    output = args.output_dir.expanduser()
    output.mkdir(parents=True, exist_ok=True)
    log = EventLog(output)
    started = datetime.now().astimezone()
    stem = output / f'session_{started:%Y%m%d_%H%M%S}'
    session = {'started': started.isoformat(timespec='seconds'), 'command': sys.argv,
               'args': vars(args), 'first_index': log.next_index}
    counts = {'triggers': 0, 'accepted': 0, 'rejected': 0, 'reconnects': 0, 'recovered': 0}
    # unread: live time of a triggered acquisition not yet saved (None if there is none).
    state = {'live': 0., 'armed': 0., 'last_status': time.monotonic(), 'last_event': None, 'unread': None}
    deadline = time.monotonic() + args.duration_h * 3600 if args.duration_h else None
    failures = 0
    connected_once = False
    resource = None

    def finished():
        return ((args.max_events is not None and counts['accepted'] >= args.max_events) or
                (deadline is not None and time.monotonic() >= deadline))

    def write_session(stopped=None):
        session.update(counts=counts, live_time_s=state['live'] + state['armed'],
                       last_index=log.next_index - 1, last_event=state['last_event'],
                       updated=datetime.now().astimezone().isoformat(timespec='seconds'))
        if stopped:
            session['stopped'] = stopped
        stem.with_suffix('.json').write_text(json.dumps(session, indent=2, default=str) + '\n')

    def status(armed=0.):
        """Heartbeat, also during hours without triggers."""
        state['armed'] = armed
        if time.monotonic() - state['last_status'] < args.status_s:
            return
        state['last_status'] = time.monotonic()
        live = state['live'] + armed
        rate = counts['accepted'] / live * 3600 if live else float('nan')
        print(f'--- {datetime.now():%Y-%m-%d %H:%M:%S} alive: {counts["triggers"]} triggers, '
              f'{counts["accepted"]} saved, live {live/3600:.2f} h, {rate:.2f} coincidences/h, '
              f'last event {state["last_event"] or "none"} ---', flush=True)
        write_session()

    def record(scope, armed):
        """Read, select, save and log the stopped acquisition; armed=None if unknown."""
        done = time.monotonic()
        stamp = datetime.now().astimezone()
        window = None if args.full_record else args.window_ns
        t, y, limits = read_event(scope, window)
        amplitudes = peak_amplitudes(t, y, args.polarity, args.baseline_end_ns)
        accepted = args.min_amplitude_mv is None or all(
            a >= cut for a, cut in zip(amplitudes, args.min_amplitude_mv))
        index, name = log.claim_name()
        relative = name if accepted else f'rejected/{name}' if args.save_rejected else ''
        if relative:
            write_keysight_csv(output / relative, t, y)
        counts['triggers'] += 1
        counts['accepted' if accepted else 'rejected'] += 1
        state['last_event'] = stamp.isoformat(timespec='seconds')
        log.append({'index': index, 'file': relative, 'accepted': int(accepted),
                    'pc_time': stamp.replace(tzinfo=None).isoformat(timespec='milliseconds'),
                    'utc_time': stamp.astimezone(timezone.utc).isoformat(timespec='milliseconds'),
                    'live_s': '' if armed is None else f'{armed:.3f}',
                    'readout_s': f'{time.monotonic()-done:.3f}',
                    'amplitude1_mV': f'{amplitudes[0]:.2f}', 'amplitude2_mV': f'{amplitudes[1]:.2f}',
                    'adc_limit1': limits[0], 'adc_limit2': limits[1],
                    'points': len(t), 'dt_ns': f'{t[1]-t[0]:.6g}'})
        state['unread'] = None  # logged: a later failure must not save it twice
        live = 'unknown (recovered)' if armed is None else f'{armed:8.1f} s'
        print(f'{stamp:%Y-%m-%d %H:%M:%S} #{index} C1 {amplitudes[0]:6.1f} mV  '
              f'C2 {amplitudes[1]:6.1f} mV  live {live}  {relative or "rejected, not saved"}'
              f'{"  ADC LIMIT" if any(limits) else ""}', flush=True)
        write_session()

    def log_error(error):
        with (output / 'errors.log').open('a') as stream:
            stream.write(f'{datetime.now().astimezone().isoformat(timespec="seconds")} {error!r}\n')

    print(f'Run directory: {output.resolve()} (next event {log.next_index})', flush=True)
    try:
        while not finished():
            try:
                with connect(args) as scope:
                    resource = getattr(scope, 'resource_name', resource)
                    # After an outage, save a trigger that arrived meanwhile before re-arming.
                    rescue = connected_once and pending_trigger(scope, state['unread'] is not None)
                    prepare(scope, args, session)
                    if not connected_once:
                        save_setup(scope, stem.with_name(stem.name + '_setup.bin'))
                    connected_once, failures = True, 0
                    write_session()
                    if rescue:
                        counts['recovered'] += 1
                        record(scope, state['unread'] or None)
                    state['unread'] = None
                    print('Waiting for triggers (Ctrl-C to stop) ...', flush=True)
                    while not finished():
                        armed = wait_for_trigger(scope, args.poll_s, deadline, on_idle=status)
                        state['armed'] = 0.
                        if armed is None:
                            continue
                        state['live'] += armed
                        state['unread'] = armed
                        record(scope, armed)
                        state['unread'] = None
                        status()
            except ConfigurationError:
                raise
            except RECOVERABLE as error:
                if not connected_once:
                    raise
                # Live time accumulated before the failure still counts.
                state['live'] += state['armed']
                state['armed'] = 0.
                failures += 1
                counts['reconnects'] += 1
                log_error(error)
                if args.max_reconnects is not None and failures > args.max_reconnects:
                    raise
                delay = min(60, 2 ** failures)
                print(f'{datetime.now():%Y-%m-%d %H:%M:%S} communication error: {error!r}. '
                      f'Reconnecting in {delay} s ...', flush=True)
                usbtmc_clear(resource)
                write_session()
                time.sleep(max(delay, 5))  # also lets the device re-enumerate after the reset
    except KeyboardInterrupt:
        print('\nStopped by user.', flush=True)
    finally:
        log.close()
        write_session(stopped=datetime.now().astimezone().isoformat(timespec='seconds'))
    live = state['live'] + state['armed']
    print(f'{counts["triggers"]} triggers, {counts["accepted"]} saved coincidences, '
          f'{counts["rejected"]} rejected, {counts["recovered"]} recovered after reconnection, '
          f'live time {live/3600:.3f} h.', flush=True)
    return counts


def check(args):
    with open_scope(args) as scope:
        session = {}
        prepare(scope, args, session)
        width = max(map(len, session['settings']))
        for query, value in session['settings'].items():
            print(f'  {query:<{width}}  {value}')


def main():
    args = parse_args()
    try:
        if args.check:
            check(args)
        else:
            acquire(args)
    except ConfigurationError as error:
        sys.exit(f'ERROR: {error}')


if __name__ == '__main__':
    main()
