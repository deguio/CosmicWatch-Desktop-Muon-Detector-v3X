#!/usr/bin/env python3
"""Offline tests of acquire_coincidences.py against a scripted InfiniiVision SCPI model."""

import csv
import sys
import tempfile
import unittest
from unittest import mock
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pyvisa

import acquire_coincidences as acq

# analyze_pulses.py è condiviso con RTO6 e sta in analysis_code/.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis_code'))
from analyze_pulses import amplitude_selection, expand_inputs, load_waveform, read_event_timestamp


class FakeScope:
    """Single acquisitions of two positive pulses; C2 is 10 ns after C1."""

    LSB_MV = .4

    def __init__(self, events, polls=3, reject=(), fail=None):
        self.events = list(events)  # (C1 mV, C2 mV) per trigger
        self.polls = polls
        self.reject = set(reject)   # commands answered with an SCPI error
        self.errors, self.commands = [], []
        self.aer = self.ter = 0
        self.remaining = 0
        self.source = 1
        self.current = None
        self.timeout = 10000
        self.fail = dict(fail or {})  # query -> calls before one simulated VISA timeout
        self.calls = {}

    def write(self, command):
        self.commands.append(command)
        if command in self.reject:
            self.errors.append('-224,"Illegal parameter value"')
        # Like the DSOX1202A: more points than the record is an error, not a clamp.
        if command.startswith(':WAVeform:POINts ') and command.split()[1].isdigit() and int(command.split()[1]) > 1000:
            self.errors.append('-222,"Data out of range"')
        if command == ':SINGle':
            if not self.events:
                raise KeyboardInterrupt  # end of the scripted run
            self.current = self.events.pop(0)
            self.aer, self.remaining = 1, self.polls
        elif command.startswith(':WAVeform:SOURce CHANnel'):
            self.source = int(command[-1])

    def timeout_once(self, query):
        self.calls[query] = self.calls.get(query, 0) + 1
        if self.fail.get(query) == self.calls[query]:
            del self.fail[query]
            if query == ':OPERegister:CONDition?':
                self.remaining, self.ter = 0, 1   # the trigger arrives during the outage
            raise pyvisa.errors.VisaIOError(pyvisa.constants.StatusCode.error_timeout)

    def query(self, query):
        self.commands.append(query)
        self.timeout_once(query)
        if query in (':AER?', ':TER?'):
            name = query[1:4].lower()
            value = getattr(self, name)
            setattr(self, name, 0)
            return f'+{value}'
        if query == ':OPERegister:CONDition?':
            if self.remaining:
                self.remaining -= 1
                if not self.remaining:
                    self.ter = 1
            return f'+{8 * bool(self.remaining) + 32}'
        if query == ':SYSTem:ERRor?':
            return self.errors.pop(0) if self.errors else '+0,"No error"'
        if query == '*IDN?':
            return 'KEYSIGHT TECHNOLOGIES,DSOX1202A,FAKE0001,0.0'
        if query == ':WAVeform:PREamble?':
            return f'+0,+0,+1000,+1,+1.000000E-09,-3.000000E-07,+0,{self.LSB_MV/1e3:+E},+5.0E-03,+128'
        if query == ':ACQuire:SRATe?':
            return '+2.0E+09'
        if query == ':TRIGger:PATTern?':
            return '"11X",NONE,POS'
        if query == ':TRIGger:SWEep?':
            return 'NORM'
        if query.endswith('?') or '? ' in query:
            return '+1'
        raise AssertionError(f'Unexpected query {query}')

    def write_binary_values(self, message, values, datatype):
        self.commands.append(message + f'<{len(values)} bytes>')

    def query_binary_values(self, query, datatype, container, expect_termination=True):
        if query == ':SYSTem:SETup?':
            return container(b'SETUP')
        assert query == ':WAVeform:DATA?' and datatype == 'B'
        self.timeout_once(query)
        t = np.arange(1000) - 300.
        u = np.maximum(t - 10. * (self.source - 1), 0)
        shape = (1 - np.exp(-u / 10)) * np.exp(-u / 190)
        mv = self.current[self.source - 1] * shape / shape.max()
        rng = np.random.default_rng(self.source)
        codes = np.clip(np.round(128 + (mv + rng.normal(0, .4, len(t))) / self.LSB_MV), 0, 255)
        return container(codes.astype(int))


def connect_to(scope):
    @contextmanager
    def connect(args):
        yield scope
    return connect


def arguments(directory, *extra):
    return acq.parse_args(['--host', 'fake', '-o', str(directory), '--poll-s', '0.001', *extra])


class AcquisitionTest(unittest.TestCase):
    def test_saved_events_are_read_by_the_analysis(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            scope = FakeScope([(40, 30), (35, 2), (55, 25)])
            counts = acq.acquire(arguments(run, '--min-amplitude-mv', '5'), connect_to(scope))
            self.assertEqual(counts, {'triggers': 3, 'accepted': 2, 'rejected': 1, 'reconnects': 0, 'recovered': 0})
            self.assertEqual([p.name for p in expand_inputs([str(run)])], ['scope_000000.csv', 'scope_000002.csv'])
            with (run / 'events.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([r['file'] for r in rows], ['scope_000000.csv', '', 'scope_000002.csv'])
            t, y, labels, meta = load_waveform(run / 'scope_000000.csv')
            self.assertEqual(labels, ['C1', 'C2'])
            self.assertAlmostEqual(meta['dt_ns'], 1.)
            self.assertAlmostEqual(t[0], -250.)
            self.assertEqual(len(t), 950)  # record ends at +699 ns, inside the window
            selection = amplitude_selection(t, y, labels, t[0] + .15 * np.ptp(t), 20, 100)
            self.assertTrue(selection['accepted'])
            self.assertAlmostEqual(float(rows[0]['amplitude1_mV']), 40, delta=2.5)
            self.assertEqual(read_event_timestamp(run / 'scope_000000.csv').isoformat(timespec='milliseconds'),
                             rows[0]['pc_time'])
            self.assertTrue(list(run.glob('session_*_setup.bin')))
            self.assertIn(':TRIGger:SWEep NORMal', scope.commands)

    def test_restart_never_overwrites_events(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            acq.acquire(arguments(run), connect_to(FakeScope([(40, 30)])))
            acq.acquire(arguments(run, '--max-events', '1'), connect_to(FakeScope([(30, 40), (20, 20)])))
            self.assertEqual(sorted(p.name for p in run.glob('scope_*.csv')), ['scope_000000.csv', 'scope_000001.csv'])
            _, y, _, _ = load_waveform(run / 'scope_000001.csv')
            self.assertGreater(np.ptp(y[:, 1]), np.ptp(y[:, 0]))  # second run: C2 is the larger pulse

    def test_rejected_events_are_kept_apart(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            acq.acquire(arguments(run, '--min-amplitude-mv', '5', '10', '--save-rejected'),
                        connect_to(FakeScope([(40, 6)])))
            self.assertTrue((run / 'rejected' / 'scope_000000.csv').exists())
            self.assertFalse(list(run.glob('scope_*.csv')))
            self.assertIsNotNone(read_event_timestamp(run / 'rejected' / 'scope_000000.csv'))

    def test_pattern_trigger_commands(self):
        with tempfile.TemporaryDirectory() as directory:
            scope = FakeScope([])
            acq.configure(scope, arguments(directory, '--trigger', 'pattern', '--trigger-level-mv', '6', '8'))
            for command in (':TRIGger:MODE PATTern', ':TRIGger:PATTern "11X"', ':WAVeform:POINts:MODE RAW',
                            ':TRIGger:EDGE:LEVel 0.006,CHANnel1', ':TRIGger:EDGE:LEVel 0.008,CHANnel2'):
                self.assertIn(command, scope.commands)

    def test_rejected_setting_aborts_before_acquiring(self):
        with tempfile.TemporaryDirectory() as directory:
            scope = FakeScope([(40, 30)], reject={':TRIGger:PATTern "00X"'})
            args = arguments(directory, '--trigger', 'pattern', '--trigger-level-mv', '6', '--polarity', 'negative')
            with self.assertRaisesRegex(acq.ConfigurationError, 'PATTern "00X".*front panel'):
                acq.acquire(args, connect_to(scope))
            self.assertNotIn(':SINGle', scope.commands)

    def test_too_many_points_names_the_command(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(acq.ConfigurationError, 'POINts 1000000.*out of range'):
                acq.configure(FakeScope([]), arguments(directory, '--points', '1000000'))

    def test_baseline_stops_before_the_rising_edge(self):
        t = np.arange(-55., 355., .5)
        y = np.where(t >= -10, 50., 0.)[:, None] * np.ones(2)
        self.assertEqual(acq.peak_amplitudes(t, y, 'positive', -20.).tolist(), [50., 50.])
        # Fallback: first 15% (-55..+6.5 ns) contains part of the step.
        self.assertLess(acq.peak_amplitudes(t, y, 'positive', -60.)[0], 40.)

    def run_with_outage(self, fail):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(acq.time, 'sleep'):
            run = Path(directory)
            scope = FakeScope([(40, 30), (35, 25)], polls=5, fail=fail)
            counts = acq.acquire(arguments(run), connect_to(scope))
            with (run / 'events.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            errors = (run / 'errors.log').read_text()
            files = sorted(p.name for p in run.glob('scope_*.csv'))
        return counts, rows, errors, files

    def test_trigger_during_outage_is_recovered(self):
        # Timeout on the 3rd poll; the oscilloscope triggers while the link is down.
        counts, rows, errors, files = self.run_with_outage({':OPERegister:CONDition?': 3})
        self.assertEqual((counts['reconnects'], counts['recovered'], counts['accepted']), (1, 1, 2))
        self.assertEqual(rows[0]['live_s'], '')               # unknown for the recovered event
        self.assertNotEqual(rows[1]['live_s'], '')
        self.assertEqual(files, ['scope_000000.csv', 'scope_000001.csv'])
        self.assertIn('VisaIOError', errors)

    def test_failed_readout_is_saved_once(self):
        counts, rows, errors, files = self.run_with_outage({':WAVeform:DATA?': 1})
        self.assertEqual((counts['reconnects'], counts['recovered'], counts['accepted']), (1, 1, 2))
        self.assertNotEqual(rows[0]['live_s'], '')            # armed time was measured before the failure
        self.assertEqual(len(rows), 2)
        self.assertEqual(files, ['scope_000000.csv', 'scope_000001.csv'])

    def test_heartbeat_while_waiting(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            scope = FakeScope([(40, 30)], polls=50)
            with mock.patch('builtins.print') as printed:
                acq.acquire(arguments(run, '--status-s', '0.001'), connect_to(scope))
            beats = [c for c in printed.call_args_list if 'alive' in str(c)]
            self.assertGreater(len(beats), 5)
            self.assertIn('last event none', str(beats[0]))

    def test_window_outside_record_is_an_error(self):
        scope = FakeScope([(40, 30)])
        scope.current = scope.events[0]
        with self.assertRaises(acq.ReadoutError):
            acq.read_event(scope, (1000., 2000.))


if __name__ == '__main__':
    unittest.main()
