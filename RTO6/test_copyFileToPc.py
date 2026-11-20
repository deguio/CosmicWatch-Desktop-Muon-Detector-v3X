"""Offline tests for connection selection and SSH tunnel lifecycle."""
import argparse
import subprocess
import unittest
import tempfile
from pathlib import Path, PureWindowsPath
from types import SimpleNamespace
from unittest.mock import patch
from copyFileToPc import instrument_connection, copy_files, filter_time_range, parse_time_bound

class ConnectionTests(unittest.TestCase):
    def test_time_filter_boundaries_and_companions(self):
        names = ['RefCurve_2026-09-23_1_150000.csv',
                 'RefCurve_2026-09-24_1_141459.csv',
                 'RefCurve_2026-09-24_2_141500.csv',
                 'RefCurve_2026-09-24_2_141500.Wfm.csv',
                 'RefCurve_2026-09-24_3_150000.csv',
                 'RefCurve_2026-09-24_4_150001.csv',
                 'RefCurve_2026-09-25_1_141500.csv',
                 'RefCurve_2026-09-24_5_256000.csv', 'unknown.csv']
        paths = [PureWindowsPath('C:/data') / name for name in names]
        start = parse_time_bound('2026-09-24 14:15')
        stop = parse_time_bound('2026-09-24 15:00:00')
        self.assertEqual(filter_time_range(paths, start, stop), paths[2:5])
        self.assertEqual(filter_time_range(paths, start), paths[2:7])
        self.assertEqual(filter_time_range(paths, until=start), paths[:4])
        self.assertEqual(filter_time_range(paths), paths)
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_time_bound('14:15')

    def test_time_filter_applied_before_transfer(self):
        fake = self.fake_instrument()
        fake.query_str = lambda command: '6,999,"RefCurve_2026-09-24_1_141459.csv,ASC,3","RefCurve_2026-09-24_2_141500.csv,ASC,3"'
        with tempfile.TemporaryDirectory() as directory, patch('copyFileToPc.RsInstrument', return_value=fake):
            copy_files('test', 'socketio', PureWindowsPath(r'C:\data\RefCurve*'), True,
                       Path(directory), argparse.ArgumentParser(),
                       since=parse_time_bound('2026-09-24 14:15'))
            self.assertEqual([p.name for p in Path(directory).iterdir()], ['RefCurve_2026-09-24_2_141500.csv'])
            self.assertEqual(fake.transfers, 1)

    def fake_instrument(self, fail=False):
        class Instrument:
            idn_string = 'Test RTO'
            events = SimpleNamespace(on_read_handler=None)
            closed = False
            transfers = 0
            def query_str(self, command):
                return '3,999,"RefCurve_47.csv,ASC,3"'
            def read_file_from_instrument_to_pc(self, source, destination):
                self.transfers += 1
                Path(destination).write_bytes(b'abc')
                if fail:
                    raise TimeoutError('simulated timeout')
                self.events.on_read_handler(SimpleNamespace(end_of_transfer=True, total_size=3, transferred_size=3))
            def close(self):
                self.closed = True
        return Instrument()

    def test_missing_exact_name_never_starts_download(self):
        fake = self.fake_instrument()
        with tempfile.TemporaryDirectory() as directory, patch('copyFileToPc.RsInstrument', return_value=fake):
            with self.assertRaises(SystemExit):
                copy_files('test', 'socketio', PureWindowsPath(r'C:\data\RefCurve_47'), False,
                           Path(directory), argparse.ArgumentParser())
            self.assertEqual(fake.transfers, 0)
            self.assertFalse(list(Path(directory).iterdir()))
            self.assertTrue(fake.closed)

    def test_success_publishes_complete_file_and_timeout_leaves_none(self):
        for fail in (False, True):
            fake = self.fake_instrument(fail)
            with tempfile.TemporaryDirectory() as directory, patch('copyFileToPc.RsInstrument', return_value=fake):
                def download():
                    copy_files('test', 'socketio', PureWindowsPath(r'C:\data\RefCurve_47.csv'), False,
                               Path(directory), argparse.ArgumentParser(), 1000)
                if fail:
                    with self.assertRaises(TimeoutError):
                        download()
                    self.assertFalse(list(Path(directory).iterdir()))
                else:
                    download()
                    self.assertEqual((Path(directory) / 'RefCurve_47.csv').read_bytes(), b'abc')
                self.assertTrue(fake.closed)
                self.assertEqual(fake.visa_timeout, 1000)

    def test_existing_match_is_preserved_and_remaining_files_are_copied(self):
        fake = self.fake_instrument()
        fake.query_str = lambda command: '6,999,"RefCurve_50.csv,ASC,3","RefCurve_51.csv,ASC,3"'
        with tempfile.TemporaryDirectory() as directory, patch('copyFileToPc.RsInstrument', return_value=fake):
            folder = Path(directory)
            (folder / 'RefCurve_50.csv').write_bytes(b'existing data')
            copy_files('test', 'socketio', PureWindowsPath(r'C:\data\RefCurve_5*'), True,
                       folder, argparse.ArgumentParser())
            self.assertEqual((folder / 'RefCurve_50.csv').read_bytes(), b'existing data')
            self.assertEqual((folder / 'RefCurve_51.csv').read_bytes(), b'abc')
            self.assertEqual(fake.transfers, 1)
            # Re-running with every match present succeeds without any downloads.
            copy_files('test', 'socketio', PureWindowsPath(r'C:\data\RefCurve_5*'), True,
                       folder, argparse.ArgumentParser())
            self.assertEqual(fake.transfers, 1)

    def args(self, **overrides):
        values = dict(ssh_bastion=None, resource=None, instrument_host='192.0.2.10', instrument_port=5025, local_port=15025)
        values.update(overrides)
        return argparse.Namespace(**values)

    def test_direct_and_manual_tunnel(self):
        with patch('copyFileToPc.subprocess.run') as run:
            with instrument_connection(self.args()) as connection:
                self.assertEqual(connection, ('TCPIP::192.0.2.10::INSTR', 'pyvisa-py'))
            with instrument_connection(self.args(resource='TCPIP::127.0.0.1::15025::SOCKET')) as connection:
                self.assertEqual(connection[1], 'socketio')
            run.assert_not_called()

    def test_owned_tunnel_closed_on_transfer_error(self):
        with patch('copyFileToPc.subprocess.run') as run:
            with self.assertRaisesRegex(ValueError, 'transfer failed'):
                with instrument_connection(self.args(ssh_bastion='user@bastion')) as connection:
                    self.assertEqual(connection, ('TCPIP::127.0.0.1::15025::SOCKET', 'socketio'))
                    raise ValueError('transfer failed')
            start, stop = [call.args[0] for call in run.call_args_list]
            self.assertIn('127.0.0.1:15025:192.0.2.10:5025', start)
            self.assertIn('ExitOnForwardFailure=yes', start)
            self.assertEqual(start[start.index('-S')+1], stop[stop.index('-S')+1])
            self.assertEqual(stop[-4:], ['-O', 'exit', '--', 'user@bastion'])

    def test_failed_ssh_never_yields_connection(self):
        with patch('copyFileToPc.subprocess.run', side_effect=[subprocess.CalledProcessError(255, 'ssh'), None]):
            with self.assertRaises(subprocess.CalledProcessError):
                with instrument_connection(self.args(ssh_bastion='university')):
                    self.fail('Should not open an instrument session after SSH failure')

if __name__ == '__main__':
    unittest.main()
