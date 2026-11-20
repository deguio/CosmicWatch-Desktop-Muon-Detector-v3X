"""Offline regression tests; no oscilloscope connection."""
import tempfile
import unittest
import argparse
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.stats import norm
from analyze_pulses import (constant_fraction, fit_pulse, load_rto, prepare, pulse,
                            expand_inputs, gaussian_statistics, build_summary, amplitude_selection,
                            read_event_timestamp, rate_histogram)


class TimingTests(unittest.TestCase):
    def test_rate_common_bins_and_partial_exposure(self):
        start = datetime(2026, 9, 21, 23)
        stamps = [start+timedelta(hours=h) for h in (0, 1, 2, 2.5)]
        result = rate_histogram(stamps, [True, False, True, True], 2.)
        self.assertEqual(result['all_counts'], [2, 2])
        self.assertEqual(result['selected_counts'], [1, 2])
        self.assertEqual(result['exposure_hours'], [2., .5])
        self.assertEqual(result['all_rate_per_hour'], [1., 4.])
        self.assertEqual(result['selected_rate_per_hour'], [.5, 4.])
        self.assertEqual(result['selected_mean_per_hour'], 3/2.5)
        self.assertIsNone(rate_histogram([start, start], [True, False], 2.))

    def test_instrument_timestamp_not_filename(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'unrelated_name.csv'
            path.write_text('Timestamp:2026:09:21 18:45 9,694.305.997 s:\n')
            self.assertEqual(read_event_timestamp(path), datetime(2026, 9, 21, 18, 45, 9, 694305))
            path.write_text('No timestamp\n')
            with self.assertRaises(ValueError):
                read_event_timestamp(path)

    def test_upper_amplitude_cut_and_combined_interval(self):
        t = np.arange(50.)
        y = np.tile([100., -50.], (50, 1))
        y[30] += [20., -10.]
        self.assertTrue(amplitude_selection(t, y, ['C1', 'C2'], 20., maximum=20.)['accepted'])
        self.assertTrue(amplitude_selection(t, y, ['C1', 'C2'], 20., 10., 20.)['accepted'])
        for columns in (y, y[:, ::-1]):
            self.assertFalse(amplitude_selection(t, columns, ['C1', 'C2'], 20., maximum=15.)['accepted'])
        result = amplitude_selection(t, y, ['C1', 'C2'], 20., 12., 18.)
        self.assertFalse(result['accepted'])
        self.assertIn('Below 12 mV: C2', result['reason'])
        self.assertIn('Above 18 mV: C1', result['reason'])

    def test_amplitude_cut_requires_both_and_handles_baseline_polarity_equality(self):
        t = np.arange(50.)
        y = np.tile([100., -50.], (50, 1))
        y[30] += [20., -10.]
        result = amplitude_selection(t, y, ['C1', 'C2'], 20., 10.)
        self.assertTrue(result['accepted'])
        self.assertEqual(result['peak_amplitudes_mV'], {'C1': 20., 'C2': 10.})
        for columns in (y, y[:, ::-1]):
            self.assertFalse(amplitude_selection(t, columns, ['C1', 'C2'], 20., 15.)['accepted'])
        self.assertFalse(amplitude_selection(t, y, ['C1', 'C2'], 20., 25.)['accepted'])

    def test_batch_discovery_deduplicates_and_ignores_unrelated_csv(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            for name in ('event45', 'event46', 'event47'):
                (folder / f'{name}.csv').write_text('SignalRecordLength:40:\nResolution:5e-10:\n')
                (folder / f'{name}.Wfm.csv').write_text('1,2\n')
            (folder / 'unrelated.csv').write_text('column1,column2\n')
            paths = expand_inputs([directory, str(folder / '*.csv'), str(folder / 'event45.Wfm.csv')])
            self.assertEqual(len(paths), 3)
            self.assertTrue(all(not path.name.endswith('.Wfm.csv') for path in paths))
            with self.assertRaises(ValueError):
                expand_inputs([str(folder / 'missing*.csv')])

    def test_gaussian_summary_statistics(self):
        stats = gaussian_statistics([1., 2., 3.])
        self.assertEqual(stats['mean'], 2.)
        self.assertAlmostEqual(stats['sigma_mle'], np.sqrt(2/3))
        self.assertAlmostEqual(stats['sem'], 1/np.sqrt(3))
        self.assertLess(stats['mean_ci95'][0], 1.)
        self.assertIsNone(gaussian_statistics([1.])['sem'])
        self.assertNotEqual(gaussian_statistics([2., 2.])['fit_status'], 'ok')
        with self.assertRaises(ValueError):
            gaussian_statistics([np.nan])

    def test_summary_keeps_superluminal_estimates_and_accounts_for_exclusions(self):
        reports = []
        for index, (time, beta) in enumerate(((6., 1.1), (-3., 2.2), (0., None))):
            timing = {'delta_ns': time+16., 'corrected_delta_ns': time}
            if beta is not None:
                timing.update(beta=beta, speed_m_per_s=beta*299792458.)
            reports.append({'input': {'header': f'event{index}.csv'},
                            'timing': {'constant_fraction': timing}, 'channels': {}})
        args = argparse.Namespace(summary_method='constant_fraction', distance_m=2.,
                                  offset_ns=16., fraction=.2, show=False)
        with tempfile.TemporaryDirectory() as directory:
            summary = build_summary(reports, [{'input': 'bad.csv', 'error': 'test'}], args, Path(directory))
            self.assertEqual(summary['velocity_events'], 2)
            self.assertEqual(len(summary['excluded_from_velocity']), 1)
            self.assertEqual(len(summary['failed_events']), 1)
            self.assertAlmostEqual(summary['gaussian_velocity_beta']['mean'], 1.65)
            self.assertIsNone(summary['speed_from_mean_time'])
            for suffix in ('png', 'json', 'csv'):
                self.assertTrue((Path(directory) / f'summary_constant_fraction.{suffix}').is_file())
            self.assertFalse((Path(directory) / 'summary_constant_fraction.pdf').exists())

    def test_gaussian_convolution_against_quadrature(self):
        for time in (-5., 0., 7., 30., 200.):
            numerical = quad(lambda u: pulse(np.array([u]), 0., 1., 0., 7., 80.)[0]
                             * norm.pdf(time-u, scale=2.), 0., max(time+20., 30.),
                             epsabs=1e-11)[0]
            analytical = pulse(np.array([time]), 0., 1., 0., 7., 80., 2.)[0]
            self.assertAlmostEqual(numerical, analytical, places=9)

    def test_fraction_invariant_under_amplitude_baseline_polarity(self):
        t = np.arange(-100., 400., .25)
        wave = pulse(t, 0., 1., 3.17, 7., 80.)
        noise = np.random.default_rng(42).normal(0, .0001, len(t))
        times = []
        for scale, baseline in ((10., 2.), (40., -3.), (-25., 5.)):
            p = prepare(t, baseline + scale*(wave+noise), -50., 5)
            times.append(constant_fraction(t, p, .2)['time_ns'])
        np.testing.assert_allclose(times, times[0], atol=1e-10)
        shifted = pulse(t, 0., 1., 3.17+13.5, 7., 80.)
        measured = constant_fraction(t, prepare(t, shifted+noise, -50., 5), .2)['time_ns']
        self.assertAlmostEqual(measured-times[0], 13.5, delta=.02)

    def test_free_t0_recovery_for_both_models(self):
        t = np.arange(-100., 600., .5)
        rng = np.random.default_rng(11)
        for model, sigma in (('causal', 0.), ('gaussian', 1.8)):
            y = pulse(t, 2., 25., 3.17, 7., 80., sigma) + rng.normal(0, .1, len(t))
            fit = fit_pulse(t, y, prepare(t, y, -50., 5), model)
            self.assertAlmostEqual(fit['parameters']['t0_ns'], 3.17, delta=.15)
            self.assertAlmostEqual(fit['parameters']['tau_decay_ns'], 80., delta=1.)
            self.assertLess(fit['chi2_ndof'], 1.5)

    def test_export_time_axis_and_companion_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            header = Path(directory) / 'test.csv'
            data = Path(directory) / 'test.Wfm.csv'
            header.write_text('SignalRecordLength:40:\nResolution:5e-10:\nXStart:-1e-8:\n'
                              'XStop:1e-8:\nMultiChannelExport:On:\n'
                              'MultiChannelExportState:4:On:On:Off:Off:\n'
                              'MultiChannelSource:4:C1W1:C2W1:None:None:\n'
                              'MultiChannelViewUnit:4:V:V:::\n')
            np.savetxt(data, np.ones((40, 2))*.003, delimiter=',')
            t, y, labels, _ = load_rto(data)
            self.assertEqual(labels, ['C1W1', 'C2W1'])
            self.assertAlmostEqual(t[-1], 9.5)
            np.testing.assert_allclose(y, 3.)
            header.write_text(header.read_text()+'NumberOfAcquisitions:2:\n')
            with self.assertRaises(ValueError):
                load_rto(header)


if __name__ == '__main__':
    unittest.main()
