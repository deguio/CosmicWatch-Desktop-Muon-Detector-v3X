#!/usr/bin/env python3
"""Write synthetic two-channel muon events in the Keysight InfiniiVision CSV format.

Each file mimics one oscilloscope export: 1 GSa/s, trigger on C1 at t=0,
negative SiPM-like pulses, DC offsets, white noise and ADC quantization.
The true parameters of every event are saved in truth.json next to the CSVs.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import moyal

C = 299792458.


def pulse(t, t0, amplitude, tau_rise, tau_decay):
    """Negative pulse with unit-normalized peak, same shape as the analysis model."""
    u = np.maximum(t - t0, 0.)
    peak = tau_rise * np.log1p(tau_decay / tau_rise)
    norm = -np.expm1(-peak / tau_rise) * np.exp(-peak / tau_decay)
    return -amplitude * -np.expm1(-u / tau_rise) * np.exp(-u / tau_decay) / norm


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).parent / 'synthetic')
    parser.add_argument('-n', '--events', type=int, default=100)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--distance-m', type=float, default=3., help='Path length C1 -> C2')
    parser.add_argument('--beta', type=float, default=1., help='Muon speed v/c')
    parser.add_argument('--offset-ns', type=float, default=2., help='Instrumental (cable) delay added to C2')
    parser.add_argument('--jitter-ns', type=float, default=.5, help='Gaussian time jitter per channel')
    parser.add_argument('--noise-mv', type=float, default=1., help='RMS white noise per channel')
    parser.add_argument('--lsb-mv', type=float, default=.5, help='ADC quantization step; 0 disables')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    t = np.arange(-200, 800) * 1.  # ns, 1000 samples at 1 GSa/s
    tof = args.distance_m / (args.beta * C) * 1e9
    channels = [  # baseline, Landau MPV, Landau width (mV), tau_rise, tau_decay (ns)
        {'baseline_mV': 5.4, 'mpv_mV': 50., 'width_mV': 10., 'tau_rise_ns': 3., 'tau_decay_ns': 40.},
        {'baseline_mV': -4.6, 'mpv_mV': 40., 'width_mV': 8., 'tau_rise_ns': 3.5, 'tau_decay_ns': 45.},
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    truth = {'distance_m': args.distance_m, 'beta': args.beta, 'tof_ns': tof,
             'offset_ns': args.offset_ns, 'expected_delta_ns': tof + args.offset_ns,
             'jitter_ns_per_channel': args.jitter_ns, 'noise_mV': args.noise_mv,
             'lsb_mV': args.lsb_mv, 'seed': args.seed, 'channels': channels, 'events': []}
    for k in range(args.events):
        t0 = [rng.normal(0., args.jitter_ns), rng.normal(tof + args.offset_ns, args.jitter_ns)]
        event = {'file': f'scope_{k:03d}.csv', 't0_ns': t0, 'amplitude_mV': []}
        traces = []
        for ch, start in zip(channels, t0):
            amplitude = float(np.clip(moyal.rvs(ch['mpv_mV'], ch['width_mV'], random_state=rng), 10, 400))
            event['amplitude_mV'].append(amplitude)
            y = ch['baseline_mV'] + pulse(t, start, amplitude, ch['tau_rise_ns'], ch['tau_decay_ns'])
            y += rng.normal(0., args.noise_mv, len(t))
            if args.lsb_mv > 0:
                y = np.round(y / args.lsb_mv) * args.lsb_mv
            traces.append(y / 1000)
        event['delta_ns'] = t0[1] - t0[0]
        truth['events'].append(event)
        with (args.output_dir / event['file']).open('w') as stream:
            stream.write('x-axis,1,2\nsecond,Volt,Volt\n')
            stream.write(f'{(t[0] - 1) * 1e-9:+.6E},,\n')  # empty first row, as the scope writes it
            for time, a, b in zip(t * 1e-9, *traces):
                stream.write(f'{time:+.6E},{a:+.9E},{b:+.9E}\n')
    (args.output_dir / 'truth.json').write_text(json.dumps(truth, indent=2) + '\n')
    print(f'{args.events} events in {args.output_dir.resolve()}; '
          f'expected deltaT = {tof:.3f} (TOF) + {args.offset_ns:g} (offset) = {tof + args.offset_ns:.3f} ns')


if __name__ == '__main__':
    main()
