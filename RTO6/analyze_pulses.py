#!/usr/bin/env python3
"""Offline timing of a two-channel RTO6 CSV export; all internal units are ns/mV."""

import argparse
import csv
import glob
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
from scipy.special import log_ndtr
from scipy.stats import norm, t as student_t


PARAMETERS = ['baseline_mV', 'scale_mV', 't0_ns', 'tau_rise_ns', 'tau_decay_ns']


class AmplitudeRejected(Exception):
    def __init__(self, selection):
        self.selection = selection
        super().__init__(selection['reason'])


def amplitude_selection(t, y, labels, baseline_end, threshold=None, maximum=None):
    """Cut on raw peak magnitude above each channel's own pre-pulse baseline."""
    mask = t <= baseline_end
    if mask.sum() < 20:
        raise ValueError('Select at least 20 pre-pulse samples for the baseline.')
    amplitudes = np.max(np.abs(y - np.mean(y[mask], axis=0)), axis=0)
    below = [label for label, amplitude in zip(labels, amplitudes)
             if threshold is not None and amplitude < threshold]
    above = [label for label, amplitude in zip(labels, amplitudes)
             if maximum is not None and amplitude > maximum]
    reasons = []
    if below:
        reasons.append(f'Below {threshold:g} mV: {", ".join(below)}')
    if above:
        reasons.append(f'Above {maximum:g} mV: {", ".join(above)}')
    return {'accepted': not reasons, 'min_amplitude_mV': threshold,
            'max_amplitude_mV': maximum,
            'definition': 'Maximum absolute baseline-subtracted raw sample, before smoothing or fitting',
            'peak_amplitudes_mV': dict(zip(labels, map(float, amplitudes))),
            'reason': '; '.join(reasons) if reasons else 'Both channels pass'}


def load_rto(path):
    """Read the metadata CSV and companion .Wfm.csv, without inventing a time axis."""
    path = Path(path)
    header = path.with_name(path.name[:-8] + '.csv') if path.name.endswith('.Wfm.csv') else path
    waveform = header.with_name(header.stem + '.Wfm.csv')
    meta = {}
    for line in header.read_text(encoding='utf-8-sig').splitlines():
        fields = line.rstrip(':').split(':')
        meta[fields[0]] = fields[1:]
    y = np.loadtxt(waveform, delimiter=',', ndmin=2)
    n = int(meta['SignalRecordLength'][0])
    if len(y) != n or n < 30 or not np.isfinite(y).all():
        raise ValueError('Invalid sample count or nonfinite waveform data.')
    if int(meta.get('NumberOfAcquisitions', ['1'])[0]) != 1:
        raise ValueError('This script expects one simultaneous acquisition, not history records.')
    if meta.get('XAxisTDRDomain', ['time'])[0] != 'time':
        raise ValueError('A time-domain export is required.')
    if meta.get('MultiChannelExport', ['Off'])[0] == 'On':
        states = meta['MultiChannelExportState'][1:]
        sources = meta['MultiChannelSource'][1:]
        units = meta['MultiChannelViewUnit'][1:]
        labels = [source for source, state in zip(sources, states) if state == 'On']
        if any(unit != 'V' for unit, state in zip(units, states) if state == 'On'):
            raise ValueError('Only voltage waveforms in V are supported.')
    else:
        labels = [meta['Source'][0]]
        if meta.get('BaseUnit', [''])[0] != 'V':
            raise ValueError('Only voltage waveforms in V are supported.')
    if len(labels) != y.shape[1] or y.shape[1] != 2:
        raise ValueError('Export exactly two simultaneous voltage channels.')
    dt = float(meta['Resolution'][0]) * 1e9
    start = float(meta['XStart'][0]) * 1e9
    if dt <= 0:
        raise ValueError('Invalid time resolution.')
    t = start + dt * np.arange(n)  # XStop is the exclusive end, not the last sample.
    return t, y * 1000, labels, {'header': str(header), 'waveform': str(waveform),
                                'dt_ns': dt, 'samples': n}


def pulse(t, baseline, scale, t0, tau_rise, tau_decay, sigma=0.):
    """B + S * H(u)*(1-exp(-u/tr))*exp(-u/td), optionally convolved with N(0,sigma).

    The causal unblurred pulse is normalized to unit peak. With Gaussian blur,
    t0 is the onset of the underlying causal pulse, not a hard edge of the output.
    """
    u = np.asarray(t) - t0
    fast = tau_rise * tau_decay / (tau_rise + tau_decay)
    peak_time = tau_rise * np.log1p(tau_decay / tau_rise)
    normalization = -np.expm1(-peak_time / tau_rise) * np.exp(-peak_time / tau_decay)
    if sigma == 0:
        positive = np.maximum(u, 0.)
        shape = -np.expm1(-positive / tau_rise) * np.exp(-positive / tau_decay)
    else:
        def exp_gauss(tau):
            return np.exp(-u / tau + .5 * (sigma / tau)**2
                          + log_ndtr(u / sigma - sigma / tau))
        shape = exp_gauss(tau_decay) - exp_gauss(fast)
    return baseline + scale * shape / normalization


def prepare(t, y, baseline_end, window):
    base_mask = t <= baseline_end
    if base_mask.sum() < 20:
        raise ValueError('Select at least 20 pre-pulse samples for the baseline.')
    baseline = float(np.mean(y[base_mask]))
    noise = float(np.std(y[base_mask], ddof=1))
    if noise <= 0:
        raise ValueError('Zero baseline noise: cannot estimate fit uncertainty.')
    polarity = 1 if np.max(y - baseline) >= np.max(baseline - y) else -1
    signal = polarity * (y - baseline)
    smooth = savgol_filter(signal, window, 2) if window > 1 else signal.copy()
    peak = int(np.argmax(smooth))
    amplitude = float(smooth[peak])
    if t[peak] <= baseline_end or peak >= len(t) - window or amplitude < 5 * noise:
        raise ValueError('Pulse missing, truncated, low SNR, or inside the baseline interval.')
    if np.max(signal[base_mask]) > .3 * amplitude:
        raise ValueError('Baseline may contain the pulse; choose an earlier --baseline-end-ns.')
    return {'baseline': baseline, 'noise': noise, 'polarity': polarity,
            'smooth': smooth, 'peak': peak, 'amplitude': amplitude, 'baseline_mask': base_mask,
            'baseline_lag1_correlation': float(np.corrcoef(y[base_mask][:-1], y[base_mask][1:])[0, 1])}


def constant_fraction(t, prepared, fraction):
    """Offline fraction-of-peak crossing; not a delayed analog CFD emulation."""
    y = prepared['smooth']
    threshold = fraction * prepared['amplitude']
    peak = prepared['peak']
    crossings = np.flatnonzero((y[:peak] < threshold) & (y[1:peak+1] >= threshold))
    crossings = crossings[t[crossings] > t[prepared['baseline_mask']][-1]]
    if not len(crossings):
        raise ValueError(f'No rising-edge crossing at fraction {fraction}.')
    i = int(crossings[-1])
    slope = (y[i+1] - y[i]) / (t[i+1] - t[i])
    time = t[i] + (threshold-y[i])/slope
    return {'time_ns': float(time), 'fraction': float(fraction),
            'threshold_mV': float(threshold), 'slope_mV_per_ns': float(slope),
            'noise_over_slope_ns': float(prepared['noise']/slope),
            'threshold_snr': float(threshold/prepared['noise'])}


def fit_pulse(t, y, prepared, model):
    polarity = prepared['polarity']
    z = polarity * y
    baseline = polarity * prepared['baseline']
    amplitude = prepared['amplitude']
    noise = prepared['noise']
    dt = float(np.median(np.diff(t)))
    span = float(np.ptp(t))
    t20 = constant_fraction(t, prepared, .2)['time_ns']
    t80 = constant_fraction(t, prepared, .8)['time_ns']
    rise = max(t80 - t20, dt)
    peak_time = t[prepared['peak']]
    post = np.flatnonzero((t > peak_time) & (prepared['smooth'] < amplitude / np.e))
    decay = max(t[post[0]] - peak_time, rise * 2) if len(post) else span / 4
    low = [baseline - 10*noise, amplitude*.1, max(t[0], t20-10*rise), dt*.02, dt*.1]
    high = [baseline + 10*noise, amplitude*10, t80, span/2, span*10]
    if model == 'gaussian':
        low.append(dt*.02)
        high.append(max(rise*5, dt))
    def residual(p):
        return (pulse(t, *p) - z) / noise
    candidates = []
    for factor in (.2, 1., 3.):
        initial = [baseline, amplitude, t20-rise*.3, rise*factor, decay]
        if model == 'gaussian':
            initial.append(max(dt*.1, rise*.2/factor))
        initial = np.clip(initial, np.array(low)+1e-7, np.array(high)-1e-7)
        result = least_squares(residual, initial, bounds=(low, high), x_scale='jac', max_nfev=3000)
        if result.success:
            candidates.append(result)
    if not candidates:
        raise ValueError(f'{model} fit failed to converge.')
    result = min(candidates, key=lambda r: np.sum(r.fun**2))
    names = PARAMETERS + (['sigma_ns'] if model == 'gaussian' else [])
    chi2 = float(np.sum(result.fun**2))
    ndof = len(t)-len(result.x)
    _, singular, vt = np.linalg.svd(result.jac, full_matrices=False)
    rank_ok = singular[-1] > singular[0] * 1e-12
    covariance = (vt.T / singular**2) @ vt if rank_ok else None
    errors = np.sqrt(np.diag(covariance)) if covariance is not None else np.full(len(names), np.nan)
    correlations = covariance / np.outer(errors, errors) if covariance is not None else None
    flags = []
    if chi2/ndof > 2:
        flags.append('Model residuals exceed baseline noise: formal errors are not total timing errors.')
    near_bound = (np.minimum(result.x-np.array(low), np.array(high)-result.x)
                  < 1e-5*(np.array(high)-np.array(low)))
    if np.any(result.active_mask) or np.any(near_bound):
        flags.append('One or more fit parameters are at/near a bound; covariance may be unreliable.')
    if not rank_ok:
        flags.append('Singular covariance: parameters are not identifiable.')
    elif np.max(np.abs(np.delete(correlations[2], 2))) > .95:
        flags.append('t0 strongly correlated with pulse-shape parameters.')
    residual_mV = z - pulse(t, *result.x)
    lag1 = float(np.corrcoef(residual_mV[:-1], residual_mV[1:])[0, 1])
    if abs(lag1) > .2:
        flags.append('Correlated residuals: independent-noise covariance can underestimate uncertainty.')
    return {'model': model, 'parameters': dict(zip(names, map(float, result.x))),
            'formal_iid_errors': dict(zip(names, map(float, errors))),
            'covariance': covariance.tolist() if covariance is not None else None,
            'correlation': correlations.tolist() if correlations is not None else None,
            'parameter_order': names, 'chi2': chi2, 'ndof': ndof, 'chi2_ndof': chi2/ndof,
            'residual_rms_mV': float(np.sqrt(np.mean(residual_mV**2))),
            'residual_lag1_correlation': lag1, 'flags': flags, 'polarity': polarity}


def evaluate_fit(t, fit):
    return fit['polarity'] * pulse(t, *[fit['parameters'][name] for name in fit['parameter_order']])


def clean_json(value):
    if isinstance(value, dict):
        return {k: clean_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean_json(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', nargs='+', help='One or more CSV files, directories, or quoted glob patterns')
    parser.add_argument('--fraction', type=float, default=.2)
    parser.add_argument('--min-amplitude-mv', type=float,
                        help='Require BOTH raw peak magnitudes above their baselines to be >= this value in mV; default: no amplitude cut')
    parser.add_argument('--max-amplitude-mv', type=float,
                        help='Require BOTH raw peak magnitudes above their baselines to be <= this value in mV; default: no upper cut')
    parser.add_argument('--smooth-samples', type=int, default=5, help='Odd Savitzky-Golay window; 1 disables smoothing')
    parser.add_argument('--baseline-end-ns', type=float, help='End of pre-pulse baseline; default first 15%% of record')
    parser.add_argument('--fit-range-ns', type=float, nargs=2, metavar=('START', 'STOP'))
    parser.add_argument('--output-dir', type=Path, help='Default: timing_results next to input')
    parser.add_argument('--show', action='store_true', help='Also open matplotlib plots')
    parser.add_argument('--distance-m', type=float, help='Muon path length between detectors, not just vertical separation')
    parser.add_argument('--offset-ns', type=float, help='Calibrated instrumental t(C2)-t(C1) offset to subtract')
    parser.add_argument('--summary-method', choices=['constant_fraction', 'causal', 'gaussian'],
                        default='constant_fraction', help='Timing estimator used for the batch velocity histogram')
    args = parser.parse_args()
    if args.min_amplitude_mv is not None and (not np.isfinite(args.min_amplitude_mv) or args.min_amplitude_mv < 0):
        parser.error('--min-amplitude-mv must be finite and nonnegative.')
    if args.max_amplitude_mv is not None and (not np.isfinite(args.max_amplitude_mv) or args.max_amplitude_mv < 0):
        parser.error('--max-amplitude-mv must be finite and nonnegative.')
    if (args.min_amplitude_mv is not None and args.max_amplitude_mv is not None
            and args.min_amplitude_mv > args.max_amplitude_mv):
        parser.error('--min-amplitude-mv must not exceed --max-amplitude-mv.')
    if not 0 < args.fraction < 1:
        parser.error('--fraction must be between 0 and 1.')
    if args.smooth_samples != 1 and (args.smooth_samples < 3 or args.smooth_samples % 2 == 0):
        parser.error('--smooth-samples must be 1 or an odd integer >= 3.')
    if args.distance_m is not None and (args.distance_m <= 0 or args.offset_ns is None):
        parser.error('A positive --distance-m requires an explicit calibrated --offset-ns (0 if calibrated zero).')
    return args, parser


def analyze_event(args):
    t, y, labels, metadata = load_rto(args.input)
    if args.smooth_samples >= len(t):
        raise ValueError('Smoothing window exceeds record length.')
    baseline_end = args.baseline_end_ns if args.baseline_end_ns is not None else t[0]+.15*np.ptp(t)
    selection = None
    if args.min_amplitude_mv is not None or args.max_amplitude_mv is not None:
        selection = amplitude_selection(t, y, labels, baseline_end,
                                        args.min_amplitude_mv, args.max_amplitude_mv)
        if not selection['accepted']:
            raise AmplitudeRejected(selection)
    mask = np.ones(len(t), dtype=bool)
    if args.fit_range_ns:
        mask = (t >= args.fit_range_ns[0]) & (t <= args.fit_range_ns[1])
        if mask.sum() < 30:
            raise ValueError('Fit interval needs at least 30 samples.')
    prepared = [prepare(t, y[:, j], baseline_end, args.smooth_samples) for j in range(2)]
    fit_prepared = prepared if mask.all() else [
        prepare(t[mask], y[mask, j], baseline_end, args.smooth_samples) for j in range(2)]
    fits = [{model: fit_pulse(t[mask], y[mask, j], p, model)
             for model in ('causal', 'gaussian')} for j, p in enumerate(fit_prepared)]
    cfd = [constant_fraction(t, p, args.fraction) for p in prepared]
    scan = []
    for f in np.linspace(.1, .8, 15):
        pair = [constant_fraction(t, p, float(f)) for p in prepared]
        scan.append({'fraction': float(f), 'delta_ns': pair[1]['time_ns']-pair[0]['time_ns'],
                     'min_threshold_snr': min(r['threshold_snr'] for r in pair)})
    deltas = {'constant_fraction': {'delta_ns': cfd[1]['time_ns']-cfd[0]['time_ns']}}
    for model in ('causal', 'gaussian'):
        deltas[model] = {
            'delta_ns': fits[1][model]['parameters']['t0_ns']-fits[0][model]['parameters']['t0_ns'],
            'formal_iid_error_ns': float(np.hypot(*[f[model]['formal_iid_errors']['t0_ns'] for f in fits]))}
    if args.offset_ns is not None:
        for result in deltas.values():
            result['corrected_delta_ns'] = result['delta_ns'] - args.offset_ns
            if args.distance_m is not None and result['corrected_delta_ns'] != 0:
                result['speed_m_per_s'] = args.distance_m / abs(result['corrected_delta_ns']) * 1e9
                result['beta'] = result['speed_m_per_s'] / 299792458.
    smoothing_scan = []
    for window in (1, 3, 5, 7, 9):
        pair = [constant_fraction(t, prepare(t, y[:, j], baseline_end, window), args.fraction) for j in range(2)]
        smoothing_scan.append({'samples': window, 'delta_ns': pair[1]['time_ns']-pair[0]['time_ns']})
    report = {'input': metadata, 'amplitude_selection': selection,
              'convention': f'deltaT = t({labels[1]}) - t({labels[0]})',
              'baseline_end_ns': baseline_end, 'smooth_samples': args.smooth_samples,
              'fraction': args.fraction, 'fit_range_ns': args.fit_range_ns,
              'offset_ns': args.offset_ns, 'distance_m': args.distance_m,
              'channels': {label: {'baseline_mV': p['baseline'], 'noise_rms_mV': p['noise'],
                                  'baseline_lag1_correlation': p['baseline_lag1_correlation'],
                                  'polarity': p['polarity'], 'amplitude_mV': p['amplitude'],
                                  'constant_fraction': cf, 'fits': fit}
                           for label, p, cf, fit in zip(labels, prepared, cfd, fits)},
              'timing': deltas, 'fraction_scan': scan, 'smoothing_scan': smoothing_scan,
              'uncertainty_note': 'Fit covariance assumes correct model and independent baseline-level noise. '
                                  'No detector, cable, timebase, calibration, or path-length uncertainties included. '
                                  'CFD noise/slope is a diagnostic, not a complete error estimate.'}
    output = args.output_dir or args.input.parent / 'timing_results'
    output.mkdir(parents=True, exist_ok=True)
    stem = args.input.name.removesuffix('.Wfm.csv').removesuffix('.csv')
    (output / f'{stem}_timing.json').write_text(json.dumps(clean_json(report), indent=2)+'\n')
    with (output / f'{stem}_fraction_scan.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(scan[0]))
        writer.writeheader()
        writer.writerows(scan)
    import matplotlib
    if getattr(args, 'headless', not args.show):
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10,
                         'axes.labelsize': 9, 'legend.fontsize': 8})
    fig = plt.figure(figsize=(12, 7.5), dpi=110, layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=[1.25, 1])
    colors = {'causal': 'tab:orange', 'gaussian': 'tab:green'}
    channel_colors = ['#1674b1', '#c8443c']
    for j, (label, p) in enumerate(zip(labels, prepared)):
        pair_grid = grid[0, j].subgridspec(2, 1, height_ratios=[3.5, 1], hspace=0.02)
        ax = fig.add_subplot(pair_grid[0])
        residual_ax = fig.add_subplot(pair_grid[1], sharex=ax)
        ax.plot(t, y[:, j], '.', ms=2, color='0.5', alpha=.65, label='Samples')
        for model, fit in fits[j].items():
            ax.plot(t[mask], evaluate_fit(t[mask], fit), color=colors[model],
                    lw=1.2, ls='-' if model == 'causal' else '--', label=model)
            residual_ax.plot(t[mask], y[mask, j]-evaluate_fit(t[mask], fit),
                             color=colors[model], lw=.55, alpha=.8)
        ax.set(title=label, ylabel='Voltage [mV]')
        ax.tick_params(labelbottom=False)
        ax.legend(loc='upper left', fontsize=7, framealpha=.9)
        fit = fits[j]['causal']
        values, errors = fit['parameters'], fit['formal_iid_errors']
        parameter_lines = ['Causal fit (formal errors)']
        for name, display, unit in [
            ('baseline_mV', 'B', 'mV'), ('scale_mV', 'S', 'mV'),
            ('t0_ns', r'$t_0$', 'ns'), ('tau_rise_ns', r'$\tau_r$', 'ns'),
            ('tau_decay_ns', r'$\tau_d$', 'ns'),
        ]:
            value = values[name] * (fit['polarity'] if name == 'baseline_mV' else 1)
            parameter_lines.append(f'{display} = {value:.2f} ± {errors[name]:.2f} {unit}')
        parameter_lines.append(r'$\chi^2$/ndf' + f' = {fit["chi2_ndof"]:.2f}')
        ax.text(.98, .97, '\n'.join(parameter_lines), transform=ax.transAxes,
                va='top', ha='right', fontsize=8,
                bbox=dict(boxstyle='round,pad=.4', facecolor='white', edgecolor='0.85', alpha=.92))
        residual_ax.axhline(0, color='black', lw=.6)
        residual_ax.set(xlabel='Time [ns]', ylabel='Res. [mV]')
        residual_ax.locator_params(axis='y', nbins=3)
        ax.grid(alpha=.12)
        residual_ax.grid(alpha=.15)

    edge = fig.add_subplot(grid[1, 0])
    for j, (label, p) in enumerate(zip(labels, prepared)):
        edge.plot(t, p['polarity']*(y[:, j]-p['baseline'])/p['amplitude'],
                  '.', ms=2, alpha=.3, color=channel_colors[j])
        edge.plot(t, p['smooth']/p['amplitude'], color=channel_colors[j], lw=1.3, label=label)
        edge.axvline(cfd[j]['time_ns'], color=channel_colors[j], ls=':', lw=1)
    edge.axhline(args.fraction, color='0.3', ls='--', lw=.7)
    timing = deltas['constant_fraction']
    timing_lines = [f'Observed Δt = {timing["delta_ns"]:.2f} ns']
    if args.offset_ns is not None:
        timing_lines.append(f'Shift C2 = {args.offset_ns:g} ns  →  Δt corrected = {timing["corrected_delta_ns"]:.2f} ns')
    else:
        timing_lines.append('Shift not supplied: no correction')
    if args.distance_m is not None:
        if 'beta' in timing:
            timing_lines.append(f'L = {args.distance_m:g} m  |  vμ = {timing["speed_m_per_s"]/1e8:.3f} × 10⁸ m/s = {timing["beta"]:.3f} c')
        else:
            timing_lines.append(f'L = {args.distance_m:g} m  |  speed undefined (zero corrected Δt)')
    else:
        timing_lines.append('L not supplied: use --distance-m to estimate vμ')
    edge.text(.025, .98, '\n'.join(timing_lines), transform=edge.transAxes,
              va='top', fontsize=8, linespacing=1.4)
    edge.set(xlim=(min(c['time_ns'] for c in cfd)-15, max(t[p['peak']] for p in prepared)+10),
             ylim=(-.1, 1.65), xlabel='Time [ns]', ylabel='Normalized signal',
             title=f'Constant fraction {args.fraction:.0%} · Δt = t(C2) − t(C1)')
    edge.legend(loc='lower right', framealpha=.9)
    edge.grid(alpha=.12)
    ax = fig.add_subplot(grid[1, 1])
    ax.plot([s['fraction'] for s in scan], [s['delta_ns'] for s in scan],
            'o-', ms=3, lw=1, label='Fraction crossing')
    for model in ('causal', 'gaussian'):
        ax.axhline(deltas[model]['delta_ns'], color=colors[model], ls='--',
                   label=f'{model}: {deltas[model]["delta_ns"]:.2f} ns')
    ax.set(xlabel='Fraction of pulse amplitude', ylabel='Observed Δt [ns]',
           title='Sensitivity to fraction / fit model')
    ax.set_title('Sensitivity to fraction / fit model', pad=30)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.01), ncol=3,
              borderaxespad=0, frameon=False, fontsize=7)
    ax.grid(alpha=.15)
    fig.suptitle(f'{stem} — timing analysis', fontsize=11)
    fig.savefig(output / f'{stem}_timing.png', dpi=110)
    fig.savefig(output / f'{stem}_timing.pdf')

    # Separate oscilloscope-like view: raw voltages and common acquisition time axis.
    # The cable correction is used in timing calculations, not to move these traces.
    scope, scope_ax = plt.subplots(figsize=(11, 5.5), dpi=110, layout='constrained')
    scope.patch.set_facecolor('#101820')
    scope_ax.set_facecolor('#101820')
    scope_colors = ['#ffd166', '#53d8fb']
    for j, label in enumerate(labels):
        color = scope_colors[j]
        scope_ax.plot(t, y[:, j], '.', ms=2.5, alpha=.55, color=color, label=f'{label} samples')
        scope_ax.plot(t[mask], evaluate_fit(t[mask], fits[j]['causal']),
                      color=color, lw=1.7, label=f'{label} causal fit')
        scope_ax.plot(t[mask], evaluate_fit(t[mask], fits[j]['gaussian']),
                      color=color, lw=1, ls='--', alpha=.8, label=f'{label} Gaussian fit')
    scope_ax.set(xlabel='Acquisition time [ns]', ylabel='Voltage [mV]',
                 title='Both channels · full pulse fits · original amplitudes and timing')
    scope_ax.tick_params(colors='#e5edf4')
    scope_ax.xaxis.label.set_color('#e5edf4')
    scope_ax.yaxis.label.set_color('#e5edf4')
    scope_ax.title.set_color('#e5edf4')
    for spine in scope_ax.spines.values():
        spine.set_color('#667785')
    scope_ax.minorticks_on()
    scope_ax.grid(which='major', color='#506270', alpha=.55)
    scope_ax.grid(which='minor', color='#506270', alpha=.2, ls=':')
    scope_ax.legend(loc='upper right', facecolor='#182631', edgecolor='#506270',
                    labelcolor='#e5edf4', fontsize=8)
    scope.savefig(output / f'{stem}_scope.png', dpi=110, facecolor=scope.get_facecolor())
    scope.savefig(output / f'{stem}_scope.pdf', facecolor=scope.get_facecolor())
    for method, result in deltas.items():
        print(f'{method}: deltaT = {result["delta_ns"]:.4f} ns')
        if 'corrected_delta_ns' in result:
            print(f'  After subtracting {args.offset_ns:g} ns: {result["corrected_delta_ns"]:.4f} ns')
        if 'beta' in result:
            print(f'  v/c = {result["beta"]:.4f} (conditional on supplied path and offset; no total uncertainty)')
    for label, fit in zip(labels, fits):
        for model, result in fit.items():
            print(f'{label} {model}: t0={result["parameters"]["t0_ns"]:.4f} ns; chi2/ndf={result["chi2_ndof"]:.2f}')
            for flag in result['flags']:
                print(f'  {flag}')
    print(f'Results: {output.resolve()}')
    if args.show:
        plt.show()
    plt.close(fig)
    plt.close(scope)
    return report


def expand_inputs(inputs):
    """Resolve headers and companions exactly once; directories are nonrecursive."""
    headers = set()
    for expression in inputs:
        candidate = Path(expression).expanduser()
        matches = sorted(candidate.glob('*.csv')) if candidate.is_dir() else [
            Path(value) for value in sorted(glob.glob(str(candidate)))]
        if not matches:
            raise ValueError(f'No input matches: {expression}')
        added = False
        for path in matches:
            if path.name.endswith('.Wfm.csv'):
                path = path.with_name(path.name[:-8] + '.csv')
            if not path.is_file():
                raise ValueError(f'Missing metadata file: {path}')
            # Skip unrelated CSV files when scanning a directory or wildcard.
            content = path.read_text(encoding='utf-8-sig')
            if 'SignalRecordLength:' not in content or 'Resolution:' not in content:
                if len(matches) == 1 and not candidate.is_dir() and not glob.has_magic(str(candidate)):
                    raise ValueError(f'Not an RTO metadata file: {path}')
                continue
            headers.add(path.resolve())
            added = True
        if not added:
            raise ValueError(f'No RTO exports found in: {expression}')
    paths = sorted(headers)
    if len({path.stem for path in paths}) != len(paths):
        raise ValueError('Different input directories contain identical export names; use distinct names to avoid overwriting outputs.')
    return paths


def gaussian_statistics(values):
    """Unbinned normal MLE; SEM and Student interval assume independent normal data."""
    values = np.asarray(values, dtype=float)
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError('Gaussian summary requires finite, nonempty data.')
    mu, sigma = norm.fit(values)
    result = {'n': len(values), 'mean': float(mu), 'sigma_mle': float(sigma),
              'sem': None, 'mean_ci95': None,
              'fit_status': 'ok' if len(values) >= 2 and sigma > 0 else 'insufficient_spread_or_count'}
    if len(values) >= 2 and sigma > 0:
        sem = float(np.std(values, ddof=1) / np.sqrt(len(values)))
        width = float(student_t.ppf(.975, len(values)-1) * sem)
        result.update(sem=sem, mean_ci95=[float(mu-width), float(mu+width)])
    return result


def build_summary(reports, failures, args, output, rejected=None):
    rejected = rejected or []
    method = args.summary_method
    rows, exclusions = [], []
    for report in reports:
        event = Path(report['input']['header']).stem
        timing = report['timing'][method]
        flags = [f'{label}: {flag}' for label, channel in report['channels'].items()
                 for flag in channel['fits']['causal' if method == 'constant_fraction' else method]['flags']]
        row = {'event': event, 'method': method, 'delta_ns': timing['delta_ns'],
               'corrected_delta_ns': timing['corrected_delta_ns'],
               'speed_m_per_s': timing.get('speed_m_per_s'), 'beta': timing.get('beta'),
               'distance_m': args.distance_m, 'offset_ns': args.offset_ns,
               'fraction': args.fraction, 'flags': ' | '.join(flags)}
        rows.append(row)
        if row['beta'] is None or not np.isfinite(row['beta']):
            exclusions.append({'event': event, 'reason': 'Zero corrected time or nonfinite velocity'})
    usable = [row for row in rows if row['beta'] is not None and np.isfinite(row['beta'])]
    statistics = gaussian_statistics([row['beta'] for row in usable]) if usable else None
    tof_statistics = gaussian_statistics([row['corrected_delta_ns'] for row in usable]) if usable else None
    messages = ['Gaussian MLE uses individual events, not histogram bins. SEM/CI assume independent normal events.',
                'Common distance, cable and detector calibration uncertainties are not included.',
                'No cut on v/c > 1 or pulse-fit quality; inspect individual fits and flagged events.']
    if len(usable) < 10:
        messages.append('Very small sample: Gaussian shape and precision cannot be established; preliminary result.')
    tof_velocity = None
    if usable:
        times = np.array([row['corrected_delta_ns'] for row in usable])
        if np.all(times > 0) or np.all(times < 0):
            mean_t = tof_statistics['mean']
            speed = args.distance_m / abs(mean_t) * 1e9
            interval = tof_statistics['mean_ci95']
            beta_ci = None
            if interval is not None and interval[0] * interval[1] > 0:
                beta_ci = sorted(args.distance_m / abs(bound) * 1e9 / 299792458. for bound in interval)
            tof_velocity = {'speed_m_per_s': speed, 'beta': speed/299792458., 'beta_ci95': beta_ci,
                            'assumption': 'Common flight path and direction; Gaussian timing errors.'}
        else:
            messages.append('Mixed signs of corrected times: no common-direction speed from mean time reported.')
    summary = {'method': method, 'distance_m': args.distance_m, 'offset_ns': args.offset_ns,
               'min_amplitude_mV': getattr(args, 'min_amplitude_mv', None),
               'max_amplitude_mV': getattr(args, 'max_amplitude_mv', None),
               'amplitude_rejected_events': rejected,
               'fraction': args.fraction, 'successful_events': len(reports), 'velocity_events': len(usable),
               'failed_events': failures, 'excluded_from_velocity': exclusions,
               'gaussian_velocity_beta': statistics, 'gaussian_corrected_time_ns': tof_statistics,
               'speed_from_mean_time': tof_velocity, 'events': rows, 'notes': messages}
    if statistics is not None:
        summary['gaussian_velocity_m_per_s'] = {
            'mean': statistics['mean']*299792458., 'sigma_mle': statistics['sigma_mle']*299792458.,
            'sem': statistics['sem']*299792458. if statistics['sem'] is not None else None}
    prefix = output / f'summary_{method}'
    prefix.with_suffix('.json').write_text(json.dumps(clean_json(summary), indent=2)+'\n')
    with prefix.with_suffix('.csv').open('w', newline='') as stream:
        columns = ['event', 'method', 'delta_ns', 'corrected_delta_ns', 'speed_m_per_s', 'beta',
                   'distance_m', 'offset_ns', 'fraction', 'flags']
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    # Separate selection log includes both accepted and rejected exports.
    selection_rows = [{'event': Path(report['input']['header']).stem,
                       **(report.get('amplitude_selection') or {'accepted': True, 'reason': 'No amplitude cut'})}
                      for report in reports]
    selection_rows += rejected
    with (output / f'summary_{method}_selection.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['event', 'accepted', 'min_amplitude_mV', 'max_amplitude_mV',
                                                   'peak_amplitudes_mV', 'reason'], extrasaction='ignore')
        writer.writeheader()
        for row in selection_rows:
            writer.writerow({**row, 'peak_amplitudes_mV': json.dumps(row.get('peak_amplitudes_mV', {}))})
    import matplotlib
    if not args.show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.3), dpi=110, layout='constrained')
    for ax, values, stats, label, reference in [
        (axes[0], [row['beta'] for row in usable], statistics, 'Muon speed v / c', 1.),
        (axes[1], [row['corrected_delta_ns'] for row in usable], tof_statistics,
         'Corrected Δt [ns]', args.distance_m / 299792458. * 1e9),
    ]:
        ax.set(xlabel=label, ylabel='Probability density')
        ax.grid(alpha=.15)
        if not values:
            ax.text(.5, .5, 'No usable velocities', ha='center', transform=ax.transAxes)
            continue
        values = np.asarray(values)
        mean, sigma = stats['mean'], stats['sigma_mle']
        width = max(float(np.ptp(values)), abs(mean)*.02, .01)
        # Recompute after all selections; no fixed minimum/maximum bin count.
        # Equal-width bins are only a display choice; the fit is unbinned.
        bin_count = int(np.ceil(np.sqrt(len(values))))
        ax.hist(values, bins=np.linspace(values.min()-width*.12, values.max()+width*.12,
                                        bin_count+1),
                density=True, color='#6aaed6', edgecolor='white', alpha=.65,
                label=f'Events (N={len(values)}, {bin_count} bins)')
        ax.plot(values, np.zeros_like(values), '|', ms=12, color='#184b70', clip_on=False)
        if stats['fit_status'] == 'ok':
            x = np.linspace(min(values.min()-.2*width, mean-3*sigma),
                            max(values.max()+.2*width, mean+3*sigma), 400)
            ax.plot(x, norm.pdf(x, mean, sigma), color='#c8443c', lw=1.5, label='Gaussian MLE (unbinned)')
        ax.axvline(reference, color='0.35', ls=':', label='c reference')
        lines = [f'μ = {mean:.4f}', f'σ (MLE) = {sigma:.4f}']
        if stats['sem'] is not None:
            lines[0] += f' ± {stats["sem"]:.4f} (SEM)'
            lines.append(f'95% CI μ: [{stats["mean_ci95"][0]:.4f}, {stats["mean_ci95"][1]:.4f}]')
        else:
            lines.append('Gaussian width / mean uncertainty not estimable')
        if ax is axes[0]:
            lines.append(f'Mean v = {mean*299792458./1e8:.3f} × 10⁸ m/s')
        elif tof_velocity is not None:
            lines.append(f'L / |mean Δt| = {tof_velocity["beta"]:.4f} c')
        ax.set_title('\n'.join(lines), fontsize=9, pad=10)
        ax.legend(loc='upper center', bbox_to_anchor=(.5, -.16), fontsize=7, ncol=2, frameon=False)
    label = f'CFD {args.fraction:.0%}' if method == 'constant_fraction' else f'{method} t0'
    threshold = getattr(args, 'min_amplitude_mv', None)
    maximum = getattr(args, 'max_amplitude_mv', None)
    if threshold is not None and maximum is not None:
        cut_label = f' · {threshold:g} ≤ both peaks ≤ {maximum:g} mV'
    elif maximum is not None:
        cut_label = f' · both peaks ≤ {maximum:g} mV'
    else:
        cut_label = f' · both peaks ≥ {threshold:g} mV' if threshold is not None else ''
    fig.suptitle(f'Muon velocity summary · {label} · L = {args.distance_m:g} m · shift = {args.offset_ns:g} ns\n'
                 f'{len(usable)} velocities · {len(rejected)} amplitude rejects · {len(failures)} failed exports · {len(exclusions)} undefined velocities'
                 + cut_label
                 + (' · PRELIMINARY: very small sample' if len(usable) < 10 else ''), fontsize=10)
    fig.supxlabel('Statistical errors only; distance/offset systematics excluded. Gaussian speed and Gaussian time are different models.', fontsize=8)
    fig.savefig(prefix.with_suffix('.png'), dpi=110)
    fig.savefig(prefix.with_suffix('.pdf'))
    if args.show:
        plt.show()
    plt.close(fig)
    if statistics is not None:
        print(f'Summary ({method}): N={statistics["n"]}, mean v/c={statistics["mean"]:.5f}, SEM={statistics["sem"]}')
    print(f'Summary files: {prefix}')
    return summary


def main():
    args, parser = parse_args()
    try:
        paths = expand_inputs(args.inputs)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    batch = len(paths) > 1
    if batch and (args.distance_m is None or args.offset_ns is None):
        parser.error('Multiple exports require --distance-m and --offset-ns to create the velocity summary.')
    output = args.output_dir or paths[0].parent / 'timing_results'
    output.mkdir(parents=True, exist_ok=True)
    reports, failures, rejected = [], [], []
    for index, path in enumerate(paths, 1):
        print(f'\n[{index}/{len(paths)}] {path.name}')
        event_args = argparse.Namespace(**vars(args))
        event_args.input = path
        event_args.output_dir = output
        event_args.headless = not args.show
        event_args.show = args.show and not batch
        try:
            report = analyze_event(event_args)
            reports.append(report)
            selection = report['amplitude_selection'] or {'accepted': True, 'reason': 'No amplitude cut'}
        except AmplitudeRejected as error:
            selection = error.selection
            rejected.append({'event': path.stem, **selection})
            print(f'SKIPPED: {path.name}: {error}; peaks={selection["peak_amplitudes_mV"]}')
        except (ValueError, OSError, RuntimeError) as error:
            if not batch:
                raise
            failures.append({'input': str(path), 'error': str(error)})
            print(f'FAILED: {path.name}: {error}')
            continue
        (output / f'{path.stem}_selection.json').write_text(json.dumps(selection, indent=2)+'\n')
    if batch:
        build_summary(reports, failures, args, output, rejected)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
