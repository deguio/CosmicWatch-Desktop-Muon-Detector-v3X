#!/usr/bin/env python3
"""Offline timing of a two-channel RTO6 CSV export; all internal units are ns/mV."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
from scipy.special import log_ndtr


PARAMETERS = ['baseline_mV', 'scale_mV', 't0_ns', 'tau_rise_ns', 'tau_decay_ns']


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='RTO metadata .csv or companion .Wfm.csv')
    parser.add_argument('--fraction', type=float, default=.2)
    parser.add_argument('--smooth-samples', type=int, default=5, help='Odd Savitzky-Golay window; 1 disables smoothing')
    parser.add_argument('--baseline-end-ns', type=float, help='End of pre-pulse baseline; default first 15%% of record')
    parser.add_argument('--fit-range-ns', type=float, nargs=2, metavar=('START', 'STOP'))
    parser.add_argument('--output-dir', type=Path, help='Default: timing_results next to input')
    parser.add_argument('--show', action='store_true', help='Also open matplotlib plots')
    parser.add_argument('--distance-m', type=float, help='Muon path length between detectors, not just vertical separation')
    parser.add_argument('--offset-ns', type=float, help='Calibrated instrumental t(C2)-t(C1) offset to subtract')
    args = parser.parse_args()
    if not 0 < args.fraction < 1:
        parser.error('--fraction must be between 0 and 1.')
    if args.smooth_samples != 1 and (args.smooth_samples < 3 or args.smooth_samples % 2 == 0):
        parser.error('--smooth-samples must be 1 or an odd integer >= 3.')
    if args.distance_m is not None and (args.distance_m <= 0 or args.offset_ns is None):
        parser.error('A positive --distance-m requires an explicit calibrated --offset-ns (0 if calibrated zero).')
    t, y, labels, metadata = load_rto(args.input)
    if args.smooth_samples >= len(t):
        parser.error('Smoothing window exceeds record length.')
    baseline_end = args.baseline_end_ns if args.baseline_end_ns is not None else t[0]+.15*np.ptp(t)
    mask = np.ones(len(t), dtype=bool)
    if args.fit_range_ns:
        mask = (t >= args.fit_range_ns[0]) & (t <= args.fit_range_ns[1])
        if mask.sum() < 30:
            parser.error('Fit interval needs at least 30 samples.')
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
    report = {'input': metadata, 'convention': f'deltaT = t({labels[1]}) - t({labels[0]})',
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
    if not args.show:
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


if __name__ == '__main__':
    main()
