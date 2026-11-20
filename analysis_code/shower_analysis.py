#!/usr/bin/env python3
"""Ricerca di sciami estesi con N CosmicWatch affiancati (non impilati).

Ogni rivelatore è letto via USB da import_data.py. Il tempo PC (colonna Time)
ha una latenza USB di ~0.5-1 ms e jitter di ~ms; il timestamp interno
(Timestamp[s]) ha risoluzione 1 µs ma ogni rivelatore ha il proprio oscillatore
(offset diverso, deriva di qualche ppm). Si allinea quindi il clock interno di
ogni rivelatore al tempo PC:

    t_al = t_int + offset(t_int),   offset = mediana(t_PC - t_int) in bin di
                                    --align-bin-s secondi, interpolata.

La latenza USB comune a tutti i rivelatori si cancella nelle differenze.
La precisione dell'allineamento è stimata dividendo gli eventi pari/dispari.

Poi:
  1. distribuzione dei Δt fra coppie di rivelatori (tempo allineato), con fondo
     piatto stimato dalle bande laterali;
  2. cluster con span max(t)-min(t) <= W, molteplicità = rivelatori distinti;
  3. accidentali per ogni molteplicità k, in due modi:
       - analitico: sum_{|S|=k} k W^(k-1) prod_{i in S} R_i (Poisson, R_i W << 1)
       - off-time: si ripete la ricerca dopo aver traslato ogni rivelatore di
         ritardi diversi >> W (include non stazionarietà dei rate);
  4. spettro d'ampiezza delle particelle di sciame (coincidenze con un "tag"
     sopra --tag-mv, sottratte le off-time) ed efficienza della soglia;
  5. attese dallo spettro di densità al livello del mare
         N(>ρ) = K ρ^-γ  [h^-1],  K ≈ 540, γ ≈ 1.39  (ρ in m^-2)
     con P(hit | ρ) = 1 - exp(-ε ρ A) uguale per tutti (distanze di pochi m,
     molto minori del raggio di Molière ~ 80 m);
  6. con --positions: stesse attese con le posizioni reali e profilo laterale
     NKG (MC), ed eccesso per coppia in funzione della distanza. Senza
     posizioni si confronta comunque l'eccesso per coppia con un modello
     fattorizzato a_i·a_j (differenze di efficienza fra rivelatori).

Uso (dalla radice del repository, venv attivo):
    python analysis_code/shower_analysis.py \
        Data/showers_analysis/CW_seven_parallel.txt Data/showers_analysis/CW_seven_parallel_2.txt \
        --positions Data/showers_analysis/positions_IMG_6986.txt
"""

from __future__ import annotations

import argparse
import itertools
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
from scipy import integrate
from scipy.stats import binom

mycolors = ['#c70039', '#ff5733', '#ff8d1a', '#ffc300', '#eddd53', '#add45c',
            '#57c785', '#00baad', '#2a7b9b', '#3d3d6b', '#511849', '#900c3f']

# Riferimento fisso per il tempo PC: total_seconds() non dipende dall'unità
# interna del datetime (pandas >= 3 usa µs con questo formato, non ns).
EPOCH = pd.Timestamp('2000-01-01')

COLUMNS = ['ev', 'tint', 'coinc', 'adc', 'mv', 'dead', 'T', 'P', 'acc', 'gyro',
           'name', 'time', 'date']


# ---------------------------------------------------------------------------
# Lettura e allineamento
# ---------------------------------------------------------------------------
def load_file(path):
    """Legge un file USB di import_data.py (11 sensori + Name/Time/Date)."""
    df = pd.read_csv(path, sep='\t', comment='#', header=None, names=COLUMNS,
                     usecols=['ev', 'tint', 'adc', 'mv', 'dead', 'name', 'time', 'date'],
                     on_bad_lines='skip', engine='c')
    df = df.dropna(subset=['tint', 'mv', 'name', 'time', 'date'])
    # µs di risoluzione sono più che sufficienti per il tempo PC (jitter ~ms).
    stamp = pd.to_datetime(df['date'] + ' ' + df['time'].str.slice(0, 15),
                           format='%d/%m/%Y %H:%M:%S.%f', errors='coerce')
    df = df[stamp.notna()].copy()
    df['tpc'] = (stamp[stamp.notna()] - EPOCH).dt.total_seconds().to_numpy()
    df['name'] = df['name'].astype(str)
    df['file'] = os.path.basename(path)
    return df.drop(columns=['time', 'date'])


def _median_offset(x, d, bin_s):
    b = np.floor((x - x[0]) / bin_s).astype(np.int64)
    s = pd.DataFrame({'b': b, 'x': x, 'd': d}).groupby('b').agg(
        x=('x', 'mean'), d=('d', 'median'), n=('d', 'size'))
    s = s[s.n >= 10]
    return s.x.to_numpy(), s.d.to_numpy()


def align_clocks(df, bin_s, t_ref):
    """Aggiunge la colonna 'tal' (tempo interno allineato, s da t_ref).

    Ritorna un dict con, per rivelatore, curva di offset e precisione stimata.
    """
    info = {}
    df['tal'] = np.nan
    for name, g in df.groupby('name'):
        g = g.sort_values('tint')
        x = g.tint.to_numpy()
        d = g.tpc.to_numpy() - t_ref - x
        xc, dc = _median_offset(x, d, bin_s)
        xe, de = _median_offset(x[0::2], d[0::2], bin_s)
        xo, do = _median_offset(x[1::2], d[1::2], bin_s)
        grid = x[::20]
        # Ogni metà ha metà statistica: sigma(full) ~ rms(diff) / 2.
        sigma_us = np.std(np.interp(grid, xe, de) - np.interp(grid, xo, do)) / 2 * 1e6
        df.loc[g.index, 'tal'] = x + np.interp(x, xc, dc)
        slope = np.polyfit(xc, dc, 1)[0]
        info[name] = dict(x=xc, d=dc, sigma_us=sigma_us, drift_ppm=slope * 1e6,
                          resid_us=(d - np.interp(x, xc, dc)) * 1e6)
    return info


# ---------------------------------------------------------------------------
# Coincidenze
# ---------------------------------------------------------------------------
def pair_dt(times_by_det, max_abs_s):
    """Δt = t_B - t_A (s) per ogni coppia A<B, entro ±max_abs_s (vettoriale)."""
    out = {}
    names = sorted(times_by_det)
    for a, b in itertools.combinations(names, 2):
        ta, tb = times_by_det[a], times_by_det[b]
        lo = np.searchsorted(tb, ta - max_abs_s, side='left')
        hi = np.searchsorted(tb, ta + max_abs_s, side='right')
        cnt = hi - lo
        if cnt.sum() == 0:
            out[(a, b)] = np.zeros(0)
            continue
        ia = np.repeat(np.arange(len(ta)), cnt)
        ib = np.concatenate([np.arange(l, h) for l, h in zip(lo[cnt > 0], hi[cnt > 0])])
        out[(a, b)] = tb[ib] - ta[ia]
    return out


def find_clusters(t, det, window_s):
    """Cluster greedy: si apre al primo hit libero, si includono gli hit con
    t - t_start <= W. Ritorna (start_index, stop_index, molteplicità)."""
    order = np.argsort(t, kind='mergesort')
    t = t[order]
    det = det[order]
    cand = np.nonzero(np.diff(t) <= window_s)[0]
    starts, stops, mult = [], [], []
    consumed = -1
    n = len(t)
    for s in cand:
        if s <= consumed:
            continue
        j = s + 1
        while j < n and t[j] - t[s] <= window_s:
            j += 1
        m = len(np.unique(det[s:j]))
        if m >= 2:
            starts.append(s)
            stops.append(j)
            mult.append(m)
        consumed = j - 1
    return order, np.array(starts, int), np.array(stops, int), np.array(mult, int)


def multiplicity_counts(t, det, window_s, n_det):
    _, _, _, mult = find_clusters(t, det, window_s)
    return np.bincount(mult, minlength=n_det + 1)[:n_det + 1]


def offtime_counts(t, det_idx, window_s, n_det, n_shifts, base_shift_s=0.0371):
    """Accidentali: ogni rivelatore i traslato di s*i*base_shift (s=1..n_shifts)."""
    acc = np.zeros((n_shifts, n_det + 1))
    for s in range(1, n_shifts + 1):
        ts = t + det_idx * s * base_shift_s
        acc[s - 1] = multiplicity_counts(ts, det_idx, window_s, n_det)
    return acc.mean(axis=0), acc.std(axis=0, ddof=1) / np.sqrt(n_shifts)


def analytic_accidentals(rates, window_s, n_det):
    """Rate (Hz) di k-ple accidentali: sum_{|S|=k} k W^(k-1) prod R_i."""
    rates = np.asarray(rates)
    out = np.zeros(n_det + 1)
    for k in range(2, n_det + 1):
        out[k] = sum(k * window_s ** (k - 1) * np.prod(rates[list(c)])
                     for c in itertools.combinations(range(len(rates)), k))
    return out


# ---------------------------------------------------------------------------
# Attese dallo spettro di densità
# ---------------------------------------------------------------------------
def expected_shower_rates(n_det, area_m2, eff, K=540.0, gamma=1.39,
                          rho_break=1000.0, gamma_high=2.0):
    """Rate (h^-1) di eventi con esattamente k rivelatori colpiti su n_det.

    dN/dρ da N(>ρ) = K ρ^-γ (γ_high sopra rho_break). Densità uguale su tutti
    i rivelatori, Poisson per rivelatore con media ε ρ A.
    """
    from math import comb

    def dndrho(rho):
        if rho < rho_break:
            return K * gamma * rho ** (-gamma - 1)
        return K * rho_break ** (-gamma) * gamma_high * rho_break ** gamma_high \
            * rho ** (-gamma_high - 1)

    out = np.zeros(n_det + 1)
    for k in range(2, n_det + 1):
        def f(lnrho):
            rho = np.exp(lnrho)
            p = -np.expm1(-eff * rho * area_m2)
            return comb(n_det, k) * p ** k * (1 - p) ** (n_det - k) * dndrho(rho) * rho
        val = 0.0
        edges = np.log([1e-3, 1.0, 10.0, 100.0, rho_break, 1e4, 1e7])
        for lo, hi in zip(edges[:-1], edges[1:]):
            val += integrate.quad(f, lo, hi, limit=200)[0]
        out[k] = val
    return out


def save_figure(fig, outdir, name):
    """Salva la figura in PDF e in PNG (il PNG serve alla documentazione .md)."""
    fig.savefig(os.path.join(outdir, name + '.pdf'))
    fig.savefig(os.path.join(outdir, name + '.png'), dpi=150)


def read_positions(path):
    """File di testo 'label x_cm y_cm' (commenti con #). Ritorna (labels, xy_m)."""
    labels, xy = [], []
    with open(path) as fh:
        for line in fh:
            line = line.split('#', 1)[0].split()
            if len(line) >= 3:
                labels.append(line[0])
                xy.append((float(line[1]) / 100.0, float(line[2]) / 100.0))
    return labels, np.array(xy)


def nkg_shape(r, age, r_moliere):
    """Forma radiale NKG (non normalizzata): (r/rM)^(s-2) (1+r/rM)^(s-4.5)."""
    x = r / r_moliere
    return x ** (age - 2.0) * (1.0 + x) ** (age - 4.5)


def geometric_shower_rates(xy_m, area_m2, eff, K=540.0, gamma=1.39,
                           rho_break=1000.0, gamma_high=2.0, age=1.3,
                           r_moliere=80.0, n_samples=2_000_000, seed=1,
                           r_min_m=0.03):
    """Rate (h^-1) per molteplicità e per coppia con le posizioni reali.

    Si campiona la densità ρ_c al centro dell'array dallo spettro N(>ρ) (stessa
    normalizzazione del calcolo uniforme) e la distanza del core r dalla sua
    distribuzione condizionata. Per uno spettro di dimensioni a legge di
    potenza con indice γ e forma laterale fissa, p(r | ρ_c) ∝ r f(r)^γ,
    indipendente da ρ_c. Ogni rivelatore vede ρ_i = ρ_c f(r_i)/f(r).
    Le probabilità di molteplicità si calcolano esattamente (Poisson-binomiale),
    senza campionare gli hit. Ritorna anche il caso uniforme (ρ_i = ρ_c) con gli
    stessi campioni, per isolare l'effetto della geometria.
    """
    rng = np.random.default_rng(seed)
    n_det = len(xy_m)
    center = xy_m.mean(axis=0)
    # distribuzione di r: pdf in ln r ∝ r^2 f(r)^γ
    lnr = np.linspace(np.log(1e-3), np.log(20 * r_moliere), 20000)
    r_grid = np.exp(lnr)
    pdf = r_grid ** 2 * nkg_shape(r_grid, age, r_moliere) ** gamma
    cdf = np.concatenate([[0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(lnr))])
    cdf /= cdf[-1]
    ln_lo, ln_hi = np.log(1e-2), np.log(1e7)

    def dn_dlnrho(rho):
        low = K * gamma * rho ** (-gamma)
        high = K * rho_break ** (gamma_high - gamma) * gamma_high * rho ** (-gamma_high)
        return np.where(rho < rho_break, low, high)

    mult = np.zeros(n_det + 1)
    mult_uni = np.zeros(n_det + 1)
    pairs = np.zeros((n_det, n_det))
    chunk = 250_000
    done = 0
    while done < n_samples:
        n = min(chunk, n_samples - done)
        done += n
        rho_c = np.exp(rng.uniform(ln_lo, ln_hi, n))
        w = dn_dlnrho(rho_c) * (ln_hi - ln_lo) / n_samples
        r = np.exp(np.interp(rng.uniform(0, 1, n), cdf, lnr))
        phi = rng.uniform(0, 2 * np.pi, n)
        core = center + np.column_stack([r * np.cos(phi), r * np.sin(phi)])
        ri = np.hypot(core[:, None, 0] - xy_m[None, :, 0], core[:, None, 1] - xy_m[None, :, 1])
        ri = np.maximum(ri, r_min_m)
        rho_i = rho_c[:, None] * nkg_shape(ri, age, r_moliere) / nkg_shape(r, age, r_moliere)[:, None]
        p = -np.expm1(-eff * area_m2 * rho_i)
        P = np.zeros((n, n_det + 1))
        P[:, 0] = 1.0
        for i in range(n_det):
            P[:, 1:] = P[:, 1:] * (1 - p[:, i:i + 1]) + P[:, :-1] * p[:, i:i + 1]
            P[:, 0] *= 1 - p[:, i]
        mult += w @ P
        pc = -np.expm1(-eff * area_m2 * rho_c)
        mult_uni += w @ binom.pmf(np.arange(n_det + 1)[None, :], n_det, pc[:, None])
        pairs += (p * w[:, None]).T @ p
    np.fill_diagonal(pairs, np.nan)
    return mult, mult_uni, pairs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('files', nargs='+')
    ap.add_argument('-t', '--threshold-mv', type=float, default=30.0,
                    help='soglia SiPM per gli hit usati nelle coincidenze (default 30 mV)')
    ap.add_argument('-w', '--window-us', type=float, default=300.0,
                    help='span massimo del cluster W (default 300 µs)')
    ap.add_argument('--align-bin-s', type=float, default=30.0)
    ap.add_argument('--n-shifts', type=int, default=20)
    ap.add_argument('--tag-mv', type=float, default=50.0)
    ap.add_argument('--area-cm2', type=float, default=25.0)
    ap.add_argument('--K', type=float, default=540.0, help='N(>ρ)=K ρ^-γ [h^-1]')
    ap.add_argument('--gamma', type=float, default=1.39)
    ap.add_argument('--scan', type=str, default='0,15,20,30,40,50,75',
                    help='soglie (mV) per la scansione')
    ap.add_argument('-o', '--outdir', default=None,
                    help='cartella figure (default: <cartella del primo file>/Figures)')
    ap.add_argument('--positions', default=None,
                    help="file 'label x_cm y_cm' con le posizioni dei rivelatori; se le "
                         "etichette sono i nomi dei rivelatori si confronta coppia per coppia")
    ap.add_argument('--nkg-age', type=float, default=1.3)
    ap.add_argument('--moliere-m', type=float, default=80.0)
    ap.add_argument('--mc-samples', type=int, default=2_000_000)
    ap.add_argument('--show', action='store_true')
    args = ap.parse_args()
    if not args.show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'

    outdir = args.outdir or os.path.join(os.path.dirname(os.path.abspath(args.files[0])), 'Figures')
    os.makedirs(outdir, exist_ok=True)
    W = args.window_us * 1e-6

    # --- lettura e allineamento, file per file --------------------------------
    frames, align_info, live = [], {}, []
    t_ref = None
    for path in args.files:
        print('Lettura', path)
        df = load_file(path)
        if t_ref is None:
            t_ref = df.tpc.min()
        # Dopo l'allineamento, 'tal' è sulla scala PC (s da t_ref): i file in
        # date diverse si concatenano senza sovrapposizioni.
        align_info[os.path.basename(path)] = align_clocks(df, args.align_bin_s, t_ref)
        t0 = df.groupby('name').tal.min().max()
        t1 = df.groupby('name').tal.max().min()
        df = df[(df.tal >= t0) & (df.tal <= t1)]
        live.append((os.path.basename(path), t0, t1))
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    names = sorted(df.name.unique())
    n_det = len(names)
    det_index = {n: i for i, n in enumerate(names)}
    df['di'] = df.name.map(det_index).astype(int)
    T_tot = sum(t1 - t0 for _, t0, t1 in live)

    print('\n=== Dati ===')
    for f, t0, t1 in live:
        print('  %-28s  durata comune %.2f h' % (f, (t1 - t0) / 3600))
    print('  Totale: %.2f h, %d rivelatori: %s' % (T_tot / 3600, n_det, ', '.join(names)))
    print('\n  %-8s %9s %8s %10s %10s %9s' % ('det', 'N', 'R[Hz]', 'R>thr[Hz]', 'dead[%]', 'sync[µs]'))
    for n in names:
        g = df[df.name == n]
        dead = sum((h.dead.max() - h.dead.min()) / (h.tint.max() - h.tint.min())
                   for _, h in g.groupby('file')) / g.file.nunique()
        sig = np.mean([align_info[f][n]['sigma_us'] for f in align_info])
        print('  %-8s %9d %8.3f %10.3f %10.2f %9.1f' % (
            n, len(g), len(g) / T_tot, (g.mv > args.threshold_mv).sum() / T_tot,
            100 * dead, sig))

    # --- figura: offset dei clock ---------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    first = list(align_info)[0]
    for i, n in enumerate(names):
        a = align_info[first][n]
        lin = np.polyval(np.polyfit(a['x'], a['d'], 1), a['x'])
        axs[0].plot(a['x'] / 3600, (a['d'] - lin) * 1e3, lw=1, color=mycolors[i % 12],
                    label='%s (%.2f ppm)' % (n, a['drift_ppm']))
        axs[1].hist(a['resid_us'] / 1e3, bins=np.arange(-2, 6, 0.05), histtype='step',
                    color=mycolors[i % 12])
    axs[0].set_xlabel('t interno [h]')
    axs[0].set_ylabel('offset - fit lineare [ms]')
    axs[0].set_title(first, fontsize=10)
    axs[0].legend(fontsize=8)
    axs[1].set_yscale('log')
    axs[1].set_xlabel(r'$t_{PC}-t_{al}$ [ms] (latenza USB rispetto alla mediana)')
    for ax in axs:
        ax.grid(ls='--', alpha=0.5)
    fig.tight_layout()
    save_figure(fig, outdir, 'shower_clock_alignment')

    # --- Δt fra coppie --------------------------------------------------------
    def times_above(thr):
        s = df[df.mv > thr]
        return {n: np.sort(s.tal.to_numpy()[s.name.to_numpy() == n]) for n in names}

    max_dt = 3e-3
    dts = pair_dt(times_above(args.threshold_mv), max_dt)
    all_dt = np.concatenate(list(dts.values())) * 1e6
    side = (np.abs(all_dt) > 1000) & (np.abs(all_dt) < 3000)
    flat_per_us = side.sum() / 4000.0
    in_win = np.abs(all_dt) <= args.window_us
    exc_pairs = in_win.sum() - flat_per_us * 2 * args.window_us
    print('\n=== Δt fra coppie (SiPM > %.0f mV) ===' % args.threshold_mv)
    print('  |Δt| <= %.0f µs: %d coppie, fondo piatto %.1f, eccesso %.1f ± %.1f' % (
        args.window_us, in_win.sum(), flat_per_us * 2 * args.window_us, exc_pairs,
        np.sqrt(in_win.sum())))
    sig_dt = all_dt[in_win]
    print('  RMS dei Δt nella finestra (fondo incluso): %.0f µs' % sig_dt.std())

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, rng, bw in [(axs[0], 3000, 50), (axs[1], 800, 20)]:
        bins = np.arange(-rng, rng + bw / 2, bw)
        ax.hist(all_dt, bins=bins, histtype='stepfilled', color=mycolors[8], alpha=0.6,
                label='tutte le coppie, SiPM > %.0f mV' % args.threshold_mv)
        ax.axhline(flat_per_us * bw, color='k', ls=':', label='fondo dalle bande laterali')
        ax.axvline(-args.window_us, color='k', ls='--')
        ax.axvline(args.window_us, color='k', ls='--', label=r'$\pm W$')
        ax.set_xlabel(r'$\Delta t = t_B - t_A$ [µs]')
        ax.set_ylabel('coppie / %d µs' % bw)
        ax.grid(ls='--', alpha=0.5)
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, outdir, 'shower_pair_dt')

    # eccesso per coppia
    pair_exc = np.full((n_det, n_det), np.nan)
    pair_err = np.full((n_det, n_det), np.nan)
    print('\n  Eccesso per coppia (|Δt|<=W, fondo dalle bande laterali):')
    for (a, b), v in dts.items():
        v = v * 1e6
        bkg = ((np.abs(v) > 1000) & (np.abs(v) < 3000)).sum() / 4000.0 * 2 * args.window_us
        e = (np.abs(v) <= args.window_us).sum() - bkg
        pair_exc[det_index[a], det_index[b]] = pair_exc[det_index[b], det_index[a]] = e
        # errore: Poisson sul conteggio in finestra (almeno 1, così una coppia senza
        # conteggi non ha errore nullo) + fondo dalle bande (4000 µs)
        n_side = ((np.abs(v) > 1000) & (np.abs(v) < 3000)).sum()
        pair_err[det_index[a], det_index[b]] = pair_err[det_index[b], det_index[a]] = np.sqrt(
            max((np.abs(v) <= args.window_us).sum(), 1) + n_side * (2 * args.window_us / 4000.0) ** 2)
    for i, a in enumerate(names):
        print('   %-8s ' % a + ' '.join('%6.1f' % pair_exc[i, j] if i != j else '     -'
                                         for j in range(n_det)))
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(pair_exc, cmap='viridis')
    ax.set_xticks(range(n_det), names, rotation=45, fontsize=8)
    ax.set_yticks(range(n_det), names, fontsize=8)
    for i in range(n_det):
        for j in range(n_det):
            if i != j:
                ax.text(j, i, '%.0f' % pair_exc[i, j], ha='center', va='center',
                        color='w', fontsize=7)
    fig.colorbar(im, label='coppie in eccesso')
    fig.tight_layout()
    save_figure(fig, outdir, 'shower_pair_excess')

    # --- cluster e accidentali per molteplicità -------------------------------
    def analyse_threshold(thr, n_shifts):
        s = df[df.mv > thr]
        t, di = s.tal.to_numpy(), s.di.to_numpy()
        obs = multiplicity_counts(t, di, W, n_det)
        off, off_err = offtime_counts(t, di, W, n_det, n_shifts)
        ana = np.zeros(n_det + 1)
        for f, t0, t1 in live:
            sf = s[s.file == f]
            rates = [(sf.di == i).sum() / (t1 - t0) for i in range(n_det)]
            ana += analytic_accidentals(rates, W, n_det) * (t1 - t0)
        return obs, off, off_err, ana, s

    obs, off, off_err, ana, sel = analyse_threshold(args.threshold_mv, args.n_shifts)
    print('\n=== Cluster con span <= %.0f µs, SiPM > %.0f mV, %.2f h ===' % (
        args.window_us, args.threshold_mv, T_tot / 3600))
    print('  %3s %8s %12s %12s %16s %12s' % ('k', 'osservati', 'acc.off-time', 'acc.analit.',
                                           'eccesso', 'ecc. [1/h]'))
    excess = obs - off
    excess_err = np.sqrt(obs + off_err ** 2)
    for k in range(2, n_det + 1):
        print('  %3d %8d %12.2f %12.2f %9.1f ± %5.1f %12.3f' % (
            k, obs[k], off[k], ana[k], excess[k], excess_err[k], excess[k] / T_tot * 3600))

    # --- spettro delle particelle di sciame (tag) -----------------------------
    tag = df[df.mv > args.tag_mv]
    print('\n=== Spettro dei partner di un tag (SiPM_tag > %.0f mV, |Δt| <= W) ===' % args.tag_mv)
    bins_mv = np.concatenate([np.arange(0, 100, 5), np.arange(100, 400, 25)])
    h_on = np.zeros(len(bins_mv) - 1)
    h_off = np.zeros(len(bins_mv) - 1)
    off_delays = np.array([k * 0.0137 for k in range(-10, 11) if k != 0])
    for n in names:
        tt = np.sort(tag.tal.to_numpy()[tag.name.to_numpy() == n])
        for m in names:
            if m == n:
                continue
            o = df[df.name == m].sort_values('tal')
            tp, mvp = o.tal.to_numpy(), o.mv.to_numpy()
            for delays, h, wgt in [(np.zeros(1), h_on, 1.0), (off_delays, h_off, 1.0 / len(off_delays))]:
                for dl in delays:
                    lo = np.searchsorted(tp, tt + dl - W, side='left')
                    hi = np.searchsorted(tp, tt + dl + W, side='right')
                    idx = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi) if b > a] or [np.zeros(0, int)])
                    h += np.histogram(mvp[idx], bins=bins_mv)[0] * wgt
    sig_h = h_on - h_off
    tot = sig_h.sum()
    centers = 0.5 * (bins_mv[1:] + bins_mv[:-1])
    print('  partner on-time %.0f, off-time %.1f, eccesso %.1f' % (h_on.sum(), h_off.sum(), tot))
    eff_thr = {}
    for thr in [float(x) for x in args.scan.split(',')]:
        eff_thr[thr] = sig_h[bins_mv[:-1] >= thr].sum() / tot if tot > 0 else np.nan
        print('  frazione dei partner di sciame sopra %4.0f mV: %.2f' % (thr, eff_thr[thr]))
    eff = eff_thr.get(args.threshold_mv,
                      sig_h[bins_mv[:-1] >= args.threshold_mv].sum() / tot)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    widths = np.diff(bins_mv)
    ax.step(centers, h_on / widths, where='mid', color=mycolors[0], label='on-time')
    ax.step(centers, h_off / widths, where='mid', color='gray', label='off-time (accidentali)')
    ax.errorbar(centers, sig_h / widths, yerr=np.sqrt(h_on + h_off / len(off_delays)) / widths,
                fmt='o', ms=3, color=mycolors[8], label='on - off (particelle di sciame)')
    ax.axvline(args.threshold_mv, color='k', ls='--', label='soglia %.0f mV' % args.threshold_mv)
    ax.set_xscale('log')
    ax.set_xlim(3, 400)
    ax.set_xlabel('SiPM partner [mV]')
    ax.set_ylabel('hit / mV')
    ax.legend(fontsize=8)
    ax.grid(ls='--', alpha=0.5)
    fig.tight_layout()
    save_figure(fig, outdir, 'shower_partner_spectrum')

    # --- attese ----------------------------------------------------------------
    A = args.area_cm2 * 1e-4
    print('\n=== Attese da N(>ρ) = %.0f ρ^-%.2f h^-1, A = %.0f cm², %d rivelatori ===' % (
        args.K, args.gamma, args.area_cm2, n_det))
    exp_rows = {}
    for label, e_, g_ in [('ε=1, γ=%.2f' % args.gamma, 1.0, args.gamma),
                          ('ε=%.2f, γ=%.2f' % (eff, args.gamma), eff, args.gamma),
                          ('ε=%.2f, γ=1.30' % eff, eff, 1.30),
                          ('ε=%.2f, γ=1.50' % eff, eff, 1.50)]:
        exp_rows[label] = expected_shower_rates(n_det, A, e_, args.K, g_)
    print('  %3s %12s ' % ('k', 'oss. [1/h]') + ' '.join('%16s' % l for l in exp_rows))
    for k in range(2, n_det + 1):
        print('  %3d %5.3f±%5.3f ' % (k, excess[k] / T_tot * 3600, excess_err[k] / T_tot * 3600)
              + ' '.join('%16.4f' % r[k] for r in exp_rows.values()))
    ob2 = excess[2]
    if ob2 > 0:
        print('\n  Rapporti indipendenti dalla normalizzazione K (CREDO):')
        for k in range(3, min(n_det, 5) + 1):
            preds = ', '.join('%.3f' % (r[k] / r[2]) for r in exp_rows.values())
            print('   R%d/R2 oss. = %.3f ± %.3f   attesi: %s' % (
                k, excess[k] / ob2, excess[k] / ob2 * np.sqrt(
                    (excess_err[k] / max(excess[k], 1e-9)) ** 2 + (excess_err[2] / ob2) ** 2), preds))

    # --- coppie: modello fattorizzato (non richiede le posizioni) --------------
    # Per sciami con densità uniforme l'eccesso nella coppia (i,j) è ∝ ε_i ε_j:
    # le differenze di guadagno fra rivelatori spiegano già una parte della
    # dispersione. Si confronta un eccesso costante con a_i a_j.
    from scipy.optimize import least_squares
    iu = np.triu_indices(n_det, 1)
    y, sy = pair_exc[iu], pair_err[iu]
    w_mean = np.sum(y / sy ** 2) / np.sum(1 / sy ** 2)
    chi2_const = np.sum(((y - w_mean) / sy) ** 2)
    fit = least_squares(lambda a: (y - (a[:, None] * a[None, :])[iu]) / sy,
                        np.full(n_det, np.sqrt(max(w_mean, 1.0))))
    a_fit = fit.x
    model_fact = (a_fit[:, None] * a_fit[None, :])[iu]
    chi2_fact = np.sum(((y - model_fact) / sy) ** 2)
    from scipy.stats import chi2 as chi2_dist
    print('\n=== Eccesso per coppia: costante vs fattorizzato a_i·a_j ===')
    print('  costante:     χ²/ndf = %.1f/%d  (p = %.3f)' % (
        chi2_const, len(y) - 1, chi2_dist.sf(chi2_const, len(y) - 1)))
    print('  fattorizzato: χ²/ndf = %.1f/%d  (p = %.3f)' % (
        chi2_fact, len(y) - n_det, chi2_dist.sf(chi2_fact, len(y) - n_det)))
    rate_thr = np.array([(df[df.name == n].mv > args.threshold_mv).sum() / T_tot for n in names])
    print('  %-8s %8s %12s' % ('det', 'a_i', 'R>thr [Hz]'))
    for n, a_, r_ in zip(names, a_fit, rate_thr):
        print('  %-8s %8.2f %12.3f' % (n, a_, r_))

    # --- geometria: posizioni reali e profilo laterale NKG ---------------------
    geo = None
    if args.positions:
        labels, xy = read_positions(args.positions)
        if len(labels) != n_det:
            print('\n!! %s contiene %d posizioni, i rivelatori sono %d: salto la geometria'
                  % (args.positions, len(labels), n_det))
        else:
            geo, geo_uni, geo_pairs = geometric_shower_rates(
                xy, A, eff, args.K, args.gamma, age=args.nkg_age,
                r_moliere=args.moliere_m, n_samples=args.mc_samples)
            dist = np.hypot(*(xy[:, None, :] - xy[None, :, :]).transpose(2, 0, 1))
            print('\n=== Attese con le posizioni reali (%s) ===' % os.path.basename(args.positions))
            print('  NKG s = %.2f, r_M = %.0f m, ε = %.2f, γ = %.2f; distanze %.0f-%.0f cm'
                  % (args.nkg_age, args.moliere_m, eff, args.gamma,
                     dist[iu].min() * 100, dist[iu].max() * 100))
            print('  %3s %14s %12s %12s %10s %10s' % ('k', 'oss. [1/h]', 'uniforme', 'geometria',
                                                    'geo/uni', 'oss/geo'))
            for k in range(2, n_det + 1):
                print('  %3d %6.3f±%5.3f %12.4f %12.4f %10.3f %10.2f' % (
                    k, excess[k] / T_tot * 3600, excess_err[k] / T_tot * 3600, geo_uni[k],
                    geo[k], geo[k] / geo_uni[k], excess[k] / T_tot * 3600 / geo[k]))
            k3 = slice(3, None)
            print('  k>=3: oss. %.3f ± %.3f /h, atteso %.3f /h  (rapporto %.1f)' % (
                excess[k3].sum() / T_tot * 3600, np.sqrt(np.sum(excess_err[k3] ** 2)) / T_tot * 3600,
                geo[k3].sum(), excess[k3].sum() / T_tot * 3600 / geo[k3].sum()))
            gp = geo_pairs[iu] * T_tot / 3600
            print('  coppie (inclusive): attese %.1f-%.1f per coppia in %.1f h (variazione '
                  'con la distanza %.1f%%); osservate in media %.1f'
                  % (gp.min(), gp.max(), T_tot / 3600, 100 * (gp.max() / gp.min() - 1),
                     y.mean()))

            mapped = set(labels) == set(names)
            fig, ax = plt.subplots(figsize=(7, 4.5))
            if mapped:
                order_lab = [labels.index(n) for n in names]
                xy_names = xy[order_lab]

                # Fit congiunto y_ij = a_i a_j (1 + β (d_ij - <d>)): efficienze e
                # dipendenza dalla distanza insieme (un fit in due passi farebbe
                # assorbire parte dell'effetto agli a_i).
                def joint_fit(xy_):
                    d_ = np.hypot(*(xy_[iu[0]] - xy_[iu[1]]).T)
                    dc = d_ - d_.mean()
                    r = least_squares(
                        lambda p: (y - p[iu[0]] * p[iu[1]] * (1 + p[-1] * dc)) / sy,
                        np.r_[a_fit, 0.0])
                    cov = np.linalg.inv(r.jac.T @ r.jac)
                    return r.x, np.sqrt(cov[-1, -1]), np.sum(r.fun ** 2), d_
                p_j, beta_err, chi2_j, d_m = joint_fit(xy_names)
                beta = p_j[-1]
                model_j = p_j[iu[0]] * p_j[iu[1]] * (1 + beta * (d_m - d_m.mean()))
                print('\n  Fit congiunto a_i·a_j·(1 + β(d - %.2f m)):' % d_m.mean())
                print('   β = %.2f ± %.2f per m (atteso da EAS: ≈ %.3f per m)'
                      % (beta, beta_err, (gp.max() / gp.min() - 1) / -np.ptp(dist[iu])))
                print('   χ²/ndf = %.1f/%d (senza β: %.1f/%d), Δχ² = %.1f'
                      % (chi2_j, len(y) - n_det - 1, chi2_fact, len(y) - n_det,
                         chi2_fact - chi2_j))
                # Test a permutazioni: tutte le assegnazioni nome -> posizione.
                # Non dipende dal modello degli errori.
                betas = np.array([joint_fit(xy_names[list(perm)])[0][-1]
                                  for perm in itertools.permutations(range(n_det))])
                print('   permutazioni nome->posizione: frazione con β <= osservato = %.4f '
                      '(%d permutazioni, β medio %.2f, dispersione %.2f)'
                      % (np.mean(betas <= beta), len(betas), betas.mean(), betas.std()))
                wts = 1 / sy ** 2
                d_edges = np.quantile(d_m, [0, 1 / 3, 2 / 3, 1]) + np.array([-1e-6, 0, 0, 1e-6])
                print('   eccesso medio per coppia in classi di distanza (EAS atteso %.1f):'
                      % gp.mean())
                bin_c, bin_y, bin_e = [], [], []
                for lo, hi in zip(d_edges[:-1], d_edges[1:]):
                    m_ = (d_m > lo) & (d_m <= hi)
                    bin_c.append(np.mean(d_m[m_]) * 100)
                    bin_y.append(np.sum(wts[m_] * y[m_]) / wts[m_].sum())
                    bin_e.append(1 / np.sqrt(wts[m_].sum()))
                    print('    %4.0f-%4.0f cm (%d coppie): %5.1f ± %4.1f'
                          % (lo * 100, hi * 100, m_.sum(), bin_y[-1], bin_e[-1]))
                ax.errorbar(d_m * 100, y, yerr=sy, fmt='o', color='k', alpha=0.5,
                            label='osservato (coppia)')
                ax.errorbar(bin_c, bin_y, yerr=bin_e, fmt='D', ms=7, color=mycolors[8],
                            label='osservato, media pesata in classi di distanza')
                ax.plot(d_m * 100, model_j, 's', mfc='none', color=mycolors[6],
                        label=r'fit $a_i a_j\,(1+\beta(d-\langle d\rangle))$, '
                              r'$\beta$=%.2f$\pm$%.2f/m' % (beta, beta_err))
                for (i, j), d_, y_, s_, m_, mj in zip(zip(*iu), d_m, y, sy, model_fact, model_j):
                    print('   %-8s-%-8s d=%5.0f cm  oss %5.1f ± %4.1f  fatt. %5.1f  '
                          'fatt.+β %5.1f' % (names[i], names[j], d_ * 100, y_, s_, m_, mj))
            else:
                ax.axhspan(y.mean() - y.std(), y.mean() + y.std(), color='gray', alpha=0.3,
                           label='osservato: media ± dispersione (assegnazione ignota)')
                print('  Etichette %s diverse dai nomi dei rivelatori: per il confronto '
                      'coppia per coppia scrivere i nomi nel file delle posizioni.' % labels)
            ax.plot(dist[iu] * 100, gp, '^', color=mycolors[0],
                    label='atteso EAS (NKG, posizioni reali)')
            ax.set_xlabel('distanza fra i rivelatori [cm]')
            ax.set_ylabel('coppie in eccesso in %.1f h' % (T_tot / 3600))
            ax.set_ylim(0, None)
            ax.legend(fontsize=8)
            ax.grid(ls='--', alpha=0.5)
            fig.tight_layout()
            save_figure(fig, outdir, 'shower_pair_vs_distance')

    # --- scansione in soglia ---------------------------------------------------
    print('\n=== Scansione in soglia (W = %.0f µs, accidentali off-time) ===' % args.window_us)
    print('  %6s ' % 'thr' + ' '.join('  k=%d: oss / acc  ' % k for k in range(2, 6)))
    scan_rows = []
    for thr in [float(x) for x in args.scan.split(',')]:
        o_, f_, fe_, a_, _ = analyse_threshold(thr, max(5, args.n_shifts // 4))
        scan_rows.append((thr, o_, f_))
        print('  %6.0f ' % thr + ' '.join('%7d / %8.1f ' % (o_[k], f_[k]) for k in range(2, 6)))

    # --- figura molteplicità ---------------------------------------------------
    ks = np.arange(2, n_det + 1)
    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.errorbar(ks, obs[2:] / T_tot * 3600, yerr=np.sqrt(obs[2:]) / T_tot * 3600, fmt='o',
                color='k', label='osservati')
    ax.step(ks, off[2:] / T_tot * 3600, where='mid', color='gray', label='accidentali (off-time)')
    ax.plot(ks, ana[2:] / T_tot * 3600, 'x', color='gray', label='accidentali (analitici)')
    ax.errorbar(ks + 0.08, np.clip(excess[2:], 1e-9, None) / T_tot * 3600,
                yerr=excess_err[2:] / T_tot * 3600, fmt='s', color=mycolors[0],
                label='osservati - accidentali')
    for i, (label, r) in enumerate(exp_rows.items()):
        ax.plot(ks, r[2:], '-', color=mycolors[(2 * i + 6) % 12], label='atteso: ' + label)
    if geo is not None:
        ax.plot(ks, geo[2:], '--', color='k', label='atteso: ε=%.2f, posizioni reali (NKG)' % eff)
    ax.set_yscale('log')
    ax.set_ylim(1e-4, None)
    ax.set_xlabel('molteplicità k (rivelatori colpiti)')
    ax.set_ylabel('eventi / h')
    ax.set_title('SiPM > %.0f mV, W = %.0f µs, %.1f h' % (args.threshold_mv, args.window_us,
                                                          T_tot / 3600), fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(ls='--', alpha=0.5, which='both')
    fig.tight_layout()
    save_figure(fig, outdir, 'shower_multiplicity')

    # --- elenco eventi ad alta molteplicità e uniformità temporale -------------
    t, di, mv = sel.tal.to_numpy(), sel.di.to_numpy(), sel.mv.to_numpy()
    order, st, sp, mult = find_clusters(t, di, W)
    print('\n=== Eventi con k >= 4 ===')
    stamps = EPOCH + pd.to_timedelta(t_ref + t[order][st], unit='s')
    n_printed = 0
    for s0, s1, m, ts in zip(st, sp, mult, stamps):
        if m >= 4 and n_printed < 100:
            n_printed += 1
            ii = order[s0:s1]
            span = (t[ii].max() - t[ii].min()) * 1e6
            hits = ', '.join('%s:%.0f' % (names[d], v) for d, v in zip(di[ii], mv[ii]))
            print('  %s  k=%d  span=%4.0f µs  [%s]' % (ts.strftime('%d/%m %H:%M:%S'), m, span, hits))
    # test di uniformità nel tempo dei cluster k>=2 rispetto agli accidentali
    from scipy.stats import kstest
    tc = np.sort(t[order][st[mult >= 3]])
    if len(tc) > 5:
        # mappa su tempo "vivo" cumulativo
        u = np.zeros_like(tc)
        acc = 0.0
        for (f, t0, t1) in live:
            m_ = (tc >= t0) & (tc <= t1)
            u[m_] = (tc[m_] - t0 + acc) / T_tot
            acc += t1 - t0
        print('\n  KS uniformità temporale dei cluster k>=3: p = %.2f (N=%d)' % (
            kstest(u, 'uniform').pvalue, len(tc)))

    print('\nFigure salvate in', outdir)
    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
