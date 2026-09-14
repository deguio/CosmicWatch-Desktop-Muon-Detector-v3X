#!/usr/bin/env python3
"""Event display 2D per la torre CosmicWatch a nove scintillatori.

Le coincidenze sono ricostruite usando le colonne Time/Date assegnate dal PC,
perche' nel file ``coincidence_nine.txt`` il flag Coincident vale sempre zero.
Non essendoci informazione sulla posizione del passaggio nello scintillatore,
il grafico non attribuisce una singola pendenza alla traccia: mostra invece le
rette campionate che intersecano tutti e soli i rivelatori accesi. Per i sette
piani interni viene inoltre stimata l'efficienza con il
rapporto N(sotto + centro + sopra) / N(sotto + sopra), usando i vicini
immediati come telescopio di riferimento.

Esempi
-------
    python display.py
    python display.py --window-ms 1 --min-devices 3
    python display.py --event 12 --save evento_12.png --no-show

Nel display interattivo usare freccia destra/sinistra (oppure n/p) per cambiare
evento, Home/End per andare al primo/ultimo evento e q per chiudere.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import date
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import numpy as np


DEVICE_NAMES = (
    "UNO", "DUE", "TRE", "QUATTRO", "CINQUE",
    "SEI", "SETTE", "OTTO", "NOVE",
)

# Vista laterale: larghezza x spessore. La profondita' (5 cm) non e' visibile.
SCINTILLATOR_WIDTH_CM = 6.0
SCINTILLATOR_DEPTH_CM = 5.0
SCINTILLATOR_THICKNESS_CM = 1.1
DEFAULT_TOWER_HEIGHT_CM = 38.5


@dataclass(frozen=True)
class Hit:
    time_ns: int
    device: str
    event_number: int
    adc: int
    sipm_mv: float
    coincident_flag: bool


@dataclass(frozen=True)
class Coincidence:
    hits: tuple[Hit, ...]

    @property
    def start_ns(self) -> int:
        return min(hit.time_ns for hit in self.hits)

    @property
    def span_ms(self) -> float:
        return (max(hit.time_ns for hit in self.hits) - self.start_ns) / 1e6

    @property
    def multiplicity(self) -> int:
        return len(self.hits)


@dataclass(frozen=True)
class EfficiencyEstimate:
    device: str
    lower_device: str
    upper_device: str
    detected: int
    reference_events: int

    @property
    def efficiency(self) -> float:
        return self.detected / self.reference_events

    @property
    def standard_error(self) -> float:
        probability = self.efficiency
        return (probability * (1.0 - probability) / self.reference_events) ** 0.5


def _computer_timestamp_ns(time_text: str, date_text: str) -> int:
    """Converte Date e Time in nanosecondi, conservando le nove cifre del file."""
    try:
        day, month, year = (int(value) for value in date_text.split("/"))
        hour_text, minute_text, second_text = time_text.split(":")
        if "." in second_text:
            seconds_text, fraction = second_text.split(".", 1)
        else:
            seconds_text, fraction = second_text, ""
        fraction_ns = int((fraction + "000000000")[:9])
        seconds = (
            int(hour_text) * 3600 + int(minute_text) * 60 + int(seconds_text)
        )
        return (
            date(year, month, day).toordinal() * 86_400_000_000_000
            + seconds * 1_000_000_000
            + fraction_ns
        )
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"timestamp PC non valido: {date_text!r} {time_text!r}"
        ) from error


def load_hits(path: Path, min_sipm_mv: float = 0.0) -> tuple[list[Hit], int]:
    """Legge il file CosmicWatch; restituisce hit validi e righe scartate."""
    hits: list[Hit] = []
    skipped = 0

    with path.open("r", encoding="utf-8", errors="replace") as data_file:
        for line_number, line in enumerate(data_file, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            try:
                # Le prime sei colonne sono fisse; Name, Time e Date sono le
                # ultime tre anche se il firmware cambia il blocco dei sensori.
                if len(fields) < 9:
                    raise ValueError("numero di colonne insufficiente")
                device = fields[-3].upper()
                if device not in DEVICE_NAMES:
                    raise ValueError(f"device sconosciuto {device!r}")
                hit = Hit(
                    time_ns=_computer_timestamp_ns(fields[-2], fields[-1]),
                    device=device,
                    event_number=int(fields[0]),
                    adc=int(fields[3]),
                    sipm_mv=float(fields[4]),
                    coincident_flag=bool(int(fields[2])),
                )
            except (ValueError, IndexError) as error:
                skipped += 1
                print(
                    f"Attenzione: riga {line_number} ignorata ({error})",
                    file=sys.stderr,
                )
                continue
            if hit.sipm_mv >= min_sipm_mv:
                hits.append(hit)

    hits.sort(key=lambda hit: hit.time_ns)
    return hits, skipped


def find_coincidences(
    hits: list[Hit], window_ms: float = 1.0, min_devices: int = 2
) -> list[Coincidence]:
    """Trova gruppi disgiunti di device diversi entro una finestra stretta.

    La finestra parte dal primo hit, quindi per ogni gruppo vale rigorosamente
    max(t)-min(t) <= window. Se un device produce piu' hit nella stessa
    finestra viene tenuto quello con ampiezza SiPM maggiore. Quando una finestra
    forma una coincidenza, tutti i suoi hit vengono consumati: ogni misura puo'
    appartenere al massimo a un evento mostrato.
    """
    if window_ms <= 0:
        raise ValueError("--window-ms deve essere maggiore di zero")
    if not 2 <= min_devices <= len(DEVICE_NAMES):
        raise ValueError("--min-devices deve essere compreso tra 2 e 9")

    window_ns = round(window_ms * 1e6)
    coincidences: list[Coincidence] = []
    first = 0

    while first < len(hits):
        last = first + 1
        limit_ns = hits[first].time_ns + window_ns
        while last < len(hits) and hits[last].time_ns <= limit_ns:
            last += 1

        best_by_device: dict[str, Hit] = {}
        for hit in hits[first:last]:
            previous = best_by_device.get(hit.device)
            if previous is None or hit.sipm_mv > previous.sipm_mv:
                best_by_device[hit.device] = hit

        if len(best_by_device) >= min_devices:
            ordered_hits = tuple(
                sorted(best_by_device.values(), key=lambda hit: DEVICE_NAMES.index(hit.device))
            )
            coincidences.append(Coincidence(ordered_hits))
            first = last
        else:
            first += 1

    return coincidences


def estimate_detector_efficiencies(
    coincidences: list[Coincidence],
) -> list[EfficiencyEstimate]:
    """Stima l'efficienza dei sette scintillatori interni con il sandwich.

    Per il device i-esimo il campione di riferimento contiene gli eventi nei
    quali sono presenti i suoi vicini immediati i-1 e i+1. Il riferimento
    seleziona quindi una particella compatibile con l'attraversamento del piano
    centrale senza richiedere a priori che quest'ultimo abbia risposto.
    """
    estimates: list[EfficiencyEstimate] = []
    device_sets = [{hit.device for hit in event.hits} for event in coincidences]

    for level in range(1, len(DEVICE_NAMES) - 1):
        lower = DEVICE_NAMES[level - 1]
        device = DEVICE_NAMES[level]
        upper = DEVICE_NAMES[level + 1]
        reference_events = 0
        detected = 0
        for active_devices in device_sets:
            if lower in active_devices and upper in active_devices:
                reference_events += 1
                detected += device in active_devices
        if reference_events:
            estimates.append(
                EfficiencyEstimate(
                    device=device,
                    lower_device=lower,
                    upper_device=upper,
                    detected=detected,
                    reference_events=reference_events,
                )
            )
    return estimates


def print_summary(
    path: Path,
    hits: list[Hit],
    coincidences: list[Coincidence],
    efficiency_estimates: list[EfficiencyEstimate],
    window_ms: float,
    min_devices: int,
    skipped: int,
) -> None:
    detector_counts = Counter(hit.device for hit in hits)
    multiplicities = Counter(event.multiplicity for event in coincidences)
    firmware_coincident = sum(hit.coincident_flag for hit in hits)

    print(f"File: {path}")
    print(f"Hit selezionati: {len(hits)} (righe non valide: {skipped})")
    print(
        "Hit per device: "
        + ", ".join(f"{name}={detector_counts[name]}" for name in DEVICE_NAMES)
    )
    print(f"Flag Coincident=1 nel file: {firmware_coincident}")
    print(
        f"Coincidenze temporali: {len(coincidences)} "
        f"(finestra {window_ms:g} ms, almeno {min_devices} device)"
    )
    if multiplicities:
        print(
            "Molteplicita': "
            + ", ".join(
                f"M={multiplicity}: {multiplicities[multiplicity]}"
                for multiplicity in sorted(multiplicities)
            )
        )
    print("\nEfficienza di rivelazione (metodo dei vicini immediati):")
    print("  device     riferimento      rilevati/totali       efficienza")
    for estimate in efficiency_estimates:
        print(
            f"  {estimate.device:<9} "
            f"{estimate.lower_device:>7}+{estimate.upper_device:<7} "
            f"{estimate.detected:>7}/{estimate.reference_events:<7} "
            f"{100.0 * estimate.efficiency:8.2f} +/- "
            f"{100.0 * estimate.standard_error:.2f} %"
        )
    if not efficiency_estimates:
        print("  nessuna coppia di riferimento disponibile")


def _draw_possible_trajectories(
    ax: plt.Axes,
    active_levels: list[int],
    level_centres_cm: list[float],
    tower_height_cm: float,
) -> tuple[str, int]:
    """Disegna rette che colpiscono esattamente i rivelatori accesi.

    Il controllo di intersezione usa l'intero rettangolo 6 x 1.1 cm di ciascun
    piano, non soltanto il suo centro. Le rette sono poi prolungate al di sopra
    e al di sotto della torre per rendere evidente da dove entrano e dove escono.
    """
    first_level = min(active_levels)
    last_level = max(active_levels)
    if set(active_levels) != set(range(first_level, last_level + 1)):
        return "hole", 0

    z_lower = level_centres_cm[first_level]
    z_upper = level_centres_cm[last_level]
    half_width = SCINTILLATOR_WIDTH_CM / 2
    active_set = set(active_levels)
    plot_z_min = -2.5
    plot_z_max = tower_height_cm + 2.5
    segments = []

    # Campioniamo le due coordinate di attraversamento nei piani accesi estremi.
    # Ogni coppia identifica una retta; viene conservata soltanto se l'elenco
    # dei rettangoli intersecati coincide esattamente con quello osservato.
    endpoints = np.linspace(-half_width, half_width, 61)
    for x_lower in endpoints:
        for x_upper in endpoints:
            slope = (x_upper - x_lower) / (z_upper - z_lower)
            compatible = True
            for level, centre_z in enumerate(level_centres_cm):
                detector_bottom = centre_z - SCINTILLATOR_THICKNESS_CM / 2
                detector_top = centre_z + SCINTILLATOR_THICKNESS_CM / 2
                x_bottom = x_lower + slope * (detector_bottom - z_lower)
                x_top = x_lower + slope * (detector_top - z_lower)
                line_min_x = min(x_bottom, x_top)
                line_max_x = max(x_bottom, x_top)
                intersects = (
                    line_max_x >= -half_width and line_min_x <= half_width
                )
                if intersects != (level in active_set):
                    compatible = False
                    break
            if compatible:
                segments.append(
                    (
                        (x_lower + slope * (plot_z_min - z_lower), plot_z_min),
                        (x_lower + slope * (plot_z_max - z_lower), plot_z_max),
                    )
                )

    if not segments:
        return "none", 0

    # Mantiene leggibile e veloce il display anche quando le soluzioni sono molte.
    if len(segments) > 700:
        selected = np.linspace(0, len(segments) - 1, 700, dtype=int)
        segments = [segments[index] for index in selected]
    ax.add_collection(
        LineCollection(
            segments,
            colors="#18a999",
            linewidths=0.9,
            alpha=max(0.025, min(0.10, 12.0 / len(segments))),
            # Sopra ai rettangoli: si vede direttamente quali scintillatori
            # vengono intersecati e quali vengono evitati.
            zorder=4,
        )
    )
    return "drawn", len(segments)


def draw_event(
    ax: plt.Axes,
    info_ax: plt.Axes,
    event: Coincidence,
    event_index: int,
    event_count: int,
    tower_height_cm: float,
    normalizer: colors.Normalize,
    color_map,
) -> None:
    ax.clear()
    info_ax.clear()
    info_ax.axis("off")

    pitch_cm = (
        tower_height_cm - SCINTILLATOR_THICKNESS_CM
    ) / (len(DEVICE_NAMES) - 1)
    level_centres_cm = [
        level * pitch_cm + SCINTILLATOR_THICKNESS_CM / 2
        for level in range(len(DEVICE_NAMES))
    ]
    hit_by_device = {hit.device: hit for hit in event.hits}
    active_levels = [DEVICE_NAMES.index(hit.device) for hit in event.hits]
    trajectory_status, trajectory_count = _draw_possible_trajectories(
        ax, active_levels, level_centres_cm, tower_height_cm
    )

    for level, device in enumerate(DEVICE_NAMES):
        z_cm = level * pitch_cm
        hit = hit_by_device.get(device)
        facecolor = color_map(normalizer(hit.sipm_mv)) if hit else "#e8edf2"
        edgecolor = "#8b1e3f" if hit else "#65727c"
        ax.add_patch(
            Rectangle(
                (-SCINTILLATOR_WIDTH_CM / 2, z_cm),
                SCINTILLATOR_WIDTH_CM,
                SCINTILLATOR_THICKNESS_CM,
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=2.2 if hit else 1.1,
                zorder=3,
            )
        )
        ax.text(
            0,
            level_centres_cm[level],
            device,
            ha="center",
            va="center",
            fontsize=8.5,
            color="white" if hit else "#34414a",
            fontweight="bold",
            zorder=5,
        )

    ax.set_xlim(-4.2, 4.2)
    ax.set_ylim(-2.5, tower_height_cm + 2.5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("z [cm]")
    ax.set_xticks([-3, 0, 3])
    ax.set_yticks(np.arange(0, tower_height_cm + 0.1, 5.0))
    ax.grid(color="#d8dde1", linewidth=0.55, alpha=0.65)
    ax.set_axisbelow(True)
    ax.set_title("Vista laterale della torre", fontsize=12, pad=10)

    info_ax.text(
        0.0, 0.98, f"EVENTO {event_index + 1} / {event_count}",
        transform=info_ax.transAxes, ha="left", va="top",
        fontsize=14, fontweight="bold", color="#24323b",
    )
    info_ax.text(
        0.0, 0.925,
        f"Molteplicita': {event.multiplicity}\nSpan temporale: {event.span_ms:.3f} ms",
        transform=info_ax.transAxes, ha="left", va="top",
        fontsize=10.5, linespacing=1.5, color="#46545d",
    )
    info_ax.text(
        0.0, 0.82, "HIT", transform=info_ax.transAxes,
        ha="left", va="top", fontsize=10, fontweight="bold", color="#24323b",
    )
    y_position = 0.78
    for hit in reversed(event.hits):
        delta_us = (hit.time_ns - event.start_ns) / 1e3
        info_ax.text(
            0.0, y_position, hit.device,
            transform=info_ax.transAxes, ha="left", va="top",
            fontsize=10, fontweight="bold", color="#71162f",
        )
        info_ax.text(
            0.30, y_position,
            f"{hit.sipm_mv:6.1f} mV    Δt {delta_us:7.1f} µs",
            transform=info_ax.transAxes, ha="left", va="top",
            fontsize=9.5, color="#46545d", family="monospace",
        )
        y_position -= 0.047

    gap_cm = pitch_cm - SCINTILLATOR_THICKNESS_CM
    info_ax.text(
        0.0, 0.27,
        "GEOMETRIA\n"
        f"Torre: {tower_height_cm:g} cm\n"
        f"Scintillatore: {SCINTILLATOR_WIDTH_CM:g} × "
        f"{SCINTILLATOR_DEPTH_CM:g} × {SCINTILLATOR_THICKNESS_CM:g} cm³\n"
        f"Gap calcolato: {gap_cm:.3f} cm",
        transform=info_ax.transAxes, ha="left", va="top",
        fontsize=9.5, linespacing=1.45, color="#46545d",
    )
    if trajectory_status == "drawn":
        trajectory_note = (
            f"Banda verde: {trajectory_count} traiettorie campionate.\n"
            "Ognuna colpisce tutti e soli i rivelatori accesi\n"
            "ed e' prolungata oltre la torre."
        )
        note_color = "#087f73"
    elif trajectory_status == "hole":
        trajectory_note = (
            "Banda non disegnata: almeno un rivelatore\n"
            "intermedio non si e' acceso (possibile inefficienza)."
        )
        note_color = "#9a5b13"
    else:
        trajectory_note = (
            "Nessuna retta attraversa tutti e soli i rivelatori\n"
            "accesi con questa geometria."
        )
        note_color = "#a13737"
    info_ax.text(
        0.0, 0.10, trajectory_note,
        transform=info_ax.transAxes, ha="left", va="top",
        fontsize=9, linespacing=1.35, color=note_color,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": note_color,
              "edgecolor": "none", "alpha": 0.08},
    )


def show_events(
    coincidences: list[Coincidence],
    initial_event: int,
    tower_height_cm: float,
    save_path: Path | None,
    show: bool,
) -> None:
    if not coincidences:
        raise ValueError("nessuna coincidenza soddisfa i criteri richiesti")
    if not 0 <= initial_event < len(coincidences):
        raise ValueError(
            f"--event deve essere compreso tra 1 e {len(coincidences)}"
        )
    minimum_height = len(DEVICE_NAMES) * SCINTILLATOR_THICKNESS_CM
    if tower_height_cm < minimum_height:
        raise ValueError(
            f"--tower-height-cm deve essere almeno {minimum_height:g} cm"
        )

    all_amplitudes = np.asarray(
        [hit.sipm_mv for event in coincidences for hit in event.hits], dtype=float
    )
    vmax = max(float(np.percentile(all_amplitudes, 98)), 1.0)
    normalizer = colors.Normalize(vmin=0.0, vmax=vmax, clip=True)
    color_map = plt.get_cmap("plasma")

    fig, (ax, info_ax) = plt.subplots(
        1, 2, figsize=(10.6, 9.2),
        gridspec_kw={"width_ratios": (0.85, 1.15)},
    )
    fig.subplots_adjust(left=0.09, right=0.96, top=0.88, bottom=0.11, wspace=0.20)
    fig.suptitle(
        "CosmicWatch · Event display 2D",
        fontsize=17, fontweight="bold", color="#24323b",
    )
    colorbar = fig.colorbar(
        ScalarMappable(norm=normalizer, cmap=color_map),
        ax=info_ax, orientation="horizontal", fraction=0.035,
        pad=0.035, location="bottom", aspect=30,
    )
    colorbar.set_label("Ampiezza SiPM [mV]")
    state = {"index": initial_event}

    def redraw() -> None:
        draw_event(
            ax,
            info_ax,
            coincidences[state["index"]],
            state["index"],
            len(coincidences),
            tower_height_cm,
            normalizer,
            color_map,
        )
        fig.canvas.draw_idle()

    def on_key(event) -> None:
        if event.key in ("right", "n", " "):
            state["index"] = (state["index"] + 1) % len(coincidences)
        elif event.key in ("left", "p", "backspace"):
            state["index"] = (state["index"] - 1) % len(coincidences)
        elif event.key == "home":
            state["index"] = 0
        elif event.key == "end":
            state["index"] = len(coincidences) - 1
        elif event.key in ("q", "escape"):
            plt.close(fig)
            return
        else:
            return
        redraw()

    fig.canvas.mpl_connect("key_press_event", on_key)
    redraw()

    if save_path is not None:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=180, bbox_inches="tight")
        print(f"Display salvato in: {save_path}")
    if show:
        print("Comandi: ←/p precedente, →/n successivo, Home/End, q per uscire")
        plt.show()
    else:
        plt.close(fig)


def parse_arguments() -> argparse.Namespace:
    default_input = Path(__file__).with_name("coincidence_nine.txt")
    parser = argparse.ArgumentParser(
        description="Ricostruisce coincidenze temporali e mostra la torre 2D."
    )
    parser.add_argument(
        "input", nargs="?", type=Path, default=default_input,
        help=f"file dati (default: {default_input.name})",
    )
    parser.add_argument(
        "-w", "--window-ms", type=float, default=1.0,
        help="ampiezza massima del gruppo temporale in ms (default: 1)",
    )
    parser.add_argument(
        "-m", "--min-devices", type=int, default=2,
        help="numero minimo di device attraversati (default: 2)",
    )
    parser.add_argument(
        "--min-sipm-mv", type=float, default=0.0,
        help="scarta hit sotto questa ampiezza SiPM (default: 0)",
    )
    parser.add_argument(
        "--tower-height-cm", type=float, default=DEFAULT_TOWER_HEIGHT_CM,
        help="altezza totale della torre in cm (default: 38.5)",
    )
    parser.add_argument(
        "-e", "--event", type=int, default=1,
        help="evento iniziale, numerato da 1 (default: 1)",
    )
    parser.add_argument("--save", type=Path, help="salva l'evento iniziale in PNG/PDF")
    parser.add_argument(
        "--no-show", action="store_true",
        help="non apre la finestra (utile insieme a --save)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    try:
        hits, skipped = load_hits(args.input, args.min_sipm_mv)
        # L'efficienza usa il campione minimo a due piani. --min-devices filtra
        # soltanto gli eventi presentati nel display e non cambia la stima.
        all_coincidences = find_coincidences(hits, args.window_ms, min_devices=2)
        efficiency_estimates = estimate_detector_efficiencies(all_coincidences)
        coincidences = [
            event
            for event in all_coincidences
            if event.multiplicity >= args.min_devices
        ]
        print_summary(
            args.input,
            hits,
            coincidences,
            efficiency_estimates,
            args.window_ms,
            args.min_devices,
            skipped,
        )
        show_events(
            coincidences,
            initial_event=args.event - 1,
            tower_height_cm=args.tower_height_cm,
            save_path=args.save,
            show=not args.no_show,
        )
    except (OSError, ValueError) as error:
        print(f"Errore: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
