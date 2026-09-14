#!/usr/bin/env python3
"""Event display 2D per la torre CosmicWatch a nove scintillatori.

Le coincidenze sono ricostruite usando le colonne Time/Date assegnate dal PC,
perche' nel file ``coincidence_nine.txt`` il flag Coincident vale sempre zero.
Non essendoci informazione sulla posizione del passaggio nello scintillatore,
il grafico mostra i piani attraversati ma non attribuisce una pendenza alla
traccia. Per i sette piani interni viene inoltre stimata l'efficienza con il
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
DEFAULT_GAP_CM = 5.0


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


def draw_event(
    ax: plt.Axes,
    event: Coincidence,
    event_index: int,
    event_count: int,
    gap_cm: float,
    normalizer: colors.Normalize,
    color_map,
) -> None:
    ax.clear()
    pitch_cm = SCINTILLATOR_THICKNESS_CM + gap_cm
    hit_by_device = {hit.device: hit for hit in event.hits}
    active_levels = [DEVICE_NAMES.index(hit.device) for hit in event.hits]

    # La fascia indica tutte le coordinate x compatibili, non una traiettoria
    # ricostruita. I nove piani forniscono infatti solo una risposta binaria.
    lower_z = min(active_levels) * pitch_cm
    upper_z = max(active_levels) * pitch_cm + SCINTILLATOR_THICKNESS_CM
    ax.add_patch(
        Rectangle(
            (-SCINTILLATOR_WIDTH_CM / 2, lower_z),
            SCINTILLATOR_WIDTH_CM,
            upper_z - lower_z,
            facecolor="#f6c453",
            edgecolor="none",
            alpha=0.10,
            zorder=0,
        )
    )

    for level, device in enumerate(DEVICE_NAMES):
        z_cm = level * pitch_cm
        hit = hit_by_device.get(device)
        facecolor = color_map(normalizer(hit.sipm_mv)) if hit else "#e4e9ee"
        linewidth = 2.0 if hit else 1.0
        edgecolor = "#8b1e3f" if hit else "#5d6973"
        ax.add_patch(
            Rectangle(
                (-SCINTILLATOR_WIDTH_CM / 2, z_cm),
                SCINTILLATOR_WIDTH_CM,
                SCINTILLATOR_THICKNESS_CM,
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=linewidth,
                zorder=2,
            )
        )
        ax.text(
            -SCINTILLATOR_WIDTH_CM / 2 - 0.35,
            z_cm + SCINTILLATOR_THICKNESS_CM / 2,
            device,
            ha="right",
            va="center",
            fontsize=10,
            fontweight="bold" if hit else "normal",
        )
        if hit:
            delta_us = (hit.time_ns - event.start_ns) / 1e3
            ax.text(
                SCINTILLATOR_WIDTH_CM / 2 + 0.35,
                z_cm + SCINTILLATOR_THICKNESS_CM / 2,
                f"{hit.sipm_mv:.1f} mV   Δt={delta_us:.1f} µs",
                ha="left",
                va="center",
                fontsize=9,
                color="#71162f",
            )

    total_height = (len(DEVICE_NAMES) - 1) * pitch_cm + SCINTILLATOR_THICKNESS_CM
    # Il margine destro ospita ampiezza e ritardo senza sovrapporli alla torre.
    ax.set_xlim(-6.5, 12.5)
    ax.set_ylim(-1.2, total_height + 2.5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [cm] (posizione del passaggio non misurata)")
    ax.set_ylabel("z [cm]")
    ax.set_xticks([-3, 0, 3])
    ax.grid(axis="y", color="#c8cdd2", linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)
    ax.set_title(
        f"Evento {event_index + 1}/{event_count}  |  "
        f"molteplicita' {event.multiplicity}  |  span {event.span_ms:.3f} ms\n"
        f"Device: {', '.join(hit.device for hit in event.hits)}",
        pad=12,
    )
    ax.text(
        0,
        total_height + 0.75,
        "fascia gialla = regione compatibile; angolo non ricostruibile",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#76520d",
    )


def show_events(
    coincidences: list[Coincidence],
    initial_event: int,
    gap_cm: float,
    save_path: Path | None,
    show: bool,
) -> None:
    if not coincidences:
        raise ValueError("nessuna coincidenza soddisfa i criteri richiesti")
    if not 0 <= initial_event < len(coincidences):
        raise ValueError(
            f"--event deve essere compreso tra 1 e {len(coincidences)}"
        )
    if gap_cm < 0:
        raise ValueError("--gap-cm non puo' essere negativo")

    all_amplitudes = np.asarray(
        [hit.sipm_mv for event in coincidences for hit in event.hits], dtype=float
    )
    vmax = max(float(np.percentile(all_amplitudes, 98)), 1.0)
    normalizer = colors.Normalize(vmin=0.0, vmax=vmax, clip=True)
    color_map = plt.get_cmap("plasma")

    fig, ax = plt.subplots(figsize=(8.2, 9.2))
    fig.subplots_adjust(left=0.16, right=0.84, top=0.90, bottom=0.10)
    colorbar = fig.colorbar(
        ScalarMappable(norm=normalizer, cmap=color_map), ax=ax,
        fraction=0.045, pad=0.04,
    )
    colorbar.set_label("Ampiezza SiPM [mV]")
    state = {"index": initial_event}

    def redraw() -> None:
        draw_event(
            ax,
            coincidences[state["index"]],
            state["index"],
            len(coincidences),
            gap_cm,
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
        "--gap-cm", type=float, default=DEFAULT_GAP_CM,
        help="spazio libero fra due scintillatori in cm (default: 5)",
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
            gap_cm=args.gap_cm,
            save_path=args.save,
            show=not args.no_show,
        )
    except (OSError, ValueError) as error:
        print(f"Errore: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
