#!/usr/bin/env python3
"""Algoritmi condivisi per coincidenze e gruppi temporali CosmicWatch.

Tutte le finestre espresse in millisecondi rappresentano lo span massimo del
gruppo: ``max(t) - min(t) <= window``. Gli indici restituiti si riferiscono
sempre agli array di input originali.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class TemporalGroup:
    """Un gruppo di hit appartenenti a device differenti."""

    indices: tuple[int, ...]
    center_s: float
    span_s: float
    score: float


def _validated_inputs(
    timestamps_s: Sequence[float], detector_names: Sequence[str], window_ms: float
) -> tuple[np.ndarray, np.ndarray, float]:
    timestamps = np.asarray(timestamps_s, dtype=float)
    names = np.asarray(detector_names, dtype=str)
    if timestamps.ndim != 1 or names.ndim != 1:
        raise ValueError("Timestamps and detector names must be one-dimensional")
    if len(timestamps) != len(names):
        raise ValueError("Timestamps and detector names must have the same length")
    if window_ms <= 0:
        raise ValueError("The coincidence window must be greater than zero")
    if len(timestamps) and not np.all(np.isfinite(timestamps)):
        raise ValueError("Timestamps must be finite")
    return timestamps, names, window_ms / 1000.0


def find_compact_groups(
    timestamps_s: Sequence[float],
    detector_names: Sequence[str],
    window_ms: float = 1.0,
    min_devices: int = 2,
    required_devices: Sequence[str] | None = None,
) -> list[TemporalGroup]:
    """Trova gruppi disgiunti premiando molteplicità e compattezza.

    Sono generati tutti i candidati contigui con span entro la finestra. Se un
    intervallo contiene più hit dello stesso device, viene conservato quello
    più vicino alla mediana. I candidati sono ordinati con un punteggio che
    favorisce gruppi numerosi e temporalmente compatti; ogni hit può essere
    assegnato a un solo gruppo.
    """
    timestamps, names, window_s = _validated_inputs(
        timestamps_s, detector_names, window_ms
    )
    if min_devices < 2:
        raise ValueError("min_devices must be at least two")
    if not len(timestamps):
        return []

    required = set(required_devices or ())
    if required and min_devices > len(required):
        raise ValueError("min_devices cannot exceed the number of required devices")

    order = np.argsort(timestamps, kind="mergesort")
    sorted_times = timestamps[order]
    sorted_names = names[order]
    relative_times = sorted_times - sorted_times[0]
    timing_scale_s = window_s / 4.0
    candidate_by_indices: dict[
        frozenset[int], tuple[float, TemporalGroup]
    ] = {}

    for first in range(len(order)):
        last = int(np.searchsorted(
            relative_times,
            relative_times[first] + window_s,
            side="right",
        ))
        for stop in range(first + 2, last + 1):
            positions = range(first, stop)
            interval_median = float(np.median(relative_times[first:stop]))
            best_by_device: dict[str, int] = {}
            for position in positions:
                device = str(sorted_names[position])
                previous = best_by_device.get(device)
                if (
                    previous is None
                    or abs(relative_times[position] - interval_median)
                    < abs(relative_times[previous] - interval_median)
                ):
                    best_by_device[device] = position

            if required:
                if not required.issubset(best_by_device):
                    continue
                selected_positions = [best_by_device[name] for name in required]
            else:
                if len(best_by_device) < min_devices:
                    continue
                selected_positions = list(best_by_device.values())

            selected_positions.sort(key=lambda position: relative_times[position])
            selected_times = relative_times[selected_positions]
            span_s = float(selected_times[-1] - selected_times[0])
            if span_s > window_s:
                continue
            center_s = float(np.median(selected_times))
            normalized_mean_square = float(np.mean(
                ((selected_times - center_s) / timing_scale_s) ** 2
            ))
            score = len(selected_positions) / (1.0 + normalized_mean_square)
            original_indices = tuple(
                int(order[position]) for position in selected_positions
            )
            key = frozenset(original_indices)
            group = TemporalGroup(
                indices=original_indices,
                center_s=center_s + float(sorted_times[0]),
                span_s=span_s,
                score=score,
            )
            previous_candidate = candidate_by_indices.get(key)
            if previous_candidate is None or score > previous_candidate[0]:
                candidate_by_indices[key] = (score, group)

    ranked = sorted(
        candidate_by_indices.values(),
        key=lambda item: (-item[0], -len(item[1].indices), item[1].span_s),
    )
    used_indices: set[int] = set()
    groups: list[TemporalGroup] = []
    for _, group in ranked:
        if any(index in used_indices for index in group.indices):
            continue
        groups.append(group)
        used_indices.update(group.indices)
    groups.sort(key=lambda group: group.center_s)
    return groups


def time_coincidence_analysis(
    timestamps_s: Sequence[float],
    detector_names: Sequence[str],
    window_ms: float = 10.0,
) -> tuple[np.ndarray, int]:
    """Restituisce la mask e tutte le associazioni pairwise nella finestra."""
    timestamps, names, window_s = _validated_inputs(
        timestamps_s, detector_names, window_ms
    )
    if len(set(names)) < 2:
        raise ValueError(
            "Time-based coincidences require at least two different detector names"
        )

    coincident = np.zeros(len(timestamps), dtype=bool)
    pair_count = 0
    order = np.argsort(timestamps, kind="mergesort")
    sorted_times = timestamps[order]
    sorted_names = names[order]
    left = 0
    for right in range(len(order)):
        while sorted_times[right] - sorted_times[left] > window_s:
            left += 1
        for candidate in range(left, right):
            if sorted_names[candidate] != sorted_names[right]:
                coincident[order[candidate]] = True
                coincident[order[right]] = True
                pair_count += 1
    return coincident, pair_count


def coincidences_from_computer_time(
    timestamps_s: Sequence[float],
    detector_names: Sequence[str],
    window_ms: float = 10.0,
) -> np.ndarray:
    return time_coincidence_analysis(
        timestamps_s, detector_names, window_ms
    )[0]


def unique_time_coincidence_pairs(
    timestamps_s: Sequence[float],
    detector_names: Sequence[str],
    window_ms: float = 10.0,
) -> list[tuple[int, int, float]]:
    """Associa gli hit più vicini uno-a-uno per ogni coppia di device."""
    timestamps, names, window_s = _validated_inputs(
        timestamps_s, detector_names, window_ms
    )
    unique_pairs: list[tuple[int, int, float]] = []
    unique_names = sorted(set(names))

    for name_index, name_a in enumerate(unique_names):
        indices_a = np.flatnonzero(names == name_a)
        indices_a = indices_a[np.argsort(timestamps[indices_a], kind="mergesort")]
        for name_b in unique_names[name_index + 1:]:
            indices_b = np.flatnonzero(names == name_b)
            indices_b = indices_b[
                np.argsort(timestamps[indices_b], kind="mergesort")
            ]
            times_b = timestamps[indices_b]
            candidates: list[tuple[float, int, int]] = []
            for index_a in indices_a:
                time_a = timestamps[index_a]
                first = int(np.searchsorted(times_b, time_a - window_s, side="left"))
                last = int(np.searchsorted(times_b, time_a + window_s, side="right"))
                for position_b in range(first, last):
                    index_b = int(indices_b[position_b])
                    candidates.append((
                        abs(timestamps[index_b] - time_a), int(index_a), index_b
                    ))

            used_a: set[int] = set()
            used_b: set[int] = set()
            for delta_t, index_a, index_b in sorted(candidates):
                if index_a in used_a or index_b in used_b:
                    continue
                used_a.add(index_a)
                used_b.add(index_b)
                unique_pairs.append((index_a, index_b, delta_t))
    return unique_pairs


def strict_all_detector_coincidence_groups(
    timestamps_s: Sequence[float],
    detector_names: Sequence[str],
    window_ms: float = 10.0,
) -> list[tuple[tuple[int, ...], float]]:
    """Trova gruppi compatti e disgiunti con un hit da ogni device."""
    timestamps, names, _ = _validated_inputs(
        timestamps_s, detector_names, window_ms
    )
    required_devices = tuple(sorted(set(names)))
    if len(required_devices) < 2:
        raise ValueError(
            "Strict time coincidences require at least two different detector names"
        )
    groups = find_compact_groups(
        timestamps,
        names,
        window_ms=window_ms,
        min_devices=len(required_devices),
        required_devices=required_devices,
    )
    return [(group.indices, group.span_s) for group in groups]


def build_reference_groups(
    timestamps_s: Sequence[float],
    detector_names: Sequence[str],
    reference_devices: Sequence[str],
    window_ms: float,
) -> list[TemporalGroup]:
    """Costruisce combinazioni compatte usando soltanto i riferimenti indicati."""
    timestamps, names, window_s = _validated_inputs(
        timestamps_s, detector_names, window_ms
    )
    references = tuple(str(device) for device in reference_devices)
    if len(references) < 2 or len(set(references)) != len(references):
        raise ValueError("Reference devices must be at least two and all different")

    indices_by_device = {
        device: np.flatnonzero(names == device)[
            np.argsort(timestamps[names == device], kind="mergesort")
        ]
        for device in references
    }
    times_by_device = {
        device: timestamps[indices]
        for device, indices in indices_by_device.items()
    }
    anchor_device = min(
        references, key=lambda device: len(indices_by_device[device])
    )
    other_devices = tuple(device for device in references if device != anchor_device)
    candidates: dict[frozenset[int], tuple[float, float, TemporalGroup]] = {}

    for anchor_index in indices_by_device[anchor_device]:
        anchor_time = timestamps[anchor_index]
        nearby_by_device: list[list[int]] = []
        for device in other_devices:
            device_indices = indices_by_device[device]
            device_times = times_by_device[device]
            begin = int(np.searchsorted(
                device_times, anchor_time - window_s, side="left"
            ))
            end = int(np.searchsorted(
                device_times, anchor_time + window_s, side="right"
            ))
            if begin == end:
                break
            nearby_by_device.append([
                int(index) for index in device_indices[begin:end]
            ])
        else:
            for combination in product(*nearby_by_device):
                group_indices = (int(anchor_index),) + combination
                group_times = timestamps[list(group_indices)]
                span_s = float(np.max(group_times) - np.min(group_times))
                if span_s > window_s:
                    continue
                center_s = float(np.median(group_times))
                dispersion = float(np.mean((group_times - center_s) ** 2))
                ordered_indices = tuple(sorted(
                    group_indices, key=lambda index: timestamps[index]
                ))
                group = TemporalGroup(
                    indices=ordered_indices,
                    center_s=center_s,
                    span_s=span_s,
                    score=1.0 / (1.0 + dispersion / (window_s / 4.0) ** 2),
                )
                candidates[frozenset(group_indices)] = (
                    dispersion, span_s, group
                )

    used_indices: set[int] = set()
    groups: list[TemporalGroup] = []
    for _, _, group in sorted(candidates.values(), key=lambda item: item[:2]):
        if any(index in used_indices for index in group.indices):
            continue
        groups.append(group)
        used_indices.update(group.indices)
    groups.sort(key=lambda group: group.center_s)
    return groups


def match_probe_to_groups(
    groups: Sequence[TemporalGroup],
    timestamps_s: Sequence[float],
    detector_names: Sequence[str],
    probe_device: str,
    window_ms: float,
) -> list[tuple[int, int]]:
    """Associa uno-a-uno il probe più vicino alla mediana dei riferimenti."""
    timestamps, names, window_s = _validated_inputs(
        timestamps_s, detector_names, window_ms
    )
    probe_indices = np.flatnonzero(names == probe_device)
    probe_indices = probe_indices[
        np.argsort(timestamps[probe_indices], kind="mergesort")
    ]
    probe_times = timestamps[probe_indices]
    candidates: list[tuple[float, int, int]] = []

    for group_index, group in enumerate(groups):
        group_times = timestamps[list(group.indices)]
        earliest_probe = float(np.max(group_times) - window_s)
        latest_probe = float(np.min(group_times) + window_s)
        begin = int(np.searchsorted(probe_times, earliest_probe, side="left"))
        end = int(np.searchsorted(probe_times, latest_probe, side="right"))
        for position in range(begin, end):
            candidates.append((
                abs(float(probe_times[position]) - group.center_s),
                group_index,
                int(probe_indices[position]),
            ))

    used_groups: set[int] = set()
    used_probes: set[int] = set()
    matches: list[tuple[int, int]] = []
    for _, group_index, probe_index in sorted(candidates):
        if group_index in used_groups or probe_index in used_probes:
            continue
        used_groups.add(group_index)
        used_probes.add(probe_index)
        matches.append((group_index, probe_index))
    return matches
