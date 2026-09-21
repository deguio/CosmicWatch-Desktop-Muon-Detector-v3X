#!/usr/bin/env python3

"""Copy existing files from the oscilloscope, optionally matching a wildcard."""

import argparse
import csv
from fnmatch import fnmatchcase
from pathlib import Path, PureWindowsPath

from RsInstrument import RsInstrument


def matching_files(rto, source):
    # The catalog has two size fields followed by quoted "name,type,size" entries.
    directory = str(source.parent).replace("'", "''")
    response = rto.query_str(f"MMEMory:CATalog? '{directory}'")
    fields = next(csv.reader([response.strip()], skipinitialspace=True, strict=True))
    if len(fields) < 2:
        raise ValueError(f"Unexpected catalog response: {response!r}")
    matches = []
    for entry in fields[2:]:
        if not entry.strip():
            continue
        parts = entry.rsplit(',', 2)
        if len(parts) != 3:
            raise ValueError(f"Unexpected catalog entry: {entry!r}")
        name, kind, size = parts
        if kind.strip().upper() == 'DIR':
            continue
        # Accept only file names, never paths outside the requested directory.
        if not name or name in ('.', '..') or any(c in name for c in '/\\:'):
            raise ValueError(f"Unexpected file name in catalog: {name!r}")
        if fnmatchcase(name.casefold(), source.name.casefold()):
            matches.append(source.parent / name)
    return sorted(set(matches), key=lambda path: path.name.casefold())


def main():
    parser = argparse.ArgumentParser(
        description="Copy existing files from the RTO6 disk (no acquisition)."
    )
    parser.add_argument(
        "-f", "--fileName", required=True,
        help="File name/pattern, or full Windows path on the oscilloscope. "
             "Bare names use RefWaveforms. Quote wildcards, e.g. 'RefCurve_2026-09-21*'.",
    )
    parser.add_argument(
        "-o", "--output", type=Path,
        help="Local destination: an existing directory for wildcards; a file path or "
             "existing directory for a single file. Default: current directory.",
    )
    args = parser.parse_args()

    source = PureWindowsPath(args.fileName)
    if not source.is_absolute():
        if source.drive or source.root or len(source.parts) != 1:
            parser.error("Use a bare file name or an absolute Windows path for --fileName.")
        source = PureWindowsPath(
            r"C:\Users\Public\Documents\Rohde-Schwarz\RTx\RefWaveforms"
        ) / source

    if any(c in str(source.parent) for c in '*?['):
        parser.error("Wildcards are supported only in the file name, not in directories.")
    wildcard = any(c in source.name for c in '*?[')
    output = (args.output or Path.cwd()).expanduser()
    if wildcard and not output.is_dir():
        parser.error(f"With a wildcard, --output must be an existing directory: {output}")
    if not output.is_dir() and not output.parent.is_dir():
        parser.error(f"Destination directory does not exist: {output.parent}")

    RsInstrument.assert_minimum_version("1.53.0")
    rto = RsInstrument(
        "TCPIP::193.206.156.231::INSTR", True, False,
        options="SelectVisa=pyvisa-py",
    )
    try:
        rto.visa_timeout = 3000000
        rto.instrument_status_checking = True
        sources = matching_files(rto, source) if wildcard else [source]
        if not sources:
            parser.error(f"No files match: {source}")
        transfers = [
            (item, output / item.name if output.is_dir() else output)
            for item in sources
        ]
        for _, destination in transfers:
            if destination.exists() or destination.is_symlink():
                parser.error(f"Destination already exists: {destination}")
        print(f"Copying {len(transfers)} file(s)")
        for item, destination in transfers:
            print(f"Copying {item} to {destination}")
            rto.read_file_from_instrument_to_pc(str(item), str(destination))
            print(f"File saved to {destination.resolve()}")
    finally:
        rto.close()


if __name__ == "__main__":
    main()
