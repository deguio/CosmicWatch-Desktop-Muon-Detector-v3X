#!/usr/bin/env python3

"""Copy existing files from the oscilloscope, optionally matching a wildcard."""

import argparse
import csv
import os
from contextlib import contextmanager
from fnmatch import fnmatchcase
from pathlib import Path, PureWindowsPath
import subprocess
import tempfile
from time import monotonic

from RsInstrument import RsInstrument


@contextmanager
def instrument_connection(args):
    """Yield the resource/backend; own and close only the SSH tunnel created here."""
    if not args.ssh_bastion:
        resource = args.resource or f"TCPIP::{args.instrument_host}::INSTR"
        backend = 'socketio' if resource.upper().endswith('::SOCKET') else 'pyvisa-py'
        yield resource, backend
        return
    # OpenSSH handles keys, ssh-agent, passwords/MFA and ~/.ssh/config itself.
    # A private control socket allows cleanup without touching other SSH sessions.
    with tempfile.TemporaryDirectory(prefix='cwssh-', dir='/tmp') as directory:
        control = str(Path(directory) / 'control')
        destination = args.ssh_bastion
        forward = f'127.0.0.1:{args.local_port}:{args.instrument_host}:{args.instrument_port}'
        command = [
            'ssh', '-M', '-S', control, '-fNT',
            '-o', 'ControlPersist=no', '-o', 'ExitOnForwardFailure=yes',
            '-o', 'ConnectTimeout=15', '-o', 'ServerAliveInterval=30',
            '-o', 'ServerAliveCountMax=3', '-L', forward, '--', destination,
        ]
        print(f'Opening SSH tunnel via {destination}: {forward}', flush=True)
        try:
            subprocess.run(command, check=True)
            yield f'TCPIP::127.0.0.1::{args.local_port}::SOCKET', 'socketio'
        finally:
            # Also attempt cleanup if startup or the instrument session fails.
            subprocess.run(['ssh', '-S', control, '-O', 'exit', '--', destination],
                           check=False, capture_output=True, timeout=10)


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
    connection = parser.add_mutually_exclusive_group()
    connection.add_argument('--ssh-bastion', help='SSH destination: user@host or an alias in ~/.ssh/config')
    connection.add_argument('--resource', help='Explicit VISA resource, e.g. TCPIP::127.0.0.1::15025::SOCKET for a manual tunnel')
    parser.add_argument('--instrument-host', default='193.206.156.231', help='Oscilloscope address reachable from the bastion (or directly on campus)')
    parser.add_argument('--instrument-port', type=int, default=5025, help='Remote SCPI socket port for SSH mode (default: 5025)')
    parser.add_argument('--local-port', type=int, default=15025, help='Loopback port for the SSH tunnel (default: 15025)')
    parser.add_argument('--timeout-ms', type=int, default=30000, help='Timeout for each VISA I/O operation, in milliseconds (default: 30000)')
    args = parser.parse_args()
    if args.timeout_ms <= 0:
        parser.error('--timeout-ms must be positive.')
    if not all(1 <= port <= 65535 for port in (args.local_port, args.instrument_port)):
        parser.error('TCP ports must be between 1 and 65535.')
    if args.ssh_bastion and (args.ssh_bastion.startswith('-') or any(c.isspace() for c in args.ssh_bastion)):
        parser.error('--ssh-bastion must be user@host or an SSH alias, not a shell command.')

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
    with instrument_connection(args) as (resource, backend):
        copy_files(resource, backend, source, wildcard, output, parser, args.timeout_ms)


def copy_files(resource, backend, source, wildcard, output, parser, timeout_ms=30000):
    print(f'Connecting to {resource} ...', flush=True)
    rto = RsInstrument(resource, True, False,
                       options=f"SelectVisa={backend}, VisaTimeout={timeout_ms}")
    try:
        rto.visa_timeout = timeout_ms
        rto.instrument_status_checking = True
        rto.data_chunk_size = 64 * 1024
        print(f'Connected: {rto.idn_string}', flush=True)
        print(f'Checking files in {source.parent} ...', flush=True)
        # Validate exact names too: a missing MMEM:DATA? file may produce no block,
        # leaving the socket reader waiting for a response until its timeout.
        sources = matching_files(rto, source)
        if not sources:
            parser.error(f"No files match: {source}. Include the full name and extension, "
                         "or quote a prefix followed by *, e.g. -f 'RefCurve_2026-09-21_47*'.")
        transfers = [
            (item, output / item.name if output.is_dir() else output)
            for item in sources
        ]
        copied = skipped = 0
        print(f"Found {len(transfers)} matching file(s)", flush=True)
        for item, destination in transfers:
            if destination.exists() or destination.is_symlink():
                print(f"Skipping existing destination: {destination}", flush=True)
                skipped += 1
                continue
            print(f"Copying {item} to {destination}", flush=True)
            last_update = 0.
            def progress(event):
                nonlocal last_update
                now = monotonic()
                if event.end_of_transfer or now-last_update >= 1:
                    total = event.total_size if event.total_size is not None else '?'
                    print(f'  Received {event.transferred_size}/{total} bytes', flush=True)
                    last_update = now
            rto.events.on_read_handler = progress
            try:
                # Do not leave a zero-length or incomplete file under the final name.
                with tempfile.TemporaryDirectory(prefix='.copy-rto-', dir=destination.parent) as directory:
                    temporary = Path(directory) / 'download'
                    rto.read_file_from_instrument_to_pc(str(item), str(temporary))
                    # Publish without overwriting a file that appeared during transfer.
                    try:
                        os.link(temporary, destination)
                    except FileExistsError:
                        print(f"Skipping destination created during transfer: {destination}", flush=True)
                        skipped += 1
                        continue
            except BaseException:
                print(f'Transfer interrupted or failed: {item}. No final file created by this attempt.', flush=True)
                raise
            finally:
                rto.events.on_read_handler = None
            print(f"File saved to {destination.resolve()} ({destination.stat().st_size} bytes)", flush=True)
            copied += 1
        print(f"Done: {copied} copied, {skipped} skipped (already present).", flush=True)
    finally:
        rto.close()


if __name__ == "__main__":
    main()
