"""Feed firmware USB lines to the real parser of Data/import_data.py.

import_data.py runs its acquisition at import time, so only its pure helper
functions are extracted (via ast) and executed here.
"""
import ast
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[3]
src = (ROOT / 'Data' / 'import_data.py').read_text()
wanted = {'current_time_ns', 'format_pc_timestamp', 'build_output_row', 'output_header'}
tree = ast.parse(src)
module = ast.Module(
    body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted],
    type_ignores=[])
ns = {}
exec('import time\nfrom datetime import datetime\n', ns)
exec(compile(module, 'import_data_helpers', 'exec'), ns)

lines = pathlib.Path(sys.argv[1]).read_bytes().split(b'\n')
lines = [l + b'\n' for l in lines if l]
expected_schemas = [(), ('temperature', 'pressure'), ('accel', 'gyro'),
                    ('temperature', 'pressure', 'accel', 'gyro')]
ok = 0
for i, raw in enumerate(lines):
    parsed = ns['build_output_row'](raw, 1_700_000_000_123_456_789)
    assert parsed is not None, raw
    row, schema, name = parsed
    assert name == 'TOP', name
    assert schema == expected_schemas[i // 3], (schema, raw)
    assert len(row) == 6 + len(schema) + 3
    float(row[1]); int(row[3]); float(row[4]); float(row[5])
    # GUI live mode: the first six fields must parse and the flag is '0'/'1'
    fields = raw.decode().rstrip('\r\n\t').split('\t')
    assert fields[2] in ('0', '1')
    ok += 1
hk = (b'# HK t=60.001 T=23.41 P=100512 HV=30.120 Thr_set=80.00 Thr_meas=80.05 base=55.30 '
      b'base_sd=1.20 base_n=60 corr=1.00000 events=171 coinc=24 dead=0.005321 lost=0/0/0 '
      b'bmp_fail=0 bmp_rej=0 src=periodic\r\n')
assert ns['build_output_row'](hk, 0) is None, 'HK line must be ignored'
print(f'import_data.py parser accepted {ok}/{len(lines)} lines; header: '
      + ns['output_header'](expected_schemas[3]))
