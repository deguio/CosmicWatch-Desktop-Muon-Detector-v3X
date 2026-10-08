#!/usr/bin/env bash
# Host tests of the hardware-independent code and of the data format against
# the parsers used by Data/import_data.py and GUI/GUI.py.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT

cc -std=c11 -Wall -Wextra -Werror -O1 -fsanitize=address,undefined -I"$HERE/../src" \
    "$HERE/test_cwcore.c" "$HERE/../src/cwcore.c" -lm -o "$OUT/test_cwcore"
"$OUT/test_cwcore"

cc -std=c11 -Wall -O1 -I"$HERE/../src" "$HERE/gen_lines.c" "$HERE/../src/cwcore.c" -lm -o "$OUT/gen_lines"
"$OUT/gen_lines" > "$OUT/lines.txt"
PY="${PYTHON:-$HERE/../../../../venv_cosmicwatch.nosync/bin/python}"
"$PY" "$HERE/check_format.py" "$OUT/lines.txt"
