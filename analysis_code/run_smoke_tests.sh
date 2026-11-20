#!/usr/bin/env bash
# Prova rapida che tutti gli script girino sui dati di Data/TestData.
# Uso (venv attivo, da qualunque cartella):  bash analysis_code/run_smoke_tests.sh
# I risultati vanno in una cartella temporanea, che viene stampata e non cancellata.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
TEST="$ROOT/Data/TestData"
OUT="$(mktemp -d)"
export MPLBACKEND=Agg
PY="${PYTHON:-python}"
FAILED=0

step() {
    local name="$1"; shift
    if "$@" > "$OUT/$name.log" 2>&1; then
        echo "  OK    $name"
    else
        echo "  FAIL  $name  (log: $OUT/$name.log)"
        FAILED=$((FAILED + 1))
    fi
}

echo "Risultati in $OUT"
cp -R "$TEST"/. "$OUT/data"

echo "Script di analisi:"
step plot "$PY" -B "$HERE/plot.py" -i "$OUT/data/METRO_C_014.txt"
step plot_time "$PY" -B "$HERE/plot.py" -i "$OUT/data/seven_parallel_10min.txt" -c time -w 1
step display "$PY" -B "$HERE/display.py" "$OUT/data/tower_nine_10min.txt" --no-show --save "$OUT/display_evento1.png"
step shower_analysis "$PY" -B "$HERE/shower_analysis.py" "$OUT/data/seven_parallel_10min.txt" \
    --positions "$ROOT/Data/showers_analysis/positions_IMG_6986.txt" \
    --n-shifts 5 --mc-samples 200000 -o "$OUT/shower"
step pulses_rto6_calibrazione "$PY" -B "$HERE/analyze_pulses.py" "$OUT/data/rto6_calibrazione_4cm" \
    --offset-ns 0 --distance-m 0.04 --fraction 0.3 --output-dir "$OUT/rto6_calibrazione"
step pulses_rto6_tof "$PY" -B "$HERE/analyze_pulses.py" "$OUT/data/rto6_tof_187cm" \
    --offset-ns 15.72 --distance-m 1.87 --fraction 0.3 --summary-method causal --output-dir "$OUT/rto6_tof"
step pulses_keysight "$PY" -B "$HERE/analyze_pulses.py" "$OUT/data/keysight_tof_3m" \
    --offset-ns 0 --distance-m 3 --fraction 0.3 --output-dir "$OUT/keysight_tof"
step synthetic "$PY" -B "$ROOT/Keysight-1000X/make_synthetic.py" --output-dir "$OUT/synthetic" -n 20

echo "Test unitari:"
step unittest_analyze_pulses bash -c "cd '$HERE' && '$PY' -B -m unittest -q test_analyze_pulses"
step unittest_rto6 bash -c "cd '$ROOT/RTO6' && '$PY' -B -m unittest -q"
step unittest_keysight bash -c "cd '$ROOT/Keysight-1000X' && '$PY' -B -m unittest -q"

if [ "$FAILED" -eq 0 ]; then
    echo "Tutto OK."
else
    echo "$FAILED passi falliti."
    exit 1
fi
