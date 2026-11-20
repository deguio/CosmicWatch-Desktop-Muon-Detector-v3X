# Dati di prova

Estratti piccoli (~11 MB in tutto) delle prese dati vere, per verificare che
gli script funzionino. Non servono per i risultati: le prese complete restano
fuori da git.

| file o cartella | da dove viene | formato | script |
|---|---|---|---|
| `seven_parallel_10min.txt` | primi 10 min di `showers_analysis/CW_seven_parallel_2.txt` (7 rivelatori affiancati, 19/11/2026) | USB, `import_data.py` | `shower_analysis.py`, `plot.py -c time` |
| `tower_nine_10min.txt` | primi 10 min di `coincidence_nine.txt` (torre di 9 rivelatori, 11/09/2026) | USB, `import_data.py` | `display.py` |
| `METRO_C_014.txt` | coppia in coincidenza hardware in metropolitana (rivelatore QUATTRO, 1.4 h) | microSD | `plot.py` |
| `rto6_calibrazione_4cm/` | 10 eventi da `RTO6/0ps` (rivelatori a contatto) | export RTO6 | `analyze_pulses.py` |
| `rto6_tof_187cm/` | 10 eventi da `RTO6/16ps` (rivelatori a 1.87 m) | export RTO6 | `analyze_pulses.py` |
| `keysight_tof_3m/` | primi 10 eventi di `Keysight-1000X/run_3m`, con `events.csv` | CSV Keysight | `analyze_pulses.py` |

Per provare tutto in un colpo (venv attivo, dalla radice del repository):

```bash
bash analysis_code/run_smoke_tests.sh
```

Lo script copia questi file in una cartella temporanea, lancia ogni programma
e i test unitari, e stampa `OK` o `FAIL` per ciascun passo con il percorso del
log. Con 10 minuti di dati i numeri non hanno valore statistico: si controlla
solo che la catena giri dall'inizio alla fine.
