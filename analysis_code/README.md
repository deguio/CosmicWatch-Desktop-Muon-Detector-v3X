# Codice di acquisizione e analisi CosmicWatch

Questa cartella raccoglie gli script Python usati nel progetto con il liceo
(autunno 2026) per acquisire e analizzare i dati dei CosmicWatch v3X e delle
misure di tempo di volo all'oscilloscopio. I dati stanno in `../Data/` (non in
git, tranne `ExampleData/` e la documentazione dei risultati); le misure con gli
oscilloscopi sono descritte in [`../RTO6/README.md`](../RTO6/README.md) e
[`../Keysight-1000X/README.md`](../Keysight-1000X/README.md).

| script | a cosa serve |
|---|---|
| [`import_data.py`](import_data.py) | acquisizione via USB da uno o più CosmicWatch, con l'ora del PC su ogni evento |
| [`plot.py`](plot.py) | analisi standard di un file: spettri, rate, coincidenze, sensori |
| [`temporal_groups.py`](temporal_groups.py) | funzioni per le coincidenze fra rivelatori (usate da `plot.py` e `display.py`) |
| [`display.py`](display.py) | event display di una torre di rivelatori impilati ed efficienza dei piani |
| [`shower_analysis.py`](shower_analysis.py) | ricerca di sciami con rivelatori affiancati: allineamento dei clock, accidentali, modello |
| [`analyze_pulses.py`](analyze_pulses.py) | analisi delle forme d'onda all'oscilloscopio (RTO6 e Keysight) e tempo di volo |
| [`test_analyze_pulses.py`](test_analyze_pulses.py) | test di `analyze_pulses.py` |
| [`run_smoke_tests.sh`](run_smoke_tests.sh) | prova rapida di tutti gli script sui dati di [`../Data/TestData/`](../Data/TestData/) |

---

## 1. Prima di cominciare

**Ambiente.** Tutto gira nel venv del progetto (`venv_cosmicwatch.nosync`, una
cartella sopra la radice del repository):

```bash
source ../../venv_cosmicwatch.nosync/bin/activate   # dalla radice del repository
```

I pacchetti necessari sono elencati in [`../requirements.txt`](../requirements.txt),
divisi per uso: analisi, acquisizione USB, Keysight, RTO6, notebook, GUI.
Per un venv nuovo:

```bash
python -m pip install -r requirements.txt            # dalla radice del repository
```

**Flusso dei dati.**

```
CosmicWatch ──USB──► import_data.py ──► file .txt (uno per presa dati, tutti i rivelatori)
     │                                        │
     └──microSD──► file NOME_C_NNN.txt        ├──► plot.py            (un file, analisi standard)
                   (un rivelatore)            ├──► display.py         (torre impilata)
                                              └──► shower_analysis.py (rivelatori affiancati)

Oscilloscopio ──► CSV per evento ──► analyze_pulses.py (tempo di volo fra due rivelatori)
```

### 1.1 Formato dei file

Ogni riga è un evento, colonne separate da tab, commenti con `#`:

```
Event  Timestamp[s]  Coincident[bool]  ADC[12b]  SiPM[mV]  Deadtime[s]  Temp[C]  Pressure[Pa]  Accel(X:Y:Z)[g]  Gyro(X:Y:Z)[deg/sec]  [Name  Time  Date]
```

- **Timestamp[s]**: clock interno del rivelatore, risoluzione 1 µs, da zero
  all'accensione. Ogni rivelatore ha il suo quarzo (differenze di 1–3 ppm).
- **Coincident**: 1 se il rivelatore collegato col cavo di coincidenza ha
  visto un impulso entro ~3 µs. È una coincidenza **hardware**, fra due soli
  rivelatori.
- **ADC / SiPM[mV]**: ampiezza dell'impulso. Le scale del firmware originale e
  del firmware open non sono confrontabili in valore assoluto.
- **Deadtime[s]**: tempo morto cumulativo. Live time = Timestamp − Deadtime.
- **Name, Time, Date**: solo nei file scritti da `import_data.py`. Nome del
  rivelatore e **ora del PC** all'arrivo del primo byte della riga.

I file su microSD (`NOME_C_NNN.txt` o `NOME_M_NNN.txt`) non hanno le ultime
tre colonne: contengono un solo rivelatore e solo il clock interno.

### 1.2 Firmware originale e firmware open

Le prese dati di questo progetto usano due firmware
(vedi [`../Firmware/CosmicWatchOpen/`](../Firmware/CosmicWatchOpen/)):

| | originale (v3X.28) | open (`CosmicWatchOpen` 1.0.0) |
|---|---|---|
| tempo morto per evento | ~0.6–0.9 ms | ~13 µs |
| timestamp interno | con jitter di ~100 µs rispetto al trigger | entro ~0.3 µs dal trigger |
| latenza al PC, 99° percentile | 0.5 ms | 2.7 ms (scritture su SD) |

Il formato dei file è lo stesso. Tutte le prese descritte qui (sciami, torre,
profondità, pressione) sono state fatte col firmware originale. Si riconosce
dall'header: l'originale scrive `Welcome to the CosmicWatchDAQ (Version: v3X.28)`,
il firmware open `CosmicWatch v3X open firmware (Version: CosmicWatch-open 1.0.0)`.

### 1.3 Coincidenze e accidentali

Si possono cercare coincidenze in tre modi, con risoluzioni molto diverse:

| metodo | risoluzione | quando usarlo |
|---|---|---|
| flag `Coincident` (hardware) | ~3 µs | due rivelatori collegati col cavo |
| ora del PC (`Time`) | ~1 ms (jitter USB) | più rivelatori via USB, analisi rapida |
| clock interni allineati al PC | ~10 µs + jitter del firmware | più rivelatori via USB, analisi fine (`shower_analysis.py`) |

**Accidentali.** Due rivelatori indipendenti con rate R₁ e R₂ danno per caso
coincidenze entro ±W con rate 2·R₁·R₂·W. Per k rivelatori entro uno span W la
formula è k·W^(k−1)·∏R_i. Con rivelatori a ~7 Hz:

| finestra | doppie casuali per coppia | triple casuali per terna |
|---|---|---|
| 3 µs (hardware) | ~1/h | trascurabili |
| 1 ms (ora PC) | ~350/h | ~4/h |
| 10 ms (default di `plot.py -c time`) | ~3500/h | ~370/h |

Morale: con le coincidenze via ora del PC le doppie sono dominate dal caso;
servono una soglia in ampiezza, finestre strette, molteplicità ≥ 3 o
l'allineamento dei clock. In ogni analisi conviene stimare gli accidentali
anche col metodo **off-time** (si spostano i tempi di un rivelatore di
decine di ms e si ricontano le coincidenze).

---

## 2. Gli script

### 2.1 `import_data.py`: acquisizione

Legge in parallelo (un thread per porta) i rivelatori collegati via USB e
scrive un unico file con tutti gli eventi. A ogni riga aggiunge nome del
rivelatore e ora del PC, presa **all'arrivo del primo byte** della riga.

```bash
cd ../Data                       # il file viene scritto nella cartella corrente
python ../analysis_code/import_data.py        # -p stampa anche gli eventi a schermo
```

Lo script chiede quali porte leggere e il nome del file (default `CW_data.txt`).

- **Nomi dei rivelatori:** vengono dal `config.txt` della microSD (o dalla
  flash, col firmware open). Devono essere **tutti diversi**: lo script scarta
  un rivelatore con un nome già usato.
- **Colonne:** devono essere uguali per tutti (stessi sensori); le righe con
  colonne diverse vengono scartate con un avviso.
- **Prese lunghe:** su Mac `caffeinate -i python ...` impedisce lo standby.

Il test del firmware [`../Firmware/CosmicWatchOpen/tests/check_format.py`](../Firmware/CosmicWatchOpen/tests/check_format.py)
usa le funzioni di parsing di questo script: se le si modifica, rilanciare
`../Firmware/CosmicWatchOpen/tests/run_tests.sh`.

### 2.2 `plot.py`: analisi standard

```bash
python analysis_code/plot.py -i Data/ExampleData/AxLab_C_001.txt
python analysis_code/plot.py -i Data/coincidence_six.txt -c time -w 1           # coincidenze col tempo PC
python analysis_code/plot.py -i Data/coincidence_six.txt -c time -w 1 --time-coincidence-mode all
```

| opzione | significato |
|---|---|
| `-i` | file di ingresso |
| `-b` | larghezza dei bin del rate nel tempo, in s (default 60) |
| `-c device\|time` | coincidenze dal flag hardware (default) o dall'ora del PC |
| `-w` | finestra in ms per `-c time` (default 10: conviene molto meno, ~1 ms) |
| `--time-coincidence-mode pairwise\|all` | coppie indipendenti, oppure un evento da **ogni** rivelatore nella finestra |
| `--delta-t-range-ms` | semiampiezza dell'istogramma dei Δt fra rivelatori (default 20) |

Produce in `Figures/` accanto al file di ingresso, in PDF:

- spettri ADC e SiPM, separati per eventi in coincidenza e non, con fit
  Landau (convoluta con una gaussiana) per gli eventi in coincidenza;
- rate nel tempo e tempo morto;
- temperatura, pressione, accelerometro e giroscopio nel tempo;
- con `-c time`: istogramma dei Δt fra rivelatori, rate di coppie (o di
  gruppi di N) osservato, accidentale atteso e differenza.

Stampa a terminale rate totale, rate in coincidenza e non, con errori.

### 2.3 `temporal_groups.py`: funzioni per le coincidenze

Libreria, non si lancia da sola. Tutte le finestre sono **span**:
max(t) − min(t) ≤ W.

| funzione | cosa fa |
|---|---|
| `find_compact_groups` | gruppi disgiunti di hit di rivelatori diversi, preferendo i più numerosi e compatti |
| `time_coincidence_analysis`, `coincidences_from_computer_time` | flag di coincidenza ricostruito dall'ora del PC |
| `unique_time_coincidence_pairs` | coppie indipendenti (un hit non viene riusato) |
| `strict_all_detector_coincidence_groups` | gruppi con un hit da ogni rivelatore |
| `build_reference_groups`, `match_probe_to_groups` | telescopio di riferimento più rivelatore sotto test (efficienze in `display.py`) |

### 2.4 `display.py`: torre impilata

Per una torre di (fino a) nove rivelatori impilati, `UNO` in basso e `NOVE` in
alto. Ricostruisce le coincidenze con l'ora del PC e:

- **mostra gli eventi** uno per uno: piani accesi e fascio di rette compatibili
  (senza posizione nel piano non si può dare una sola traccia). Frecce o `n`/`p`
  per scorrere, `q` per uscire;
- **misura l'efficienza** di ogni piano interno: N(sotto + centro + sopra) /
  N(sotto + sopra), sottraendo gli accidentali attesi. Una seconda stima usa due
  piani sotto e due sopra, con fondo molto minore.

```bash
python analysis_code/display.py                      # default: Data/coincidence_nine.txt
python analysis_code/display.py Data/coincidence_sette.txt --window-ms 1 --min-devices 3
python analysis_code/display.py --event 12 --save evento_12.png --no-show
```

Opzioni: `-w` finestra in ms (default 1), `-m` molteplicità minima,
`--min-sipm-mv` soglia per i piani di riferimento, `--tower-height-cm`
(default 38.5), `-e` evento iniziale.

### 2.5 `shower_analysis.py`: sciami con rivelatori affiancati

Allinea i clock interni dei rivelatori al tempo del PC (precisione ~10 µs),
cerca gruppi di hit in rivelatori diversi entro W, stima gli accidentali
(analitici e off-time), misura lo spettro delle particelle di sciame e
confronta il numero di eventi per molteplicità con il modello dello spettro di
densità degli sciami, anche con le posizioni reali dei rivelatori.

```bash
python analysis_code/shower_analysis.py \
    Data/showers_analysis/CW_seven_parallel.txt Data/showers_analysis/CW_seven_parallel_2.txt \
    --positions Data/showers_analysis/positions_IMG_6986.txt
```

Tutto (modello, scelte, risultati, figure) è spiegato in
[`../Data/showers_analysis/ANALISI_SCIAMI.md`](../Data/showers_analysis/ANALISI_SCIAMI.md).
L'allineamento dei clock è riutilizzabile in qualunque analisi multi-rivelatore
via USB (funzione `align_clocks`).

### 2.6 `analyze_pulses.py`: forme d'onda all'oscilloscopio

Legge gli export a due canali dell'RTO6 (`.csv` + `.Wfm.csv`) e del Keysight
InfiniiVision (CSV con colonna dei tempi). Per ogni evento:

- trova il tempo di arrivo con un **constant fraction** digitale (default 20%)
  e con due fit della forma dell'impulso (differenza di esponenziali, con e
  senza convoluzione gaussiana);
- calcola ΔT = t(C2) − t(C1), sottrae l'offset dei cavi (`--offset-ns`) e,
  con `--distance-m`, la velocità v = L/ΔT.

Su più eventi produce riepiloghi di ΔT, v/c e rate nel tempo, con tagli in
ampiezza (`--min-amplitude-mv`, `--max-amplitude-mv`). Modello, opzioni e
verifiche sono documentati in [`../RTO6/README_timing.md`](../RTO6/README_timing.md);
l'acquisizione col Keysight in [`../Keysight-1000X/README_acquisition.md`](../Keysight-1000X/README_acquisition.md).

```bash
cd RTO6
python ../analysis_code/analyze_pulses.py '16ps/RefCurve_2026-09-2*' --offset-ns 15.72 --distance-m 1.87 \
    --min-amplitude-mv 25 --max-amplitude-mv 170 --fraction 0.3 --summary-method causal
```

### 2.7 Notebook degli studenti

[`../Data/pressure_analysis/`](../Data/pressure_analysis/) contiene i notebook
`13novembre.ipynb` e `13novembre_locale.ipynb`, scritti con gli studenti e poi
rivisti. Spiegano passo passo come caricare un file CosmicWatch con pandas,
selezionare le coincidenze, fare spettri e rate nel tempo, e misurare la
**correlazione fra rate dei muoni e pressione atmosferica**. I dati vanno nella
sottocartella `esperienzeTracceCosmo/` (non in git), indicata da
`percorso_cartella` nella seconda cella.

---

## 3. Risultati ottenuti con questi script

| misura | dati | risultato | dettagli |
|---|---|---|---|
| sciami con 7 rivelatori affiancati | 26.3 h, firmware originale | 38 eventi con ≥3 rivelatori contro 0.35 casuali; rate ~3× il modello; doppie con componente locale | [`ANALISI_SCIAMI.md`](../Data/showers_analysis/ANALISI_SCIAMI.md) |
| coefficiente barometrico | coppia in coincidenza, ~14 giorni | β = −0.178 ± 0.040 %/hPa (documentazione CosmicWatch: −0.141 ± 0.007) | notebook `13novembre.ipynb` |
| superficie contro metropolitana | coppia in coincidenza hardware | 13.65 ± 0.11 /min contro 2.73 ± 0.18 /min: rapporto 0.200 ± 0.013 | sez. 3.1 |
| efficienza dei piani della torre | 9 rivelatori impilati, 1.2 h | 55–75% (sandwich a 3), 40–57% (a 5), con grandi correzioni per gli accidentali | sez. 3.2 |
| tempo di volo con RTO6 | 1.87 m, 78 eventi | ΔT = 5.99 ± 0.10 ns contro 6.24 ns attesi a c | [`../RTO6/README.md`](../RTO6/README.md) |
| tempo di volo con Keysight | 3 m, 147 eventi | ΔT = 10.30 ± 0.09 ns contro 10.01 ns attesi a c (offset non calibrato) | [`../Keysight-1000X/README.md`](../Keysight-1000X/README.md) |

### 3.1 Superficie e metropolitana

File in `../Data/misure_profondita/`: due rivelatori impilati e collegati col
cavo di coincidenza (soglia 80 mV, firmware originale). In superficie NOVE (sopra) e
OTTO (sotto), 20.4 h; in metropolitana QUATTRO e TRE, 1.36 h.

| | coincidenze | durata | rate |
|---|---|---|---|
| superficie | 16 616 | 20.38 h | 13.65 ± 0.11 /min |
| metropolitana | 222 | 1.36 h | 2.73 ± 0.18 /min |

- **Accidentali:** con la finestra hardware di ~3 µs e singoli a 5–7 Hz valgono
  ~0.02/min. Sono lo 0.1% in superficie e meno dell'1% in metropolitana:
  trascurabili.
- **Interpretazione:** in metropolitana passa il **20%** dei muoni. È
  l'attenuazione attesa per qualche decina di metri equivalenti d'acqua di
  materiale sopra la stazione. La pressione media è ~5 hPa più alta in
  metropolitana, ma le due misure sono state fatte in giorni diversi e la
  differenza include anche il meteo: non va usata per ricavare la profondità.
- Le figure di `plot.py` sono in `../Data/misure_profondita/Figures/`.

### 3.2 Efficienza dei piani della torre

`display.py` su `Data/coincidence_nine.txt` (9 rivelatori, firmware
originale, 1.19 h, finestra di 1 ms sul tempo del PC):

| piano | sandwich a 3 | sandwich a 5 |
|---|---|---|
| DUE | 56.8 ± 4.7 % | — |
| TRE | 60.3 ± 4.7 % | 57.2 ± 10.8 % |
| QUATTRO | 58.3 ± 5.5 % | 43.8 ± 12.4 % |
| CINQUE | 65.7 ± 5.6 % | 42.2 ± 11.3 % |
| SEI | 59.9 ± 5.1 % | 50.1 ± 13.4 % |
| SETTE | 55.0 ± 4.8 % | 40.0 ± 11.0 % |
| OTTO | 74.5 ± 5.8 % | — |

Le efficienze sono più basse di quanto ci si aspetta per uno scintillatore
attraversato da un muone. Nel sandwich a 3 circa metà delle coppie di
riferimento sono casuali (per esempio B2 = 254 su N2 = 513 per DUE), quindi il
risultato dipende molto dalla sottrazione. **Da riprendere** con
l'allineamento dei clock di `shower_analysis.py`, che permette finestre di
~300 µs col firmware originale e di pochi µs col firmware open, e con una presa
dati più lunga.

---

## 4. Test

**Prova rapida di tutta la catena** sui dati di
[`../Data/TestData/`](../Data/TestData/): estratti di ~10 minuti o 10 eventi
per ogni tipo di misura, ~11 MB.

```bash
bash analysis_code/run_smoke_tests.sh      # dalla radice del repository, venv attivo
```

Lancia ogni script sui dati di prova in una cartella temporanea, poi i test
unitari, e stampa `OK` o `FAIL` per ogni passo. Va rilanciata dopo ogni
modifica al codice.

**Test unitari e del firmware**, uno per uno:

```bash
cd analysis_code && python -m unittest -v test_analyze_pulses.py
cd RTO6 && python -m unittest -v                      # copia dei file dall'RTO6
cd Keysight-1000X && python -m unittest -v            # acquisizione Keysight (oscilloscopio simulato)
bash Firmware/CosmicWatchOpen/tests/run_tests.sh      # firmware + parser di import_data.py
```

I messaggi `communication error: VisaIOError` durante i test del Keysight
vengono dall'oscilloscopio simulato e sono previsti.

---

## 5. Consigli per le prossime edizioni

1. **Usare il firmware open** per tutte le prese multi-rivelatore: timestamp
   preciso, tempo morto quasi nullo. Con l'allineamento dei clock le
   coincidenze via USB diventano precise a pochi µs.
2. **Nomi unici** sulle microSD prima di collegare i rivelatori allo stesso PC.
3. **Annotare la geometria** di ogni presa dati: posizioni, distanze, foto
   (come `Data/showers_analysis/positions_IMG_6986.txt`). Senza, molte analisi
   non si possono interpretare.
4. **Calibrare l'offset dei canali** dell'oscilloscopio con i rivelatori
   sovrapposti prima di ogni misura di tempo di volo, con gli stessi cavi.
5. **Stimare sempre gli accidentali**, analitici e off-time, prima di
   interpretare un rate di coincidenze.
