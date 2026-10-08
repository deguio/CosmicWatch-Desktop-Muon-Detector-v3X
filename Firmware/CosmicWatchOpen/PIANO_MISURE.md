# Piano: calibrazione in temperatura, confronto firmware, sviluppi

Per tutto ciò che segue valgono le stesse definizioni:
- **firmware nuovo** = `cosmicwatch_open.uf2`; **firmware vecchio** = `CosmicWatch_v3X.1.1.52.uf2`;
- **A = ADC − b**, dove b è la baseline dell'ADC0 a riposo;
- **MPV** = moda della distribuzione di A per gli eventi coincidenti, stimata con una
  Moyal troncata.

Per ogni punto distinguo fra ciò che si misura, ciò che ci si aspetta (stime, da
verificare) e i criteri per decidere.

---

## 0. Aggiunte al firmware prima di iniziare

**Stato: F1-F5 implementate** (firmware 1.0.0 attuale). Le chiavi e il formato
della riga HK sono descritti nel README. Restano da provare sul detector.

| # | aggiunta | perché |
|---|---|---|
| F1 | Riga di housekeeping ogni 60 s, su file e su USB: `# HK t=… T=… P=… HV=… Thr_set=… Thr_meas=… base=… corr=… lost=…` | HV e soglia misurata non sono nelle colonne. Senza, non si separa il drift del SiPM da quello del booster e della soglia. |
| F2 | Misura periodica della baseline b durante il run (peak detector resettato e PIO disarmato per circa 100 µs, contati come tempo morto) | b può derivare con T (offset degli operazionali) e va sottratta per ricavare α (vedi 1.4). |
| F3 | Filtro sulla T del BMP280: mediana su 5 letture e rifiuto dei salti > 2 °C/s | Le letture spurie (tipo −5 °C) non devono entrare né nelle colonne né nella correzione. |
| F4 | Correzione empirica `c = exp(−α·(T − T_ref))`, con `GAIN_TEMP_COEFF` (α) letto da `config.txt` e alternativa al modello OV | Per applicare l'α misurato nella fase 1. |
| F5 | Fattore di correzione registrato (nella riga HK oppure come colonna opzionale) | Rende la correzione tracciabile e reversibile in analisi. |

Il firmware vecchio non ha F1-F5. Per i confronti, nel suo `config.txt` va
spento `USE_REALTIME_GAIN_CORRECTION`, che il file `TOP_C_002.txt` mostra acceso.

---

## 1. Calibrazione della dipendenza dalla temperatura

### 1.1 Obiettivo
Misurare per ciascun detector

```
α = d ln(MPV) / dT        [%/°C]
```

Va confrontato con il modello del solo guadagno, `α_G = −(21.5 mV/°C) / OV ≈ −0.37 %/°C`
(con OV ≈ 5.75 V). Per lo stesso intervallo di T servono anche il rate di
singole a soglia fissa, b(T), V_HV(T) e V_thr(T).

### 1.2 Setup
- **Detector.** Due detector impilati in coincidenza, entrambi con il firmware nuovo
  (con F1-F3). Dove possibile, la stessa calibrazione andrebbe ripetuta su altri
  detector, perché α dipende dalla V_br del singolo SiPM.
- **Ambiente.** Una camera climatica, oppure una scatola isolata con riscaldatore e
  termostato, o un frigorifero. Raffreddando sotto il punto di rugiada si forma
  condensa: usare un sacchetto essiccante e risalire lentamente.
- **Termometro di riferimento.** Una termocoppia o PT100 sul case di alluminio in
  corrispondenza del SiPM, registrata dal PC con timestamp NTP. Il BMP280 sta
  sulla scheda principale, vicino a Pico e booster (0.5 W in tutto): misura la
  scheda, non il SiPM.
- **Configurazione.** Soglia fissa (`trigger_voltage_mV` uguale per tutto il run),
  `USE_REALTIME_GAIN_CORRECTION = false`, dati salvati sia su SD sia con
  `import_data.py`, così ogni evento ha anche il tempo del PC.

### 1.3 Profilo termico
- **Gradini:** da 5 a 40 °C a passi di 5 °C (8 punti), in salita e poi in discesa.
- **Permanenza:** 2.5 h per gradino. I primi 45 min si scartano come assestamento:
  il criterio è che la T di riferimento sia stabile entro 0.2 °C su 15 min.
- **Durata totale:** circa 40 h.
- **Statistica attesa** (stima): ~0.4 Hz di coincidenze, quindi ~2500 eventi utili
  per gradino. Con σ(MPV) ≈ 1.4·s/√N, dove s ≈ 5 mV è la larghezza della Moyal, si
  ottiene ≈ 0.7 % sull'MPV.
- **Effetto atteso:** su 35 °C una variazione dell'MPV del 13-17 %. È quindi
  misurabile anche la non linearità.

### 1.4 Analisi
1. **Selezione:** eventi con `Flag = 1`, nelle finestre stabili di ciascun gradino.
2. **Baseline:** A = ADC − b(t), con b presa dalle righe HK (F2). Se b non viene
   sottratta, α risulta diluito di un fattore (MPV − b)/MPV ≈ 0.8.
3. **Fit per gradino.** Fit Moyal di massima verosimiglianza non binnato, con la pdf
   normalizzata sull'intervallo `[μ − 1.5 s, μ + 3 s]` e iterato.
   - pdf: `f(x) = exp(−(λ + e^(−λ))/2) / (s·√(2π))`, con `λ = (x − μ)/s`; è
     `scipy.stats.moyal`, dove `loc` coincide con μ e quindi con l'MPV.
   - Non binnato perché l'ADC dell'RP2040 ha errori di non linearità differenziale
     ai codici 512, 1536, …, che distorcerebbero un istogramma.
4. **Fit di α.** Fit lineare pesato `ln MPV_i = ln MPV_0 + α·(T_i − T_0)`, fatto
   separatamente per:
   - T di riferimento e T del BMP280: il confronto stima l'offset e il ritardo
     termico fra i due;
   - rampa in salita e rampa in discesa: la differenza è l'isteresi termica.
5. **Fit alternativo in funzione dell'OV.** `ln MPV = a + β·ln(OV)`, con
   `OV = V_HV(T) − V_br(T)`. Separa il contributo del SiPM da quello del booster.
   Un β > 1 indica il contributo di PDE e cross-talk oltre al guadagno.
6. **Altre dipendenze da misurare:** rate di singole contro T a soglia fissa, che
   riflette gamma e dark noise sopra soglia; b(T); V_thr misurata (ADC1) contro T,
   che deve restare piatta grazie al trim in anello chiuso.

### 1.5 Sistematiche da quotare

| effetto | come si stima |
|---|---|
| Intervallo di fit | si variano gli estremi di ±0.5 s |
| Modello della pdf | si confronta con la mediana di A sopra 0.5·MPV, un taglio relativo ricalcolato per ogni gradino |
| T del SiPM diversa da quella del BMP280 | differenza fra α ottenuto con la T di riferimento e con quella del BMP280, e isteresi fra salita e discesa |
| Fondo di coincidenze a bassa ampiezza | si sposta l'estremo sinistro del fit; l'MPV deve stare almeno 3 volte sopra la soglia a tutte le T |
| Pressione | l'MPV non ne risente; serve solo per correggere i rate (1.6) |

### 1.6 Validazione
- **Run di controllo:** si ripete il profilo termico (anche solo 3 gradini) con
  `USE_REALTIME_GAIN_CORRECTION = true`,
  `GAIN_CORRECTION_MODEL = empirical` e `GAIN_TEMP_COEFF = α` (F4).
- **Criteri di accettazione:**
  - MPV piatto entro ±1 % su tutto l'intervallo di temperatura;
  - rate di coincidenze, corretto per la pressione con
    `R = R_0·exp(β_p·(P − P_0))` e β_p fittato (tipicamente −0.1…−0.2 %/hPa), piatto
    entro la statistica.
- **Decisione sul modello:** se α è compatibile con il modello OV entro le
  incertezze, si tiene il modello fisico; altrimenti si usa α come default
  per-detector, salvato in `config.txt`.

---

## 2. Confronto fra firmware nuovo e vecchio

### 2.1 Regole generali
- Stesso detector, stessa posizione, stessa soglia (`trigger_voltage_mV`
  ha lo stesso significato nei due firmware).
- `USE_REALTIME_GAIN_CORRECTION = false` in entrambi.
- **Alternanza:** run alternati di 12 h vecchio / 12 h nuovo, ripetuti per almeno
  3 cicli. In questo modo le variazioni di T e P si mediano invece di finire
  tutte su uno dei due firmware.
- **Coincidenze:** il protocollo di coincidenza è diverso, quindi i test 2.6 vanno
  fatti con due detector entrambi vecchi oppure entrambi nuovi.
- **Iniezione di impulsi:** in diversi test si iniettano impulsi dal BNC, che
  arriva sul nodo del segnale SiPM attraverso C13 (200 pF). La forma differisce da
  quella del SiPM, ma per test relativi va bene. Serve un generatore a due canali
  con fronti veloci e ritardo regolabile.

### 2.2 Tempo morto

| metodo | cosa misura | atteso nuovo | atteso vecchio |
|---|---|---|---|
| Incrementi della colonna `Deadtime[s]` per evento | τ dichiarato | ≈ 12 µs | ≈ 375-400 µs (dai dati esistenti) |
| Distribuzione degli intervalli Δt fra eventi consecutivi | τ reale: esponenziale che va a zero per Δt < τ | taglio a circa 12 µs | taglio a circa 0.4 ms |
| Impulsi periodici a frequenza f (10 Hz … 20 kHz) | frazione registrata, contro `1/(1 + f·τ)` (non paralizzabile) | perdite circa f·τ: 1.2 % a 1 kHz | perdite evidenti già a 100 Hz |
| Oscilloscopio su GPIO22 (trigger), GPIO21 (reset) e GPIO0/1 (coincidenza) | sequenza temporale per evento | coincidenza tenuta 3.1 µs, reset dopo circa 5 µs | da documentare |

Note:
- **Statistica del test Δt** sui dati cosmici: a 3 Hz la frazione di eventi con
  Δt < 12 µs è circa 3.6·10⁻⁵, cioè circa 9 al giorno. Per il firmware nuovo servono
  alcuni giorni, per il vecchio basta un giorno.
- **Massimo rate sostenuto:** si aumenta f finché la coda di 2048 eventi o l'USB
  perdono dati; il comando `status` riporta `lost: queue/USB/SD`. Va confrontato
  con i 700 Hz dichiarati per il firmware vecchio.

### 2.3 Risposta in ampiezza
- **Spettri cosmici:** stesso detector, stessa soglia. Si confrontano lo spettro ADC
  di tutti gli eventi (test di Kolmogorov-Smirnov), l'MPV delle coincidenze e la
  baseline.
  - Atteso: compatibili. Possibile piccola differenza dovuta all'istante di
    campionamento dell'ADC, circa 3.1 µs dopo il trigger nel nuovo e circa 2.3 µs
    nel vecchio, combinata con la scarica lenta di C5.
- **Scansione in ampiezza con impulsi iniettati:** ADC in funzione dell'ampiezza
  iniettata, dalla soglia alla saturazione. Si misurano linearità, punto di
  saturazione e risoluzione (σ dell'ADC a ampiezza fissa).
- **Scansione del ritardo di campionamento** (solo firmware nuovo, con
  `reset_us` e la finestra resi configurabili): ADC in funzione dell'istante di
  lettura. Misura la velocità di scarica del peak detector e verifica che 3 µs sia
  sul plateau.
- **Calibrazione della costante k (mV/LSB):** si confrontano le ampiezze del SiPM
  viste all'oscilloscopio dal BNC con l'ADC. Il valore attuale, 0.0706 mV/LSB, è
  stato ricavato dai dati del firmware vecchio.

### 2.4 Soglia
- **Accuratezza e stabilità:** soglia misurata (ADC1) contro quella impostata, in
  funzione del tempo e della temperatura. Atteso: nuovo entro 0.3 mV grazie al trim
  in anello chiuso; vecchio entro circa 13 mV, che è il passo del suo PWM.
- **Curva di soglia:** rate in funzione della soglia. Nel nuovo, a passi di 1 mV
  da 30 a 200 mV; nel vecchio solo a passi di 13 mV. La curva del nuovo permette di
  vedere i plateau a singolo e multiplo fotoelettrone e di scegliere la soglia di
  lavoro.

### 2.5 Dipendenza dalla temperatura (con i dati della fase 1)
Si confrontano tre configurazioni:
1. firmware vecchio con la sua correzione accesa;
2. firmware nuovo senza correzione;
3. firmware nuovo con α (F4).

L'indicatore è la pendenza residua di MPV e rate di coincidenze in funzione di T.
Basta una rampa ridotta a 3 gradini per ciascuna.

### 2.6 Coincidenze (due detector con lo stesso firmware)

| test | come | atteso nuovo |
|---|---|---|
| Finestra | stesso impulso nei due BNC con ritardo Δ fra −6 e +6 µs, a passi di 0.25 µs; frazione con flag 1 in funzione di Δ | gradino a ±3.0 µs (lo dice la simulazione del PIO) |
| Accidentali | detector lontani o affiancati, singole alzate con una sorgente gamma fino a circa 100 Hz; rate di coincidenze contro `2·W·R1·R2` | con W = 3.0 µs e 100 Hz ciascuno: 0.06 Hz, misurabile in poche ore |
| Efficienza sui muoni | detector impilati; frazione di coincidenze ricostruite offline dai timestamp del PC che hanno il flag 1 | vicino al 100 % |
| Avvio | 20 accensioni con un solo cavo di alimentazione, 20 reset sfalsati di 0-2 s | modo C in tutti i casi |

### 2.7 Affidabilità
- **Run lungo:** almeno 1 settimana. Si controlla che `Event` sia continuo, che i
  timestamp siano monotoni, che non ci siano reset del watchdog (lo segnala la
  diagnostica) e che i file sulla SD siano integri.
- **Situazioni anomale:** USB staccato e riattaccato, programma che tiene la porta
  aperta senza leggere, microSD tolta durante il run, mancanza di alimentazione.
  Dopo una mancanza di alimentazione si perdono al massimo 10 s di dati, cioè
  l'intervallo fra due `f_sync`.
- **Consumo:** con un misuratore USB. Atteso circa 0.5 W in entrambi.
- **Deriva dell'orologio:** timestamp del Pico contro tempo NTP del PC
  (`import_data.py`). La deriva viene dal quarzo, quindi dovrebbe essere uguale per
  i due firmware.

### 2.8 Ordine consigliato e durata

| passo | durata |
|---|---|
| Accensione e checklist del README | ½ giornata |
| 2.2 e 2.6 con generatore | 1 giornata |
| 2.4 curva di soglia | ½ giornata |
| 2.3 spettri a confronto | 3 giorni con run alternati |
| Fase 1, calibrazione in temperatura | circa 2 giorni per coppia di detector |
| 2.5 e validazione 1.6 | 1-2 giorni |
| 2.7 run lungo | 1 settimana, in parallelo al resto dell'analisi |

---

## 3. Sviluppi possibili

In ordine di rapporto fra utilità e sforzo.

### Firmware: misure e qualità dei dati
1. **F1-F5** della sezione 0.
2. **Timestamp del trigger preso in hardware.** Il PIO può memorizzare l'istante
   del trigger, con un contatore di cicli oppure facendo copiare dal DMA il
   registro del timer al momento del marker. Si elimina l'incertezza di al massimo
   2 µs dovuta alle letture di housekeeping.
3. **Δt fra i due detector con passo di 8 ns.** Il PIO può contare i cicli fra il
   proprio trigger e il fronte di discesa della linea del partner, e il risultato
   può diventare una colonna in più. È una misura di tempo di volo fra detector a
   costo hardware zero, interessante per l'analisi della velocità dei muoni. La
   risoluzione sarebbe limitata dal time walk del comparatore e dalla
   sincronizzazione degli ingressi (circa 16 ns), non dai clock indipendenti delle
   due schede.
4. **Finestra di coincidenza e `reset_us` configurabili** da `config.txt`:
   servono ai test 2.3 e 2.6.
5. **Soglia relativa alla baseline**, per esempio `trigger_above_baseline_mV`:
   b viene misurata periodicamente (F2) e la soglia la segue.
6. **Correzione della non linearità differenziale dell'ADC dell'RP2040**, con una
   tabella ai codici noti. Utile per gli spettri a bassa ampiezza.
7. **Più di due detector:** oggi la coincidenza è solo a coppie. Si possono
   portare le linee SPI0 del connettore RJ45 su uno schema a bus, oppure
   aggiungere un modulo master.

### Firmware: tempo assoluto e gestione
8. **Ora assoluta:** un comando USB `time <epoch>` imposta l'orologio. Così i file
   avrebbero date di creazione corrette (oggi FatFs usa una data fissa) e si
   potrebbe scrivere un timestamp UTC nell'intestazione.
9. **Ingresso PPS da GPS** su un GPIO libero (16, 17 o 18): disciplina il clock e
   dà un tempo assoluto con precisione dell'ordine del µs. Permette coincidenze
   fra detector lontani analizzate offline.
10. **Rotazione dei file** (per esempio uno al giorno) e comandi `start`, `stop`,
    `newfile` da USB.
11. **Interfaccia OLED:** etichette sui rate, pagine alternate (tempo morto, eventi
    persi, baseline, α), avviso quando la correzione è attiva.
12. **Buzzer:** tono o sequenza al boot quando viene trovato il partner.

### Hardware: modifiche piccole
13. **Stabilizzazione dell'HV:** PWM filtrato da GPIO7, attraverso circa 330 kΩ,
    verso il nodo FB del MAX5026, con anello di controllo sull'ADC2. Tiene
    costanti guadagno, PDE e cross-talk invece di correggerli a posteriori.
    Da sviluppare dopo la fase 1, che dice quanto il modello software residuo
    sia già sufficiente.
14. **Sensore di temperatura vicino al SiPM:** un TMP117 sull'I²C esistente
    oppure un DS18B20 su un GPIO libero. Elimina la principale sistematica della
    fase 1.

### Infrastruttura
15. **Prove su banco automatizzate:** generatore pilotato da PyVISA (già nel venv)
    più uno script che esegue 2.2, 2.4 e 2.6 e confronta i risultati con quelli di
    riferimento a ogni nuova versione del firmware.
16. **Integrazione continua:** build e `tests/run_tests.sh` a ogni commit. La
    toolchain è già scriptata in `build.sh`.
