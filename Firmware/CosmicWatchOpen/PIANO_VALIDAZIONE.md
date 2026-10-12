# Piano di validazione del firmware open (dopo il confronto su DIECI)

Questo piano chiude i punti rimasti aperti dopo il primo confronto con il
firmware originale ([confronto_fw/CONFRONTO_FW_DIECI.md](confronto_fw/CONFRONTO_FW_DIECI.md))
e porta il firmware a uno stato affidabile per misure di fisica.
[PIANO_MISURE.md](PIANO_MISURE.md) resta il riferimento per la calibrazione in
temperatura (§1); qui viene richiamato, non ripetuto.

Per ogni misura sono indicati: punto aperto, obiettivo, strumenti, procedura,
criterio di accettazione e durata stimata.

---

## 0. Punti aperti

| # | punto aperto | dove è emerso | effetto se non risolto |
|---|---|---|---|
| P1 | Istante di lettura del peak detector: 3.1 µs è sul plateau? C'è overshoot? | spettro ADC diverso dall'originale | scala d'ampiezza e linearità incerte |
| P2 | Scala k (mV/LSB) del firmware nuovo | k è stato ricavato dall'originale | `SiPM[mV]` non è assoluto |
| P3 | Curva di saturazione sopra ~3490 LSB | ricavata dai dati dell'originale | valori sbagliati sugli impulsi grandi |
| P4 | Lettura dell'HV: 28.63 V con l'originale, 29.63 V con il nuovo | stessa scheda | OV e correzione in temperatura incerti |
| P5 | Accuratezza assoluta di ADC1 (soglia) e ADC0 | fit 1577 contro 1540 mV/duty | soglia vera incerta di qualche % |
| P6 | Riscaldamento della scheda: +0.12 contro +0.05 °C/min | BMP280 | elettronica e SiPM più caldi; deriva iniziale |
| P7 | IMU: valori grezzi (nuovo) contro calibrati all'avvio (originale) | accelerometro a 0.784 g | dati IMU non confrontabili |
| P8 | Durata del reset di C5 (5 µs): basta, oppure resta un residuo? | non verificato | sovrastima degli impulsi che seguono da vicino |
| P9 | Tempo morto e perdite misurati in modo indipendente | solo tramite la colonna e il Δt minimo | correzione del live time non validata |
| P10 | Finestra di coincidenza reale e accidentali | solo simulazione del PIO | purezza delle coincidenze non quantificata |
| P11 | Rate e stabilità contro l'originale su un detector sano | DIECI difettoso | equivalenza dei conteggi non dimostrata |
| P12 | Latenza USB: coda al 99° percentile di 2.7 ms (originale 0.5 ms) | scritture su SD | analisi temporali offline fra detector |
| P13 | Timestamp: latenza ≤ 2 µs non contata durante le conversioni di housekeeping | FIRMWARE.md §4.2 | timing fine fra detector |
| P14 | Affidabilità a lungo termine (giorni) e guasti (SD, USB, alimentazione) | mai provata | perdita di dati in run lunghi |
| P15 | Molti detector su USB (9) con il firmware nuovo | non provato | presa dati di gruppo |
| P16 | Dipendenza dalla temperatura (α) | PIANO_MISURE §1 | correzione in temperatura non validata |
| P17 | Guasto hardware di DIECI | rate instabile con entrambi i firmware | detector inutilizzabile |
| P18 | Non linearità differenziale dell'ADC dell'RP2040 | nota per il chip | artefatti negli spettri |

---

## 1. Lavoro sul firmware e sugli strumenti, prima delle misure

Tutto fattibile senza detector. Ogni cambiamento della lettura dell'ADC comporta
un **aumento di versione** e il **ritardo di lettura scritto nell'intestazione**
dei file.

| # | aggiunta | serve a |
|---|---|---|
| S1 | **Modalità diagnostica di lettura** (`ADC_DIAG = true`): 8 conversioni per evento a ritardi fissi (0.5, 1, 2, 3, 5, 10, 30, 100 µs) in colonne aggiuntive; tempo morto ~110 µs per evento, solo in questa modalità | P1, P8 |
| S2 | Chiavi `ADC_DELAY_US` (istante di lettura, ammesso anche sotto i 3 µs perché il PIO gestisce la finestra da solo) e `RESET_US` | P1, P8 |
| S3 | Chiave `COINC_WINDOW_US` (ritardi del PIO calcolati al boot) | P10 |
| S4 | **Attesa a basso consumo:** il core 0 dorme fra un giro e l'altro (con una sveglia periodica); il core 1 dorme finché l'interrupt del PIO non lo risveglia | P6 |
| S5 | **Comando `selftest`:** legge ADC0, ADC1 e ADC2 con e senza i pull-down interni attivi, e stampa i valori grezzi in LSB e in mV | P4, P5 |
| S6 | Opzione `IMU_CALIBRATE_AT_BOOT` (offset sottratti all'avvio, come l'originale) e stampa degli offset nel log | P7 |
| S7 | Scrittura su SD a blocchi più piccoli e uscita USB servita fra un blocco e l'altro | P12 |
| S8 | Timestamp del trigger preso in hardware: il PIO conta i cicli dal marker | P13 |
| S9 | Strumenti nel repo (`tools/`): `longcap.py` (presa con timestamp del PC), `analyze_run.py` (metriche standard), `compare_runs.py` (report A/B con grafici) | tutti |
| S10 | Test automatici estesi: simulazione del PIO (pioemu) nei test, a condizione che tu approvi l'installazione nel venv | regressioni |

---

## 2. Misure al banco (strumenti, detector aperto)

**Strumenti:**
- multimetro;
- oscilloscopio: l'RTO6 o il Keysight 1000X, con gli script già presenti in
  `RTO6/` e `Keysight-1000X/`;
- generatore di impulsi a due canali con fronti di pochi ns e ritardo regolabile.
  Il WaveGen del Keysight 1000X può bastare per un canale;
- tester di corrente USB.

### B1. Tensioni di riferimento e accuratezza dell'ADC (P4, P5)
- **Procedura.** Con il multimetro si misurano HV (P1, pin HV), VRef (C25, atteso
  2.500 V), TriggerThreshold (test point o R19) a 3–4 soglie diverse, e la
  tensione Signal a riposo. Lo stesso si legge con `selftest` (S5) nel firmware
  nuovo, e dalla diagnostica di boot nell'originale.
- **Accettazione.**
  - ADC del firmware nuovo entro ±0.5 % dal multimetro su ciascun canale;
  - la differenza con/senza pull-down spiega lo scarto dell'originale (atteso
    −2/−3 %).
- **Durata:** 1 ora.

### B2. Ritardo di lettura e linearità del peak detector (P1, P2, P3)
- **Procedura.** Si iniettano impulsi dal BNC: forma simile al SiPM (salita di
  pochi ns, discesa ~100 ns), 20 ampiezze da appena sopra soglia fino alla
  saturazione, 1000 impulsi per ampiezza, a 10–50 Hz.
  1. Con S1 si ricava, per ogni ampiezza, la curva ADC(ritardo): salita, plateau,
     droop.
  2. Con S2 si conferma l'ADC medio al ritardo scelto.
  3. In parallelo, all'oscilloscopio: ampiezza sul BNC e su TP1 (Amplifier).
- **Risultati attesi:**
  - ritardo ottimale: inizio del plateau, con margine;
  - droop sul plateau, in %/µs;
  - curva ADC contro ampiezza al BNC: linearità e valore di k;
  - punto di saturazione e curva sopra soglia (LUT nuova).
- **Accettazione.**
  - plateau piatto entro lo 0.5 % su almeno 5 µs al ritardo scelto;
  - non linearità integrale entro l'1 % sotto la saturazione;
  - nessun valore sotto il minimo garantito dalla soglia.
- **Durata:** mezza giornata più l'analisi.

### B3. Reset del peak detector (P8)
- **Procedura.** Coppie di impulsi (uno grande seguito da uno piccolo) con
  intervallo da 15 a 200 µs; in alternativa, con S1, l'ADC letto subito dopo il
  reset. Si fa variare `RESET_US` da 1 a 20 µs.
- **Accettazione.** L'impulso piccolo non risente del grande per qualunque
  intervallo maggiore del tempo morto.
- **Durata:** 2 ore.

### B4. Tempo morto e perdite (P9)
- **Procedura.** Impulsi periodici a 10 Hz, 100 Hz, 1, 5, 10 e 20 kHz, per 60 s
  ciascuno. Si confronta la frazione registrata con `1/(1 + f·τ)` e la si
  verifica sui contatori `lost` e `stuck`. Poi si misura il rate massimo sostenuto
  su USB e su SD con il firmware nuovo, e lo stesso con l'originale.
- **Accettazione.**
  - τ ricavato dal fit uguale entro il 10 % alla mediana della colonna `Deadtime[s]`;
  - perdite pari a zero fino al rate massimo dichiarato;
  - nessun evento perso in silenzio.
- **Durata:** 2 ore.

### B5. Finestra di coincidenza (P10)
- **Procedura.** Due detector con il firmware nuovo, collegati con il cavo RJ45,
  e lo stesso impulso nei due BNC con ritardo Δ da −6 a +6 µs a passi di 0.25 µs.
  Si misura la frazione di eventi con flag 1 in funzione di Δ, ripetendo con
  `COINC_WINDOW_US` = 1, 3 e 10 µs.
- **Accettazione.** Fronti della finestra a ±W, entro 0.1 µs; finestra
  simmetrica; efficienza al centro del 100 %.
- **Durata:** 2 ore.

### B6. Consumo e riscaldamento (P6)
- **Procedura.** Corrente con il tester USB per: firmware originale, firmware
  nuovo 1.0.0, firmware nuovo con S4. Poi un'ora di riscaldamento per ciascuno,
  dall'accensione a freddo, nella stessa posizione e con la scheda chiusa nel
  case: si registra T(t) dal BMP280 e con un termometro a contatto sul case,
  vicino al SiPM.
- **Accettazione.** Con S4 consumo e ΔT di equilibrio uguali o inferiori a quelli
  dell'originale, a parità di latenza (verificata come in C4).
- **Durata:** mezza giornata, in parallelo ad altro.

### B7. IMU (P7)
- **Procedura.** Calibrazione a 6 posizioni: il detector appoggiato su ciascuna
  faccia per 60 s. Si ricavano offset e scala di ciascun asse, con il firmware
  nuovo in modo grezzo e poi con S6.
- **Accettazione.** |a| = 1.000 ± 0.02 g in tutte le orientazioni dopo la
  calibrazione; giroscopio fermo entro 0.5 °/s.
- **Durata:** 30 minuti.

### B8. Diagnosi di DIECI (P17)
- **Procedura.** All'oscilloscopio, sul BNC (segnale SiPM) e sull'HV a valle del
  connettore (lato SiPM, P1), durante un episodio di rate basso; si muovono
  delicatamente il connettore e la schedina SiPM. Si controllano saldature,
  connettore, gel ottico e la nastratura (perdite di luce).
- **Accettazione.** Causa individuata e rate stabile (χ²/ndf ≈ 1 su 1 h).
- **Durata:** 1–2 ore.

### B9. ADC dell'RP2040: non linearità differenziale (P18)
- **Procedura.** Rampa lenta su ADC0 tramite impulsi a ampiezza crescente
  (oppure con un generatore DC sul nodo Signal, staccando il buffer, solo se
  sicuro). Si fa l'istogramma della densità dei codici.
- **Risultato:** tabella di correzione, da applicare offline o nel firmware se gli
  artefatti superano il 5 % per bin.
- **Durata:** 2 ore.

---

## 3. Misure con raggi cosmici (detector chiuso, nessuno strumento)

Per i confronti si usano **detector sani**, non DIECI. In ogni run vengono
registrati la versione del firmware, il `config.txt` e la diagnostica di boot.
Per l'originale la correzione in temperatura va **spenta**.

### C1. Confronto A/B/A su un detector sano (P11)
- **Procedura.** Firmware nuovo 30 minuti, originale 30 minuti, nuovo 30 minuti;
  ripetuto 2 volte su due detector diversi. Stessa soglia, stessa posizione,
  nessun intervento fra i run.
- **Metriche:** rate sul live time, stabilità (χ² su bin da 60 s), Δt minimo e
  distribuzione degli intervalli, spettro (solo in forma, perché la scala ADC è
  diversa), sensori, clock.
- **Accettazione.** Rate del nuovo uguale alla media dei due run dell'originale
  entro 2σ statistici (≈ 2 % con 30 minuti a 3–4 Hz). Per il nuovo il χ²/ndf
  deve essere compatibile con 1.
- **Durata:** 3 ore per detector.

### C2. Curva rate-soglia (P5, P11)
- **Procedura.** Soglia da 30 a 200 mV: a passi di 1 mV fino a 100 mV, poi di
  5 mV, con 5 minuti per punto. Con l'originale bastano pochi punti, perché ha
  passi di 13 mV.
- **Risultati attesi:** plateau a 1 e 2 fotoelettroni; scelta motivata della
  soglia di lavoro; confronto della posizione dei gradini fra i due firmware,
  che verifica l'accuratezza della soglia.
- **Durata:** 1 notte (automatica, con il comando `threshold` inviato dallo script).

### C3. Coincidenze con detector impilati (P10)
- **Procedura.** Due detector impilati, 12 ore con il firmware nuovo e 12 ore con
  l'originale (entrambi i detector con lo stesso firmware).
- **Metriche:** rate di coincidenza e spettro degli eventi coincidenti (MPV,
  confrontabile solo all'interno dello stesso firmware); accidentali stimate dai
  rate di singola con `2·W·R1·R2`; coincidenze ricostruite offline dai timestamp
  del PC che hanno il flag (efficienza del flag).
- **Accettazione.** Rate di coincidenza uguali fra i due firmware entro 2σ.
  Efficienza del flag superiore al 99 % sulle coincidenze offline entro 3 µs, se
  la latenza USB lo permette (altrimenti si usa C4).
- **Durata:** 2 notti.

### C4. Latenza e timestamp (P12, P13)
- **Procedura.** Presa con `longcap.py` per 1 h, con le varianti: SD sì/no, OLED
  sì/no, versione con S7 e S8. Con due detector che ricevono lo stesso impulso dal
  generatore si misura la differenza dei timestamp del Pico (risoluzione 1 µs) e
  la sua dispersione.
- **Accettazione.** Latenza al PC sotto 1 ms al 99° percentile; dispersione dei
  timestamp fra i due Pico entro 2 µs più la deriva lineare dei quarzi.
- **Durata:** mezza giornata.

### C5. Calibrazione in temperatura (P16)
Come in [PIANO_MISURE.md §1](PIANO_MISURE.md), **dopo** B2 (scala stabile) e S4
(riscaldamento ridotto). Usare la T del termometro a contatto oltre a quella del
BMP280.

---

## 4. Affidabilità a lungo termine

### L1. Run lungo (P14)
- **Procedura.** Una settimana, 2 detector in coincidenza, firmware nuovo, su SD e
  USB contemporaneamente.
- **Verifiche:**
  - numerazione continua;
  - nessun reset del watchdog;
  - SD integra; righe HK presenti ogni 60 s;
  - deriva della soglia (trim) e della baseline;
  - confronto dei conteggi fra file SD e flusso USB (devono coincidere).
- **Durata:** 1 settimana, senza presidio.

### L2. Prove di guasto (P14)
Ognuna ripetuta 5 volte; il comportamento atteso è quello di FIRMWARE.md §13:

| prova | atteso |
|---|---|
| microSD tolta e rimessa durante il run | logging su SD fermato con un avviso, USB attivo, nessun blocco |
| USB staccato e riattaccato | log di boot e intestazione reinviati alla riconnessione |
| programma che tiene la porta aperta senza leggere per 10 minuti | acquisizione e SD intatte; righe USB perse contate |
| interruzione dell'alimentazione durante la scrittura | file leggibile, al massimo 10 s persi |
| 20 accensioni con il partner (stesso cavo), 20 reset sfalsati di 0–2 s | modo C in tutti i casi |
| alimentazione da power bank, senza PC | acquisizione normale su SD |

**Durata:** mezza giornata.

### L3. Molti detector su USB (P15)
- **Procedura.** 9 detector su un hub USB con `import_data.py`, per 12 ore.
- **Accettazione.** Nessun buco nella numerazione su alcun detector; latenza al
  99° percentile sotto 2 ms; nomi tutti distinti senza microSD (flash o `CW-xxxx`).
- **Durata:** 1 notte.

---

## 5. Ordine consigliato

| fase | contenuto | prerequisiti | durata |
|---|---|---|---|
| 1 | S1, S2, S5, S9; poi B1, B2, B3, B8 | detector aperto, oscilloscopio, generatore | 2 giorni |
| 2 | Versione 1.1.0 con il ritardo e il k scelti, LUT nuova; B4, B5, B9 | fase 1 | 1 giorno |
| 3 | S4, S6, S7, S8 (versione 1.2.0); B6, B7, C4 | fase 2 | 2 giorni |
| 4 | C1, C2, C3 su detector sani | fase 3 | 1 settimana |
| 5 | C5 (temperatura), L1, L2, L3 | fase 4 | 2 settimane |

Al termine di ogni fase si aggiornano FIRMWARE.md (comportamento) e un registro
delle versioni, con le differenze di scala ADC fra le versioni, in modo che ogni
file di dati sia riconducibile alla sua calibrazione.

## 6. Criteri di "firmware validato"

Il firmware si considera pronto per misure di fisica quando:
1. il rate coincide con l'originale entro il 2 % su almeno due detector sani (C1);
2. tempo morto e perdite sono verificati con il generatore (B4);
3. l'istante di lettura è sul plateau e k è misurato con un'incertezza sotto il
   2 % (B2);
4. HV e soglia lette concordano con il multimetro entro lo 0.5 % (B1);
5. la finestra di coincidenza è misurata (B5) e il rate di coincidenza è
   compatibile con l'originale (C3);
6. una settimana di run e le prove di guasto passano senza perdite non dichiarate
   (L1, L2).
