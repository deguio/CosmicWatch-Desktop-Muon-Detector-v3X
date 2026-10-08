# CosmicWatch v3X open firmware: descrizione della logica

Descrive il firmware 1.0.0 (`src/`), con clock di sistema a 125 MHz. Tempi e
periodi sono quelli impostati nel codice. Le prestazioni misurate sul detector
reale vanno verificate con [PIANO_MISURE.md](PIANO_MISURE.md): il firmware non è
ancora stato provato su hardware.

Indice:
1. [Architettura](#1-architettura)
2. [Sequenza di boot](#2-sequenza-di-boot)
3. [Acquisizione di un evento](#3-acquisizione-di-un-evento)
4. [Tempo morto](#4-tempo-morto)
5. [Coincidenza fra due detector](#5-coincidenza-fra-due-detector)
6. [Soglia di trigger](#6-soglia-di-trigger)
7. [Da ADC a SiPM[mV] e correzione in temperatura](#7-da-adc-a-sipmmv-e-correzione-in-temperatura)
8. [Sensori e housekeeping: cosa viene letto e quando](#8-sensori-e-housekeeping-cosa-viene-letto-e-quando)
9. [Uscita dei dati](#9-uscita-dei-dati)
10. [config.txt](#10-configtxt)
11. [Comandi da seriale](#11-comandi-da-seriale)
12. [OLED, LED e buzzer](#12-oled-led-e-buzzer)
13. [Robustezza e gestione degli errori](#13-robustezza-e-gestione-degli-errori)
14. [Risorse dell'RP2040](#14-risorse-dellrp2040)
15. [Differenze rispetto al firmware originale](#15-differenze-rispetto-al-firmware-originale)

---

## 1. Architettura

Il lavoro è diviso fra quattro esecutori hardware. Solo il PIO e il core 1
partecipano all'acquisizione; tutto il resto avviene sul core 0, quindi non
aggiunge tempo morto.

| esecutore | compiti | vincoli |
|---|---|---|
| **PIO0, macchina a stati 0** (`daq.pio`) | aspetta il fronte del comparatore, abbassa la propria linea di coincidenza, campiona quella del partner per 3.07 µs, segnala alla CPU | 16 istruzioni, tempi deterministici a 8 ns per ciclo |
| **Core 1** (`daq.c`) | timestamp, conversione ADC del picco, reset del peak detector, riarmo del PIO, coda degli eventi; fra un evento e l'altro esegue le conversioni ADC chieste dal core 0 e i campioni di baseline | gira da RAM con gli interrupt disabilitati; non chiama funzioni e non accede alla flash (verificato sul disassemblato) |
| **DMA, un canale** (`threshold.c`) | aggiorna il compare del PWM della soglia a 1 MHz (dithering) | nessun intervento della CPU, salvo un interrupt di riavvio ogni ~71 min |
| **Core 0** (`main.c` e moduli) | boot e diagnostica, formattazione delle righe, USB, microSD, OLED, BMP280, MPU-6050, LED, buzzer, soglia, correzione in temperatura, comandi, watchdog | gestisce tutti gli interrupt (USB, DMA) |

I due core comunicano attraverso:
- **la coda degli eventi:** un buffer circolare in RAM di 2048 eventi da 24 byte
  (timestamp, tempo morto cumulativo, ADC, flag). Il core 1 è l'unico che scrive,
  il core 0 l'unico che legge; non servono lock, bastano barriere di memoria;
- **le richieste di servizio** dal core 0 al core 1: media di N conversioni su un
  canale ADC, un campione di baseline, pausa (per le scritture in flash).

Il core 1 è l'unico che usa l'ADC durante la presa dati. Così si evitano
conflitti sul multiplexer e una conversione non può mai interrompere la lettura
di un picco.

**File sorgente:**

| file | contenuto |
|---|---|
| `board.h` | mappa dei pin, costanti analogiche (VREF, partitori, bias), costanti del SiPM |
| `main.c` | boot, ciclo principale del core 0, formattazione, USB, OLED, LED, buzzer, comandi |
| `daq.c`, `daq.pio` | acquisizione (core 1 e PIO) |
| `coinc.c` | handshake di coincidenza al boot |
| `threshold.c` | soglia (PWM, DMA, calibrazione, trim) |
| `cwcore.c` | logica pura e testata su host: conversioni, BMP280, filtro T/P, `config.txt`, formato delle righe |
| `sensors.c` | I²C, BMP280, MPU-6050 |
| `oled.c`, `font5x7.h` | display SSD1306 |
| `sdcard.c`, `logger.c` | microSD in SPI per FatFs; file di dati e buffer |
| `flashcfg.c` | nome e soglia salvati nell'ultimo settore della flash |

---

## 2. Sequenza di boot

I tempi sono indicativi: dipendono soprattutto dalla microSD.

| # | passo | durata | note |
|---|---|---|---|
| 1 | Inizializzazione USB CDC; GPIO21 alto (C5 tenuto scarico); LED spenti; GPIO6 basso | ms | l'USB si enumera in background sugli interrupt |
| 2 | **Handshake di coincidenza** (§5) | 0.05–3 s | è il primo passo, così due detector accesi insieme lo raggiungono insieme. Con partner trovato si accendono i due LED |
| 3 | Lettura di nome e soglia dalla flash; se mancano, nome `CW-xxxx` dall'ID unico della scheda | µs | |
| 4 | Bus I²C a 400 kHz e scansione degli indirizzi; OLED e schermata iniziale | ~30 ms | |
| 5 | microSD: inizializzazione SPI a 400 kHz, mount (exFAT/FAT32), lettura o creazione di `config.txt`, salvataggio in flash se nome o soglia sono cambiati, apertura del nuovo file `<nome>_<M|C>_<NNN>.txt`, spazio libero | 0.1–1 s | senza scheda si va avanti con il solo USB |
| 6 | Buzzer configurato in DC o in PWM a seconda di `BUZZER_ACTIVE` | | |
| 7 | BMP280: reset, calibrazione, modo normal; prima lettura attraverso il filtro. MPU-6050: reset e configurazione | ~170 ms | |
| 8 | Schermata **SENSOR CHECK** se un sensore abilitato non risponde | 4 s | solo in caso di problemi |
| 9 | ADC e **calibrazione della soglia**: scansione del duty, fit lineare, impostazione, trim, verifica entro il 20 % | ~0.4 s | §6 |
| 10 | Controllo dell'HV (ADC2) entro il 20 % di 30.2 V | | in caso di errore, messaggio su seriale e OLED |
| 11 | Tensione del segnale, **baseline** (256 cicli di reset, 50 µs di attesa, conversione: media ± deviazione standard), controllo che il trigger non sia bloccato alto | ~15 ms | la baseline iniziale serve alla conversione in mV (§7) |
| 12 | Avvio del core 1 e del PIO; stampa della finestra di coincidenza; header | | |
| 13 | Scrittura su SD di tutto il log di boot e dell'header; header su USB | | |
| 14 | `run_start`, primo riarmo del PIO, riga `# HK … src=start`, watchdog a 8 s | | da qui partono i timestamp |

Tutte le righe di diagnostica iniziano con `#` e il loro formato riprende quello
del firmware originale (`# (1234ms) ...`). Vengono salvate in un buffer di 8 kB,
scritte all'inizio del file su SD e reinviate su USB ogni volta che un programma
apre la porta.

---

## 3. Acquisizione di un evento

### 3.1 Programma PIO (`daq.pio`)

```
pull block          ; aspetta il riarmo dalla CPU
wait 0 pin 0        ; il comparatore deve essere basso
armed:
wait 1 pin 0        ; fronte del trigger (GPIO22)
set pindirs, 1      ; abbassa la propria linea di coincidenza
mov isr, ~null
push                ; marker 0xFFFFFFFF -> la CPU prende il timestamp
set y, 0
set x, 31
sample:             ; 32 campioni ogni 12 cicli = 384 cicli = 3.07 us
jmp pin, high [5]   ; linea del partner alta?
set y, 1            ; no: è bassa, quindi coincidenza
high:
jmp x--, sample [5]
set pins, 1         ; spinge brevemente la linea alta (fronte di salita veloce)
set pindirs, 0      ; la rilascia ai pull-up
set pins, 0
mov isr, y
push                ; flag di coincidenza 0/1
```

In modo singolo (M) la linea di uscita non è collegata al PIO: le istruzioni
`set` non hanno effetto sui pin e la CPU ignora il flag.

### 3.2 Sequenza temporale (core 1, `read_event`)

| t (dal fronte del trigger) | azione |
|---|---|
| 0 | il comparatore U5B passa alto: l'impulso amplificato supera la soglia |
| ~30 ns | il PIO abbassa la propria linea di coincidenza (2 cicli di sincronizzazione dell'ingresso più 2 istruzioni) |
| ~50 ns – 0.3 µs | marker nella FIFO; il core 1, che la controlla continuamente, prende **t0 = timestamp** con risoluzione di 1 µs |
| 3.1 µs | fine della finestra: il PIO rilascia la linea e mette il flag nella FIFO |
| 3.1 → 5.1 µs | **una conversione ADC0** (2 µs) del picco tenuto su C5, letto attraverso il buffer U5A e il partitore 1.3/2.3 |
| 5.1 → 10.1 µs | **reset del peak detector**: GPIO21 alto per 5 µs (Q2 scarica C5 attraverso 1 kΩ, τ ≈ 0.2 µs) |
| | attesa che il comparatore sia tornato basso (di solito lo è già; con un timeout di 1 ms aumenta il contatore `stuck`) |
| ~11–12 µs | GPIO21 basso, 1 µs di assestamento, **riarmo** del PIO = t1 |
| | `dead += t1 − t0`; accensione dei LED (scrittura di un registro PWM); evento nella coda |

Il core 0 prende poi l'evento dalla coda (fino a 64 per giro del ciclo
principale), calcola `SiPM[mV]`, aggiunge gli ultimi valori dei sensori,
formatta la riga e la scrive su SD e USB.

**Eventi persi:**
- **coda piena** (core 0 in ritardo di più di 2048 eventi): l'evento viene
  contato in `lost` (primo campo) e non scritto. A 3 Hz servirebbero più di
  10 minuti di blocco del core 0, che il watchdog interromperebbe prima;
- **trigger durante il tempo morto:** l'evento non viene visto, e il tempo è già
  contato in `Deadtime[s]`.

---

## 4. Tempo morto

### 4.1 Cosa viene contato

`Deadtime[s]` è la somma cumulativa, dall'inizio del run, di tre contributi:

| contributo | durata | frequenza | frazione del tempo a 3 Hz |
|---|---|---|---|
| lettura di un evento, `t1 − t0` | ≈ 12 µs (3.1 finestra + 2 ADC + 5 reset + 2 di assestamento e overhead) | ogni evento | 3.6·10⁻⁵ |
| campione di baseline (§8.3) | ≈ 63 µs (5 + 50 + 2 + 5 + 1) | ogni 5 s | 1.3·10⁻⁵ |
| pausa del core 1 durante una scrittura in flash | ~50 ms | solo al comando `name`/`threshold` o se cambia la config al boot | — |

Totale a 3 Hz: circa 5·10⁻⁵ (0.005 %), contro circa 0.12 % del firmware
originale (≈ 400 µs per evento).

Il live time si calcola come `live = Timestamp − Deadtime`, e il rate come
`R = N / live`. È lo stesso calcolo che fanno l'OLED e la GUI.

### 4.2 Cosa non viene contato (limiti noti)

- **Latenza fra il fronte del trigger e t0**, cioè fra il trigger e il momento in
  cui il core 1 legge il marker. Di solito è sotto 0.3 µs. Se però il core 1 sta
  facendo una conversione di housekeeping (§8.2) arriva a ~2 µs, e in quel caso
  si spostano anche il timestamp e l'istante di lettura dell'ADC. Le conversioni
  di housekeeping occupano ~134 × 2 µs ≈ 0.27 ms al secondo, quindi riguardano
  ~3·10⁻⁴ degli eventi. Il flag di coincidenza resta comunque corretto, perché lo
  determina il PIO.
- **Il periodo fra l'avvio del core 1 e `run_start`** (pochi ms), perché i
  timestamp partono da `run_start`.

### 4.3 Perché è così basso rispetto all'originale

Nell'originale la stampa su seriale e la gestione dei LED stanno nel percorso
critico dell'evento. Qui il core 1 fa solo le operazioni analogiche
indispensabili; formattazione, USB, SD e display sono sul core 0 e non
bloccano il trigger. Il passo che limita è la finestra di coincidenza (3.1 µs),
fissata per compatibilità con il manuale.

---

## 5. Coincidenza fra due detector

### 5.1 Linee fisiche
- **Collegamento:** GPIO0 (linea A) e GPIO1 (linea B) sono collegati ai pin 4 e 3
  del connettore RJ45, con un cavo dritto.
- **Pull-up:** interno su entrambe le estremità. A riposo le linee sono alte.
- **Pilotaggio:** una linea viene solo tirata bassa o rilasciata. Il ritorno
  alto è spinto attivamente per un ciclo, poi la linea viene rilasciata. Due
  uscite non possono quindi trovarsi a livelli opposti.

### 5.2 Handshake al boot (`coinc.c`)

1. Entrambe le linee sono configurate come ingressi con pull-up.
2. **Ascolto** per un tempo casuale fra 50 e 300 ms (numero casuale dell'RP2040).
   Se in questo intervallo la linea A resta bassa per almeno 5 ms, il partner si
   sta annunciando. Il detector allora tira bassa la linea B, aspetta che A torni
   alta (al massimo 3 s), attende 10 ms e rilascia B. Diventa **follower**:
   segnala gli eventi su B e ascolta A.
3. Altrimenti **si annuncia**: tiene bassa la linea A per al massimo 2.5 s e
   aspetta che B resti bassa per almeno 5 ms. Se succede, rilascia A, aspetta che
   B torni alta e diventa **leader**: segnala gli eventi su A e ascolta B.
4. Nessuna risposta: **modo M** (singolo).

Conseguenze:
- i due detector devono avviarsi entro circa 2 s l'uno dall'altro;
- quando il partner viene trovato, entrambi i detector accendono i due LED fino
  all'avvio dell'acquisizione;
- se un detector si riavvia durante il run, non ritrova il partner e passa in
  modo M. Vanno riavviati entrambi.

Il protocollo non è compatibile con il firmware originale: entrambi i detector
devono avere questo firmware.

### 5.3 Finestra
Ogni detector tiene bassa la propria linea per tutta la durata della propria
finestra di campionamento, 3.1 µs. Due eventi risultano coincidenti se
|t_A − t_B| ≤ 3.01 µs. La simulazione del programma PIO (emulatore pioemu) dà
una finestra simmetrica e contigua da −376 a +376 cicli.

Il rate di coincidenze accidentali atteso è `2 · W · R_A · R_B`, con W = 3.0 µs:
con 3 Hz su entrambi i detector fa 5.4·10⁻⁵ Hz.

---

## 6. Soglia di trigger

### 6.1 Circuito
- **Generazione:** il PWM su GPIO20 passa per il partitore R21/R19 (1 kΩ / 1 kΩ)
  e il condensatore C11 (1 µF), con costante di tempo di 0.5 ms, e arriva
  all'ingresso invertente del comparatore U5B. Il valore va da 0 a 1.65 V.
- **Rilettura:** l'ADC1 legge lo stesso nodo (TriggerThreshold).

### 6.2 PWM con dithering
- **Il problema:** il PWM è a 1 MHz come nell'originale, quindi a 125 MHz il
  contatore ha 125 passi. Ogni passo vale 3.3 V/2/125 = 13.2 mV.
- **La soluzione:** un canale DMA, sincronizzato con il ciclo del PWM, scrive nel
  registro di compare una tabella di 256 valori. La tabella alterna i livelli L e
  L+1 secondo una sequenza sigma-delta del primo ordine, e l'RC ne fa la media.
- **Risultato:** risoluzione di 13.2/256 ≈ 0.05 mV. Il ripple residuo è dominato
  dal PWM a 1 MHz ed è circa 0.16 mV picco-picco a 80 mV.

### 6.3 Calibrazione al boot
1. **Scansione:** livelli interi 1, 2, … fino a quando la soglia supera 1.2 V.
   Per ciascuno si attendono 4 ms (8τ) e si fa una media di 64 conversioni
   dell'ADC1.
2. **Fit:** `V_thr[mV] = a · duty + b`. È accettato se 1000 < a < 2500 mV (atteso
   ≈ 1650) e R² > 0.99; altrimenti si usano i valori nominali (a = 1650, b = 0).
3. **Impostazione e trim:** il valore richiesto viene impostato, misurato e
   corretto una volta. Poi viene verificato entro il 20 %.

### 6.4 Durante il run
- **Lettura:** ogni 1 s la soglia viene misurata (64 conversioni) e mostrata
  sull'OLED.
- **Trim in anello chiuso:** ogni 10 s, se |misurata − impostata| > 0.3 mV, il
  duty viene corretto. Compensa le derive dei 3.3 V che alimentano il PWM. L'ADC
  è riferito all'LM4040, più stabile. Il trim viene saltato nel giro in cui la
  correzione in temperatura ha appena spostato la soglia.

---

## 7. Da ADC a SiPM[mV] e correzione in temperatura

### 7.1 Conversione

```
lin(ADC) = ADC                                   per ADC <= 3490
lin(ADC) = interpolazione lineare della tabella  per 3490 < ADC < 3818
lin(ADC) = 19430                                 per ADC >= 3818

SiPM[mV] = k · ( b + (lin(ADC) − b) · c )
```

- **k = `SiPM_mV_per_LSB`**, 0.0706 di default. È ricavato dai file prodotti dal
  firmware originale sui detector v3X.
- **La tabella** (20 nodi, `cwcore.c`) riproduce la curva di saturazione del
  firmware originale, ricostruita dai dati. Un evento saturo vale ≈ 1370 mV. Si
  disattiva con `SIPM_SATURATION_LUT = false`, e in quel caso lin(ADC) = ADC.
- **b** è la baseline in LSB: all'inizio quella misurata al boot, poi la media
  dei campioni di ogni periodo HK, se sono almeno 5.
- **c** è il fattore di correzione in temperatura (1 se la correzione è spenta).
  Con c = 1 la formula si riduce a k · lin(ADC), identica all'originale.

La colonna `ADC[12b]` contiene sempre il valore grezzo.

### 7.2 Correzione in temperatura (`USE_REALTIME_GAIN_CORRECTION = true`)

Ogni 10 s, se la temperatura filtrata (§8.4) ha meno di 30 s, viene calcolato c
con il modello scelto:

```
ov:         c = (HV − Vbr(Tref)) / (HV − Vbr(T))
            Vbr(T) = 24.45 V + 0.0215 V/°C · (T − 21 °C)        (datasheet MicroFC-60035)
empirical:  c = exp( −GAIN_TEMP_COEFF/100 · (T − Tref) )        (GAIN_TEMP_COEFF in %/°C, negativo)
```

Il valore viene accettato se 0.8 < c < 1.25 e cambia di più di 10⁻⁴. In quel caso:
- **ampiezza:** `SiPM[mV]` di tutti gli eventi successivi usa il nuovo c;
- **soglia:** `V_thr' = 24.75 mV + (trigger_voltage_mV − 24.75 mV) / c`, dove
  24.75 mV è il bias dell'amplificatore. La soglia resta così costante in unità
  di ampiezza del segnale;
- **registrazione:** se c è cambiato di più dello 0.2 % dall'ultima volta, viene
  scritta una riga HK con `src=corr`.

Nel modello OV un errore assoluto sull'HV pesa solo al secondo ordine, perché
compare sia al numeratore sia al denominatore.

---

## 8. Sensori e housekeeping: cosa viene letto e quando

### 8.1 Riepilogo dei periodi (ciclo principale del core 0)

| cosa | periodo | dove | note |
|---|---|---|---|
| MPU-6050 (accelerometro, giroscopio) | 100 ms | core 0, I²C | 14 byte; sensore campionato a 100 Hz con filtro passa-basso a 44 Hz, ±2 g, ±250 °/s |
| BMP280 (T, P) | 1 s | core 0, I²C | il sensore misura ogni 0.5 s con filtro IIR ×16; il firmware applica anche la mediana su 5 letture |
| HV (ADC2) | 1 s | richiesta al core 1 | media di 64 conversioni |
| Soglia (ADC1) | 1 s | richiesta al core 1 | media di 64 conversioni |
| Baseline (ADC0, PIO fermo) | 5 s | core 1 | 1 conversione, ~63 µs di tempo morto |
| Correzione in temperatura e trim della soglia | 10 s | core 0 | il trim fa altre 64 conversioni ADC1 |
| Riga HK | `HK_PERIOD_S`, 60 s di default | core 0 | anche all'avvio e quando c cambia |
| OLED | 1 s | core 0, I²C | circa 25 ms per il frame completo |
| Scrittura del buffer SD | 1 s, oppure con 8 kB in attesa | core 0, SPI | `f_sync` ogni 10 s |
| Watchdog | 0.5 s | core 0 | aggiornato solo se anche il core 1 è vivo |
| LED, buzzer | 20 ms dopo l'evento | core 0 | spegnimento |

**I sensori non vengono letti a ogni evento.** Ogni riga contiene l'ultimo valore
disponibile: al massimo 100 ms di ritardo per accelerometro e giroscopio, 1 s per
T e P. Le letture I²C non toccano mai il core 1.

### 8.2 Conversioni ADC di housekeeping
Il core 0 chiede "N conversioni sul canale X". Il core 1 le esegue **una alla
volta** (~2 µs ciascuna), solo quando la FIFO del PIO è vuota, e dopo ognuna
riporta il multiplexer sul canale 0. In totale sono circa 134 conversioni al
secondo. L'effetto sui tempi è descritto in §4.2.

### 8.3 Baseline durante il run
Ogni 5 s il core 0 chiede un campione di baseline. Il core 1:
1. ferma la macchina a stati del PIO; se non è ferma sull'istruzione `armed`
   (`wait 1 pin`) o se c'è un evento nella FIFO, la riavvia e riprova più tardi;
2. scarica C5 per 5 µs, attende 50 µs come nella misura al boot, fa una
   conversione;
3. se nel frattempo il comparatore è scattato, aspetta che torni basso e
   scarica di nuovo C5, per non lasciare un picco parziale;
4. riavvia il PIO e aggiunge la durata al tempo morto.

Media e deviazione standard dei campioni compaiono nella riga HK (`base`,
`base_sd`, `base_n`). Se i campioni sono almeno 5, la media diventa la nuova b
della conversione (§7.1).

### 8.4 Filtro su temperatura e pressione (`cw_envfilter_*`)
- **Letture non fisiche:** scartate se T è fuori da (−40, 85) °C o P è fuori da
  (1, 115) kPa. Il pattern T = −5 °C, P = 0, che il firmware originale scrive
  quando la lettura fallisce, ricade qui.
- **Salti:** scartati se |T − ultima lettura accettata| supera 2 °C/s · max(Δt, 1 s)
  più 0.5 °C.
- **Valori in uscita:** mediana delle ultime 5 letture accettate, separatamente
  per T e P.
- **Ripartenza:** dopo 5 scarti consecutivi per salto il filtro riparte dalla
  nuova lettura, per accettare un cambiamento persistente.
- **Contatori:** le letture fallite sul bus (`bmp_fail`) e quelle scartate
  (`bmp_rej`) compaiono nella riga HK e con il comando `status`.
- **Uso:** i valori filtrati vanno nelle colonne `Temp[C]` e `Press[Pa]`,
  sull'OLED e nella correzione.

---

## 9. Uscita dei dati

### 9.1 Riga di evento
Colonne separate da tab, terminatore `\n` su SD e `\r\n` su USB:

```
Event  Timestamp[s]  Flag  ADC[12b]  SiPM[mV]  Deadtime[s]  [Temp[C]  Press[Pa]]  [Accel(X:Y:Z)[g]  Gyro(X:Y:Z)[deg/sec]]  [Name]
```

| colonna | formato | note |
|---|---|---|
| Event | intero | da 1 |
| Timestamp[s] | `%.6f`, calcolato con interi | dall'inizio del run (`run_start`), risoluzione 1 µs, nessun arrotondamento in virgola mobile |
| Flag | 0/1 | sempre 0 in modo M |
| ADC[12b] | 0–4095 | grezzo |
| SiPM[mV] | `%.1f` | §7 |
| Deadtime[s] | `%.6f` | cumulativo (§4) |
| Temp[C], Press[Pa] | `%.1f`, `%.0f` | solo se il BMP280 è abilitato e risponde al boot |
| Accel, Gyro | `%.3f:%.3f:%.3f`, `%.1f:%.1f:%.1f` | solo se l'MPU-6050 è abilitato e risponde al boot |
| Name | testo | **solo su USB** (lo legge `import_data.py`) |

Le colonne dei sensori vengono decise al boot e non cambiano durante il run: se
un sensore smette di rispondere, rimane l'ultimo valore valido.

### 9.2 Intestazione
È compatibile con l'originale. La larghezza delle righe di `#` segue quella
della riga delle colonne:

```
#######################################################################################################################
#                                     CosmicWatch: The Desktop Muon Detector v3X
#                                          Firmware: CosmicWatch-open 1.0.0
#                                                 Detector Name: TOP
# Event  Timestamp[s]  Flag  ADC[12b]  SiPM[mV]  Deadtime[s]  Temp[C]  Press[Pa]  Accel(X:Y:Z)[g]  Gyro(X:Y:Z)[deg/sec]
#######################################################################################################################
```

### 9.3 Riga HK (`HK_PERIOD_S`)

```
# HK t=60.001 T=23.41 P=100512 HV=30.120 Thr_set=80.00 Thr_meas=80.05 base=55.30 base_sd=1.20 base_n=12 corr=1.00000 events=171 coinc=24 dead=0.005321 lost=0/0/0 bmp_fail=0 bmp_rej=0 src=periodic
```

| campo | significato |
|---|---|
| t | secondi dall'inizio del run |
| T, P | valori filtrati, oppure `NA` |
| HV | V, dall'ADC2 |
| Thr_set, Thr_meas | soglia impostata (eventualmente corretta) e misurata, in mV |
| base, base_sd, base_n | statistica dei campioni di baseline dalla riga precedente, in LSB |
| corr | fattore c in uso |
| events, coinc, dead | contatori cumulativi |
| lost | eventi persi per coda piena / righe USB perse / byte SD persi |
| bmp_fail, bmp_rej | letture BMP280 fallite / scartate dal filtro |
| src | `start`, `periodic` oppure `corr` |

### 9.4 microSD (`logger.c`)
- **Nome del file:** `<detector_name>_<M|C>_<NNN>.txt`. NNN è il massimo dei
  numeri già presenti per quel nome più 1, indipendentemente dal modo.
- **Contenuto:** log di boot, header, eventi, righe HK e messaggi.
- **Buffer:** 32 kB in RAM, scritti con `f_write` ogni secondo o quando superano
  8 kB, con `f_sync` ogni 10 s. In caso di mancanza di alimentazione si perdono
  al massimo circa 10 s di dati.
- **Buffer pieno:** succede se la scheda resta bloccata più a lungo di quanto il
  buffer riesca a coprire, cioè circa 0.45 s a 700 Hz (~70 kB/s) o diversi minuti
  a 3 Hz. Le righe in eccesso vengono scartate e contate in `lost` (terzo campo).
- **Errore di scrittura:** un tentativo di ripresa; se fallisce, il logging su SD
  si ferma con un avviso e l'USB continua.
- **Scheda e filesystem:** SPI1 a 12.5 MHz; FatFs R0.16 con exFAT, FAT32 e nomi
  lunghi. Data dei file fissa (1/1/2026): non c'è un orologio reale.

### 9.5 USB
- **Nome:** CDC seriale, che si presenta come "CosmicWatch v3X". Il baud rate
  non conta.
- **Apertura della porta:** quando un programma apre la porta (DTR), il
  firmware reinvia il log di boot e l'header con la colonna `Name`.
- **Host che non legge:** se una scrittura non riesce entro 3 ms (buffer USB da
  1 kB pieno), l'uscita USB viene sospesa per 1 s e le righe saltate sono
  contate in `lost` (secondo campo). La SD non ne risente.
- **Nessun host:** non viene scritto nulla e non si perde tempo.

---

## 10. config.txt

Il file si chiama **`config.txt`** (non `config.cfg`) e sta nella radice della
microSD, come nel firmware originale.

### 10.1 Da dove vengono i parametri

1. Valori di default (tabella sotto).
2. **Flash:** nome e soglia, se presenti; altrimenti il nome è `CW-xxxx`, dalle
   ultime 4 cifre esadecimali dell'ID unico della scheda.
3. **`config.txt`**, se c'è la microSD; ha la precedenza.
4. Se `config.txt` esiste e nome o soglia sono diversi da quelli in flash,
   vengono salvati in flash. Così il detector mantiene la propria identità
   anche senza scheda, e `import_data.py` non trova nomi duplicati.

Se `config.txt` manca, viene creato con i valori correnti, commentati. Un
`config.txt` già presente **non viene mai riscritto**: le chiavi nuove vanno
aggiunte a mano, e quelle assenti valgono come default.

### 10.2 Sintassi
- una riga `chiave = valore` per volta; `#` inizia un commento;
- le chiavi e i valori booleani non distinguono maiuscole e minuscole
  (`true/false`, `1/0`, `yes/no`, `on/off`);
- nei valori numerici sono ammesse unità di misura in coda (`80 mV`, `10%`);
- **chiave sconosciuta:** viene segnalata nel log e ignorata;
- **valore non valido:** viene segnalato e resta il valore precedente;
- **valore fuori intervallo:** viene segnalato e si usa il default;
- `detector_version`, scritta dal firmware originale, viene ignorata.

### 10.3 Chiavi

| chiave | tipo | default | intervallo / valori | effetto |
|---|---|---|---|---|
| `detector_name` | testo | `CW-xxxx` | 1–20 caratteri; tutto ciò che non è lettera, cifra o `-` diventa `-` | nome nei file, nell'intestazione e nell'ultima colonna USB |
| `trigger_voltage_mV` | mV | 80 | 30–1200 | soglia al comparatore; la baseline del segnale è circa 25 mV |
| `LED_BRIGHTNESS` | % | 10 | 0–100 | intensità dei LED sugli eventi |
| `USE_LEDs` | bool | true | | LED sugli eventi |
| `USE_OLED` | bool | true | | display (con false viene spento) |
| `USE_SERIAL` | bool | true | | righe di evento su USB (i commenti vengono inviati comunque) |
| `USE_BMP280` | bool | true | | colonne Temp/Press e correzione |
| `USE_MPU6050` | bool | true | | colonne Accel/Gyro |
| `USE_BUZZER` | bool | false | | buzzer |
| `BUZZER_ACTIVE` | bool | true | | true: GPIO6 in DC (buzzer attivo); false: onda quadra a 2.7 kHz (passivo) |
| `BUZZER_EVENTS` | testo | coincidence | `coincidence` / `all` | quando suona |
| `USE_REALTIME_GAIN_CORRECTION` | bool | false | | correzione in temperatura (§7.2) |
| `GAIN_CORRECTION_MODEL` | testo | ov | `ov` / `empirical` | modello |
| `GAIN_TEMP_COEFF` | %/°C | −0.37 | da −5 a 5 | α del modello empirico |
| `GAIN_REF_TEMP_C` | °C | 25 | −40…85 | temperatura di riferimento |
| `HK_PERIOD_S` | s | 60 | 0 (disattivato) oppure 5–3600 | riga HK |
| `SiPM_mV_per_LSB` | mV/LSB | 0.0706 | 0–10 | k della conversione |
| `SIPM_SATURATION_LUT` | bool | true | | curva di saturazione sopra 3490 LSB |

### 10.4 File creato automaticamente (default)

```
# CosmicWatch v3X configuration (CosmicWatch-open 1.0.0)
# Lines are 'key = value'; '#' starts a comment.

# Up to 20 characters: letters, digits and '-'. Used in file names and in
# the last column of the USB data stream.
detector_name = CW-1A2B

# Comparator threshold in mV, measured at the TriggerThreshold node.
# The amplified signal baseline sits at about 25 mV.
trigger_voltage_mV = 80.0

LED_BRIGHTNESS = 10
USE_LEDs = true
USE_OLED = true
USE_SERIAL = true
USE_BMP280 = true
USE_MPU6050 = true
USE_BUZZER = false
# true: buzzer with internal oscillator (sounds with DC), false: passive buzzer
BUZZER_ACTIVE = true
# coincidence: beep only on coincident events; all: beep on every event
BUZZER_EVENTS = coincidence

# Normalise SiPM[mV] to the gain at GAIN_REF_TEMP_C (needs the BMP280) and move
# the threshold accordingly. Models: ov = SiPM overvoltage from datasheet Vbr(T)
# and measured HV; empirical = exp(-GAIN_TEMP_COEFF/100 * (T - GAIN_REF_TEMP_C)),
# GAIN_TEMP_COEFF = d ln(MPV)/dT in %/C from a temperature calibration.
USE_REALTIME_GAIN_CORRECTION = false
GAIN_CORRECTION_MODEL = ov
GAIN_TEMP_COEFF = -0.370
GAIN_REF_TEMP_C = 25.0

# Housekeeping comment line (T, P, HV, threshold, baseline, correction) every
# HK_PERIOD_S seconds in the data stream; 0 = off.
HK_PERIOD_S = 60

# SiPM[mV] = SiPM_mV_per_LSB * ADC. With SIPM_SATURATION_LUT the values above
# ~3490 ADC follow the saturation curve of the original firmware.
SiPM_mV_per_LSB = 0.0706
SIPM_SATURATION_LUT = true
```

Un `config.txt` scritto dal firmware originale viene letto senza problemi:
contiene solo un sottoinsieme di queste chiavi.

---

## 11. Comandi da seriale

Si inviano come righe di testo terminate da `\r` o `\n`. Le risposte iniziano
con `#`, vanno su USB e, durante il run, anche nel file su SD.

| comando | effetto |
|---|---|
| `help` | elenco dei comandi |
| `status` | modo e ruolo, file, tempo di run e live time, conteggi, soglia impostata e misurata, HV, c, baseline, stato del BMP280, eventi persi, trigger bloccati |
| `name <NOME>` | salva il nome in flash; vale dal boot successivo, ma `config.txt` sulla SD ha la precedenza |
| `threshold <mV>` | imposta subito la soglia (30–1200 mV) e la salva in flash; il core 1 si ferma per circa 50 ms |
| `reboot` | riavvio |
| `bootsel` | riavvio nel bootloader USB, per caricare un nuovo `.uf2` senza premere BOOTSEL |

---

## 12. OLED, LED e buzzer

**OLED** (SSD1306 128×64, 8 righe da 21 caratteri), aggiornato ogni secondo:

```
TOP        C 01:23:45      nome, modo M/C, tempo di run
TOP_C_005.txt              file, oppure "No microSD logging"
Events: 14203
2.851+/-0.024Hz            rate sul live time, errore sqrt(N)/live
Coinc.: 1987               solo in modo C
0.399+/-0.009Hz            solo in modo C
Thr 80.0mV HV 30.1V        soglia misurata e HV
23.4°C 1013.2hPa           oppure "Dead x.xxx%" senza BMP280
```

Ogni minuto il testo si sposta di 1 pixel (0–2 px) per ridurre il burn-in.
Schermate speciali: logo al boot, SENSOR CHECK (§2), `Threshold Error!`,
`HV Error!`.

**LED:**

| LED | quando si accende |
|---|---|
| 5 mm (GPIO12) | a ogni evento |
| 3 mm (GPIO13) | solo sugli eventi coincidenti |

- Accensione: dal core 1, subito dopo il riarmo, con un PWM a 1 kHz.
- Spegnimento: dal core 0 dopo 20 ms.
- Intensità: `LED_BRIGHTNESS`.
- Al boot, con partner trovato, entrambi restano accesi al 100 % fino all'avvio
  dell'acquisizione.

**Buzzer:** 20 ms per evento.
- Pilotaggio: in DC oppure a 2.7 kHz, secondo `BUZZER_ACTIVE`.
- Eventi: con `BUZZER_EVENTS = coincidence` solo quelli coincidenti (in modo M
  resta muto); con `all` tutti.

---

## 13. Robustezza e gestione degli errori

| situazione | comportamento |
|---|---|
| Core 0 o core 1 bloccati | watchdog a 8 s. Il core 0 lo aggiorna solo se il contatore di attività del core 1 avanza, poi riavvia; al boot successivo compare `WARNING: the previous run ended with a watchdog reset` |
| Comparatore bloccato alto (soglia sotto la baseline) | avviso al boot. Durante il run il PIO non si riarma finché il comparatore non torna basso; ogni timeout di 1 ms incrementa `stuck triggers` (`status`) |
| HV o soglia fuori del 20 % al boot | messaggio su seriale e OLED con i componenti da controllare; l'acquisizione parte comunque |
| microSD assente o non risponde | solo USB; l'OLED mostra "No microSD logging" |
| microSD rimossa o guasta durante il run | un tentativo di ripresa, poi il logging su SD si ferma con un avviso; USB e OLED continuano |
| Host USB che non legge | uscita USB sospesa a intervalli di 1 s, righe contate; SD e acquisizione non ne risentono |
| BMP280 che fallisce | resta l'ultimo valore valido; `bmp_fail` sale; correzione sospesa dopo 30 s senza letture buone |
| Scrittura in flash | core 1 in pausa nella RAM, interrupt del core 0 disabilitati per ~50 ms; la pausa viene contata come tempo morto |
| Coda degli eventi piena | evento scartato e contato in `lost` |

---

## 14. Risorse dell'RP2040

| risorsa | uso |
|---|---|
| Flash | 116 kB di programma; l'ultimo settore da 4 kB contiene nome e soglia |
| RAM | 108 kB di dati statici: coda di eventi 48 kB, buffer SD 32 kB, log di boot 8 kB, framebuffer OLED 1 kB, FatFs pochi kB |
| Core | 0: I/O e controllo; 1: acquisizione |
| PIO | PIO0, macchina a stati 0, 16 istruzioni su 32 |
| DMA | 1 canale (soglia), interrupt DMA_IRQ_0 sul core 0 |
| PWM | slice 2A (soglia, GPIO20), 3A (buzzer, GPIO6, se passivo), 6A/6B (LED, GPIO12/13) |
| ADC | canale 0 (picco), 1 (soglia), 2 (HV); riferimento esterno 2.5 V (LM4040) |
| SPI1 | microSD (GPIO8–11) |
| I²C1 | OLED, BMP280, MPU-6050 (GPIO14/15), 400 kHz |
| GPIO liberi | 2–5 (SPI0 sull'RJ45), 7, 16, 17, 18 |

---

## 15. Differenze rispetto al firmware originale

| aspetto | originale (1.1.52) | questo firmware |
|---|---|---|
| Tempo morto per evento | ≈ 400 µs | ≈ 12 µs |
| Trigger e coincidenza | software, 32 letture × 72 ns | PIO hardware, finestra simmetrica ±3.0 µs |
| Protocollo di coincidenza | non documentato | §5; i due firmware non sono compatibili fra loro |
| Risoluzione della soglia | 13.2 mV (PWM a 125 passi) | ~0.05 mV con dithering via DMA, trim in anello chiuso |
| BMP280 non letto | scrive T = −5, P = 0 | ultimo valore valido, filtro, contatori |
| Correzione in temperatura | presente, modello non documentato | modello OV o empirico, configurabile e registrato nelle righe HK |
| Baseline | solo al boot | anche durante il run, ogni 5 s |
| Housekeeping nei dati | assente | riga `# HK` periodica |
| Nome senza microSD | default comune | salvato in flash oppure `CW-xxxx` unico |
| Comandi da seriale | nessuno | `status`, `name`, `threshold`, `reboot`, `bootsel` |
| Formato dei dati, `config.txt`, nomi dei file | — | compatibili (solo chiavi aggiunte) |
| EEPROM CAT24C64 (scheda v3C), melodie del buzzer | presenti | non implementate |
