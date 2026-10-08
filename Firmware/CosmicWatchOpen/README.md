# Firmware open source per CosmicWatch v3X (Raspberry Pi Pico / RP2040)

Firmware scritto da zero a partire dallo schema elettrico (manuale, Fig. B.2, e
`PCB/Circuit_Diagram.pdf`) e dalla descrizione del funzionamento nei documenti del
progetto. Sostituisce `CosmicWatch_v3X.1.1.52.uf2` mantenendo lo stesso formato
dei dati: `import_data.py`, la GUI e gli script in `Data/` funzionano senza modifiche.

La logica del firmware è descritta in dettaglio in [FIRMWARE.md](FIRMWARE.md); il piano
di verifica e calibrazione è in [PIANO_MISURE.md](PIANO_MISURE.md).

> **Stato:** compila senza warning e la logica è verificata con test su host e
> simulazione del PIO, ma **non è ancora stato provato su un detector reale**.
> Prima di usarlo per una misura, fare le verifiche della sezione
> [Prima accensione](#prima-accensione). Per tornare al firmware originale basta
> ricaricare `../CosmicWatch_v3X.1.1.52.uf2`.

## Caricare il firmware

1. Tenere premuto BOOTSEL sul Pico mentre si collega l'USB (oppure, con questo
   firmware già caricato, inviare il comando `bootsel` sulla seriale).
2. Copiare `build/cosmicwatch_open.uf2` nel disco `RPI-RP2` che compare.

## Compilare

```bash
./build.sh            # produce build/cosmicwatch_open.uf2
tests/run_tests.sh    # test su host (logica pura + formato dati)
```

`build.sh` usa solo strumenti dentro `venv_cosmicwatch.nosync`:

| strumento | dove | come è stato installato |
|---|---|---|
| cmake, ninja | `bin/` | `pip install cmake ninja` |
| Arm GNU Toolchain 14.3 | `opt/arm-gnu-toolchain-14.3.rel1-darwin-arm64-arm-none-eabi` | tarball da developer.arm.com |
| pico-sdk 2.3.1 + tinyusb | `opt/pico-sdk` | `git clone -b 2.3.1 --depth 1` + submodule `lib/tinyusb` |

Al primo build il pico-sdk compila due tool per il Mac (pioasm, picotool, senza
supporto USB). Su questo Mac i Command Line Tools hanno una cartella
`usr/include/c++/v1` incompleta: lo script lo rileva e punta il compilatore host
agli header dell'SDK di macOS.

FatFs R0.16 (ChaN, licenza BSD a 1 clausola) è incluso in `lib/fatfs` con exFAT
e nomi lunghi abilitati (`ffconf.h`).

## Come funziona

**Core 1: solo acquisizione.** Gira da RAM con gli interrupt disabilitati e non
accede mai alla flash (verificato sul disassemblato: nessuna chiamata a funzione,
solo accessi a RAM e periferiche). In questo modo USB, SD e OLED sul core 0 non
introducono latenza. Per ogni evento:

1. Un programma PIO (`src/daq.pio`) aspetta il fronte del comparatore (GPIO22).
   Entro 16 ns abbassa la propria linea di coincidenza e per 3.1 µs campiona
   quella del partner. Coincidenza vuol dire |Δt| ≤ 3.0 µs (finestra simmetrica,
   verificata con l'emulatore PIO).
2. La CPU prende il timestamp (µs) e legge l'ADC0, cioè il picco tenuto su C5
   attraverso R24/R27, circa 3 µs dopo il trigger come nel firmware originale.
3. Porta alto GPIO21 (Q2 scarica C5) per 5 µs, aspetta che il comparatore torni
   basso e riarma il PIO.

Il tempo morto per evento è di circa 12 µs, contro i circa 400 µs dell'originale.

**Baseline durante il run.** Ogni 5 s il core 1 misura la baseline
del peak detector. Lo fa solo se il PIO è armato e in attesa, e mai durante un
evento: ferma la macchina a stati, scarica C5, attende 50 µs, fa una conversione
e riparte. Il costo è di circa 63 µs di tempo morto per campione (0.0013 % del tempo), incluso
nella colonna `Deadtime[s]`.
Viene misurato evento per evento e accumulato nella colonna `Deadtime[s]`.

**Core 0: tutto il resto.** Formatta gli eventi presi da una coda di 2048 eventi
e li scrive su USB e microSD (buffer da 32 kB, scrittura ogni secondo, `f_sync`
ogni 10 s). Gestisce poi OLED, BMP280, MPU-6050, LED, buzzer, i comandi da
seriale e il watchdog (8 s, alimentato solo se anche il core 1 è vivo).

**Soglia.** È data dal PWM a 1 MHz su GPIO20 filtrato da R21/R19/C11, come
nell'originale: 0-1.65 V al comparatore. A 125 MHz il PWM ha solo 125 passi
(13 mV ciascuno). Un canale DMA alterna due passi adiacenti secondo una sequenza
sigma-delta di 256 valori, così la risoluzione diventa circa 0.05 mV e il ripple
resta trascurabile. Al boot la soglia viene calibrata rileggendola con l'ADC1
(fit lineare, come l'originale). Durante la presa dati viene ricontrollata ogni
10 s e corretta se si discosta di più di 0.3 mV.

**Coincidenza.** Il cavo RJ45 è dritto, quindi i due detector sono simmetrici.
Al boot (prima di tutto il resto) ciascuno ascolta per un tempo casuale fra 50 e
300 ms. Se trova la linea A tenuta bassa risponde sulla linea B e diventa
"follower"; altrimenti abbassa A per 2.5 s aspettando la risposta e diventa
"leader". Le linee sono solo tirate basse o rilasciate (pull-up interni), quindi
due uscite non possono mai entrare in conflitto. Durante la presa dati il leader
segnala su GPIO0 e ascolta su GPIO1, il follower il contrario. Se il partner
viene trovato, entrambi i LED si accendono al boot.
- I due detector vanno resettati entro circa 2 s l'uno dall'altro (alimentarne
  uno solo, che alimenta l'altro via cavo, va bene).
- **Entrambi devono avere questo firmware**: il protocollo dell'originale non è
  documentato.

## Pin (schema rev. v2.1)

| GPIO | funzione | GPIO | funzione |
|---|---|---|---|
| 0, 1 | linee di coincidenza (RJ45 4, 3) | 19 | card-detect microSD (solo diagnostica) |
| 6 | buzzer (PWM) | 20 | PWM soglia |
| 8-11 | microSD su SPI1 (RX, CS, SCK, TX) | 21 | reset del peak detector (Q2) |
| 12 | LED 5 mm, tutti gli eventi | 22 | uscita comparatore (trigger) |
| 13 | LED 3 mm, eventi coincidenti | 26 / ADC0 | peak detector ×0.565 |
| 14, 15 | I2C1: OLED, BMP280, MPU-6050 | 27 / ADC1 | lettura della soglia |
| 2-5 | SPI0 su RJ45, non usati | 28 / ADC2 | HV ÷ 23.1 |

L'ADC ha come riferimento i 2.5 V dell'LM4040. Il GPIO7 ("TriggerCommand" nella
revisione dello schema stampata nel manuale) non è collegato e non viene usato.

## Formato dei dati

È identico all'originale: colonne separate da tab, commenti con `#`.

```
Event  Timestamp[s]  Flag  ADC[12b]  SiPM[mV]  Deadtime[s]  [Temp[C]  Press[Pa]]  [Accel(X:Y:Z)[g]  Gyro(X:Y:Z)[deg/sec]]
```

- **Colonne dei sensori.** Ci sono solo se il sensore è abilitato e risponde al
  boot, come nell'originale. Su USB si aggiunge in fondo il nome del detector
  (letto da `import_data.py`).
- **Nome e intestazione del file su microSD.** Il file si chiama
  `<nome>_<M|C>_<NNN>.txt`, con il numero che continua da quelli già presenti.
  L'intestazione contiene tutta la diagnostica di boot.
- **Nuova connessione USB.** Quando un programma apre la porta, il firmware
  ri-invia la diagnostica di boot e l'header.
- **Conversione in `SiPM[mV]`.** Vale `SiPM_mV_per_LSB × ADC`, con 0.0706 di
  default (ricavato dai file prodotti dal firmware originale su questi
  detector). Sopra ~3490 LSB si applica la stessa curva di saturazione
  dell'originale, ricostruita dai dati: gli eventi saturi valgono circa 1370 mV.
  Si disattiva con `SIPM_SATURATION_LUT = false`.
- **Tempo morto.** Le differenze rispetto all'originale vengono solo dal tempo
  morto molto più piccolo e dalla riga `Firmware:` nell'intestazione.

## config.txt (microSD)

Le chiavi sono quelle dell'originale: `detector_name`, `trigger_voltage_mV`,
`LED_BRIGHTNESS`, `USE_LEDs`, `USE_OLED`, `USE_SERIAL`, `USE_BMP280`,
`USE_MPU6050`, `USE_BUZZER`, `USE_REALTIME_GAIN_CORRECTION`. In più ci sono
`SiPM_mV_per_LSB`, `SIPM_SATURATION_LUT`, `BUZZER_ACTIVE`, `BUZZER_EVENTS`,
`GAIN_CORRECTION_MODEL`, `GAIN_TEMP_COEFF`, `GAIN_REF_TEMP_C` e `HK_PERIOD_S`. Se `config.txt` manca, viene creato
con valori di default commentati. Un `config.txt` esistente non viene mai
sovrascritto.

- **Nome del detector.** Al massimo 20 caratteri fra lettere, cifre e `-`; gli
  altri caratteri diventano `-`. Il `_` non è ammesso perché separa i campi nel
  nome del file.
- **`USE_REALTIME_GAIN_CORRECTION`** (spenta di default). Ogni 10 s calcola un
  fattore c a partire dalla temperatura filtrata, solo se l'ultima lettura valida
  del BMP280 ha meno di 30 s. Poi:
  - riscala la parte di segnale sopra la baseline b:
    `SiPM[mV] = k · (b + (ADC − b) · c)`;
  - sposta la soglia: `Vthr' = 24.75 mV + (Vthr − 24.75 mV) / c`.

  I modelli disponibili (`GAIN_CORRECTION_MODEL`) sono due:
  - `ov` (default): `c = (HV − Vbr(Tref)) / (HV − Vbr(T))`, con
    `Vbr(T) = 24.45 V + 21.5 mV/°C · (T − 21 °C)` (datasheet MicroFC-60035) e HV
    misurata;
  - `empirical`: `c = exp(−GAIN_TEMP_COEFF/100 · (T − Tref))`, dove
    `GAIN_TEMP_COEFF` è d ln(MPV)/dT in %/°C ottenuto dalla calibrazione
    (`PIANO_MISURE.md`, sezione 1).

  In entrambi `Tref = GAIN_REF_TEMP_C`, 25 °C di default.
- **Filtro su temperatura e pressione.** Si usa la mediana delle ultime 5 letture
  del BMP280, letto una volta al secondo. Vengono scartate le letture non fisiche
  e i salti di T superiori a 2 °C/s rispetto all'ultima lettura accettata; dopo 5
  scarti consecutivi il filtro riparte. I valori filtrati vanno nelle colonne,
  sull'OLED e nella correzione. Il firmware originale, quando la lettura
  fallisce, scrive T = −5 °C e P = 0; questo firmware tiene l'ultimo valore valido
  e conta le letture fallite (`bmp_fail` nella riga HK).
- **`HK_PERIOD_S`** (default 60, 0 = disattivata). Ogni `HK_PERIOD_S` secondi
  scrive una riga di commento nel file e su USB; i parser la ignorano perché
  inizia con `#`:
  ```
  # HK t=… T=… P=… HV=… Thr_set=… Thr_meas=… base=… base_sd=… base_n=… corr=… events=… coinc=… dead=… lost=q/usb/sd bmp_fail=… bmp_rej=… src=…
  ```
  `base`, `base_sd` e `base_n` riassumono i campioni di baseline presi dalla riga
  precedente. `src` vale `start`, `periodic`, oppure `corr` quando il fattore di
  correzione cambia di più dello 0.2 %. Il fattore c applicato a ogni evento si
  ricostruisce da `SiPM[mV]`, ADC e `base`.
- **Buzzer** (attivo solo con `USE_BUZZER = true`):
  - `BUZZER_ACTIVE` decide il tipo di pilotaggio. Con `true` (default) il GPIO6
    viene alzato in DC per 20 ms: è il caso di un buzzer attivo, che suona da solo
    in continua come quello della lista acquisti. Con `false` il GPIO6 genera
    un'onda quadra a 2.7 kHz, per un buzzer passivo.
  - `BUZZER_EVENTS` decide quando suona. Con `coincidence` (default) suona solo
    sugli eventi coincidenti, quindi in modo M resta muto; con `all` suona a ogni
    evento.
- **Memoria in flash.** Nome e soglia vengono salvati anche nell'ultimo settore
  della flash. Senza microSD il detector mantiene il proprio nome; senza alcuna
  configurazione usa `CW-xxxx`, dove xxxx viene dall'ID unico della scheda.
  Questo evita i nomi duplicati che `import_data.py` rifiuta.

## Comandi da seriale

Si scrivono una riga per volta; tutte le risposte iniziano con `#`, quindi i
parser le ignorano.

| comando | effetto |
|---|---|
| `status` | tempo di run e live time, conteggi, soglia impostata e misurata, HV, eventi persi |
| `name TOP` | salva il nome in flash (attivo dal boot successivo; `config.txt` sulla SD ha la precedenza) |
| `threshold 80` | imposta subito la soglia in mV e la salva in flash |
| `reboot`, `bootsel` | riavvio normale o nel bootloader USB |

## Prima accensione

Da fare in quest'ordine, guardando la diagnostica sulla seriale (115200, anche
se su USB CDC la velocità non conta):

1. **HV.** `Checking: HV to SiPM voltage` deve dare circa 30 V; il valore teorico
   è 30.2 V e il file `TOP_C_002.txt` del firmware originale riporta 28.6 V.
2. **Fit della soglia.** La pendenza attesa è circa 1650 mV/duty (3.3 V/2) con
   R² > 0.999; l'originale riportava 1541-1655. La riga
   `Trigger Threshold (Measured, Expected)` deve essere entro qualche decimo di mV.
3. **Baseline dell'ADC0.** Nei dati dell'originale l'ADC minimo è circa 50-60 LSB.
4. **Rate.** Confrontarlo con quello del firmware originale alla stessa soglia; il
   manuale indica circa 2.5 Hz, di cui circa 0.4 Hz di muoni.
5. **Oscilloscopio.** Sul fronte di GPIO22 il GPIO0/1 deve scendere per circa
   3.1 µs, e il reset su GPIO21 deve durare circa 5 µs a partire da circa 5 µs
   dopo il trigger. Con il comando `status` controllare che `stuck triggers` resti 0.
6. **Coincidenza.** Con due detector e il cavo RJ45 devono accendersi entrambi i
   LED al boot e il file deve chiamarsi `_C_`. Sovrapponendo i detector il rate
   di coincidenza deve risultare dell'ordine del rate di muoni attraverso
   entrambi.
7. **Confronto con l'originale.** Spettro ADC e rate con lo stesso detector e la
   stessa soglia. La differenza attesa è solo il tempo morto, circa 0.003 %
   invece di circa 0.1 %.

## Struttura

```
src/main.c        boot, diagnostica, uscita dati, OLED, sensori, comandi (core 0)
src/daq.c/.pio    acquisizione (core 1 + PIO)
src/threshold.c   DAC della soglia (PWM + DMA) e calibrazione
src/coinc.c       handshake di coincidenza
src/cwcore.c      logica pura: conversioni, BMP280, config.txt, formato righe
src/sensors.c     BMP280, MPU-6050      src/oled.c  SSD1306 + font 5x7
src/sdcard.c      microSD in SPI (diskio per FatFs)   src/logger.c  file e buffer
src/flashcfg.c    nome/soglia in flash
tests/            test su host (ASan/UBSan) e verifica con il parser di import_data.py
```
