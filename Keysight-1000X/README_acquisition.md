# Acquisizione delle coincidenze dal Keysight DSOX1202A

`acquire_coincidences.py` arma l'oscilloscopio in single, aspetta il trigger,
legge in binario C1 e C2 della stessa acquisizione e salva un file per evento
nel formato CSV dell'export Keysight (`x-axis,1,2` / `second,Volt,Volt`).
`analyze_pulses.py` legge questi file senza conversioni. Il ciclo si ripete
finché non si interrompe con Ctrl-C, oppure al raggiungimento di `--max-events`
o `--duration-h`.

Usa PyVISA e PyVISA-py del venv `venv_cosmicwatch.nosync`. Per l'USB il venv
contiene anche `pyusb` e `libusb-package` (libreria libusb inclusa nel pacchetto,
nessuna installazione di sistema).

## Connessione

USB (configurazione attuale, verificata il 9 ottobre 2026 con il DSOX1202A
seriale CN61482183): `--usb` cerca da solo il primo oscilloscopio Keysight.

```bash
python acquire_coincidences.py --usb -o run_3m --check
```

La risorsa trovata è `USB0::10893::903::CN61482183::0::INSTR` e si può passare
anche esplicitamente con `--resource`.

### USB su Linux

Nel venv usato su Linux servono anche lì `pyusb` e `libusb-package`:

```bash
pip install pyvisa pyvisa-py pyusb libusb-package numpy scipy matplotlib
```

Da utente normale Linux non permette di aprire il dispositivo USB: pyvisa-py
allora non lo elenca e `--usb` risponde che non trova l'oscilloscopio, indicando
la causa. Regola udev (su openSUSE Aeon `/etc` è scrivibile senza
`transactional-update`), poi staccare e riattaccare il cavo:

```bash
echo 'SUBSYSTEM=="usb", ATTR{idVendor}=="2a8d", ATTR{idProduct}=="0387", MODE="0666"' | sudo tee /etc/udev/rules.d/99-keysight.rules
sudo udevadm control --reload-rules && sudo udevadm trigger
```

`MODE="0666"` e non `TAG+="uaccess"`: lo script deve poter girare anche come
servizio, senza un utente connesso. Se il dispositivo risulta occupato e
`lsmod | grep usbtmc` mostra il driver del kernel:

```bash
echo "blacklist usbtmc" | sudo tee /etc/modprobe.d/no-usbtmc.conf
sudo modprobe -r usbtmc
```

Con Python dentro distrobox/toolbox, la regola udev va comunque sul sistema host.

LAN diretta (VXI-11), con il PC sulla stessa rete dell'oscilloscopio:

```bash
python acquire_coincidences.py --host 192.168.1.20 -o run_3m --check
```

L'IP si legge sull'oscilloscopio in Utility → I/O. Altre risorse VISA si passano con
`--resource`, per esempio `TCPIP::192.168.1.20::5025::SOCKET` o una risorsa `USB0::...::INSTR`.
Attraverso il bastion SSH, come per l'RTO6, si usa il socket SCPI (porta 5025):

```bash
python acquire_coincidences.py --ssh-bastion deguio@cerbero.mib.infn.it --host <IP oscilloscopio> -o run_3m --check
```

Per una presa dati di giorni conviene un PC in laboratorio connesso direttamente.
Su Mac impedire lo standby: `caffeinate -i python acquire_coincidences.py ...`.

## Procedura consigliata

Scale e finestra sono **fisse nello script** (default), scelte sugli impulsi
misurati con questo oscilloscopio il 9 ottobre 2026: positivi, mediana circa
25–40 mV, coda di decadimento circa 190 ns, picchi fino a oltre 170 mV.

- **25 mV/div con 87.5 mV al centro**, su entrambi i canali: schermo da −12.5 a
  +187.5 mV, baseline mezza divisione sopra il fondo. L'ADC resta valido circa
  5.1 divisioni sopra il centro, cioè fino a circa +215 mV. Con 20 mV/div e 70 mV
  al centro un evento su 26 di `run_3m` è saturato a 173 mV. Il gradino dell'ADC
  è circa 1.0 mV.
- **100 ns/div con il centro a +300 ns**: finestra da −200 a +800 ns, 2000 punti
  a 0.5 ns (2 GSa/s, il massimo con due canali). Circa 190 ns di baseline prima
  del fronte e coda sotto il 5% del picco.
- Il trigger resta quello del pannello (`--trigger keep`, default): pattern
  `"11X"`, C1 e C2 sopra 10 mV. Si salvano solo coincidenze.

Le opzioni `--scale-mv`, `--offset-mv`, `--timebase-ns`, `--delay-ns` cambiano
questi valori. Lo script avvisa se il campionamento scende sotto 2 GSa/s; ogni
campione all'estremo dell'ADC è segnalato come `ADC LIMIT` e contato in `events.csv`.

**1. Controllo** (`--check`): applica le impostazioni, stampa quelle lette e
esce senza acquisire.

```bash
python acquire_coincidences.py --usb -o run_3m --check
```

**2. Prova breve** con i rivelatori vicini, dove il rate è alto:

```bash
python acquire_coincidences.py --usb -o test_vicini --max-events 5
python analyze_pulses.py test_vicini/scope_000000.csv --show
```

**3. Presa dati a 3 m**, in una cartella dedicata:

```bash
caffeinate -is python acquire_coincidences.py --usb -o run_3m
```

Ordine di grandezza atteso, assumendo scintillatori di 5×5 cm² e intensità
verticale di circa 70 m⁻² s⁻¹ sr⁻¹: I·A²/d² ≈ 5·10⁻⁵ Hz, cioè circa 0.2
coincidenze/ora (pochi eventi al giorno). Il rate misurato è riportato
periodicamente sul terminale e in `session_*.json`.

### Trigger alternativo: edge e coincidenza software

```bash
python acquire_coincidences.py --host 192.168.1.20 -o run_3m_edge \
    --trigger edge --edge-source 1 --trigger-level-mv 5 --min-amplitude-mv 5
```

Trigger su un solo canale. L'evento viene salvato solo se **entrambi** i picchi
superano `--min-amplitude-mv` rispetto alla baseline (primo 15% della finestra,
come nell'analisi). Utile se il pattern trigger non funziona come previsto.
Si trasferiscono però tutti i singoli, e ogni trasferimento è tempo morto.
`--save-rejected` salva anche gli scarti nella sottocartella `rejected/`.

## File prodotti

| File | Contenuto |
|---|---|
| `scope_NNNNNN.csv` | un evento: tempo [s], C1 e C2 [V], finestra `--window-ns` (default −250…800 ns) |
| `events.csv` | una riga per trigger, anche se scartato: ora PC (locale e UTC), tempo vivo, tempo di lettura, picchi, campioni al limite dell'ADC, passo |
| `session_<data>.json` | comando, IDN, impostazioni lette dallo strumento, avvisi, conteggi, tempo vivo |
| `session_<data>_setup.bin` | setup binario dello strumento (`:SYSTem:SETup?`), se supportato |

Rilanciare con la stessa cartella riprende la numerazione: nessun file viene
sovrascritto. I file compaiono con il nome definitivo solo a scrittura completata.
`adc_limit1/2 > 0` indica campioni a fondo scala (probabile clipping): in quel
caso aumentare la scala o usare `--max-amplitude-mv` nell'analisi.

`live_s` è il tempo fra l'armamento e la fine dell'acquisizione, a meno del
passo di polling (`--poll-s`, default 50 ms). La somma dei `live_s` è il tempo
vivo per il calcolo del rate. Il tempo morto è dato dalle letture e dalle pause
di riconnessione.

In caso di errore di comunicazione durante la presa dati, lo script:

1. scrive l'errore, con data e ora, in `errors.log` nella cartella della presa dati;
2. via USB esegue il device clear USBTMC (INITIATE_CLEAR, CHECK_CLEAR_STATUS,
   sblocco degli endpoint), che pyvisa-py non implementa; l'acquisizione in corso
   sull'oscilloscopio non viene toccata;
3. riapre la sessione con attesa crescente fino a 60 s, senza limite salvo
   `--max-reconnects`;
4. **prima di riarmare**, se l'oscilloscopio è fermo con un trigger arrivato
   durante l'interruzione, o con un evento la cui lettura era fallita, legge e
   salva quell'evento. In `events.csv` il suo `live_s` è vuoto se il tempo vivo
   non è noto; il conteggio `recovered` del JSON di sessione li riporta.

Se le impostazioni lette dopo una riconnessione differiscono da quelle iniziali,
lo script lo segnala e le registra nel JSON della sessione. Un errore alla prima
connessione, o un'impostazione rifiutata dallo strumento, interrompe subito.

### Prese dati lunghe con pochi eventi

Via USB lo script esegue il clear USBTMC anche **a ogni apertura** della
connessione. Il 9 ottobre 2026 un trasferimento rimasto a metà aveva lasciato il
DSOX1202A in uno stato in cui ogni comando successivo andava in timeout in
scrittura (`[Errno 60] Operation timed out` su macOS, timeout USB su Linux),
anche dopo aver rilanciato lo script; il reset della porta USB non bastava, il
clear USBTMC sì. Nello stesso giorno un'attesa di 296 s senza trigger, con
interrogazione ogni 0.25 s, si è conclusa regolarmente con il salvataggio
dell'evento.

Ogni `--status-s` (default 300 s) lo script stampa una riga `alive` e aggiorna
il JSON di sessione (`updated`, `last_event`, conteggi, tempo vivo) **anche
mentre aspetta il trigger**: l'ora di `updated` distingue "nessun muone" da
"acquisizione ferma". Lo stato viene interrogato ogni `--poll-s` (default 0.25 s).

Sul PC Linux disattivare la sospensione automatica (GNOME sospende di default
dopo circa 15 minuti di inattività, e al risveglio la comunicazione USB va in
timeout):

```bash
sudo systemctl mask sleep.target suspend.target hibernate.target hybrid-sleep.target
```

## Analisi

```bash
python analyze_pulses.py run_3m --distance-m 3 --offset-ns <offset calibrato> --min-amplitude-mv 10 --max-amplitude-mv 130
```

Per i file prodotti da questo script, il grafico del rate usa l'ora del PC
registrata in `events.csv`, perché il CSV Keysight non contiene un timestamp.
Gli export fatti dal pannello frontale restano senza timestamp. L'offset
strumentale va calibrato con **questo** oscilloscopio e questi cavi: per esempio
rivelatori sovrapposti, oppure scambiando i rivelatori sui canali. Il valore
dell'RTO6 non è trasferibile.

## Note sull'hardware

- Il DSOX1202A ha solo ingressi da 1 MΩ. Se il segnale arriva su cavo coassiale
  e con l'RTO6 era terminato a 50 Ω, usare un terminatore passante da 50 Ω: forma
  e ampiezza dipendono dalla terminazione.
- Il fattore sonda deve essere 1 per un cavo diretto: lo script avvisa se non lo è,
  perché le tensioni lette includono quel fattore.
- Con entrambi i canali attivi il campionamento massimo misurato è 2 GSa/s per canale (0.5 ns).
  `--points-mode raw` (default) trasferisce il record completo; chiedere meno punti
  del record lo decima.

## Prova sullo strumento (9 ottobre 2026, USB)

Con le impostazioni del pannello (`--trigger keep`): pattern trigger C1 e C2
sopra 10 mV, 20 mV/div con 70 mV al centro, 41 ns/div con il centro a +150 ns,
2 GSa/s. Sono state acquisite 5 coincidenze in `prova_usb/`, analizzate con
`analyze_pulses.py`. Verificato:

- armamento con `:AER?` (circa 25 ms) e fine acquisizione con il bit di run di
  `:OPERegister:CONDition?`;
- record `raw` di 820 punti, passo 0.5 ns, finestra da −55 a +355 ns (lo schermo);
  lettura di entrambi i canali e scrittura in circa 70–80 ms per evento;
- `:WAVeform:POINts` con un numero superiore al record è un **errore** (−222),
  non viene ridotto al massimo: per questo di default lo script non lo invia;
- il pattern del modello a 2 canali ha tre caratteri, C1 C2 EXT (`"11X"`, letto
  in formato ASCII): lo script costruisce la stringa sulla lunghezza letta;
- PyVISA-py non supporta il device clear via USB: all'apertura lo script scarta
  eventuali risposte rimaste da una sessione interrotta.

Questa prova usava la finestra del pannello (−55…+355 ns), poi sostituita
dai valori fissi descritti sopra. Con soli 55 ns prima del trigger, il primo 15%
del record arriva a +6 ns, dentro il fronte di salita. Lo script calcola le ampiezze con la baseline fino a −20 ns
(`--baseline-end-ns`); nell'analisi usare lo stesso accorgimento:

```bash
python analyze_pulses.py prova_usb --distance-m 3 --offset-ns 0 --baseline-end-ns -15
```

Gli eventi 0–2 di `prova_usb` hanno in `events.csv` ampiezze registrate con la
vecchia baseline (sottostimate); le forme d'onda salvate sono corrette.
Più spazio prima del trigger (`--delay-ns` minore) e una coda più lunga
(`--timebase-ns` maggiore) migliorerebbero baseline e fit di τ di decadimento.

Non ancora provati sullo strumento: l'impostazione del trigger da script
(`--trigger pattern|edge`), le riconnessioni automatiche e la modalità SSH.
I test offline (`python -m unittest -v test_acquire_coincidences.py`) usano un
modello simulato dei comandi SCPI.
