# Analisi offline delle pulse shape RTO6

`analyze_pulses.py` legge la coppia di file `.csv` (metadati) e `.Wfm.csv`
(campioni) e analizza due canali della stessa acquisizione. Non comunica con
l'oscilloscopio. Richiede numpy, scipy e matplotlib, già presenti nel venv.

Dalla cartella `RTO6`, con il venv attivo:

```bash
python analyze_pulses.py RefCurve_2026-09-21_45_175205.csv --offset-ns 16 --show
```

Con `--show` apre due figure separate: un riepilogo compatto (1320×825 pixel)
con residui sotto gli impulsi e parametri del solo fit causale in sovraimpressione,
e una vista tipo oscilloscopio (1210×605 pixel) con entrambi i canali e i fit.
I campioni e le curve smussate di ciascun canale hanno lo stesso colore.
Senza `--show` salva comunque entrambe le figure in PNG (`_timing` e `_scope`),
parametri/covarianze/risultati in JSON e scansione della frazione in CSV nella
cartella `timing_results`.
È possibile passare direttamente anche il file `.Wfm.csv`.

Opzioni utili:

```bash
python analyze_pulses.py RefCurve_2026-09-21_45_175205.csv --offset-ns 16 --fraction 0.3 --output-dir timing_f30
python analyze_pulses.py RefCurve_2026-09-21_45_175205.csv --offset-ns 16 --smooth-samples 1 --output-dir timing_raw
python analyze_pulses.py RefCurve_2026-09-21_45_175205.csv --offset-ns 16 --baseline-end-ns -100 --fit-range-ns -150 300 --output-dir timing_short
```

La finestra di fit deve includere almeno 20 campioni di fondo e l'impulso.
Il fondo predefinito occupa il primo 15% del record; se l'impulso cade lì,
selezionare esplicitamente `--baseline-end-ns`. Il file deve contenere un solo
impulso isolato per canale; pile-up e clipping richiedono una selezione diversa.
Il programma non è un selettore automatico di coincidenze o di tracce di muoni.

## Modelli e significato di t0

Il primo fit usa, per u = t − t0:

    V(t) = B + S/N · H(u) · [1 − exp(−u/τr)] · exp(−u/τd)

`B`, `S`, `t0`, `τr` e `τd` sono tutti liberi, separatamente per ciascun canale.
`N` normalizza la forma non convoluta a picco unitario; `S` è quindi il picco
sopra il fondo in questo modello. Le costanti τ sono parametri del modello,
non tempi di salita/discesa 10–90%.

La forma è una differenza di esponenziali con costanti τd e
τeff = τr·τd/(τr+τd). A normalizzazione vicino, è equivalente alla convoluzione
di due risposte esponenziali causali. È un modello fenomenologico della catena
scintillatore/SiPM/elettronica, non una ricostruzione unica dei suoi componenti.

Il secondo fit convoluisce questa forma con una gaussiana normalizzata,
aggiungendo la larghezza libera σ. La convoluzione è valutata analiticamente
con una forma numericamente stabile basata su `scipy.special.log_ndtr`.
`S` rimane la scala del picco *prima* della convoluzione. `t0` è l'inizio della
forma causale sottostante: la forma convoluta non ha un inizio netto.
I fit minimizzano i residui sull'intero record, o nella finestra richiesta,
senza smoothing dei campioni e con più inizializzazioni.

Il constant fraction digitale sottrae il fondo, riconosce la polarità e
interseca il fronte crescente con f volte il massimo. Il default usa f=0.20,
Savitzky–Golay di ordine 2 su 5 campioni per stabilizzare picco e attraversamento,
e interpolazione lineare fra i due campioni adiacenti. Si tratta di un
attraversamento a frazione del picco offline, non dell'emulazione di un CFD
analogico con ramo ritardato. Il tempo t_f non coincide con t0; la differenza
fra i due canali elimina lo spostamento t_f−t0 solo per forme sufficientemente
simili. L'amplitude walk viene rimosso idealmente per impulsi che differiscono
solo per ampiezza; rumore, forma variabile e stima del picco restano rilevanti.

Riferimenti primari:

- [ORTEC, Timing Electronics and Fast Timing Methods with Scintillation Detectors](https://www.ortec-online.com/-/media/ametekortec/technical-papers/time-spectrometry-and-timing-applications/timing-electronics-fast-timing-methods-scintillation-detectors.pdf?hash=1C19DD3A85453CD394A854BD21AA999C&la=en&revision=d3898bd7-0132-42cb-9987-4c384ef4967c)
- [SciPy least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)
- [SciPy log_ndtr](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_ndtr.html)

## Risultato sul file del 21 settembre 2026

2000 campioni per canale, passo temporale 0.5 ns; tempo costruito come
XStart + i·Resolution, senza includere XStop come ultimo campione.
C1 e C2 sono due colonne dello stesso export. La quantizzazione osservata è
circa 0.79 mV: su C2 il picco sopra il fondo occupa solo circa 8 livelli.
Il rumore RMS pre-pulse è circa 0.39 mV in entrambi i canali.

Convenzione: ΔT = t(C2) − t(C1). Si sottraggono i 16 ns di ritardo aggiuntivo
sul cavo C2 indicati dall'utente. Questo non certifica che tutti gli altri
ritardi di rivelatori/elettronica/canali siano già compensati.

| Metodo | ΔT osservato [ns] | ΔT meno 16 ns [ns] |
|---|---:|---:|
| Frazione costante 20% | 22.774 | 6.774 |
| Fit causale, t0 libero | 23.089 | 7.089 |
| Fit convoluto gaussiano, t0 libero | 23.038 | 7.038 |

Nel fit causale: t0(C1) = −51.089 ns, t0(C2) ≈ −28.000 ns;
τr ≈ 10.24/10.33 ns, τd ≈ 195.78/179.98 ns; χ²/ndf ≈ 1.28/0.86.
L'errore formale indipendente sul ΔT del fit causale è circa 0.27 ns, ma
**non è l'incertezza totale sulla misura del tempo di volo**.
La gaussiana dà σ ≈ 1.04 ns su C1, mentre su C2 raggiunge quasi il limite
inferiore (0.01 ns): su questo evento il parametro aggiuntivo non è ben
vincolato. I residui di C1 hanno correlazione al lag 1 ≈ 0.38.

Fra il 20% e il 30%, il ΔT misurato varia circa da 22.71 a 22.86 ns.
Fra il 20% e il 60% varia circa da 22.71 a 23.89 ns; oltre il 60% cresce molto.
Questi intervalli sono sensibilità alla scelta dell'estimatore su **un evento**,
non barre d'errore statistiche. La soglia del 20% su C2 è solo 3.46 volte
il rumore RMS; il 10% è troppo vicino al fondo per considerarlo robusto.
La scansione di smoothing fra 1 e 9 campioni dà circa 22.77–22.96 ns al 20%.

Un controllo del fit causale con inizio finestra −150 ns e termine variabile
150/300/500/784 ns dà rispettivamente circa 22.89/23.05/23.04/23.09 ns.
I numeri sono salvati in `timing_results/fit_window_sensitivity.json`.

**Scelta iniziale:** usare il constant fraction al 20–30% come stimatore
semplice e il fit causale come controllo della forma e di t0. Non c'è
sufficiente evidenza per preferire il fit con σ libero. Per stabilire quale
metodo abbia davvero migliore risoluzione servono molte coincidenze, una
scansione della frazione e il controllo di ΔT contro ampiezza, rapporto di
ampiezze e tempo di salita, su dati diversi da quelli usati per scegliere i
parametri. Non scegliere la frazione imponendo che la velocità risulti c.

## Da ΔT alla velocità

Un percorso di 2 m a c richiede 6.671 ns: i risultati hanno la scala temporale
attesa, ma con una sola coppia e una distanza “poco meno di 2 m” non è possibile
attribuire una precisione alla velocità. Il cammino reale del muone può essere
maggiore della separazione verticale se la traiettoria è inclinata.

Quando il cammino L e l'offset sono noti, si può aggiungere `--distance-m L`:
lo script calcola v = L / |ΔT − offset| e v/c. Il pannello in basso a sinistra
riporta shift, ΔT corretto, distanza e velocità stimata con il constant fraction.
La distanza proviene esclusivamente da `--distance-m`: se manca non viene
assunto alcun valore e non viene mostrata una velocità. Si tratta della velocità
del muone confrontata con c, non di una misura diretta della velocità della luce. Non assume la distanza di 2 m
per questo file; il report attuale non contiene una velocità numerica.
Il segno del ΔT corretto rimane nel report anche se la velocità usa il modulo.

Per una misura quantitativa occorre calibrare lo skew dei canali/cavi con
un impulso comune e caratterizzare anche l'offset relativo dei rivelatori,
per esempio mediante geometrie simmetriche o inversioni controllate.
L'incertezza sui 16 ns si trasferisce direttamente al tempo di volo.
Ottimizzare la scala verticale di C2 migliorerebbe il numero di livelli
sul fronte, evitando saturazione. Il passo di campionamento di 0.5 ns non è,
da solo, la risoluzione temporale dello stimatore.

Le covarianze JSON assumono campioni indipendenti e rumore RMS pari al fondo;
non includono correlazioni, quantizzazione modellata esplicitamente,
inadeguatezza del modello, jitter comune, calibrazione o distanza.
`noise_over_slope_ns` del CFD è solo un indicatore locale, non include tutta
l'incertezza sul picco/fondo/smoothing. Le bandiere sui parametri ai limiti e
sui residui correlati devono essere considerate insieme al χ².

## Verifiche

```bash
python -m unittest -v test_analyze_pulses.py
```

Quattro test verificano lettura/base temporale e rifiuto di history multiple,
invarianza del CFD a scala/fondo/polarità e ritardo noto, recupero di t0 libero
su impulsi sintetici per entrambi i modelli e convoluzione analitica contro
quadratura numerica indipendente. Sono stati eseguiti anche il fit sui dati
reali e il controllo visivo del grafico prodotto.

## Analisi di più acquisizioni

Lo stesso script accetta uno o più nomi di file, una cartella o pattern tra
apici. Header e `.Wfm.csv` della stessa acquisizione vengono deduplicati.
La ricerca nelle cartelle non è ricorsiva. Ogni acquisizione resta una coppia
di impulsi simultanei; non si tratta di cercare più impulsi nello stesso record.

Dalla cartella RTO6, per i tre eventi 45, 46 e 47:

```bash
python analyze_pulses.py 'RefCurve_2026-09-21_4[567]*.csv' --distance-m 2 --offset-ns 16 --show
```

Oppure, per tutti gli export RTO presenti nella cartella corrente:

```bash
python analyze_pulses.py . --distance-m 2 --offset-ns 16
```

I valori 2 m e 16 ns sono argomenti di questi esempi: non sono default.
Tutti gli eventi del batch devono condividere distanza, configurazione dei
canali e offset. Non mescolare configurazioni diverse nello stesso summary.

Per ogni acquisizione vengono salvati `_timing.png`, `_scope.png`,
`_timing.json` e `_fraction_scan.csv`. In modalità batch `--show` apre insieme,
alla fine, il summary di velocità/ΔT e quello del rate; i grafici
individuali rimangono su disco. Su un singolo file il comportamento precedente
è preservato e `--show` apre i due grafici individuali.

Il summary usa di default il constant fraction alla frazione selezionata.
Per usare invece la differenza dei t0 fitted:

```bash
python analyze_pulses.py 'RefCurve_2026-09-21_4[567]*.csv' --distance-m 2 --offset-ns 16 --summary-method causal
```

`--summary-method gaussian` usa il fit della pulse shape convoluto. Questa
opzione seleziona lo stimatore temporale, non il tipo di fit dell'istogramma:
il fit della distribuzione è sempre una normale. Per modificare la frazione
si usa `--fraction`, come nell'analisi individuale.

Il file `summary_constant_fraction.png/.pdf` (oppure `summary_causal` /
`summary_gaussian`) contiene:

- istogramma normalizzato delle velocità v/c e gaussiana fitted;
- istogramma dei tempi corretti e gaussiana fitted, con L/|media dei tempi|;
- media μ, larghezza σ, errore statistico sulla media e intervallo al 95%;
- numero di eventi, errori di lettura/analisi ed eventuali velocità indefinite.

La gaussiana è stimata mediante massima verosimiglianza sui valori individuali
(`scipy.stats.norm.fit`), **non mediante un fit dei conteggi dei bin**: con tre
eventi un fit a tre parametri sui bin sarebbe particolarmente fragile. μ è la
media aritmetica e σ_MLE è la deviazione standard con denominatore N.
L'errore sulla media usa s/√N, con s calcolata con denominatore N−1;
l'intervallo al 95% usa Student con N−1 gradi di libertà. Entrambi assumono
eventi indipendenti distribuiti normalmente. Con un solo valore o dispersione
nulla non viene disegnata una gaussiana di larghezza artificiale.

Riferimento: [SciPy, stima MLE delle distribuzioni](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html).

Le incertezze comuni di distanza, cavo e offset dei rivelatori non sono incluse
e non diminuiscono come 1/√N. Il fit gaussiano alle velocità è descrittivo:
se il rumore è additivo sui tempi e il cammino è comune, è più naturale stimare
prima il tempo medio e convertirlo in velocità. L'inversione v=L/t è non lineare:
la media delle velocità può essere spostata verso valori alti dai piccoli t.
Il pannello dei tempi permette di vedere questa differenza. Con segni dei tempi
corretti discordi non si riporta una velocità da tempo medio per direzione comune.
Una singola acquisizione non stabilisce la direzione senza ulteriori informazioni.

Il CSV del summary contiene una riga per evento, tempo misurato/corretto,
velocità, parametri comuni e segnalazioni sulla qualità del fit. Il JSON
contiene anche statistiche, ipotesi, errori ed esclusioni. Gli eventi con
v>c non vengono tagliati né limitati a c. I fit imperfetti restano inclusi e
segnalati; definire eventuali selezioni su criteri indipendenti dalla velocità.
Una velocità non definita (tempo corretto nullo) viene documentata ed esclusa
dall'istogramma; il batch continua dopo errori di lettura/analisi gestiti e
termina con codice 1 se almeno un export fallisce.

### Primo batch di tre eventi

Eseguito con L=2 m, offset=16 ns, frazione 20%:

| Evento | Tempo corretto [ns] | v/c apparente |
|---|---:|---:|
| 45 | 6.774 | 0.985 |
| 46 | 3.332 | 2.002 |
| 47 | 5.761 | 1.158 |

Gaussiana sulle velocità: μ/c ≈ 1.382, σ/c ≈ 0.444,
SEM ≈ 0.314 c. Dal tempo medio si ottiene invece circa 1.261 c.
Sono stime preliminari della distribuzione delle **misure**, non evidenza di
muoni superluminali né una misura conclusiva della velocità. Tre eventi non
consentono di verificare l'ipotesi gaussiana. L'evento 46 richiede particolare
attenzione alla coincidenza, alla forma e alla calibrazione.

I test ora includono anche discovery/deduplicazione degli export, formule del
fit gaussiano e gestione di velocità >c, tempi nulli, segni discordi ed errori.

## Soglia minima sul picco di entrambi i canali

```bash
python analyze_pulses.py '16ps/RefCurve_2026-09-21_*.csv' --distance-m 1.9 --offset-ns 16 --summary-method causal --min-amplitude-mv 20 --output-dir 16ps/timing_min20 --show
```

`--min-amplitude-mv 20` richiede che **entrambi** i canali abbiano un picco di
almeno 20 mV rispetto alla propria baseline. L'ampiezza per questa selezione è
`max(abs(V - baseline))` sui campioni originali, prima di smoothing e fit:
vale anche per impulsi negativi e non è l'ampiezza S del modello fitted.
La baseline usa la stessa regione pre-pulse dell'analisi. La soglia è inclusiva:
un picco esattamente uguale al valore richiesto passa. Il default è nessun taglio.

Gli eventi sotto soglia vengono saltati prima del fit, senza interrompere il
batch e senza essere contati come errori. Per essi non vengono prodotti nuovi
grafici; il summary usa solo gli eventi accettati e analizzati con successo.
Il titolo del summary riporta soglia e numero di scarti. Il JSON del summary
contiene `amplitude_rejected_events`, con picchi e motivo; il file
`summary_<metodo>_selection.csv` registra la selezione degli eventi letti
correttamente. Ogni evento ha inoltre un file `_selection.json`, anche quando
si analizza un solo export. Se tutti gli eventi sono scartati, il summary
indica che non ci sono velocità utilizzabili e non esegue il fit gaussiano.

Usare una cartella di output diversa per confrontare soglie: i vecchi grafici
di eventi ora scartati non vengono cancellati se erano già presenti.

Il binning di entrambi gli istogrammi del summary viene ricalcolato a ogni
esecuzione con `Nbin = ceil(2*sqrt(N))`, dove N conta solo gli eventi inclusi nel
summary dopo selezioni ed esclusioni, senza limiti fissi sul numero di bin.
Per esempio: 3 eventi → 4 bin; 17 → 9; 100 → 20; 1000 → 64.
Se la larghezza di Freedman–Diaconis (`2 × IQR / N^(1/3)`) è più fine,
si aumenta ulteriormente il numero di bin per rispettarla. Questo evita
che code estese rendano poco visibile la struttura del picco. Con IQR nullo
si mantiene la regola precedente.
Il default privilegia una visualizzazione più fine; `--summary-bins 20`
permette di scegliere esplicitamente il numero di bin per entrambi i pannelli.
L'intervallo si adatta ai valori osservati, con un margine ai bordi. La legenda
riporta N e numero di bin. La gaussiana resta stimata sui valori individuali,
quindi il risultato del fit non dipende dal binning.

### Soglia massima di ampiezza

`--max-amplitude-mv` impone anche un limite superiore sullo stesso picco dei
campioni originali, in valore assoluto rispetto alla baseline di ciascun canale.
Può essere usato da solo o insieme al minimo:

```bash
python analyze_pulses.py '16ps/RefCurve_2026-09-21_*.csv' --distance-m 1.9 --offset-ns 16 --summary-method causal --min-amplitude-mv 20 --max-amplitude-mv 100 --output-dir 16ps/timing_20_100 --show
```

In questo esempio entrambi i picchi devono essere nell'intervallo inclusivo
[20, 100] mV. Un solo canale fuori intervallo basta per scartare l'evento.
Soglia massima e motivo dello scarto sono registrati nei report di selezione;
il summary mostra l'intervallo applicato. Minimo maggiore del massimo, valori
negativi e non finiti vengono rifiutati dalla CLI.
Il massimo è un taglio in ampiezza, non un riconoscimento automatico del clipping:
per escludere saturazioni va scelto sotto il livello di saturazione osservato,
misurato rispetto alla baseline. Un picco esattamente uguale al massimo passa.

### Rate degli eventi nel tempo

In modalità batch viene generato anche `summary_rate.png/.pdf`, con due pannelli:
 tutte le acquisizioni in input e quelle che superano i tagli d'ampiezza su entrambi
 i canali. Gli eventi selezionati sono contati indipendentemente dalla riuscita
 del fit o dalla disponibilità di una velocità finita.

`--rate-bin-hours 2` imposta la larghezza dei bin in ore (default: 2).
I pannelli usano gli stessi intervalli e riportano eventi/ora, errori √N divisi
per la durata del bin, numero di eventi e rate medio. L'ultimo bin viene
normalizzato alla sua durata effettiva. `--show` apre questo grafico insieme
al summary di velocità e ΔT alla fine del processing.

I tempi provengono dal campo `Timestamp` dei metadati dell'oscilloscopio.
L'intervallo va dalla prima all'ultima acquisizione: è tempo trascorso, non
live time misurato; eventuali pause sono incluse e non si corregge il dead time.
Record senza timestamp valido vengono segnalati ed esclusi dal grafico;
se la waveform non è leggibile, il record con timestamp valido viene contato
solo tra tutte le acquisizioni. Con meno di due timestamp distinti il rate
non è definito. Conteggi, rate e stato dei record sono salvati anche in
`summary_rate.csv` e `summary_rate.json`.

### Formati dei grafici

Il formato predefinito è solo PNG. `--plot-formats png` salva solo PNG;
`--plot-formats pdf` salva solo PDF; `--plot-formats png pdf` salva entrambi.
La scelta vale per i grafici individuali (`_timing`, `_scope`) e per entrambi
 i summary (velocità/ΔT e rate). JSON e CSV vengono sempre scritti;
`--show` continua ad aprire le figure indipendentemente dal formato scelto.
I file di altri formati già presenti nella cartella non vengono cancellati
né aggiornati.

```bash
python analyze_pulses.py '16ps/RefCurve_2026-09-2*' --offset-ns 15.69 --distance-m 1.87 --min-amplitude-mv 25 --max-amplitude-mv 170 --fraction 0.3 --summary-method causal --plot-formats png --show
```
