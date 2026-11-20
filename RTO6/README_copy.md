# Copia dei file RTO6 via bastion SSH

Eseguire dal terminale del Mac, con il venv attivo:

```bash
python copyFileToPc.py --ssh-bastion <utente>@<bastion> -f 'RefCurve_2026-09-21*' -o ~/Desktop
```

La destinazione indicata è il bastion fornito dall’utente; si può usare anche un alias di `~/.ssh/config`.
Lo script usa OpenSSH di sistema: chiavi, ssh-agent, password/MFA e impostazioni
SSH rimangono gestiti da SSH, senza memorizzare credenziali nello script.
Porte SSH personalizzate, IdentityFile o ulteriori ProxyJump possono essere
configurati nell'alias SSH. Il bastion non deve avere Python o RsInstrument.

Il tunnel inoltra `127.0.0.1:15025` sul Mac verso `<IP-oscilloscopio>:5025` visto
dal bastion. RsInstrument usa `TCPIP::127.0.0.1::15025::SOCKET` con backend
`socketio`, quindi non necessita della libreria nativa rsvisa per questa modalità.
Il tunnel viene chiuso al termine e in caso di eccezioni gestite; viene chiusa
solo la sessione SSH creata dallo script. Il listener è limitato al loopback.
La modalità senza `--ssh-bastion` conserva la connessione diretta VXI-11/pyvisa-py.

Prerequisiti da verificare nell'installazione reale:

- il bastion consente il TCP forwarding e raggiunge la rete dell'oscilloscopio;
- il servizio SCPI socket dell'oscilloscopio è attivo e raggiungibile sulla
  porta indicata (5025 è la porta tipica R&S, non verificata su questo dispositivo);
- la porta locale 15025 è libera; altrimenti usare `--local-port 15026`.

L'indirizzo e la porta remoti sono configurabili con `--instrument-host` e
`--instrument-port`. L'apertura del listener SSH non garantisce che la porta
remota sia raggiungibile: questo viene verificato dalla connessione allo strumento.
Lo script continua a leggere esclusivamente file già salvati; non invia comandi
di acquisizione o esportazione.

## Tunnel manuale, utile per la diagnostica

In un terminale:

```bash
ssh -NT -o ExitOnForwardFailure=yes -L 127.0.0.1:15025:<IP-oscilloscopio>:5025 <utente>@<bastion>
```

In un secondo terminale:

```bash
python copyFileToPc.py --resource 'TCPIP::127.0.0.1::15025::SOCKET' -f 'RefCurve_2026-09-21*' -o ~/Desktop
```

Il tunnel manuale resta aperto finché non si interrompe SSH con Ctrl-C.
Non basta sostituire l'IP con localhost nella vecchia risorsa `::INSTR`:
VXI-11 usa RPC e porte aggiuntive; qui si sceglie esplicitamente il socket SCPI.
Se il socket SCPI non è disponibile, occorre valutare un protocollo alternativo
supportato dal dispositivo, una VPN universitaria o eseguire il client VISA su
un computer interno e poi trasferire i file via SSH.

Riferimenti:

- [OpenSSH: forwarding locale -L e controllo -O](https://man.openbsd.org/ssh)
- [RsInstrument: connessioni e backend socketio](https://rsinstrument.readthedocs.io/en/latest/StepByStepGuide.html)

Verificati CLI e tre test offline (`python -m unittest -v test_copyFileToPc.py`)
su selezione backend e chiusura tunnel in caso di errore. Nessuna connessione
reale al bastion è stata eseguita.

## Se rimane fermo su «Copying»

`-f` richiede un nome completo oppure un pattern esplicito. Per esempio:

```bash
python copyFileToPc.py --ssh-bastion <utente>@<bastion> -f 'RefCurve_2026-09-21_47*' -o ~/Desktop
```

`RefCurve_2026-09-21_47` senza `*` non significa «tutti i file dell'evento 47».
Lo script ora controlla il catalogo remoto anche per i nomi esatti prima di
inviare MMEM:DATA?. Il timeout di ogni operazione VISA è 30 secondi, configurabile
con `--timeout-ms`; non è il limite della durata totale della copia.
Mostra identificazione, ricerca e byte ricevuti. I file vengono scaricati in
una cartella temporanea nella destinazione e pubblicati sotto il nome definitivo
solo dopo il completamento, senza sovrascrivere file esistenti. La pubblicazione
usa un hard link locale (supportato dal filesystem APFS del Mac).
I vecchi file vuoti creati da tentativi precedenti non vengono rimossi.

I file già presenti nella destinazione vengono saltati senza interrompere il batch,
anche se tutti i match sono già presenti. Viene mostrato il conteggio di copiati
e saltati. Il controllo è sul nome: non confronta contenuto o dimensione.

## Filtro per data e ora

`--since 'YYYY-MM-DD HH:MM[:SS]'` e `--until 'YYYY-MM-DD HH:MM[:SS]'`
selezionano un intervallo inclusivo; ciascun limite è facoltativo.
Il filtro viene applicato dopo il wildcard e prima dei trasferimenti.
Usa la data/ora di esportazione nel nome `RefCurve_YYYY-MM-DD_N_HHMMSS`,
non il timestamp interno dell'acquisizione né l'orario del Mac. Non vengono
applicate conversioni di fuso orario. Se i due orologi differiscono, scegliere
il limite corrispondente nell'orologio usato dai nomi dei file.
I companion `.csv` e `.Wfm.csv` hanno lo stesso orario e passano insieme,
purché entrambi matchino il wildcard. Con un filtro attivo, nomi senza data/ora
riconoscibile vengono saltati con un messaggio. Nessun file nell'intervallo
non è un errore; i file già presenti continuano a essere saltati.

Per copiare nella cartella `0ps` le esportazioni dal 24 settembre 2026 alle
14:15 in avanti (eseguire dalla cartella RTO6):

```bash
mkdir -p 0ps
python copyFileToPc.py --ssh-bastion <utente>@<bastion> -f 'RefCurve_*' --since '2026-09-24 14:15' -o 0ps/
```

Per limitarsi alla sola giornata, usare `-f 'RefCurve_2026-09-24*'`.
Per una fascia chiusa, aggiungere ad esempio `--until '2026-09-24 18:00'`.
