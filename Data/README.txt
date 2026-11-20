Questa cartella contiene solo dati e risultati. In git ci sono solo:
  ExampleData/        file di esempio dell'autore originale
  TestData/           estratti piccoli per provare gli script (vedi TestData/README.md)
  showers_analysis/   documentazione, figure e posizioni dell'analisi sciami
  pressure_analysis/  notebook degli studenti (senza dati)

Gli script di acquisizione e analisi stanno in ../analysis_code/:
  import_data.py      acquisizione via USB da uno o più rivelatori
  plot.py             analisi di un file (rate, spettri, coincidenze, sensori)
  display.py          event display e efficienze della torre impilata
  shower_analysis.py  ricerca di sciami con rivelatori affiancati
Istruzioni complete in ../analysis_code/README.md.

Acquisizione rapida (venv attivo, rivelatori collegati via USB):
  python ../analysis_code/import_data.py
Il file viene scritto nella cartella da cui si lancia il comando.
