# Changelog

Tutte le modifiche rilevanti a questo progetto sono documentate in questo file.
Le date sono in formato AAAA-MM-GG.

## [1.1.0-beta1] - 2026-10-11

### Aggiunto

- Pagina HTML standalone `ea-non-iid-verifica-sorgenti-v1_1_0-beta1.html`: dieci stimatori SP 800-90B per sorgenti non-IID, autotest contro tre file di riferimento NIST, verdetto secondo una soglia impostabile, download della sorgente estratta e del rapporto.
- Script `ea_verifica_sorgenti.py`: stessi calcoli della pagina, sola libreria standard di Python 3, uscita con codici 0 (utilizzabile), 1 (non utilizzabile), 2 (errore).
- `README.md`, `AVVERTENZE.md`, `RISULTATI_TEST.md`.

### Verificato

- Concordanza con il programma ufficiale `ea_non_iid` (commit `87c104d0ed4cbc96103e7b8b38d6f2c7e0a6b289`) su un file di test di 1.000.000 di byte: 17 stimatori su 17 applicabili uguali a sei decimali, minimo finale uguale (6.638405).
- Concordanza tra script Python e pagina HTML sullo stesso file.
- Autotest superato con lo script Python sul Raspberry Pi, in Firefox su Windows e in Chromium headless.

### Note

- Costruzione dei file: 2026-10-08 (secondo il passaggio di consegne).
- Non testato su Android, Safari e iPhone.
- Non ancora testato su file di sorgenti reali.

### Impronte SHA-256 (versione 1.1.0-beta1)

- `ea-non-iid-verifica-sorgenti-v1_1_0-beta1.html`: `b1c981676bed8fad83242d5b92dcff528408c7d744ec98442084c4bb70dfe8eb`
- `ea_verifica_sorgenti.py`: `c19937184793c6df5ebd3c1c5bb3bf1ae619752ac561335a03529ccca7cf4c57`
