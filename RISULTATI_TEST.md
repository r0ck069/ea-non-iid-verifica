# Risultati dei test: ea-non-iid-verifica-sorgenti v1.1.0-beta1

Prove eseguite l'11/10/2026. Il file di test è preparato per verificare la concordanza tra implementazioni, non per misurare una sorgente.

## 1. Oggetto

Confronto tra la pagina HTML, lo script Python e il programma ufficiale `ea_non_iid`, sullo stesso file di dati.

Il file di test **non è una sorgente di entropia**: è un file di byte generati dal sistema del Raspberry Pi con `/dev/urandom`. Serve solo a verificare che tre implementazioni diano lo stesso risultato.

## 2. Ambiente e impronte

| Elemento | Dettaglio |
|---|---|
| Raspberry Pi | Raspberry Pi 5, Raspberry Pi OS 64 bit (Debian 13), Python 3.13.5 |
| Programma NIST | `ea_non_iid`, commit `87c104d0ed4cbc96103e7b8b38d6f2c7e0a6b289`, verificato con `git rev-parse HEAD`; compilato con `make non_iid` dal Makefile del repository |
| Script Python | SHA-256 `c19937184793c6df5ebd3c1c5bb3bf1ae619752ac561335a03529ccca7cf4c57`, verificato sul PC Windows e sul Pi |
| Pagina HTML | SHA-256 `b1c981676bed8fad83242d5b92dcff528408c7d744ec98442084c4bb70dfe8eb`, verificato sul PC Windows |
| File di test | `test_1M.bin`, 1.000.000 byte; SHA-256 riportato dal programma NIST: `f8656330a5c76576a4294354188d2af78ef26a9f2fd1ce402e28ab537e7bbd69` |

Lettura del file: byte grezzi, 8 bit per simbolo, 256 valori distinti osservati, nessuna riga scartata.

## 3. Autotest contro i riferimenti NIST

| Riferimento | Simboli | Bit | Esito (Python sul Pi) |
|---|---|---|---|
| rand1_short | 10.000 | 1 | 11/11 uguali |
| rand4_short | 10.000 | 4 | 18/18 uguali |
| rand8_short | 10.000 | 8 | 18/18 uguali |

Esito: PASS. Lo stesso autotest ha dato PASS anche sul PC Windows (script), in Firefox su Windows (pagina) e in Chromium headless (prove di sviluppo dell'8 ottobre).

## 4. Confronto sul file di test

Valori per simbolo (H per simbolo):

| Stimatore | Programma NIST | Script Python | Pagina HTML |
|---|---|---|---|
| MCV | 7.878214 | 7.878214 | 7.878214 |
| t-Tuple | 7.353758 | 7.353758 | 7.353758 |
| LRS | 7.871595 | 7.871595 | 7.871595 |
| MultiMCW | 7.925647 | 7.925647 | 7.925647 |
| Lag | 6.638405 | 6.638405 | 6.638405 |
| MultiMMC | 7.926809 | 7.926809 | 7.926809 |
| LZ78Y | 7.928938 | 7.928938 | 7.928938 |

Valori per bit (H per bit, stringa di bit):

| Stimatore | Programma NIST | Script Python | Pagina HTML |
|---|---|---|---|
| MCV | 0.997958 | 0.997958 | 0.997958 |
| Collisione | 0.957371 | 0.957371 | 0.957371 |
| Markov | 0.999360 | 0.999360 | 0.999360 |
| Compressione | 0.939843 | 0.939843 | 0.939843 |
| t-Tuple | 0.935706 | 0.935706 | 0.935706 |
| LRS | 0.997159 | 0.997159 | 0.997159 |
| MultiMCW | 0.999111 | 0.999111 | 0.999111 |
| Lag | 0.999026 | 0.999026 | 0.999026 |
| MultiMMC | 0.998747 | 0.998747 | 0.998747 |
| LZ78Y | 0.998252 | 0.998252 | 0.998252 |

Per la pagina HTML, la prova in Chromium è dello sviluppo; la prova dell'11/10 è in Firefox su Windows, con risultati identici.

Esito: su 17 stimatori applicabili, i tre strumenti danno lo stesso valore a sei decimali.

## 5. Risultato finale

- Minimo sui simboli, come nel programma NIST: **6.638405 bit per simbolo** (stimatore Lag).
- Valore per bit: il risultato finale del programma NIST è min(H_original, 8 × H_bitstring) / 8 (formula verificata nel sorgente, commit `87c104d0`, righe 484-504): **0.829801**. H_bitstring = 0.935706 (t-Tuple) è un valore intermedio. Lo script e la pagina ottengono lo stesso valore per bit, e la soglia del verdetto si applica a quello.
- Verdetto con soglia 1,00: **NON UTILIZZABILE**. Il motivo è 0.829801 < 1,00.

## 6. Cosa significa, e cosa non significa

Verificato: pagina, script e programma NIST danno gli stessi risultati su questo file.

Verificato con il secondo file di controllo (sezione 9): il minimo per simbolo di circa 6,64 si ripete su un file diverso. Nel controllo il valore deriva dal test locale del predittore MultiMCW: con N = 1.000.000 e r = 4, Plocal = 0,0100, come da tabella 3 dell'Appendice G.2 del NIST. Per il minimo Lag dell'originale non è stato verificato il valore di r.

Lo stimatore Lag usa 128 sottopredittori, ciascuno basato su un lag da 1 a 128 (NIST SP 800-90B, §6.3.8). In una prima ipotesi avevo previsto un valore vicino a 8: non era un dato, e il calcolo sopra lo chiarisce.

Non verificato: che questo comportamento sia appropriato come verdetto per una sorgente reale. Il verdetto "non utilizzabile" dipende dal test locale su byte uniformi, quindi non dice nulla sulla qualità di una sorgente fisica.

## 7. Prove non eseguite

- Edge su Windows: il browser del primo test sul PC non è stato indicato, quindi non è attribuibile.
- Android, Safari e iPhone: non provati.
- Esecuzione dello script su Windows: non fatta (lo script è stato eseguito solo sul Pi).
- Confronto dei tre riferimenti incorporati con i file `cpp/selftest/refdata` del repository NIST, file per file: non fatto. I riferimenti sono quelli dichiarati nel passaggio di consegne come presi dal repository.
- Bit per simbolo 2, 3, 5, 6 e 7: nessun riferimento disponibile per l'autotest.
- File di sorgenti reali (cucina, pioggia, acqua, radio): non ancora analizzati con questo strumento.
- Test di salute RCT e APT dello standard: non implementati nello strumento.

## 8. Prossimo passo di verifica

1. Secondo file di controllo: eseguito il 16/10/2026, vedi sezione 9.
2. Il confronto dei tre riferimenti incorporati con il repository NIST.
3. Analisi di un file reale di sorgente, con il confronto tra programma NIST e strumento.

## 9. Secondo file di controllo (16/10/2026)

- File: `test_1M_controllo.bin`, 1.000.000 byte da `/dev/urandom`, creato il 16/10/2026 alle 21:20 (ora locale del Pi) nella cartella `~/verifica-ea`. SHA-256 `0d28ea583b1b2341ed445739133695c9adcda5d16684a7557eea4c2629875f9c`.
- Il file originale `test_1M.bin` è invariato: SHA-256 `f8656330a5c76576a4294354188d2af78ef26a9f2fd1ce402e28ab537e7bbd69`.
- Programma NIST, stessi parametri dell'originale (`./ea_non_iid -i -v <file> 8`, nella cartella `cpp/` del repository): H_original 6.638383 (MultiMCW, valore letterale); H_bitstring 0.900478 (Compressione); minimo 6.638383; per bit 0.829798.
- Script Python (`--bps 8`): stessi 17 valori del programma NIST a sei decimali; autotest PASS (47/47); verdetto non utilizzabile. Rapporto: `rapporto_test_1M_controllo.txt`.
- Il minimo cambia stimatore tra i due file (Lag nell'originale, MultiMCW nel controllo), ma resta circa 6,64 in entrambi. Gli altri stimatori variano di più: per esempio LRS per simbolo da 7,871595 a 7,939199, Lag da 6,638405 a 7,910781, Compressione per bit da 0,939843 a 0,900478.
- Nel controllo il test locale di MultiMCW dà Plocal = 0,010038 per r = 4; la tabella 3 del NIST dà 0,0100 per r = 4 con N = 1.000.000.
