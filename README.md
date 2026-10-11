# ea-non-iid-verifica-sorgenti

Strumento per stimare la min-entropia di un file grezzo di una sorgente di entropia, con i dieci stimatori dello standard NIST SP 800-90B per sorgenti non-IID.

È disponibile in due forme che danno gli stessi risultati:

- una pagina HTML standalone, che funziona offline nel browser;
- uno script Python, con la sola libreria standard.

Versione: **v1.1.0-beta1**. Stato: beta. Non è una certificazione.

---

## A cosa serve

Prima di usare una sorgente di entropia (un audio, una foto, un fenomeno fisico) per costruire un seme, bisogna sapere quanta entropia contiene davvero. Non basta il suo nome: "pioggia", "radio" o "dado" non dicono nulla finché i dati non sono misurati.

Questo strumento fa la prima misura, e funziona da filtro: se una sorgente non supera la soglia scelta, si scarta.

Lo strumento non prova che una sorgente sia buona e non crea entropia. Dice quanto è prevedibile la sequenza di dati che gli si fornisce, nel caso peggiore secondo lo standard.

## Che cos'è un file non-IID

IID vuol dire "indipendenti e identicamente distribuiti": ogni valore non dipende dai precedenti, e tutti i valori seguono la stessa distribuzione. Molte stime semplici partono da questa ipotesi.

Una sorgente reale quasi mai la rispetta. Un audio ha correlazioni tra campioni vicini, perché il suono ha continuità. Una foto ha strutture spaziali. Un sensore può derivare nel tempo. Un dado o una moneta possono avere un bias. Un file è non-IID quando non si può assumere nessuna di queste proprietà.

Lo standard NIST SP 800-90B prevede due famiglie di stime: quelle per sorgenti IID e quelle per sorgenti non-IID. Le seconde non fanno l'ipotesi di indipendenza: cercano regolarità di tipi diversi (valori frequenti, ripetizioni, dipendenze tra valori, capacità di predizione). Questo strumento applica la seconda famiglia.

## Quale problema risolve

Il dubbio concreto è: questo file ha abbastanza entropia per contribuire a un seme?

La misura standard è la **min-entropia**: il numero di bit per simbolo che si può considerare garantito nel caso peggiore, cioè considerando lo stimatore che prevede meglio la sequenza. Un valore di 1 bit per bit vuol dire che, nel caso peggiore, la sequenza è difendibile come casuale per quel numero di bit.

Lo strumento dà questa misura in modo riproducibile e controllabile, con un programma pubblico. Non dipende da un giudizio a occhio né da una stima fatta con un programma che nessuno può verificare.

## Come lo fa

1. **Legge il file.** Accetta byte grezzi, testo con bit, esadecimali o decimali, e file WAV (PCM a 8, 16, 24, 32 bit, o in virgola mobile a 32 o 64 bit; viene letto il blocco dati).
2. **Estrae i simboli.** Un simbolo per byte, oppure per gruppi di bit da 1 a 8 bit per simbolo, usando i bit meno significativi.
3. **Esegue l'autotest.** Ricalcola i dieci stimatori su tre file di riferimento del NIST incorporati e confronta i valori con quelli ufficiali (tolleranza 1e-6). Se l'autotest fallisce, il verdetto è rosso.
4. **Calcola i dieci stimatori.** Sui simboli e, per simboli a più bit, sulla stringa di bit. Il valore finale è il minimo di tutti.
5. **Dà il verdetto.** È "utilizzabile" solo se sono vere tre condizioni: l'autotest è superato; ci sono almeno 1.000.000 di campioni; la min-entropia per bit è maggiore o uguale alla soglia scelta (predefinita 1,00).
6. **Produce il rapporto** in formato testo e, se richiesto, la **sorgente estratta** (un file .bin con un byte per simbolo).

I dieci stimatori sono: MCV, Collisione, Markov, Compressione, t-Tuple, LRS, MultiMCW, Lag, MultiMMC, LZ78Y. Collisione, Markov e Compressione si applicano ai dati binari, come nello standard.

## Come leggere il verdetto

- **Verde, "SORGENTE UTILIZZABILE"**: la sorgente supera la soglia scelta sul file. Non vuol dire che sia adatta a un seme: è solo il primo filtro.
- **Rosso, "SORGENTE NON UTILIZZABILE"**: la sorgente non supera la soglia, oppure il calcolo non è affidabile (campioni insufficienti, autotest fallito). Il rapporto indica il motivo.

La soglia è una scelta. Il valore predefinito, 1,00, è il più severo. Si può abbassare, ma la scelta va motivata e documentata.

## Uso

**Pagina HTML.** Aprire il file in un browser moderno. L'autotest parte da solo. Caricare il file, scegliere formato e bit per simbolo, premere Calcola. Il download della sorgente estratta e del rapporto avviene solo su richiesta.

**Script Python** (Python 3, nessuna dipendenza esterna):

```
python3 ea_verifica_sorgenti.py FILE --bps 8
python3 ea_verifica_sorgenti.py FILE --fmt bits --bps 1
python3 ea_verifica_sorgenti.py --selftest
```

Opzioni: `--fmt raw|bits|hex|dec`, `--bps 1..8`, `--offset N`, `--strict`, `--min-h 1.0`, `--out-dir DIR`, `--selftest`.

Codici di uscita: 0 per verde, 1 per rosso, 2 per errore di lettura o di calcolo.

Tempi indicativi, dalle prove di sviluppo: circa 35 secondi per 40.000 campioni; circa 9 minuti per 1.000.000 di bit. Dipendono dalla macchina.

## Posizione nella suite

- È il **primo filtro** della catena. Si applica alla sorgente grezza prima di concatenarla con le altre.
- Il riferimento ufficiale resta il programma `ea_non_iid`, compilato dal repository NIST al commit `87c104d0ed4cbc96103e7b8b38d6f2c7e0a6b289`. La pagina e lo script sono un'implementazione che deve dare gli stessi risultati. La verifica di questa concordanza è in `RISULTATI_TEST.md`.
- A valle, lo strumento di estrazione (EntropyPipeline) applica gli altri controlli e l'estrazione.

## Contenuto del repository

- `ea-non-iid-verifica-sorgenti-v1_1_0-beta1.html`: pagina standalone.
- `ea_verifica_sorgenti.py`: script Python, stessi calcoli.
- `README.md`: questo file.
- `AVVERTENZE.md`: limiti e uso corretto.
- `RISULTATI_TEST.md`: prove e concordanza con il programma NIST.
- `CHANGELOG.md`: storia delle versioni.
- `SHA256SUMS.txt`: impronte dei file, da generare alla pubblicazione.

## Stato di verifica

Verificato: concordanza con `ea_non_iid` su un file di test, per 17 stimatori su 17 applicabili; autotest superato in Chromium headless, in Firefox su Windows e con lo script Python sul Raspberry Pi.

Non verificato: Android, Safari e iPhone; Edge su Windows; il comportamento su file di sorgenti reali; i bit per simbolo 2, 3, 5, 6 e 7, per cui non ci sono riferimenti. Dettagli in `RISULTATI_TEST.md` e limiti in `AVVERTENZE.md`.

## Riferimenti

- NIST SP 800-90B, *Recommendation for the Entropy Sources Used for Random Bit Generation*.
- Repository ufficiale: https://github.com/usnistgov/SP800-90B_EntropyAssessment

## Licenza

MIT. Il testo della licenza è nel file `LICENSE` del repository.
