# Avvertenze e limiti

Leggere questo file prima di usare lo strumento o di interpretarne il risultato.

## 1. Non è una certificazione

Un risultato verde non certifica la sorgente e non dice che sia adatta a generare un seme. Lo strumento è un primo filtro: stabilisce se i dati superano una soglia, niente di più.

## 2. Misura i dati, non il fenomeno

Il calcolo riguarda il file che gli si fornisce. Un file può essere buono perché la sorgente è buona, oppure perché il sistema lo ha elaborato (filtri, compressione, normalizzazione del volume, riduzione del rumore). Lo strumento non distingue i due casi. Chi valuta una sorgente deve conoscere come è stata registrata.

## 3. Il verdetto dipende dalla soglia

Il valore predefinito, 1,00 bit per bit, è il più severo. Una soglia diversa cambia il verdetto sullo stesso file. La soglia va scelta prima di guardare il risultato e va dichiarata insieme al risultato.

## 4. Rosso non prova che la sorgente sia cattiva

Il rosso può dipendere da campioni insufficienti, da un formato di lettura sbagliato, dalla soglia o da un autotest fallito. Il rapporto indica il motivo. Va letto prima di concludere qualsiasi cosa sulla sorgente.

## 5. Servono almeno un milione di campioni

Sotto 1.000.000 di simboli il verdetto è rosso, perché la stima dello standard non è affidabile su campioni più piccoli. Con meno campioni il calcolo può comunque girare, ma il risultato non si usa.

## 6. L'autotest copre poco

L'autotest confronta il calcolo con tre file di riferimento NIST corti (10.000 simboli, da 1, 4 e 8 bit). Non copre tutti i casi. Per i bit per simbolo 2, 3, 5, 6 e 7 non ci sono riferimenti, quindi i valori vanno considerati non verificati da questo strumento.

## 7. Nessun test di salute continuo

Lo standard prevede test di salute in tempo reale (RCT e APT) per controllare una sorgente mentre lavora. Questo strumento non li implementa: misura dati già raccolti. Un guasto che compare durante la raccolta non viene segnalato.

## 8. Non valuta il condizionamento

Lo strumento non esamina né approva funzioni di hash, estrattori o altri trattamenti applicati ai dati. Il valore misurato è quello dei dati grezzi.

## 9. Non sostituisce ea_non_iid

Il riferimento ufficiale resta il programma NIST compilato dal repository al commit indicato nel README. La pagina e lo script sono un'implementazione che va verificata a ogni nuova versione, con il confronto descritto in `RISULTATI_TEST.md`.

## 10. Tempi lunghi

Il calcolo su un milione di simboli richiede minuti, e può richiedere più tempo su macchine lente. La pagina lo esegue in un Web Worker, quindi non si blocca, ma l'attesa è normale. Non interrompere il calcolo se non è necessario.

## 11. Dati e rete

La pagina HTML non invia dati in rete: la sua politica di sicurezza blocca le connessioni, e il download avviene solo su richiesta. Lo script Python lavora in locale. I file caricati restano sul dispositivo, ma chi li verifica deve decidere se sono sensibili prima di condividere il rapporto.

## 12. Il file estratto non è un seme

Il file `sorgente_estratta_*.bin` è un output intermedio: contiene un byte per simbolo e non va usato direttamente come chiave, né come seme, né come sorgente unica.

## 13. Usi che non sono previsti

Lo strumento non valuta la sicurezza di un dispositivo, non verifica l'integrità di un generatore hardware e non certifica la conformità a standard di validazione.
