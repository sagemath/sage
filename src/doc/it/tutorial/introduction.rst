************
Introduzione
************
Questo tutorial dovrebbe richiedere circa 3/4 ore per 
una lettura completa. Lo si può leggere in versione HTML o PDF.

Nonostante molto in Sage sia implementato usando Python, la conoscenza di Python
non è un prerequisito per la lettura di questo tutorial. Per chi volesse imparare
il Python (un linguaggio molto divertente!) allo stesso tempo, ci sono molte risorse 
eccellenti e libere per farlo tra le quali [PyT]_ e [PyB]_.
Se si vuole solo provare velocemente Sage, questo tutorial è il punto di partenza adatto.
Per esempio:

::

    sage: 2 + 2
    4
    sage: factor(-2007)
    -1 * 3^2 * 223
    
    sage: A = matrix(4,4, range(16)); A
    [ 0  1  2  3]
    [ 4  5  6  7]
    [ 8  9 10 11]
    [12 13 14 15]
    
    sage: factor(A.charpoly())
    x^2 * (x^2 - 30*x - 80)
    
    sage: m = matrix(ZZ,2, range(4))
    sage: m[0,0] = m[0,0] - 3
    sage: m
    [-3  1]
    [ 2  3]
    
    sage: E = EllipticCurve([1,2,3,4,5]); 
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 
    over Rational Field
    sage: E.anlist(10)
    [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
    sage: E.rank()
    1
    
    sage: k = 1/(sqrt(3)*I + 3/4 + sqrt(73)*5/9); k
    36/(20*sqrt(73) + 36*I*sqrt(3) + 27)

    sage: N(k)
    0.165495678130644 - 0.0521492082074256*I
    sage: N(k,30)      # 30 "bits"
    0.16549568 - 0.052149208*I
    sage: latex(k)
    \frac{36}{20 \, \sqrt{73} + 36 i \, \sqrt{3} + 27}

Installazione
=============

Se non si ha Sage installato su un computer e si vogliono solamente
provare alcuni comandi, si può usare online all'indirizzo https://sagecell.sagemath.org/.

Si veda la Sage Installation Guide nella sezione documentazione della homepage
di Sage [Sage]_ per istruzioni sull'installazione di Sage sul proprio computer.
Qui vengono fatti solamente due commenti.


#. Il file di download di Sage arrive con le "batterie incluse".
   In altre parole, nonostante Sage usi Python, IPython, PARI, GAP, 
   Singular, Maxima, NTL, GMP e così via, non è necessario installarli
   separatamente siccome sono incluse con la distribuzione di Sage.
   Comunque, per usare certe feature di \sage, ad es. Macaulay o KASH, 
   bisogna installare il pacchetto opzionale Sage che interessa o almeno
   avere i programmi in questioni gia installati sul proprio computer.
   Macaulay e KASH sono pacchetti di Sage (per una lista dei pacchetti 
   opzionali disponibili, digitare "sage -optional", o sfogliare la pagina
   "Download" sul sito web di Sage.

#. Le versioni binarie precompilate di Sage (che si trovano sul sito web di 
   Sage) possono essere più facili e più veloci da installare invece che la 
   versione da codice sorgente. Basta solo spacchettare il file e eseguire "sage".

Modi di usare Sage
==================

Sage si può usare in molti modi.


-  **Interfaccia grafica del notebook:** vedere la sezioni sul 
   Notebook nel manuale di riferimento e la :ref:"sezione-notebook" sotto,

-  **Linea di comando interattiva:** vedere :ref:'capitolo-shell_interattiva',

-  **Programmi:** scrivendo programmi interpretati e compilati in Sage (vedere
   :ref:'sezione-loadattach' e :ref:'sezione-compilazione'), e

-  **Scripts:** scrivendo degli script autosufficienti che usino la libreria 
   Sage (vedere :ref:'sezione-autosufficienti').


Obiettivi di lungo periodo per Sage
===================================

-  **Utile**: il pubblico per Sage il quale sage è stato pensato sono gli 
   studenti di matematica (dalla scuola superiore all'università), gli insegnanti
   e i ricercatori in matematica. Lo scopo è di fornire software che possa essere
   usato per esplorare e sperimentare le costruzioni matematiche in algebra,
   geometria, teoria dei numeri, calcolo, calcolo numerico, ecc. Sage aiuta a
   rendere più facile la sperimentazione interattiva con gli oggetti matematici.

-  **Efficiente:** essere veloce. Sage usa del software maturo e altamente
   ottimizzato come GMP, PARI, GAP e NTL e così è molto veloce con certe
   operazioni.

-  **Libero e open source:** il codice sorgente deve essere liberamente disponibile
   e leggibile, così che gli utenti possano capire cosa stia facendo veramente il 
   sistema e possano estenderlo più facilmente. Così come i matematici acquisiscono
   una comprensione più profonda di un teorema leggendo attentamente o almeno scorrendo
   velocemente la dimostrazione, le persone che fanno calcoli dovrebbero essere capaci
   di capire come funzionano i calcoli leggendo il codice sorgente documentato. Se
   si usa Sage per fare calcoli in un articolo che si pubblica, si può essere rassicurati
   dal fatto che i lettori avranno sempre libero accesso a Sage e a tutto il suo codice
   sorgente ed è persino concesso di archiviare la versione di Sage che si è utilizzata.

-  **Facile da compilare:** Sage dovrebbe essere facile da compilare dal sorgente per
   gli utenti Linux, macOS e Windows. Questo garantisce maggiore flessibilità agli utenti
   di modificare il sistema.

-  **Cooperazione:** Fornire un interfaccia robusta alla maggior parte degli altri sistemi
   di algebra computazionale, compresi: PARI, GAP, Singular, Maxima, KASH, Magma, Maple e
   Mathematica. Sage è pensato per unificare e estendere il software matematico esistente.

-  **Ben documentato:** tutorial, guida alla programmazione, manuale di riferimento e 
   how to con numerosi esempi e discussioni della matematica sottostante.

-  **Amichevole verso l'utente:** dovrebbe essere facile capire quale funzionalità è
   fornita per un dato oggetto e guardare la documentazione e il codice sorgente.
   Bisogna anche raggiungere un alto livello di supporto agli utenti.


..  [PyB] (en) The Python Beginner's Guide,
    https://wiki.python.org/moin/BeginnersGuide

..  [PyT] (en) The Python Tutorial,
    https://docs.python.org/3/tutorial/

..  [Sage] (en) Sage, https://www.sagemath.org
