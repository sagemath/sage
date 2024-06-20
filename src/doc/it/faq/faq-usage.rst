.. _chapter-faq-usage:

================
FAQ: Uso di Sage
================


Come posso partire?
"""""""""""""""""""

Puoi provare Sage senza dover scaricare nulla:

* **CoCalc™:** Vai su https://cocalc.com e crea un account gratis.

  Se accedi al sito, avrai accesso all'ultima versione di Sage e molti
  altri programmi

  Ricorda che questo sito è un servizio commerciale indipendente
  da Sage.

* **Sage cell:** Una singola versione di Sage, disponibile per
  fare una computazione alla volta. https://sagecell.sagemath.org/


Per scaricare una distribuzione binaria **precompilata** di Sage
vai al link http://www.sagemath.org/download.html e clicca sul link
relativo al file binario per il tuo sistema operativo.

Il **codice sorgente** di Sage è anche disponibile al download:
vai al link http://www.sagemath.org/download-source.html
per scaricare l'archivio TAR di qualunque rilascio di Sage.

Le sessioni di lavoro Notebook di Sage sono eseguite all'interno di
un browser web. Puoi lanciare il Notebook di sage attraverso il
seguente comando purché ``sage`` sia nella variabile ``PATH``

.. CODE-BLOCK:: shell-session

    $ sage -notebook


Quali sono i prerequisiti di Sage?
""""""""""""""""""""""""""""""""""

La maggior parte delle dipendenze sono incluse all'interno di Sage.
Nella maggior parte dei casi puoi scaricare il binario precompilato
ed usarlo senza dover installare alcun pacchetto dipendente. Se usi
Windows avrai bisogno di installare
`VirtualBox <http://www.virtualbox.org>`_,
che puoi scaricare dal link http://www.virtualbox.org/wiki/Downloads.
Dopo aver installato VirtualBox devi scaricare una distribuzione di
Sage per VirtualBox al link
http://www.sagemath.org/download-windows.html.
Segui bene le istruzioni che trovi a quella pagina. Poi puoi lanciare
la macchina virtuale Sage usando il software VirtualBox.

Puoi scaricare il codice sorgente completo di Sage per compilarlo sul
tuo sistema Linux o Mac OS X. Sage si trova in una cartella isolata e
non va ad interferire col sistema in cui si trova. Viene dato con
tutto il necessario per lo sviluppo, il codice sorgente, tutte le
dipendenze ed il changelog (cioè la lista delle modifiche operate)
completo. Sui sistemi Linux come Debian/Ubuntu puoi dover installare
il pacchetto ``build essential`` ed il processore di macro ``m4``. Il
tuo sistema deve disporre di un compilatore C funzionante se vuoi
compilare Sage da codice sorgente. Su Debian/Ubuntu puoi installare
questi prerequisiti come segue::

    sudo apt-get install build-essential m4

Se hai un sistema multiprocessore puoi scegliere una
compilazione parallela di Sage. Il comando ::

    export MAKE='make -j8'

abiliterà 8 threads per quelle parti della compilazione che supportano
il parallelismo. Al posto del numero 8 metti il numero di
processori/core del tuo sistema.


Come posso far riconoscere la mia attuale installazione di Tcl/Tk all'interprete Python di Sage?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Potresti avere la libreria Tcl/Tk installata e l'interprete Python del
tuo sistema la riconosce ma l'interprete Python di Sage no. Ad oggi (2023)
solitamente non c'è bisogno di compilare l'interprete Python di Sage, ma se ne
hai bisogno, segui queste istruzioni. Controlla di aver installato la libreria
di sviluppp Tcl/Tk. Su Ubuntu lancia, da riga di comando::

    sudo apt-get install tk8.5-dev

o qualcosa di simile. Poi reinstalla l'interprete Python di Sage con::

    make python3-clean python3-uninstall && make python3

Questo aggancerà automaticamente la libreria Tcl/Tk.
Dopo aver reinstallato correttamente l'interprete Python di Sage,
lancia i seguenti comandi dall'interfaccia a riga di comando di Sage::

    import _tkinter
    import Tkinter

Se non ti viene segnalato alcun errore di :class:`ImportError`
allora il problema è risolto.


Come faccio ad importare Sage in uno script Python?
"""""""""""""""""""""""""""""""""""""""""""""""""""

Puoi importare Sage in uno script Python come faresti con una libreria.
La cosa a cui fare attenzione è che devi lanciare quello script Python
usando l'interprete Python interno a Sage
(versione 3.7.x per Sage 9.2).
Per importare Sage metti la seguente istruzione in
cima al tuo script Python::

    from sage.all import *

Quando poi esegui il tuo script devi lanciare Sage con l'opzione
``-python`` che farà sì che venga eseguito dalla versione
dell'interprete interna a Sage. Ad esempio, se Sage è nella tua
variabile d'ambiente ``PATH``, puoi scrivere::

    sage -python /path/to/my/script.py

Un altro modo è scrivere uno script Sage e lanciarlo usando Sage
stesso. Uno script Sage ha estensione  ``.sage`` ed è simile ad uno
script Python ma utilizza funzioni e comandi specifici di Sage.
Puoi poi lanciare tale script Sage in questo modo::

    sage /path/to/my/script.sage

Questo si occuperà di caricare le variabili d'ambiente necessarie
ed eseguire gli import di default al posto tuo.


Come posso ricaricare uno script Python in una sessione di Sage?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Puoi caricare uno script Python in una sessione Sage usando il comando
**load**. Ad esempio possiamo usare Sage per importare un file di
nome "simple.py" con::

    load("simple.py")

e ripetere questo comando ogni volta che cambiamo il file.
Invece digitando::

    attach("simple.py")

ogni cambiamento al file verrà automaticamente aggiornato
anche in Sage.


Posso usare Sage con la versione 3.x di Python?
"""""""""""""""""""""""""""""""""""""""""""""""

Dalla versione 9.0 del Gennaio 2020, SageMath utilizza Python 3.


Ho usato Debian/Ubuntu per installare la versione 3.0.5 di Sage ed essa sta dando un sacco di errori. Cosa posso fare?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La versione di Sage distribuita con ``apt-get`` in Debian e Ubuntu
(tipo la 3.0.5) è molto vecchia. Nessuno ha ancora avuto tempo di
aggiornare la versione di Sage per Debian/Ubuntu.
Qualunque aiuto in tal senso sarà molto apprezzato.
Dovresti scaricare la versione più recente di Sage dal
`link di download <http://www.sagemath.org/download.html>`_ del sito
web di Sage. Se vuoi aiutarci ad aggiornare la versione di Sage per
Debian/Ubuntu manda un'email alla mailing list
`sage-devel <http://groups.google.com/group/sage-devel>`_.


Faccio meglio ad usare la versione ufficiale o quella di sviluppo?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Ti consigliamo di usare la più recente versione ufficiale di Sage.
Le versioni di sviluppo sono spesso annunciate sulle mailing list
`sage-devel <http://groups.google.com/group/sage-devel>`_ e
`sage-release <http://groups.google.com/group/sage-release>`_.
Una maniera facile di aiutare con lo sviluppo di Sage è scaricare
l'ultima versione di sviluppo, compilarla sul suo sistema,
lanciare tutti i doctest e segnalare qualunque errore di compilazione
o qualunque fallimento nei doctest.


È difficile imparare Sage?
""""""""""""""""""""""""""

Le funzionalità di base di Sage dovrebbero risultare facili da imparare
quanto le basi di Python. Molti tutorial sono disponibili in rete
per aiutarti ad imparare Sage. Per trarre il massimo da Sage
ti consigliamo di imparare qualche elemento del linguaggio Python.
Segue una lista, incompleta, di risorse su Python.
Altre risorse possono essere trovate cercando sul web.

* `Building Skills in Python <http://homepage.mac.com/s_lott/books/python.html>`_
  di Steven F. Lott
* `Dive into Python <http://www.diveintopython.net>`_
  di Mark Pilgrim
* `How to Think Like a Computer Scientist <http://www.openbookproject.net/thinkCSpy>`_
  di Jeffrey Elkner, Allen B. Downey, and Chris Meyers
* `Official Python Tutorial <http://docs.python.org/tutorial>`_
* `Python <http://www.python.org>`_ home page e
  `Python standard documentation <http://docs.python.org>`_


Posso fare X in Sage?
"""""""""""""""""""""

Ti consigliamo di usare l'autocompletamento di Sage con il tasto TAB.
Ti basta digitare qualche carattere, premere TAB e vedere se il
comando che vuoi compare nella lista di autocompletamento.
Se hai un comando che si chiama ``mycmd``,
allora digitalo e premi TAB per visualizzare la lista di funzionalità
che sono supportate da quel comando. Per leggere la documentazione di
``mycmd`` scrivi ``mycmd?`` poi premi Invio e protrai leggerla.
In modo similare, digitando ``mycmd??`` e poi Invio potrai
visualizzare il codice sorgente di tale comando.
Ti consigliamo anche di eseguire ricerche nel codice sorgente e
nella documentazione di Sage. Per eseguire ricerche nel codice
sorgente di Sage usa il comando
``search_src("<search-keyword>")``
mettendo al posto di ``<search-keyword>`` le parole chiave che
vuoi cercare.
Analogamente puoi effettuare ricerche nella documentazione di
Sage usando il comando:
``search_doc("<search-keyword>")``.


Cosa fa esattamente Sage quando digito "0.6**2"?
""""""""""""""""""""""""""""""""""""""""""""""""

Quando scrivi "0.6**2" in Python, ti viene restituito qualcosa tipo
0.35999999999999999. Ma quando fai lo stesso in Sage ti viene
restituito 0.360000000000000. Per capire perché Python si comporta in
questo modo vedi il
`Python Tutorial <http://docs.python.org/tutorial/floatingpoint.html>`_,
soprattutto il capitolo
"Aritmetica floating-point: caratteristiche e limiti".
Ciò che Sage fa è preprocessare l'input e trasformarlo come segue::

    sage: preparse("0.6**2")
    "RealNumber('0.6')**Integer(2)"

Così che ciò che viene *effettivamente* eseguito è::

    RealNumber('0.6')**Integer(2)

Gli sviluppatori Sage (in pratica Carl Witty) decisero che i numeri
floating-point di Sage dovessero, di default, stampare solo il numero
di cifre decimali corrette, quando possibile, così da evitare il
problema che ha Python. Questa decisione ha i suoi pro e contro.
Nota che ``RealNumber`` e ``Integer`` sono specifici di Sage,
quindi non puoi digitare quanto sopra nell'interprete Python ed
aspettarti che funzioni, se prima non hai eseguito delle istruzioni
di import quali::

    from sage.all import RealNumber, Integer, preparse


Perché il comando "history" di Sage è diverso da quello di Magma?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Nell'uso di Sage non disponi di una funzionalità dell'interfaccia a
riga di comando di Magma. In Magma, se immetti una linea recuperata
dalla "history" (cioè dall'elenco dei comandi digitati precedentemente
che viene automaticamente memorizzato) con il tasto "freccia in su'"
e poi premi "freccia in giu'", viene recuperata anche la linea
successiva nell'elenco. Questa funzionalità ti permette di recuperare
dalla "history" tante righe consecutive quante vuoi.
Ma Sage non ha una funzionalità simile: la riga di comando
`IPython <http://ipython.scipy.org>`_ utilizza la libreria "readline"
(via pyreadline), che evidentemente non supporta questa funzionalità.
Magma ha una sua propria libreria personalizzata simile alla
"readline" che invece supporta questa funzionalità.
(Dal momento che moltissime persone hanno richiesto questa
funzionalità, se qualcuno trovasse un modo per implementarla
sarebbe il benvenuto!)


Ho problemi di tipo nell'utilizzo da Sage di SciPy, cvxopt e NumPy.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Stai usando da Sage le librerie SciPy, cvxopt e NumPy e hai degli
errori tipo::

    TypeError: function not supported for these types, and can't coerce safely to supported types.

Quando digiti numeri in Sage, il preprocessore li converte in un
anello base, come puoi vedere facendo::

    sage: preparse("stats.uniform(0,15).ppf([0.5,0.7])")
    "stats.uniform(Integer(0),Integer(15)).ppf([RealNumber('0.5'),RealNumber('0.7')])"

Sfortunatamente il supporto che NumPy fornisce a questi tipi avanzati di
Sage, quali ``Integer`` o ``RealNumber``
(numeri reali di precisione arbitraria), non è del 100%.
Per risolvere ridefinisci ``Integer`` e/o ``RealNumber`` per cambiare
il comportamento del preprocessore di Sage così che i decimali
scritti vengano registrati come tipi float di Python anziché RealNumber
di Sage e gli interi scritti siano registrati come tipi int di Python
anziché Integer di Sage. Ad esempio::

    sage: RealNumber = float; Integer = int
    sage: from scipy import stats
    sage: stats.ttest_ind([1,2,3,4,5], [2,3,4,5,.6])
    Ttest...Result(statistic=...0.0767529..., pvalue=...0.940704...)
    sage: stats.uniform(0,15).ppf([0.5,0.7])
    array([  7.5,  10.5])

In alternativa sii esplicito circa il tipo di dato, ad esempio::

    sage: from scipy import stats
    sage: stats.uniform(int(0),int(15)).ppf([float(0.5),float(0.7)])
    array([  7.5,  10.5])

Come terza alternativa puoi usare i suffissi semplici::

    sage: from scipy import stats
    sage: stats.uniform(0r,15r).ppf([0.5r,0.7r])
    array([  7.5,  10.5])

Puoi anche disabilitare il preprocessore nel tuo codice tramite il
comando ``preparse(False)``. Puoi lanciare IPython da solo dalla riga
di comando con ``sage -ipython``, cosa che non precarica niente di
specifico di Sage. O ancora puoi cambiare il linguaggio di
sessione (Notebook language) in "Python".


Come faccio a salvare un oggetto così che non devo ridigitarlo ogni volta che apro un foglio di lavoro (worksheet)?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

I comandi ``save`` e ``load`` rispettivamente registrano e caricano un
oggetto. Nella sessione di lavoro Notebook la variabile ``DATA`` è la
collocazione dello spazio di salvataggio del foglio di lavoro
(worksheet). Per registrare l'oggetto ``my_stuff`` in un foglio di
lavoro puoi digitare::

    save(my_stuff, DATA + "my_stuff")

e, per ricaricarlo, ti basta digitare::

    my_stuff = load(DATA + "my_stuff")


Sage contiene una funzione simile alla "ToCharacterCode[]" di Mathematica?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Potresti voler convertire caratteri ASCII come "Big Mac" nel
corrispondente codice numerico per ulteriori elaborazioni. In Sage e
Python puoi usare ``ord``. Ad esempio::

    sage: list(map(ord, "abcde"))
    [97, 98, 99, 100, 101]
    sage: list(map(ord, "Big Mac"))
    [66, 105, 103, 32, 77, 97, 99]

Come posso scrivere le multiplicazioni in modo implicito come in Mathematica?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sage ha una funzione che lo rende possibile::

    sage: implicit_multiplication(True)
    sage: x 2 x  # not tested
    2*x^2
    sage: implicit_multiplication(False)

Questo viene preprocessato da Sage in codice per Python.
Potrebbe non funzionare in situazioni complicate.
Per vedere cosa il preprocessore fa::

    sage: implicit_multiplication(True)
    sage: preparse("2 x")
    'Integer(2)*x'
    sage: implicit_multiplication(False)
    sage: preparse("2 x")
    'Integer(2) x'

Vai al link https://wiki.sagemath.org/sage_mathematica per maggiori
informazioni su Mathematica vs. SageMath.

Posso far eseguire in automatico a Sage dei comandi all'accensione?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sì, ti basta creare un file ``$HOME/.sage/init.sage`` ed esso sarà
eseguito ogni volta che lanci Sage. Questo presuppone che la
variabile ambiente di Sage ``DOT_SAGE`` punti alla cartella nascosta
``$HOME/.sage``, cosa che avviene di default.


Quando compilo Sage il mio computer fa beep e si spegne o si blocca.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Compilare sage è piuttosto faticoso per il processore del computer.
Il comportamento suddetto di solito indica che il computer si è
surriscaldato. In molti casi questo si può risolvere pulendo il
ventilatore del processore del computer ed assicurando adeguata
areazione al computer. Puoi chiedere al tuo amministratore di sistema
o ad un tecnico di provvedere, qualora tu non l'abbia mai fatto.
Questa manutenzione del computer, se non fatta da persone preparate,
potrebbe anche danneggiare il computer.

Per gli utenti Linux, se pensi che la compilazione fallisca per un
problema di risorse di macchina, una soluzione potrebbe essere di
modificare il file ``/etc/inittab`` per far partire Linux al
runlevel 3. Tale file di solito contiene qualcosa del tipo::

    #   0 - halt (Do NOT set initdefault to this)
    #   1 - Single user mode
    #   2 - Multiuser, without NFS (The same as 3, if you do not have
    #   networking)
    #   3 - Full multiuser mode
    #   4 - unused
    #   5 - X11
    #   6 - reboot (Do NOT set initdefault to this)
    #
    id:5:initdefault:

Questo fa sì che la tua distribuzione Linux parta con la schermata
di login grafico. Commenta la linea ``id:5:initdefault:`` e
aggiungi la linea ``id:3:initdefault:``, così da aver qualcosa come::

    #   0 - halt (Do NOT set initdefault to this)
    #   1 - Single user mode
    #   2 - Multiuser, without NFS (The same as 3, if you do not have
    #   networking)
    #   3 - Full multiuser mode
    #   4 - unused
    #   5 - X11
    #   6 - reboot (Do NOT set initdefault to this)
    #
    # id:5:initdefault:
    id:3:initdefault:

Ora se riavvii il sistema ti troverai davanti all'interfaccia di
login testuale. Questa ti permette di accedere al sistema con una
sessione testuale all'interno di un terminale virtuale. Una tale
sessione di solito non consuma molte risorse, come farebbe invece
un'interfaccia grafica. Poi puoi compilare Sage da codice sorgente in
tale sessione testuale. Dovresti assicurarti di essere in grado di
riattivare successivamente l'interfaccia grafica, prima di tentare di
accedere tramite un'interfaccia testuale.


Quando lancio Sage, SELinux segnala che "/path/to/libpari-gmp.so.2" richiede "text-relocation" (riallocazione del testo). Come posso risolvere?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Il problema può essere risolto eseguendo il seguente comando::

    chcon -t textrel_shlib_t /path/to/libpari-gmp.so.2


L'aggiornamento di Sage è andato bene, ma adesso l'indicatore continua a mostrare la versione precedente. Come posso risolvere?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

L'indicatore (banner in inglese) è memorizzato e non ricalcolato ad
ogni esecuzione di Sage. Il fatto che non sia aggiornato non dovrebbe
impedire a Sage di funzionare regolarmente. Digita ``banner()`` in una
sessione di Sage per verificare la versione reale. Se vuoi l'indicatore
corretto allora devi ricompilare Sage digitando ``make build`` in un
terminale.


Come posso eseguire Sage come demone/servizio?
""""""""""""""""""""""""""""""""""""""""""""""

Ci sono parecchie possibilità. Puoi usare i programmi a riga di comando
``screen``, ``nohup`` o ``disown``.


Il comando show (mostra) per la visualizzazione di oggetti 3D non funziona.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La visualizzazione 3D in tempo reale per Sage dalla versione 6.4 in
avanti usa il pacchetto `Jmol/JSmol <http://jmol.sourceforge.net>`_.
Dalla linea di comando viene utilizzata l'applicazione Java Jmol,
mentre per la visualizzazione dal browser viene usato puro javascript
oppure una Java applet. In genere nei browser è usato javascript puro
per evitare problemi con quei browser che non supportano i plugin per
le applet Java (ad esempio Chrome). In ogni worksheet su browser c'è
una casella da spuntare prima di generare una vista tridimensionale
qualora l'utente voglia usare l'applet Java (essa è un po' più veloce
con viste complicate).

La ragione più probabile di un malfunzionamento è che non hai
installato l'ambiente runtime di Java (JRE) o che è più vecchio della
versione 1.7. Se le cose funzionano dalla riga di comando,
un'altra possibilità è che il tuo browser non abbia il plugin giusto
per supportare le Java applet (al momento, nel 2014, tali plugin non
lavorano con la maggior parte delle versioni di Chrome). Assicurati di
aver installato il plugin IcedTea (su Linux vedi il tuo gestore dei
pacchetti) o il plugin di Oracle Java
(vedi: `IcedTea <http://icedtea.classpath.org/wiki/IcedTea-Web>`_
e `Java <https://java.com/en/download/help/index_installing.xml>`_).

Se stai usando un server Sage sul web e anche la visualizzazione
tramite javascript non funziona, potresti avere un problema con la
funzionalità javascript del tuo browser, o potresti aver disabilitato
javascript.


Posso usare gli strumenti di Sage in un ambiente commerciale?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sì! Assolutamente! Fondamentalmente l'unico *limite* che hai è che se
fai dei delle modifiche a Sage stesso e redistribuisci pubblicamente
tale versione modificata, allora devi renderci disponibili tali
modifiche cos' che le possiamo includere nella versione standard di
Sage (se vogliamo). Altrimenti sei libero di usare quante copie di
Sage vuoi per fare soldi, ecc. senza pagare alcuna licenza.


Voglio scrivere del codice Cython che usa l'aritmetica dei campi finiti, ma l'istruzione "cimport sage.rings.finite_field_givaro" non funziona. Cosa posso fare?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Devi segnalare a Sage di usare C++ (sia Givaro che NTL sono librerie
C++) e che hai bisogno anche delle librerie GMP e STDC C++.
Ecco un piccolo esempio::

    # These comments are hints to Cython about the compiler and
    # libraries needed for the Givaro library:
    #
    # distutils: language = c++
    # distutils: libraries = givaro gmpxx gmp m
    cimport sage.rings.finite_field_givaro
    # Construct a finite field of order 11.
    cdef sage.rings.finite_field_givaro.FiniteField_givaro K
    K = sage.rings.finite_field_givaro.FiniteField_givaro(11)
    print("K is a {}".format(type(K)))
    print("K cardinality = {}".format(K.cardinality()))
    # Construct two values in the field:
    cdef sage.rings.finite_field_givaro.FiniteField_givaroElement x
    cdef sage.rings.finite_field_givaro.FiniteField_givaroElement y
    x = K(3)
    y = K(6)
    print("x is a {}".format(type(x)))
    print("x = {}".format(x))
    print("y = {}".format(y))
    print("x has multiplicative order = {}".format(x.multiplicative_order()))
    print("y has multiplicative order = {}".format(y.multiplicative_order()))
    print("x*y = {}".format(x * y))
    # Show that x behaves like a finite field element:
    for i in range(1, x.multiplicative_order() + 1):
        print("{} {}".format(i, x**i))
    assert x*(1/x) == K.one()

Per saperne di più digita quanto segue al prompt di Sage ::

    sage.rings.finite_field_givaro.FiniteField_givaro.

Poi premi TAB, ed usa ``??`` per avere più informazioni su ogni
funzione. Ad esempio::

    sage.rings.finite_field_givaro.FiniteField_givaro.one??

fornisce informazioni sull'unità moltiplicativa nel campo finito.


La compilazione su Mac OS X fallisce in maniera incomprensibile. Come posso risolvere?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Cerca il file di log della compilazione (install.log) e
controlla se c'è il seguente messaggio::

    fork: Resource temporarily unavailable.

Se è così, prova a fare questo: crea (o modifica se c'è già) il file
``/etc/launchd.conf`` includendovi quanto segue::

    limit maxproc 512 2048

Poi riavvia. Vedi
`il seguente link <http://www.macosxhints.com/article.php?story=20050709233920660>`_
per maggiori dettagli.

Come disegno la radice cubica (o altre radici dispari) di numeri negativi?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Questa è una delle domande chieste più frequentemente. Ci sono molti
metodi descritti nelle documentazione per la funzione plot, ma questo
è il più semplice::

    sage: plot(real_nth_root(x, 3), (x, -1, 1))
    Graphics object consisting of 1 graphics primitive

Tuttavia, nota che il metodo più diretto::

    sage: plot(x^(1/3), (x, -1, 1))  # not tested

produce il grafico corretto solo per valori di `x` positivi. La
*ragione* per cui ciò avviene è che Sage restituisce dei numeri
complessi per le radici dispari di numeri negativi, quando queste sono
approssimate, il che è una :wikipedia:`convenzione standard <Cube_root#Complex_numbers>`::

    sage: numerical_approx( (-1)^(1/3) )
    0.500000000000000 + 0.866025403784439*I

Come va utilizzato in Sage l'operatore XOR bitwise?
"""""""""""""""""""""""""""""""""""""""""""""""""""

L'OR esclusivo in Sage si fa con l'operatore ``^^``. C'è anche il
corrispondente "operatore inplace" ``^^=``. Ad esempio::

   sage: 3^^2
   1
   sage: a = 2
   sage: a ^^= 8
   sage: a
   10

Se definisci 2 variabili e poi confronti::

    sage: a = 5; b = 8
    sage: a.__xor__(b), 13
    (13, 13)

Puoi anche fare::

    sage: (5).__xor__(8)
    13

Le parentesi sono necessarie affinché Sage non supponga di avere a
che fare con un numero reale. Ci sono molti modi di
definire una funzione::

    sage: xor = lambda x, y: x.__xor__(y)
    sage: xor(3, 8)
    11

Un'altra possibilità, che aggira il preparser di Sage, è ::

    sage: def xor(a, b):
    ....:     return eval("%s^%s" % (a, b))
    sage: xor(3, 8)
    11

Puoi anche disattivare il preparser di Sage con il comando
``preparser(False)``, a quel punto l'operatore ``^`` funzionerà
esattamente come in Python. Puoi successivamente riattivare il
preparser con il comando ``preparser(True)``. Questo funziona solo
dalla riga di comando di Sage. Se sei in una sessione Notebook, passa
in "Python mode".


Quando provo ad usare LaTeX in una sessione Notebook, dice che non trova "fullpage.sty".
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La risposta più ampia, ma forse non la più utile, è che hai bisogno di
installare ``fullpage.sty`` in una cartella acceduta da TeX. Su Ubuntu
(e probabilmente anche su altre distribuzioni Linux) dovresti
installare il pacchetto ``texlive-latex-extra``. Se non è disponibile,
prova ad installare il pacchetto ``tetex-extra``. Se stai usando
Mac OS X dovrai usare una qualunque distribuzione TeX che hai già per
ottenere ``fullpage.sty``
(se usi MacTeX probabilmente ce l'hai già installato).
Se stai usando l'immagine VirtualBox in Windows dovrai fare login in
tale immagine ed da lì installare ``texlive-latex-extra``.


Con degli oggetti "a" e "b" ed una funzione "f" ho digitato accidentalmente "f(a)=b" anziche' "f(a)==b". Questo mi ha dato un errore "TypeError" (come mi aspettavo) ma ha anche cancellato l'oggetto "a". Perchè?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Questo è dovuto a come sono definite le funzioni in Sage con la
notazione ``f(x)=expr`` usando il preparser. Nota anche che se fai
quest'errore in un costrutto ``if``, avrai un errore :class:`SyntaxError`
prima di qualunque altro comportamento errato, quindi, in questo caso,
non hai il problema.

Come posso usare un browser internet diverso con il notebook di Sage?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Dovrai farlo dalla linea di comando. Semplicemente lancia un comando
come questo.

* Linux (assumendo che hai Sage nella cartella ``/usr/bin``):

  .. CODE-BLOCK:: shell-session

    $ env BROWSER=opera /usr/bin/sage --notebook

* Mac (assumendo che tu sia nella cartella dove hai scaricato Sage).
  Con il notebook Jupyter:

  .. CODE-BLOCK:: shell-session

    $ BROWSER='open -a Firefox %s' ./sage --notebook jupyter
    $ BROWSER='open -a Google\ Chrome %s' ./sage --notebook jupyter


Dov'è il codice sorgente di ``<function>``?
"""""""""""""""""""""""""""""""""""""""""""

Le funzioni e classi scritte in Python o Cython sono di solito
accessibili tramite la linea di comando IPython con il comando ``??``::

    sage: plot??                            # not tested
    Signature: plot(*args, **kwds)
    Source:
    ...

Tuttavia gli oggetti che sono costruiti in Python o IPython sono
compilati e non verranno visualizzati. Ci sono molte funzioni in Sage
costruite come funzioni simboliche, i.e. possono essere usate come
parte di espressioni simboliche senza dover essere calcolate.
Il loro codice sorgente potrebbe non essere accessibile dalla linea di
comando, sopratutto per le funzioni elementari, poiché sono scritte
in C++ (per ragioni di efficienza).
