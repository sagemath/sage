.. -*- coding: utf-8 -*-
.. _chapter-faq-general:

===============
FAQ: Generalità
===============


Perchè esiste questo progetto?
""""""""""""""""""""""""""""""

La missione fissata per Sage è di essere un'alternativa open-source a Magma,
Maple, Mathematica e Matlab. I predecessori di Sage, noti come HECKE e Manin,
furono creati perché William Stein ebbe bisogno di scriverli come parte della
sua ricerca sulla Teoria dei Numeri. Iniziato da William nel 2005 quando era
all'università di Harvard, Sage combina alcuni fra i migliori software open-source
per la matematica in un'unica interfaccia comune.
Da allora Sage viene utilizzato non solo da ricercatori nel campo della Teoria
dei Numeri, ma da ricercatori in tutte le scienze matematiche.

Sage si avvale ed estende le funzionalità di molti dei pacchetti inglobati.
Anche dal principio, quando Sage veniva usato principalmente per la Teoria dei
Numeri, includeva:
`Givaro <https://casys.gricad-pages.univ-grenoble-alpes.fr/givaro>`_,
`GMP <https://gmplib.org>`_,
`NTL <https://www.shoup.net/ntl>`_,
`Pari/GP <https://pari.math.u-bordeaux.fr>`_,
e molti altri troppo numerosi per essere elencati qui. Studenti, insegnanti,
professori universitari, ricercatori di tutto il mondo usano Sage perché vogliono
un pacchetto open-source comprensivo per la matematica che offra calcolo sia
simbolico che numerico. Perlopiù le persone sono contente di quanto offre Sage.

Come é comune nell'ambito del software open-source (FOSS), spesso ci sono persone
che individuano casi in cui Sage non dispone delle funzionalità richiesta da loro.
Quindi si immergono nel codice sorgente di Sage per estenderlo per il loro scopo,
o ancora per esporre delle funzionalità dei pacchetti inglobati in Sage in modo
da poter usarne le funzioni che loro preferiscono dall'interfaccia di Sage.
La squadra `Sage-Combinat <http://combinat.sagemath.org>`_ e' costituita da
ricercatori in Algebra Combinatoria. La missione che si è data tale squadra è
quella di migliorare Sage come uno strumento estendibile per la sperimentazione
a computer nell'ambito dell'Algebra Combinatoria, e favorire lo scambio di codice
sorgente fra ricercatori di questa materia.
Per informazione dettagliate sul perché esiste Sage, vedere il seguente link:
biografia matematica personale di William,
`settore software <http://sagemath.blogspot.com/2009/12/mathematical-software-and-me-very.html>`_.


Cosa vuol dire Sage e come si deve pronunciarlo?
""""""""""""""""""""""""""""""""""""""""""""""""

Nei primi anni di esistenza di Sage, il progetto era chiamato "SAGE", acronimo
inglese di "software per la sperimentazione in Algebra e Geometria".
A cominciare dal 2007 ed inizio 2008, fu largamente adottato il nome "Sage".
Considera "Sage" come il nome di un progetto software FOSS per la matematica,
esattamente come "Python" è il nome di un linguaggio FOSS di programmazione di
uso generale. Ovunque possibile, per cortesia usa "Sage" non "SAGE", che invece
è un progetto per un computer
(`SAGE <http://history.sandiego.edu/GEN/20th/sage.html>`_),
così da evitare confusioni. "Sage" si pronuncia nello stesso modo della parola
inglese "sage" che significa uomo saggio, o anche indica la pianta della salvia.
Alcune persone lo pronunciano in modo simile a "sarge", un po' come si pronuncia
`Debian <http://www.debian.org>`_ Sarge.

Comunque lo pronunci, per cortesia non confonderlo con il software di contabilità
americano che ha lo stesso nome.


Chi c'è dietro al progetto?
"""""""""""""""""""""""""""

Sage è un progetto basato sull'opera di volontari. Il suo successo è dovuto
all'opera gratuita di una grande squadra internazionale di studenti, insegnanti,
professori universitari, ricercatori, ingegneri informatici e persone che
lavorano in vari ambiti della matematica, delle scienze, dell'ingegneria, dello
sviluppo software e a tutti i livelli della scuola. Lo sviluppo di Sage ha potuto
usufruire di fondi assegnati da numerose istituzioni e ha potuto includere sia
componenti preesistenti che in corso di sviluppo da parte di numerosi autori.

Una lista di coloro che hanno dato un contributo diretto è reperibile al link
"mappa di sviluppo di Sage"
(`Sage Development Map <http://www.sagemath.org/development-map.html>`_)
e la storia delle modifiche può essere reperita al link "changelog di
alto livello" (`changelog <http://www.sagemath.org/mirror/src/changelog.txt>`_).
Fai riferimento alla
`Pagina dei riconoscimenti <http://www.sagemath.org/development-ack.html>`_
del sito web di Sage per una lista aggiornata di coloro che ci sostengono
finanziariamente o a livello di infrastruttura, a livello di siti mirror ed altri
contributi indiretti.


Perché Sage è un software libero ed open-source?
""""""""""""""""""""""""""""""""""""""""""""""""

Una regola universale nella comunità matematica è che tutto dev'essere pubblico
e aperto ad analisi. Il progetto Sage ritiene che non seguire lo stesso principio
nel software per la matematica è quantomeno scortese e maleducato, o peggio una
violazione delle pratiche comuni nella scienza. Un principio filosofico
sottostante Sage è di applicare la regola di scambio libero e confronto fra pari,
che caratterizza la comunicazione scientifica, anche allo sviluppo di software
per la matematica. Né il progetto Sage né la sua squadra di sviluppo hanno la
pretesa di essere gli originali proponenti di questo principio.

Il modello di sviluppo di Sage è largamente ispirato del movimento del software
libero di cui + stata pioniere la
`Free Software Foundation <http://www.fsf.org>`_
ed il movimento open-source. Una fonte di ispirazione all'interno della comunità
matematica è Joachim Neubüser, come espresso nell'articolo

* J. Neubüser. An invitation to computational group theory. In
  C. M. Campbell, T. C. Hurley, E. F. Robertson, S. J. Tobin, and
  J. J. Ward, editors, *Groups '93 Galway/St. Andrews, Volume 2*,
  volume 212 of London Mathematical Society Lecture Note Series, pages
  457--475. Cambridge University Press, 1995.

ed in particolare nella seguente citazione dal suo articolo::

  Puoi leggere il teorema di Sylow e la sua dimostrazione nel libro di Huppert
  in biblioteca senza nemmeno comprare il libro e poi usare questo teorema per
  il resto della tua vita senza dover pagare una tariffa, invece...per molti
  sistemi di algebra computazionale devi pagare regolarmente delle licenze
  per tutto il tempo in cui li utilizzi. Per proteggere ciò per cui devi pagare,
  non ti viene dato il codice sorgente ma soltanto l'eseguibile del programma,
  cioè un oggetto sigillato sul quale premi bottoni ed ottieni risposte così come
  ottieni belle immagini dal tuo televisore: in entrambe le situazioni non puoi
  controllare il modo in cui sono create.

  In questa situazione sono violate due delle più basilari regole di condotta in
  matematica: in essa l'informazione è passata gratuitamente e tutto può essere
  consultato e verificato. Non applicare queste regole ai sistemi di algebra
  computazionale usati per la ricerca in matematica...significa andare in una
  direzione assolutamente non desiderabile. Ancora più importante:
  possiamo aspettarci che qualcuno creda al risultato di un programma il cui
  funzionamento non gli è permesso di vedere? Inoltre: vogliamo veramente far
  pagare a colleghi in Moldavia molti anni di stipendio per un
  sistema di algebra computazionale?

Simili idee sono state anche espresse da Andrei Okounkov, come si può leggere in

* V. Muñoz and U. Persson. Interviews with three Fields
  medalists. *Notices of the American Mathematical Society*,
  54(3):405--410, 2007.

ed in particolare nella seguente citazione::

  I computer non sono una minaccia per i matematici più di quanto i robot da
  cucina lo siano per i cuochi. Poiché la matematica diviene sempre più complessa
  mentre il ritmo delle nostre vite accelera, dobbiamo delegare il più possibile
  alle macchine. Ed intendo sia il lavoro in campo numerico che in quello
  simbolico. Alcune persone possono andare avanti senza lavastoviglie, ma penso
  che le dimostrazioni vengano fuori molto più pulite quando il lavoro di
  routine è automatizzato.

  Questo porta con sè parecchie questioni. Non sono un esperto ma penso che
  abbiamo bisogno di uno standard a livello di calcolo simbolico per rendere le
  manipolazioni al computer più facili da documentare e verificare. Con tutto il
  rispetto per il libero mercato, forse in questo non dobbiamo essere dipendenti
  da un software commerciale. Un progetto open-source potrebbe, forse, trovare
  risposte migliori a problemi ovvi come la disponibilità, i bug, la
  compatibilità all'indietro, l'indipendenza dalla piattaforma, le librerie
  standard, ecc. Si può imparare dal successo di TeX e da software più
  specializzati come Macaulay2. Spero veramente che le agenzie per finanziamenti
  governativi stiano considerando questo.


Perché avete scritto Sage da zero, invece di usare software e librerie preesistenti?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sage non è stato scritto da zero. La maggior parte delle sue funzionalità sono
realizzate attraverso progetti FOSS come

* `BLAS <http://www.netlib.org/blas>`_ --- sottoprogrammi per Algebra
  Lineare di base.

* `FLINT <http://www.flintlib.org>`_ --- libreria C per Teoria dei Numeri.

* `GAP <https://www.gap-system.org>`_ --- sistema di calcolo per algebra discreta,
  con particolare enfasi sulla teoria dei gruppi.

* `Maxima <http://maxima.sourceforge.net>`_ --- sistema di calcolo
  simbolico e numerico.

* `mpmath <http://code.google.com/p/mpmath>`_ --- libreria in puro Python per
  aritmetica floating-point di precisione.

* `NumPy <http://numpy.scipy.org>`_ --- algebra lineare numerica ed altre
  funzioni di calcolo numerico per Python.

* `Pari/GP <https://pari.math.u-bordeaux.fr>`_ --- software matematico per
  calcolo veloce in Teoria dei Numeri.

* `Pynac <http://pynac.sagemath.org>`_ --- versione modificata di GiNaC che
  rimpiazza la dipendenza da CLN con Python.

* `R <http://www.r-project.org>`_ --- linguaggio ed ambiente operativo per
  calcolo statistico e grafici relativi.

* E molti altri troppo numerosi per essere elencati qui.

Una lista aggiornata può essere reperita alla seguente link:
`repository dei pacchetti standard <http://www.sagemath.org/packages/standard>`_.

I principali linguaggi di programmazione di Sage sono
`Python <http://www.python.org>`_ e `Cython <http://www.cython.org>`_.
Python è il principale linguaggio di programmazione e di interfacciamento,
mentre Cython è il principale linguaggio per ottimizzare funzionalità critiche e
per interfacciarsi con le librerie C e le estensioni C per Python. Sage integra
oltre 90 pacchetti FOSS in un'interfaccia comune. Sopra questi pacchetti sta la
libreria Sage, che consiste in oltre 700.000 righe di codice Python e Cython
scritto ex-novo. Vedi `openhub.net <https://www.openhub.net/p/sage>`_ per
l'analisi del codice sorgente dell'ultima release stabile di Sage.


Come posso ricevere aiuto?
""""""""""""""""""""""""""

Per ricevere aiuto sull'uso di Sage, ci sono due opzioni

* Il sito web di domande e risposte
  ``ask.sagemath.org`` : https://ask.sagemath.org/questions/
* La lista email ``sage-support``: http://groups.google.com/group/sage-support

Per aiuto sullo sviluppo di Sage, c'è una list email
``sage-devel`` : https://groups.google.com/forum/#!forum/sage-devel

Consulta http://www.sagemath.org/help.html per una lista con altre risorse.


Non sarebbe meglio se Sage non fosse distribuito come un gigantesco aggregato di pacchetti?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La distribuzione SageMath continua a construire le proprie versione dei sofware
richiesti ("SPKGs") che funzionano bene assieme.

Tuttavia, per ridurre i tempi di compilazione e le dimensioni dell'installazione
di Sage, fin dalle versioni 8.x c'è stato uno sforzo a livello di sviluppo che ha
reso possibile utilizzare molti pacchetti software già presenti nel sistema
operativo (or da distribuzioni di Homebrew o conda-forge) invece di costruire
delle copie solo per SageMath.

Il meccanismo nominato "spkg-configure" è utilizzato all'inizio del processo di
installazione dal codice sorgente durante la fase ``./configure``.

Per assicurasi chec SageMath si installi e funzioni correttamente su un'ampia
gamma di sistemi, noi utilizziamo dei test automatici.
Vai al capitolo
`Test di portabilità <https://doc.sagemath.org/html/en/developer/portability_testing.html>`_
nella Guida Per Sviluppatori per maggiori informazioni.

Perché ci sono così tanti bug in Sage, con centinaia di modifiche in corso, perché non producete una versione stabilizzata?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Ogni software contiene bug. In qualcosa di così complesso come Sage nessuno, né
la squadra di sviluppo di Sage né la sua comunità, ha alcuna pretesa che esso sia
libero da bug. Farlo sarebbe un atto di disonestà.

Un ciclo di rilascio di Sage di solito dura alcuni mesi, con molte versioni beta
che si susseguono a intervalli di 2-3 settimane. Ogni ciclo di rilascio è
presieduto da un singolo gestore che si occupa dell'albero di integrazione
pacchetti per tutta la durata del ciclo. In questa fase tale gestore deve spesso
dedicare del tempo equivalente ad un lavoro a tempo pieno alla gestione della
qualità e deve interagire attivamente con la comunità internazionale degli
utenti, degli sviluppatore e dei potenziali contributori a Sage.

Ci sono stati alcuni casi in cui due contributori a Sage sono stati affiancati
come gestori di rilascio per un ciclo di rilascio di Sage. Comunque in genere
poche persone hanno l'equivalente di 3 settimane di tempo libero per dedicarsi
alla gestione del rilascio. Se vuoi aiutare nella gestione del rilascio iscriviti
alla mailing list `sage-release <http://groups.google.com/group/sage-release>`_.

Fin dall'inizio del progetto Sage, i contributori hanno cercato di ascoltare e
di riflettere su cosa potesse aumentare la possibilità che altri potenziali
validi contributori dessero effettivamente un aiuto. Cosa incoraggia un
contributore può scoraggiare un altro, quindi bisogna trovare degli equilibri.
Decidere che un rilascio stabilizzato dovrebbe includere le patch di correzione
dei bug, e solo quelle, probabilmente scoraggerebbe qualcuno dal contribuire,
nel momento in cui gli fosse detto in anticipo che la sua aggiunta, anche se
giudicata positivamente, non verrebbe integrata nel rilascio.

La comunità Sage crede nel principio "pubblica subito, pubblica spesso". Il modo
in cui il progetto Sage è organizzato e gestito differisce parecchio da quello
di una azienda di software commerciale. I contributori sono tutti volontari e
questo cambia totalmente la dinamica del progetto da quella che sarebbe se Sage
fosse un'iniziativa software commerciale con sviluppatori stipendiati a
tempo pieno.


Come posso scaricare la documentazione di Sage così da poterla leggere offline?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Per scaricare la documentazione standard di Sage in formato HTML o PDF, visita
`Help and Support <http://www.sagemath.org/help.html>`_ sul sito web di Sage.
Ogni release di Sage dispone della documentazione completa che costituisce la
documentazione standard di Sage. Se hai scaricato un rilascio di Sage in formato
binario, la versione HTML della sua documentazione si trova già disponibile nella
cartella :sage_root:`src/doc/output/html/`. Nel corso della compilazione da
sorgente viene preparata anche la documentazione HTML.
Per construire la versione HTML della documentazione, lancia il seguente comando
dopo essersi posizionati in ``SAGE_ROOT``::

    $ ./sage -docbuild --no-pdf-links all html

La preparazione della documentazione in formato PDF richiede che sul tuo sistema
sia installata una versione funzionante di LaTeX. Per preparare la documentazione
in formato PDF puoi lanciare il comando,
dopo esserti posizionato in ``SAGE_ROOT``::

    $ ./sage -docbuild all pdf

Per altre maggiori opzioni disponibili da riga di comando fai riferimento alle
istruzioni stampate dai seguenti comandi::

    $ ./sage -help
    $ ./sage -advanced
