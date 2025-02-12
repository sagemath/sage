.. _chapter-faq-contribute:

========================
FAQ: Contribuire a Sage.
========================


Come posso iniziare a contribuire a Sage?
"""""""""""""""""""""""""""""""""""""""""

Il primo passo è usare Sage e incoraggiare i tuoi amici a usare Sage.
Se trovi bug o della documentazione poco chiara, per favore riporta i problemi!

Due modi comuni per contribuire a Sage sono scrivere codice sorgente e creare
documentazione o tutorial. Alcuni passi in entrambe le direzioni sono descritti
di seguito.

Voglio contribuire a Sage. Da dove inizio?
""""""""""""""""""""""""""""""""""""""""""

Dai un'occhiata alla
`guida ufficiale per lo sviluppo <https://doc.sagemath.org/html/en/developer>`_
di Sage. Come minimo, la lettura del primo capitolo è richiesta per ogni
sviluppatore di Sage. Inoltre sii molto attento alle
`linee guida per trac <https://doc.sagemath.org/html/en/developer/trac.html>`_.
Puoi entrare nella lista email
`sage-devel <https://groups.google.com/group/sage-devel>`_.
Mentre inizi a conoscere la comunità prendi una copia del codice sorgente di Sage
e familiarizza con `git <https://git-scm.com>`_, un software per il controllo
versione.

Il miglior modo per diventare familiare con il processo di sviluppo di Sage
è quello di scegliere un ticket dal
`server trac <https://trac.sagemath.org>`_
ed esaminare i cambiamenti proposti in quel ticket.
Se vuoi costruire qualcosa, è buona pratica discutere le tue idee sulla
lista email ``sage-devel``, in modo che altri sviluppatori abbiano l'opportunità
di esprimersi sulle tue idee/proposte. Sono anche aperti a nuove idee, come
tutti i matematici dovrebbero essere.

La principale lingua di programmazione per Sage è
`Python <https://www.python.org>`_.
Alcune parte di Sage potrebbero essere scritte in altre lingue,
specialmente le componenti che fanno l'elaborazione numerica più impegnativa,
ma la maggior parte delle funzionalità sono realizzate in Python,
incluso il codice di collegamento. Un aspetto valido di Python, che Sage eredita,
è che è più importante scrivere codice funzionante piuttosto che codice veloce.
Del codice veloce è prezioso, ma del codice chiaro e leggibile è importante.
Nella comunità matematica i risultati inaccurati sono inaccettabili.
La correttezza viene prima dell'ottimizzazione. Nel seguente articolo

* D. Knuth. Structured Programming with go to Statements.
  *ACM Journal Computing Surveys*, 6(4), 1974.

Don Knuth osserva che "dovremmo dimenticarci dell'efficienza di poco conto,
diciamo per il 97% del tempo: l'ottimizzazione prematura sta alla
radice di ogni male".

Se non conosci Python dovresti iniziare ad imparare il linguaggio.
Un buon posto dove iniziare è il
`Tutorial Ufficiale per Python <https://docs.python.org/3/tutorial>`_
e altra documentazione si trova in
`Documentazione standard di Python <https://docs.python.org>`_.
Un altro posto da guardare è al link
`Dive Into Python <https://diveintopython3.net>`_ di Marc Pilgrim,
che può essere veramente d'aiuto su temi specifici come
lo sviluppo guidato dai test. Il libro
`Building Skills in Python <http://itmaybeahack.com/homepage/books/python.html>`_
di Steven F. Lott è adatto a chiunque sia già a suo agio nel programmare.

Se desideri, puoi iniziare a imparare Python usando Sage.
Tuttavia, è utile sapere cosa è semplice Python e quando Sage sta usando la
sua "magia". Ci sono molte cose che funzionano in Python, ma non in Sage e
vice versa. Per di più anche quando la sintassi è identica, molti concetti
di programmazione sono spiegati più dettagliatamente in risorse focalizzate
su Python piuttosto che in risorse centrate su Sage; in quest'ultime,
la priorità viene data alla matematica.

Non sono un programmatore. C'è qualche altro modo in cui posso aiutare?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Certo. Come in ogni progetto FOSS ci sono molti modi in cui puoi dare il tuo
aiuto nella comunità di Sage e programmare è solo uno di questi.

Molte persone amano scrivere tutorial tecnici. Una delle gioie di farlo è che
impari qualcosa di nuovo nel farlo. Allo stesso tempo trasmetti delle conoscenze
ai principianti, un'abilità che è utile anche in campi estranei alla stesura
di documentazione tecnica. Un aspetto fondamentale della documentazione tecnica
è che espone un dato tecnico a dei principianti, pertanto il gergo
tecnico dev'essere ridotto al minimo. Darrell Anderson ha scritto
`alcuni suggerimenti sulla scrittura di documentazione tecnica
<http://web.archive.org/web/20130128102724/http://humanreadable.nfshost.com:80/howtos/technical_writing_tips.htm>`_,
il quale consigliamo vivamente.

Per i designer grafici o gli artisti c'è la possibilità di aiutare migliorando
l'aspetto del sito di Sage.

Se sai parlare, leggere e scrivere in un'altra lingua, ci sono molti modi in cui
il tuo contributo può essere molto prezioso all'intera comunità di Sage.
Diciamo che conosci l'italiano. Allora puoi scrivere un tutorial per Sage in
italiano, o aiutare a tradurre i tutorial ufficiali di Sage in italiano.

La suddetta è una lista molto breve.
Ci sono molti, molti più modi in cui puoi dare il tuo aiuto. Sentiti libero di
inviare una email alla mailing list
`sage-devel <https://groups.google.com/group/sage-devel>`_ per chiedere in quali
modi potresti essere d'aiuto, o per suggerire un'idea sul progetto.


Dove posso trovare risorse su Python e Cython?
""""""""""""""""""""""""""""""""""""""""""""""

Ecco una lista incompleta di risorse su Python e Cython.
Ulteriori risorse possono essere trovate cercando sul web.

**Risorse generali**

* `Cython <https://cython.org>`_
* `pep8 <https://pypi.org/project/pep8>`_
* `pydeps <https://pypi.org/project/pydeps>`_
* `pycallgraph <https://pycallgraph.readthedocs.io>`_
* `PyChecker <http://pychecker.sourceforge.net>`_
* `PyFlakes <https://pypi.org/project/pyflakes>`_
* `Pylint <https://www.logilab.org/project/pylint>`_
* `Python <https://www.python.org>`_ home page e la
  `Documentazione standard su Python <https://docs.python.org>`_
* `Snakefood <http://furius.ca/snakefood>`_
* `Sphinx <https://www.sphinx-doc.org>`_
* `XDot <https://github.com/jrfonseca/xdot.py>`_

**Tutorial e libri**

* `Cython Tutorial <http://conference.scipy.org/proceedings/SciPy2009/paper_1/>`_
  di Stefan Behnel, Robert W. Bradshaw, e Dag Sverre Seljebotn
* `Dive Into Python 3 <http://www.diveintopython3.net>`_ di Mark Pilgrim
* `Fast Numerical Computations with Cython <http://conference.scipy.org/proceedings/SciPy2009/paper_2/>`_
  di Dag Sverre Seljebotn
* `Tutorial ufficiale di Python <https://docs.python.org/3/tutorial/>`_

**Articoli e HOWTO**

* `decorator <https://pypi.org/project/decorator>`_
* `Functional Programming HOWTO <https://docs.python.org/3/howto/functional.html>`_
  di A. M. Kuchling
* `Python Functional Programming for Mathematicians <https://wiki.sagemath.org/devel/FunctionalProgramming>`_
  di Minh Van Nguyen
* `Regular Expression HOWTO <https://docs.python.org/3/howto/regex.html>`_
  di A. M. Kuchling
* `reStructuredText <https://docutils.sourceforge.io/rst.html>`_


Ci sono delle convenzioni di scrittura del codice sorgente che devo seguire?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Dovresti seguire le convenzioni standard di Python come documentato in
:pep:`8` e :pep:`257`.
Consulta anche la Guida dello Sviluppo Sage, specialmente il capitolo
`Convenzioni per scrivere codice sorgente in Sage <https://doc.sagemath.org/html/en/developer/#sage-coding-details>`_.


Ho inviato al GitHub Sage repo una correzione molte settimane fa. Perché la state ignorando?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Non stiamo cercando di ignorare la tua correzione.
Le persone che lavorano su Sage lo fanno nel loro tempo libero.
Con centinaia di ticket aperti aventi un impatto molto variabile sulla comunità
di Sage, le persone che ci lavorano devono dedicare il loro tempo e sforzo
principalmente a quei ticket che li interessano.
A volte potresti essere la sola persona che comprende la tua correzione.
In tal caso, ti invitiamo a fare uno sforzo supplementare per rendere
l'esaminazione della tua patch il più semplice possibile.
Ecco alcuni suggerimenti su come rendere la tua correzione facile da esaminare

* Hai descritto in modo chiaro il problema che la tua correzione vuole risolvere?
* Hai fornito ogni informazione di base rilevante al problema che la tua
  correzione vuole risolvere? Tali informazioni includono link a risorse online
  e ad articoli, libri, o altro materiale di riferimento.
* Hai descritto in modo chiaro come la tua correzione risolve il
  problema in oggetto?
* Hai descritto chiaramente nella tua correzione come effettuare i test
  dei cambiamenti?
* Hai elencato eventuali tickets da cui dipende la tua correzione?
* Se vi sono più correzioni, hai indicato chiaramente l'ordine in cui devono
  essere applicate ?
* La tua correzione segue le
  `convenzioni importanti <https://doc.sagemath.org/html/en/developer/#writing-code-for-sage>`_
  indicate nella "Guida dello sviluppatore"?

Se la tua correzione non ha la possibilità di essere aggiunta nell'albero dei
sorgenti di Sage, non la ignoreremo ma semplicemente chiuderemo il ticket
relativo con una spiegazione sul perché non possiamo includerla.

Come e quando posso ricordare alla comunità di Sage una correzione a cui tengo?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Ti suggeriamo di fare uno sforzo ulteriore sul come ricordare alla comunità di
Sage una correzione che vuoi venga inclusa nell'albero dei sorgenti di Sage.
Potrebbe esserci un prossimo evento "bug squash sprint" o "Sage days" che è
in relazione alla tua correzione. Tieni d'occhio le mailing list relative e
rispondi educatamente ad ogni scambio di email relativo,
spiegando chiaramente perché la tua correzione ha importanza.


Ho scritto del codice sorgente e voglio venga incluso in Sage. Però dopo aver rinominato il mio file ``a.sage`` in ``a.py`` ho degli errori di sintassi. Devo riscrivere tutto il mio codice in Python anziché in Sage?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

La risposta sostanzialmente è sì, ma riscrivere è una parola grossa per ciò che
bisogna realmente fare. C'è ben poco da fare dal momento che Sage per lo più
segue la sintassi di Python. Le 2 maggiori differenze sono la gestione degli
interi (vedi anche il link `afterword`_ per maggiori informazioni sul
preparser di Sage) e la necessità di importare quello che ti serve.

- **Gestione degli interi:** dei fare i seguenti cambiamenti:

  - Notazione per l'elevamento a potenza: In Python ``**`` significa elevamento
    a potenza e ``^`` significa “xor”.
  - Se devi restituire un intero all'utente, scrivi ``return Integer(1)``
    invece di ``return 1``. In Python, 1 è un intero Python (``int``), e
    ``Integer(1)`` è un intero Sage/Gmp. Inoltre gli ``Integer`` sono molto più
    potenti degli ``int``; ad esempio hanno collegata ad essi l'informazione di
    primalità e la fattorizzazione.
  - Dovresti anche notare che ``2 / 3`` non significa più
    ``Integer(2) / Integer(3)`` che restituisce ``2/3``, ma invece
    ``int(2) / int(3)``, e pertanto restituisce ``0`` poiché la divisione è
    intera e trascura il resto. Se stai lavorando con i tipi ``Integer``
    ma in realtà hai bisogno di eseguire una divisione intera puoi usare
    ``Integer(2) // Integer(3)``.

- **Note sull'importazione:** la seconda cosa importante da tenere presente è
  la necessità di importare tutto ciò di cui hai bisogno. Nel dettaglio, ogni
  volta che usi una funzione Sage la devi prima importare all'inizio del file.
  Ad esempio, se hai bisogno di ``PolynomialRing``, dovrai scrivere::

      from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

  Puoi chiedere a Sage dove si trove ``PolynomialRing`` usando::

      sage: PolynomialRing.__module__
      'sage.rings.polynomial.polynomial_ring_constructor'

  Questo corrisponde anche al percorso, che inizia dopo ``site-packages``,
  restituito da Sage quando richiami l'help su ``PolynomialRing``. A
  d esempio se scrivi ``PolynomialRing?`` otterrai::

      Type:    function
      [...]
      File:    /path_to_sage_root/sage/local/lib/python3.7/site-packages/sage/rings/polynomial/polynomial_ring_constructor.py
      [...]


.. _afterword: https://doc.sagemath.org/html/en/tutorial/afterword.html
