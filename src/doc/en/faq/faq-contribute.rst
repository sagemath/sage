.. _chapter-faq-contribute:

=========================
FAQ: Contributing to Sage
=========================


How can I start contributing to Sage?
"""""""""""""""""""""""""""""""""""""

The first step is to use Sage and encourage your friends to use
Sage. If you find bugs or confusing documentation along the way,
please report your problems!

Two popular ways to contribute to Sage are to write code and to
create documentation or tutorials. Some steps in each direction
are described below.


I want to contribute code to Sage. How do I get started?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Take a look at the
`official development guide <https://doc.sagemath.org/html/en/developer>`_
for Sage. At a minimum, the first chapter in that guide is required
reading for any Sage developer. Also pay special attention to the
`GitHub guidelines <https://doc.sagemath.org/html/en/developer/github.html>`_.
You can also join the
`sage-devel <https://groups.google.com/group/sage-devel>`_
mailing list or hang around on the
`SageMath Zulip chat channel <https://sagemath.zulipchat.com/>`_.
While you are getting to know the community, grab a copy of the Sage
source and familiarize yourself with the
`git <https://git-scm.com>`_ version control system.

The best way to become familiar with the Sage development process is
to choose an issue from
`the Sage repository on GitHub <https://github.com/sagemath/sage/issues>`_
and review the proposed changes contained in that issue. If you want
to implement something, it is a good practice to discuss your ideas on
the ``sage-devel`` mailing list first, so that other developers have a
chance to comment on your ideas/proposals. They are pretty open to new
ideas, too, as all mathematicians should be.

Sage's main programming language is
`Python <https://www.python.org>`_.
Some parts of Sage may be written in other languages, especially the
components that do the heavy number crunching, but most native
functionality is done using Python, including "glue code". One of the
good aspects of Python that Sage inherits is that working code is
considered more valuable than merely fast code. Fast code is valuable,
but clean, readable code is important. In the mathematics community,
inaccurate results are unacceptable. Correctness comes before
optimization. In the following paper

* D. Knuth. Structured Programming with go to Statements.
  *ACM Journal Computing Surveys*, 6(4), 1974.

Don Knuth observes that: "We should forget about small efficiencies,
say about 97% of the time: premature optimization is the root of all
evil."

If you do not know Python, you should start learning that language. A
good place to start is the
`Python Official Tutorial <https://docs.python.org/3/tutorial>`_
and other documents in the
`Python standard documentation <https://docs.python.org>`_.
Another good place to take a look at is
`Dive Into Python <https://diveintopython3.net>`_
by Mark Pilgrim, which may be pretty helpful on some specific topics
such as test-driven development. The book
`Building Skills in Python <http://itmaybeahack.com/homepage/books/python.html>`_
by Steven F. Lott is suitable for anyone who is already comfortable
with programming.

If you want, you can
try to learn Python by using Sage. However,
it is helpful to know what is pure Python and when Sage is doing its
"magic". There are many things that work in Python but not in Sage,
and vice versa. Furthermore, even when the syntax is identical, many
programming concepts are explained more thoroughly in Python-centered
resources than in Sage-centered resources; in the latter,
mathematics is usually the priority.


I am not a programmer. Is there another way I can help out?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Yes. As with any free open source software project, there are numerous
ways in which you could help out within the Sage community, and
programming is only one of many ways to contribute.

Many people like writing technical tutorials. One of the joys of doing
so is that you also learn something new in the process. At the same
time, you communicate your knowledge to beginners, a skill which is
useful in fields other than technical writing itself. A main point
about technical writing is that you communicate a technical subject to
beginners, so keep technical jargon to a minimum. Darrell Anderson
has written some
`tips on technical writing <http://web.archive.org/web/20130128102724/http://humanreadable.nfshost.com:80/howtos/technical_writing_tips.htm>`_,
which we highly recommend.

For the graphic designers or the artistically creative, you can
help out with improving the design of the Sage website.

If you can speak, read,
and write in another (natural) language, there are many ways in which
your contribution would be very valuable to the whole Sage
community. Say you know Italian. Then you can write a Sage tutorial in
Italian, or help out with translating the official Sage tutorial to
Italian.

The above is a very short list. There are many, many more ways in
which you can help out. Feel free to send an email to the
`sage-devel <https://groups.google.com/group/sage-devel>`_ mailing list
to ask about possible ways in which you could help out, or to suggest a
project idea.


Where can I find resources on Python or Cython?
"""""""""""""""""""""""""""""""""""""""""""""""

Here is an incomplete list of resources on Python and Cython. Further
resources can be found by a web search.

**General resources**

* `Cython <https://cython.org>`_
* `pep8 <https://pypi.org/project/pep8>`_
* `pydeps <https://pypi.org/project/pydeps>`_
* `pycallgraph <https://pycallgraph.readthedocs.io>`_
* `PyChecker <http://pychecker.sourceforge.net>`_
* `PyFlakes <https://pypi.org/project/pyflakes>`_
* `Pylint <https://www.logilab.org/project/pylint>`_
* `Python <https://www.python.org>`_ home page and the
  `Python standard documentation <https://docs.python.org>`_
* `Snakefood <http://furius.ca/snakefood>`_
* `Sphinx <https://www.sphinx-doc.org>`_
* `XDot <https://github.com/jrfonseca/xdot.py>`_

**Tutorials and books**

* `Cython Tutorial <http://conference.scipy.org/proceedings/SciPy2009/paper_1/>`_
  by Stefan Behnel, Robert W. Bradshaw, and Dag Sverre Seljebotn
* `Dive Into Python 3 <http://www.diveintopython3.net>`_ by Mark Pilgrim
* `Fast Numerical Computations with Cython <http://conference.scipy.org/proceedings/SciPy2009/paper_2/>`_
  by Dag Sverre Seljebotn
* `Official Python Tutorial <https://docs.python.org/3/tutorial/>`_

**Articles and HOWTOs**

* `decorator <https://pypi.org/project/decorator>`_
* `Functional Programming HOWTO <https://docs.python.org/3/howto/functional.html>`_
  by A. M. Kuchling
* `Python Functional Programming for Mathematicians <https://wiki.sagemath.org/devel/FunctionalProgramming>`_
  by Minh Van Nguyen
* `Regular Expression HOWTO <https://docs.python.org/3/howto/regex.html>`_
  by A. M. Kuchling
* `reStructuredText <https://docutils.sourceforge.io/rst.html>`_


Are there any coding conventions I need to follow?
""""""""""""""""""""""""""""""""""""""""""""""""""

You should follow the standard Python conventions as documented at
:pep:`8` and :pep:`257`.
Also consult the Sage Developer's Guide, especially the chapter
`Conventions for Coding in Sage <https://doc.sagemath.org/html/en/developer/#sage-coding-details>`_.


I submitted a bug fix to the GitHub Sage repo several weeks ago. Why are you ignoring my branch?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

We are not trying to ignore your branch. Most people who work on Sage do so
in their free time. With hundreds of open PRs of varying degrees of
impacts on the whole Sage community, people who work on PRs need
to prioritize their time and work on those PRs that interest
them. Sometimes you may be the only person who understands your
branch. In that case, you are encouraged to take extra care to make it
as easy as possible for anyone to review. Here are some
tips on making your branch easy to review:

* Have you clearly described the problem your branch is trying to
  solve?
* Have you provided any background information relevant to the problem
  your patch is trying to solve? Such information include links to
  online resources and any relevant papers, books and reference
  materials.
* Have you clearly described how your branch solves the problem under
  consideration?
* Have you clearly described how to test the changes in your branch?
* Have you listed any issues/PRs that your branch depends on?
* Is your branch based on a recent (preferably, the latest) Sage beta version?
* Does your branch
  `follow relevant conventions <https://doc.sagemath.org/html/en/developer/#writing-code-for-sage>`_
  as documented in the Developer's Guide?

If your branch stands no chance of being merged in the Sage source
tree, we will not ignore your branch but simply close the relevant
PR with an explanation why we cannot include your changes.


When and how might I remind the Sage community of a branch I care about?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You are encouraged to take extra care in how you remind the Sage
community of a branch/patch you want to get merged into the Sage source
tree. There might be an upcoming bug squash sprint or an upcoming Sage
Days workshop that relates to your patch. Monitor the relevant Sage
mailing lists and respond politely to any relevant email threads, with
clear explanation on why your patch is relevant.


I wrote some Sage code and I want it to be integrated into Sage. However, after renaming my file ``a.sage`` to ``a.py``, I got syntax errors. Do I have to rewrite all my code in Python instead of Sage?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The basic answer is yes, but rewriting is a big word for what is
really needed. There is little work to do since Sage mostly follows
Python syntax. The two main differences are handling of integer (see
also the `afterword`_ for more on the sage preparser), and the
necessity to import what you need.

- **Handling of integers:** You need to take care of the following
  changes:

  - Notation for exponentiation: In Python ``**`` means exponentiation
    and ``^`` means "xor".
  - If you need to return an integer to the user, write ``return
    Integer(1)`` instead of ``return 1``. In Python, 1 is a python
    ``int``, and ``Integer(1)`` is a Sage/Gmp integer. In addition,
    ``Integer`` are much more powerful than ``int``; for
    example, they know about being prime and factorization.
  - You should also notice that ``2 / 3`` no longer means
    ``Integer(2) / Integer(3)`` and returns ``2/3``, but rather
    ``int(2) / int(3)``, and therefore returns ``0`` due to integer
    division where it deregards any remainder. If you are dealing with
    ``Integer`` but you really need an integer division you can use
    ``Integer(2) // Integer(3)``.

- **Importing stuff:** The second big change is the necessity to
  import everything what you need. More precisely, each time you use
  some Sage function, you need to import it at the beginning of the
  file. For example, if you want ``PolynomialRing``, you need to
  write:

  .. CODE-BLOCK:: python

      from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

  You can ask Sage where to find ``PolynomialRing`` using::

      sage: PolynomialRing.__module__
      'sage.rings.polynomial.polynomial_ring_constructor'

  This also corresponds to the path starting after ``site-packages``
  given when you ask Sage for ``PolynomialRing`` help. For example,
  if you call ``PolynomialRing?``, you get:

  .. CODE-BLOCK:: text

      Type:    function
      [...]
      File:    /path_to_sage_root/sage/local/lib/python3.7/site-packages/sage/rings/polynomial/polynomial_ring_constructor.py
      [...]


.. _afterword: https://doc.sagemath.org/html/en/tutorial/afterword.html
