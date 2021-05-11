.. -*- coding: utf-8 -*-
.. _chapter-faq-general:

============
FAQ: General
============


Why does this project exist?
""""""""""""""""""""""""""""

The stated mission of Sage is to be viable free open source
alternative to Magma, Maple, Mathematica, and Matlab. Sage's
predecessors, known as HECKE and Manin, came about because William
Stein needed to write them as part of his research in number
theory. Started by William in 2005 during his time at Harvard
University, Sage combines best-of-breed free open source mathematics
software, packaging and unifying them through a common interface. Since
then Sage has become something used not just by researchers in
number theory, but throughout the mathematical sciences.

Sage builds upon and extends functionalities of many underlying
packages.  Even from early on, when Sage was primarily used for
number theory, this included
`Givaro <http://ljk.imag.fr/CASYS/LOGICIELS/givaro>`_,
`MPIR <http://www.mpir.org>`_,
`NTL <http://www.shoup.net/ntl>`_,
`Pari/GP <http://pari.math.u-bordeaux.fr>`_,
and many others too numerous to list here. Students, teachers,
professors, researchers throughout the world use Sage because they
require a comprehensive free open source mathematics package that
offers symbolic and numerical computation. Most of the time, people
are happy with what Sage has to offer.

As is common throughout the
free open source software (FOSS) world, many people often identify
cases where Sage lacks certain mathematics functionalities that they
require. And so they delve into the underlying source code that
comprises Sage in order to extend it for their purposes, or expose
functionalities of underlying packages shipped with Sage in order to
use their favourite mathematics software packages from within Sage. The
`Sage-Combinat <http://combinat.sagemath.org>`_
team is comprised of researchers in algebraic combinatorics. The
team's stated mission is to improve Sage as an extensible toolbox for
computer exploration in algebraic combinatorics, and foster code
sharing between researchers in this area.

For detailed information
about why Sage exists, see William's personal
`mathematics software biography <http://sagemath.blogspot.com/2009/12/mathematical-software-and-me-very.html>`_.


What does "Sage" mean and how do you pronounce it?
""""""""""""""""""""""""""""""""""""""""""""""""""

In the first few years of Sage's existence, the project was called
"SAGE". This acronym stood for "Software for Algebra and Geometry
Experimentation". Starting around 2007 and early 2008, the name "Sage"
was widely adopted. Think of "Sage" as a name for a free open source
mathematics software project, just as "Python" is a name for a free
open source general purpose programming language. Whenever possible,
please use the name "Sage" instead of "SAGE" to avoid confusing the
Sage project with a computer project called
`SAGE <http://history.sandiego.edu/GEN/20th/sage.html>`_.
You pronounce "Sage" similar to how you would pronounce "sage" which
refers to a wise person, or "sage" which refers to a plant. Some
people pronounce "Sage" as "sarge", similar to how you would pronounce
`Debian <http://www.debian.org>`_
Sarge.

However you pronounce "Sage", please do not confuse the Sage project
with an accounting software by the same name.


Who is behind this project?
"""""""""""""""""""""""""""

Sage is a volunteer based project. Its success is due to the voluntary
effort of a large international team of students, teachers,
professors, researchers, software engineers, and people working in
diverse areas of mathematics, science, engineering, software
development, and all levels of education. The development of Sage has
benefited from the financial support of numerous institutions, and the
previous and ongoing work of many authors of included components.

A list of (some) direct contributors can be found on the
`Sage Development Map <http://www.sagemath.org/development-map.html>`_
and the history of changes can be found in the
`changelogs <http://www.sagemath.org/changelogs/>`_. Refer
to the
`acknowledgment page <http://www.sagemath.org/development-ack.html>`_
of the Sage website for an up-to-date list of financial and
infrastructure supporters, mirror network hosting providers, and
indirect contributors.


Why is Sage free/open source?
"""""""""""""""""""""""""""""

A standard rule in the mathematics community is that everything is
laid open for inspection. The Sage project believes that not doing the
same for mathematics software is at best a gesture of impoliteness
and rudeness, and at worst a violation against standard scientific
practices. An underlying philosophical principle of Sage is to apply
the system of open exchange and peer review that characterizes
scientific communication to the development of mathematics
software. Neither the Sage project nor the Sage Development Team make
any claims to being the original proponents of this principle.

The development model of Sage is largely inspired by the free software
movement as spearheaded by the
`Free Software Foundation <http://www.fsf.org>`_,
and by the open source movement. One source of inspiration from within
the mathematics community is Joachim Neubüser as expressed in the paper

* J. Neubüser. An invitation to computational group theory. In
  C. M. Campbell, T. C. Hurley, E. F. Robertson, S. J. Tobin, and
  J. J. Ward, editors, *Groups '93 Galway/St. Andrews, Volume 2*,
  volume 212 of London Mathematical Society Lecture Note Series, pages
  457--475. Cambridge University Press, 1995.

and in particular the following quotation from his paper:

.. CODE-BLOCK:: text

    You can read Sylow's Theorem and its proof in Huppert's book in
    the library without even buying the book and then you can use
    Sylow's Theorem for the rest of your life free of charge,
    but...for many computer algebra systems license fees have to be
    paid regularly for the total time of their use. In order to
    protect what you pay for, you do not get the source, but only an
    executable, i.e. a black box. You can press buttons and you get
    answers in the same way as you get the bright pictures from your
    television set but you cannot control how they were made in either
    case.

    With this situation two of the most basic rules of conduct in
    mathematics are violated: In mathematics information is passed on
    free of charge and everything is laid open for checking. Not
    applying these rules to computer algebra systems that are made for
    mathematical research...means moving in a most undesirable
    direction. Most important: Can we expect somebody to believe a
    result of a program that he is not allowed to see? Moreover: Do we
    really want to charge colleagues in Moldava several years of their
    salary for a computer algebra system?

Similar sentiments were also expressed by Andrei Okounkov as can be
found in

* V. Muñoz and U. Persson. Interviews with three Fields
  medalists. *Notices of the American Mathematical Society*,
  54(3):405--410, 2007.

in particular the following quotation:

.. CODE-BLOCK:: text

    Computers are no more a threat to mathematicians than food
    processors are a threat to cooks. As mathematics gets more and
    more complex while the pace of our lives accelerates, we must
    delegate as much as we can to machines. And I mean both numeric
    and symbolic work. Some people can manage without dishwashers, but
    I think proofs come out a lot cleaner when routine work is
    automated.

    This brings up many issues. I am not an expert, but I think we
    need a symbolic standard to make computer manipulations easier to
    document and verify. And with all due respect to the free market,
    perhaps we should not be dependent on commercial software here. An
    open-source project could, perhaps, find better answers to the
    obvious problems such as availability, bugs, backward
    compatibility, platform independence, standard libraries, etc. One
    can learn from the success of TeX and more specialized software
    like Macaulay2. I do hope that funding agencies are looking into
    this.


Why did you write Sage from scratch, instead of using other existing software and/or libraries?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sage was not written from scratch. Most of its underlying mathematics
functionalities are made possible through FOSS projects such as

* `ATLAS <http://math-atlas.sourceforge.net>`_ --- Automatically Tuned
  Linear Algebra Software.
* `BLAS <https://www.netlib.org/blas>`_ --- Basic Linear Algebra
  Subprograms.
* `ECL <https://common-lisp.net/project/ecl>`_ --- Embeddable Common-Lisp system
* `FLINT <http://www.flintlib.org>`_ --- C library for doing number
  theory.
* `GAP <https://www.gap-system.org>`_ --- a system for computational
  discrete algebra, with particular emphasis on computational group
  theory.
* `GMP <https://gmplib.org>`_ --- GNU Multiple Precision Arithmetic Library.
* `Maxima <http://maxima.sourceforge.net>`_ --- system for symbolic
  and numerical computation.
* `mpmath <https://github.com/fredrik-johansson/mpmath>`_ --- a pure-Python
  library for multiprecision floating-point arithmetic.
* `NumPy and SciPy <https://scipy.org>`_ --- numerical linear algebra and
  other numerical computing capabilities for Python.
* `OpenBLAS <https://www.openblas.net/>`_  --- an optimized BLAS library.
* `Pari/GP <http://pari.math.u-bordeaux.fr>`_ --- a computer algebra
  system for fast computations in number theory.
* `Pynac <http://pynac.sagemath.org>`_ --- a modified version of GiNaC
  that replaces the dependency on CLN by Python.
* `R <http://www.r-project.org>`_ --- a language and environment for
  statistical computing and graphics.
* And many more too numerous to list here.

An up-to-date list can be found in the section
`External Packages <../reference/spkg/index.html>`_
in the Sage Reference Manual.

The principal programming languages of Sage are
`Python <http://www.python.org>`_
and
`Cython <http://www.cython.org>`_.
Python is the primary programming and interfacing language, while
Cython is the primary language for optimizing critical functionalities
and interfacing with C libraries and C extensions for Python. Sage
integrates over 90 FOSS packages into a common interface. On top of
these packages is the Sage library, which consists of over 700,000
lines of new Python and Cython code. See
`openhub.net <https://www.openhub.net/p/sage>`_
for source code analysis of the latest stable Sage release.


How do I get help?
""""""""""""""""""

For support about usage of Sage, there are two options:

* The question-and-answer website `ask.sagemath.org <http://ask.sagemath.org/questions/>`_
* The email list `sage-support <http://groups.google.com/group/sage-support>`_

For support about development of Sage, there is an email list
`sage-devel <http://groups.google.com/group/sage-devel>`_

See http://www.sagemath.org/help.html for a listing of other resources.


Wouldn't it be way better if Sage did not ship as a gigantic bundle?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The SageMath distribution continues to vendor versions of required
software packages ("SPKGs") that work well together.

However, in order to reduce compilation times and the size of the Sage
installation, a development effort ongoing since the 8.x release
series has made it possible to use many system packages provided by
the OS distribution (or by the Homebrew or conda-forge distributions)
instead of building SageMath's own copies.

This so-called "spkg-configure" mechanism runs at the beginning of a
build from source, during the ``./configure`` phase.

To ensure that SageMath builds and runs correctly on a wide variety of
systems, we use automated testing.  See the chapter `Portability
testing <../developer/portability_testing.html>`_ in the Developer's
Guide for details.


With so many bugs in Sage and hundreds of open tickets, why don't you produce a stabilization release?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Any software package contains bug. With something as complex as Sage,
neither the Sage community nor the Sage Development Team make any
claims that Sage is free of bugs. To do so would be an act of
dishonesty.

A Sage release cycle usually lasts for a few months, with several
betas appearing at a 2-3 week intervals.  Each
release cycle is usually chaired by a single release manager who looks
after the Sage merge tree for the duration of the release
cycle. During that time, the release manager often needs to devote the
equivalent of full-time work to quality management and actively
interacts with an international community of Sage users, developers,
and potential contributors.

There have been a number of cases where
two Sage contributors paired up to be the release managers for a Sage
release cycle. However, it is often the case that few people have the
equivalent of 3 weeks' worth of free time to devote to release
management. If you want to help out with release management, please
subscribe to the
`sage-release <http://groups.google.com/group/sage-release>`_
mailing list.

Since the beginning of the Sage project, Sage contributors have tried
to listen and think about what would increase the chances that serious
potential contributors would actually contribute. What encourages one
contributor can discourage another, so tradeoffs need to be made. To
decide that a stabilization release would merge patches with bug
fixes, and only fix bugs, would likely discourage someone from
contributing when they have been told in advance that their positively
reviewed patches will not be merged.

The Sage community believes in
the principle of "release early, release often". How the Sage project
is organized and run differ greatly from that of a commercial software
company. Contributors are all volunteers and this changes the dynamic
of the project dramatically from what it would be if Sage were a
commercial development effort with all developers being full-time
employees.


How can I download the Sage documentation to read it offline?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To download the Sage standard documentation in HTML or PDF formats,
visit the
`Help and Support <https://www.sagemath.org/help.html>`_
page on the Sage website. Each release of Sage comes with the full
documentation that makes up the Sage standard documentation. If you
have downloaded a binary Sage release, the HTML version of the
corresponding documentation comes pre-built and can be found under the
directory ``SAGE_ROOT/local/share/doc/sage/html/``.
During the compilation of Sage from source, the HTML version of the
documentation is also built in the process. To build the HTML version
of the documentation, issue the following command from ``SAGE_ROOT``:

.. CODE-BLOCK:: shell-session

    $ ./sage --docbuild --no-pdf-links all html

Building the PDF version requires that your system has a working LaTeX
installation. To build the PDF version of the documentation, issue the
following command from ``SAGE_ROOT``:

.. CODE-BLOCK:: shell-session

    $ ./sage --docbuild all pdf

For more command line options, refer to the output of any of the
following commands:

.. CODE-BLOCK:: shell-session

    $ ./sage --help
    $ ./sage --advanced
