.. _chapter-faq-usage:

===============
FAQ: Using Sage
===============


How do I get started?
"""""""""""""""""""""

You can try out Sage without downloading anything:

* **CoCalc™:** Go to https://cocalc.com and set up a free account.

  If you log in, you will gain access to the latest version of Sage and to
  many other programs.

  Note that this website is an independent commercial service.

* **Sage cell:** A "one-off" version of Sage, available for doing one
  computation at a time. https://sagecell.sagemath.org/

To download a **pre-built binary** Sage distribution, visit
http://sagemath.org/download.html and click on the link for the binary for your
operating system.

The **source code** of Sage is also available for you to download and use. Go to
http://www.sagemath.org/download-source.html to download the tar archive for any
release of Sage.

The Sage Jupyter notebook runs within a web browser. To start the notebook,
issue the following command in a terminal, if ``sage`` is in your ``PATH``

.. CODE-BLOCK:: shell-session

    $ sage -notebook

What are the prerequisites for installing a copy of Sage on my computer?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Most of the dependencies of Sage are shipped with Sage itself. In most
cases, you can download a pre-built binary and use that without
installing any dependencies. If you use Windows, you will need to
install
`VirtualBox <https://www.virtualbox.org>`_, which can be downloaded
from the page https://www.virtualbox.org/wiki/Downloads. After
installing VirtualBox, you need to download a VirtualBox distribution
of Sage available at
http://www.sagemath.org/download-windows.html. Ensure you follow the
instructions at that page, then start the Sage virtual machine
using the VirtualBox software.

You can get the complete source for Sage to compile it on your own
Linux or Mac OS X system. Sage lives in an isolated directory and does
not interfere with your surrounding system. It ships together with
everything necessary to develop Sage, the source code, all its
dependencies and the complete changelog. On Linux systems like
Debian/Ubuntu, you may have to install the ``build essential``
package and the ``m4`` macro processor. Your system
needs to have a working C compiler if you want to compile Sage
from source. On
Debian/Ubuntu, you can install these prerequisites as follows:

.. CODE-BLOCK:: shell-session

    $ sudo apt-get install build-essential m4

If you have a multi-core system, you can opt for a parallel build of
Sage. The command

.. CODE-BLOCK:: shell-session

    $ export MAKE='make -j8'

will enable 8 threads for parts of the build that support
parallelism. Change the number 8 as appropriate to suit the number of
cores on your system. Some Sage installations may have OpenMP-enabled BLAS
(and other) libraries. The amount of OpenMP parallelism is controlled by
the environment variable OMP_NUM_THREADS; however, it is known to not
play well with Python parallelism, and you might want to

.. CODE-BLOCK:: shell-session

    $ export OMP_NUM_THREADS=1

in case of crashes or hangs.


More details may be found in `Installation Manual <https://doc.sagemath.org/html/en/installation/index.html>`_.


How to get Sage's Python to recognize my system's Tcl/Tk install?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

It may be that you have Tcl/Tk installed and that your system's Python
recognizes it but Sage's Python does not. Normally speaking, there is little
need to build Sage's Python nowadays (anno 2023), but if you do, here it is.
Make sure you installed the Tcl/Tk development library.  On Ubuntu, this is the
command

.. CODE-BLOCK:: shell-session

    $ sudo apt-get install tk8.5-dev

or something along that line. Next, reinstall Sage's Python:

.. CODE-BLOCK:: shell-session

    $ make python3-clean python3-uninstall && make python3

This will pick up the Tcl/Tk library automatically. After successfully
reinstalling Sage's Python, from within the Sage command line interface,
issue these commands:

.. CODE-BLOCK:: python

    import _tkinter
    import Tkinter

If they do not raise an :class:`ImportError` then it worked.


How do I import Sage into a Python script?
""""""""""""""""""""""""""""""""""""""""""

You can import Sage as a library in a Python script. One caveat is
that you need to run that Python script using the version of Python
that is bundled with Sage (Sage 9.2 ships with Python 3.7.x).
To import Sage, put the following in your Python script:

.. CODE-BLOCK:: python

    from sage.all import *

When you want to run your script, you need to invoke Sage with the
option ``-python`` which would run your script using the version of
Python that comes with Sage. For example, if Sage is in your ``PATH``
variable then you can do this:

.. CODE-BLOCK:: shell-session

    $ sage -python /path/to/my/script.py

Another way is to write a Sage script and run that script using Sage
itself. A Sage script has the file extension ``.sage`` and is more or
less a Python script but uses Sage-specific functions and
commands. You can then run that Sage script like so:

.. CODE-BLOCK:: shell-session

    $ sage /path/to/my/script.sage

This will take care of loading the necessary environment variables and
default imports for you.

How can I reload a Python script in a Sage session?
"""""""""""""""""""""""""""""""""""""""""""""""""""

You can load a Python script in a Sage session with the command
**load**. For example, we could use Sage to import a file called
simple.py with:

.. CODE-BLOCK:: python

    load("simple.py")

and repeat this command every time that we change the file simple.py. However, if we type:

.. CODE-BLOCK:: python

    attach("simple.py")

every change applied to the file simple.py will be automatically updated in Sage.

Can I use SageMath with Python 3.x?
"""""""""""""""""""""""""""""""""""

Since release 9.0 from January 2020, SageMath is running on top of Python 3.


I used XXX to install Sage X.Y and that version is giving lots of errors. What can I do?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The version of Sage, i.e. Sage version X.Y, that is available on your XXX system
through its package manager, is very old. No one has yet
found time to update the XXX version of Sage. Any help is
greatly appreciated. You should download the latest version of Sage
from the
`download page <http://www.sagemath.org/download.html>`_.
If you would like to help with updating the XXX version of
Sage, please email the
`sage-devel <https://groups.google.com/group/sage-devel>`_
mailing list.


Should I use the official version or development version?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You are encouraged to use the latest official version of
Sage. Development versions are frequently announced on the
`sage-devel <https://groups.google.com/group/sage-devel>`_
and
`sage-release <https://groups.google.com/group/sage-release>`_
mailing lists. An easy way of helping out with Sage development is to
download the latest development release, compile it on your system,
run all doctests, and report any compilation errors or doctest
failures.


Is Sage difficult to learn?
"""""""""""""""""""""""""""

Basic features of Sage should be as easy to learn as learning the
basics of Python. Numerous tutorials are available online to help you
learn Sage. To get the most out of Sage, you are encouraged to learn
some features of the Python programming language. Here is an
incomplete list of resources on Python. Further resources can be found
by a web search.

* `Building Skills in Python <http://homepage.mac.com/s_lott/books/python.html>`_
  by Steven F. Lott
* `Dive into Python <http://www.diveintopython.net>`_ by Mark Pilgrim
* `How to Think Like a Computer Scientist <http://www.openbookproject.net/thinkCSpy>`_
  by Jeffrey Elkner, Allen B. Downey, and Chris Meyers
* `Official Python Tutorial <https://docs.python.org/tutorial>`_
* `Python <https://www.python.org>`_ home page and the
  `Python standard documentation <https://docs.python.org>`_


Can I do X in Sage?
"""""""""""""""""""

You are encouraged to use Sage's tab autocompletion. Just type a few
characters, hit the :kbd:`Tab` key, and see if the command you want appears
in the list of tab autocompletion. If you have a command called
``mycmd``, then type ``mycmd.`` and hit the :kbd:`Tab` key to get a list of
functionalities that are supported by that command. To read the
documentation of ``mycmd``, type ``mycmd?`` and press the :kbd:`Enter` key to
read the documentation for that command. Similarly, type ``mycmd??``
and hit the :kbd:`Enter` key to get the source code of that command. You are
also encouraged to search through the source code and documentation of
the Sage library. To search through the source code of the Sage
library, use the command ``search_src("<search-keyword>")`` where you
should replace ``<search-keyword>`` with the key words you are looking
for. Also, you can search through the documentation of the Sage
library using the command ``search_doc("<search-keyword>")``.


What exactly does Sage do when I type "0.6**2"?
"""""""""""""""""""""""""""""""""""""""""""""""

When you type "0.6**2" in Python, it returns something like
0.35999999999999999. But when you do the same in Sage it returns
0.360000000000000. To understand why Python behaves as it does, see
the
`Python Tutorial <https://docs.python.org/tutorial/floatingpoint.html>`_,
especially the chapter "Floating Point Arithmetic: Issues and
Limitations". What Sage does is "preparse" the input and transforms it
like this::

    sage: preparse("0.6**2")
    "RealNumber('0.6')**Integer(2)"

So what is *actually* run is:

.. CODE-BLOCK:: python

    RealNumber('0.6')**Integer(2)

The Sage developers (in fact, Carl Witty) decided that Sage floating
point numbers should by default print only the known correct decimal
digits, when possible, thus skirting the problem that Python has. This
decision has its pros and cons. Note that ``RealNumber`` and
``Integer``  are Sage specific, so you would not be able to just type
the above into Python and expect it to work without first an import
statement such as:

.. CODE-BLOCK:: python

    from sage.all import RealNumber, Integer, preparse


Why is Sage's command history different from Magma's?
"""""""""""""""""""""""""""""""""""""""""""""""""""""

Using Sage, you are missing a feature of the Magma command line
interface. In Magma, if you enter a line found in history using up
arrow key and then press down arrow key, then the next line in history
is fetched. This feature allows you to fetch as many successive lines
in history as you like. However, Sage does not have a similar
feature. The
`IPython <https://ipython.org>`_
command prompt uses the readline library (via pyreadline), which
evidently does not support this feature. Magma has its own custom
"readline-like" library, which does support this feature. (Since so
many people have requested this feature, if anybody can figure out how
to implement it, then such an implementation would certainly be
welcome!)


I have type issues using SciPy, cvxopt or NumPy from Sage.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You are using SciPy or cvxopt or NumPy from Sage and you get type
errors, e.g.

.. CODE-BLOCK:: text

    TypeError: function not supported for these types, and can't coerce safely to supported types.

When you type in numbers into Sage, the pre-processor converts them to
a base ring, which you can see by doing::

    sage: preparse("stats.uniform(0,15).ppf([0.5,0.7])")
    "stats.uniform(Integer(0),Integer(15)).ppf([RealNumber('0.5'),RealNumber('0.7')])"

Unfortunately, NumPy support of these advanced Sage types like
``Integer`` or ``RealNumber`` is not yet at 100%. As a solution,
redefine ``RealNumber`` and/or ``Integer`` to change the behavior of
the Sage preparser, so decimal literals are floats instead of Sage
arbitrary precision real numbers, and integer literals are Python
ints. For example::

    sage: RealNumber = float; Integer = int
    sage: from scipy import stats
    sage: stats.ttest_ind([1,2,3,4,5], [2,3,4,5,.6])
    Ttest...Result(statistic=...0.0767529..., pvalue=...0.940704...)
    sage: stats.uniform(0,15).ppf([0.5,0.7])
    array([  7.5,  10.5])

Alternatively, be explicit about data types, e.g. ::

    sage: from scipy import stats
    sage: stats.uniform(int(0),int(15)).ppf([float(0.5),float(0.7)])
    array([  7.5,  10.5])

As a third alternative, use the raw suffix::

    sage: from scipy import stats
    sage: stats.uniform(0r,15r).ppf([0.5r,0.7r])
    array([  7.5,  10.5])

You can also disable the preprocessor in your code via
``preparser(False)``. You can start IPython alone from the command
line ``sage -ipython`` which does not pre-load anything
Sage-specific. Or switch the Notebook language to "Python".


How do I save an object so I don't have to compute it each time I open a worksheet?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The ``save`` and ``load`` commands will save and load an object,
respectively.

Does Sage contain a function similar to Mathematica's ToCharacterCode[]?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You might want to convert ASCII characters such as "Big Mac" to ASCII
numerals for further processing. In Sage and Python, you can use ``ord``,
e.g. ::

    sage: list(map(ord, "abcde"))
    [97, 98, 99, 100, 101]
    sage: list(map(ord, "Big Mac"))
    [66, 105, 103, 32, 77, 97, 99]

How can I write multiplication implicitly as in Mathematica?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sage has a function that enables this::

    sage: implicit_multiplication(True)
    sage: x 2 x  # not tested
    2*x^2
    sage: implicit_multiplication(False)

This is preparsed by Sage into Python code. It may not work in a
complicated situation. To see what the preparser does::

    sage: implicit_multiplication(True)
    sage: preparse("2 x")
    'Integer(2)*x'
    sage: implicit_multiplication(False)
    sage: preparse("2 x")
    'Integer(2) x'

See https://wiki.sagemath.org/sage_mathematica for more information
about Mathematica vs. SageMath.

Can I make Sage automatically execute commands on startup?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Yes, just make a file ``$HOME/.sage/init.sage`` and it will be
executed any time you start Sage. This assumes that the Sage
environment variable ``DOT_SAGE`` points to the hidden directory
``$HOME/.sage``, which by default is the case.


When I compile Sage my computer beeps and shuts down or hangs.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Compiling Sage is quite taxing on the CPU. The above behavior usually
indicates that your computer has overheated. In many cases this can be
fixed by cleaning the CPU fan and assuring proper ventilation of the
system. Please ask your system administrator or a professional to do
this in case you have never done this. Such hardware maintenance, if
not performed by a skilled professional, you can potentially damage
your system.

For Linux users, if you suspect that the compilation fails because of
a resource issue, a fix might be to edit your ``/etc/inittab`` so that
Linux boots into run level 3. The file ``/etc/inittab`` usually
contains something similar to the following snippet:

.. CODE-BLOCK:: text

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

which directs your Linux distribution to boot into a graphical login
screen. Comment out the line ``id:5:initdefault:`` and add the line
``id:3:initdefault:``, so that you now have something like:

.. CODE-BLOCK:: text

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

Now if you reboot your system, you will be greeted with a text based
login screen. This allows you to log into your system with a text
based session from within a virtual terminal. A text based session
usually does not consume as much system resources as would be the case
with a graphical session. Then build your Sage source distribution
from within your text based session. You need to make sure that you
can first restore your graphical session, before you attempt to log
into a text based session.


When I start Sage, SELinux complains that "/path/to/libpari-gmp.so.2" requires text-relocation. How can I fix it?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The problem can be fixed by running the following command:

.. CODE-BLOCK:: shell-session

    $ chcon -t textrel_shlib_t /path/to/libpari-gmp.so.2


Upgrading Sage went fine, but now the banner still shows the old version. How can I fix this?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The banner is stored and not computed at every new start of Sage. If
it has not been updated, this should not prevent Sage to run
correctly. Type ``banner()`` in a Sage session to check the real
version. If you want the correct banner, you need to build Sage again
by typing ``make build`` in a terminal.


How do I run sage in daemon mode, i.e. as a service?
""""""""""""""""""""""""""""""""""""""""""""""""""""

There are several possibilities. Use ``screen``, ``nohup`` or ``disown``.


May I use Sage tools in a commercial environment?
"""""""""""""""""""""""""""""""""""""""""""""""""

Yes! Absolutely! Basically the *only* constraint is that if you make
changes to Sage itself and redistribute this changed version of Sage
publicly, then you must make these changes available to us so that we
can put them into the standard version of Sage (if we
want). Otherwise, you are free to use as many copies of Sage as you
want completely for free to make money, etc. without paying any
license fees at all.


I want to write some Cython code that uses finite field arithmetic but "cimport sage.rings.finite_field_givaro" fails. What can I do?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You need to give hints to Sage so that it uses C++ (both Givaro and
NTL are C++ libraries), and it also needs the GMP and STDC C++
libraries. Here is a small example:

.. CODE-BLOCK:: cython

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

To find out more, type

.. CODE-BLOCK:: ipython

    sage.rings.finite_field_givaro.FiniteField_givaro.

at the Sage prompt and hit tab, then use ``??`` to get more
information on each function. For example:

.. CODE-BLOCK:: ipython

    sage.rings.finite_field_givaro.FiniteField_givaro.one??

tells you more about the multiplicative unit element in the finite
field.


I'm getting weird build failures on Mac OS X. How do I fix this?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Search the build log (install.log) to see if you are getting the
following log message:

.. CODE-BLOCK:: text

    fork: Resource temporarily unavailable.

If so, try the following. Create (or edit) ``/etc/launchd.conf`` and
include the following:

.. CODE-BLOCK:: text

    limit maxproc 512 2048

then reboot. See
`this page <http://www.macosxhints.com/article.php?story=20050709233920660>`_
for more details.

How do I plot the cube root (or other odd roots) for negative input?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This is one of the most frequently asked questions.  There are several
methods mentioned in the plot documentation, but this one is easiest::

    sage: plot(real_nth_root(x, 3), (x, -1, 1))
    Graphics object consisting of 1 graphics primitive

On the other hand, note that the straightforward ::

    sage: plot(x^(1/3), (x, -1, 1))  # not tested

produces the expected plot only for positive `x`. The *reason* is that Sage
returns complex numbers for odd roots of negative numbers when numerically
approximated, which is a 
:wikipedia:`standard convention <Cube_root#Complex_numbers>`. ::

    sage: numerical_approx( (-1)^(1/3) )
    0.500000000000000 + 0.866025403784439*I

How do I use the bitwise XOR operator in Sage?
""""""""""""""""""""""""""""""""""""""""""""""

The exclusive or operator in Sage is ``^^``. This also works for
the inplace operator ``^^=``::

   sage: 3^^2
   1
   sage: a = 2
   sage: a ^^= 8
   sage: a
   10

If you define two variables and then evaluate as follows::

    sage: a = 5; b = 8
    sage: a.__xor__(b), 13
    (13, 13)

You can also do ::

    sage: (5).__xor__(8)
    13

The parentheses are necessary so that Sage does not think you have a
real number. There are several ways to define a function::

    sage: xor = lambda x, y: x.__xor__(y)
    sage: xor(3, 8)
    11

Another option, which sneaks around the Sage
preparser, is ::

    sage: def xor(a, b):
    ....:     return eval("%s^%s" % (a, b))
    sage: xor(3, 8)
    11

You can also turn off the Sage preparser with ``preparser(False)``,
then ``^`` will work just like in Python. You can later turn on the
preparser with ``preparser(True)``. That only works in command line
Sage. In a notebook, switch to Python mode.

With objects a and b and a function f, I accidentally typed f(a) = b instead of f(a) == b. This returned a TypeError (as expected), but also deleted the object a. Why?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

It is because of how functions are defined in Sage with the
``f(x) = expr`` notation using the preparser. Also notice that if you
make this mistake inside of an ``if`` statement, you will get a
:class:`SyntaxError` before anything else goes wrong. So in this case,
there is no problem.


How do I use a different browser with the Sage notebook?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You will need to do this from the command line.  Just run a command like this.

* Linux (assuming you have Sage in ``/usr/bin``):

  .. CODE-BLOCK:: shell-session

    $ env BROWSER=opera /usr/bin/sage --notebook

* Mac (assuming you are in the directory of your downloaded Sage).
  With the Jupyter notebook:

  .. CODE-BLOCK:: shell-session

    $ BROWSER='open -a Firefox %s' ./sage --notebook jupyter
    $ BROWSER='open -a Google\ Chrome %s' ./sage --notebook jupyter


Where is the source code for ``<function>``?
""""""""""""""""""""""""""""""""""""""""""""

Functions and classes written in Python or Cython are in general accessible
on the IPython command line with the ``??`` shortcut::

    sage: plot??                            # not tested
    Signature: plot(*args, **kwds)
    Source:
    ...

Objects that are built into Python or IPython are compiled and will
not show, however. There are many functions in Sage implemented as
symbolic functions, i.e., they can be used unevaluated as part of
symbolic expressions. Their source code may also not be accessible
from the command line, especially with elementary functions, because
they are coded in C++ for efficiency reasons.
