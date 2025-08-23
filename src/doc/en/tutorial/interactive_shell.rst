.. _chapter-interactive_shell:

*********************
The Interactive Shell
*********************
In most of this tutorial, we assume you start the Sage interpreter
using the ``sage`` command. This starts a customized version of the
IPython shell, and imports many functions and classes, so they are
ready to use from the command prompt. Further customization is
possible by editing the ``$SAGE_ROOT/ipythonrc`` file. Upon starting
Sage, you get output similar to the following:

.. CODE-BLOCK:: text

    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘


    sage:

To quit Sage either press Ctrl-D or type
``quit`` or ``exit``.

.. skip

::

    sage: quit
    Exiting Sage (CPU time 0m0.00s, Wall time 0m0.89s)

The wall time is the time that elapsed on the clock hanging from
your wall. This is relevant, since CPU time does not track time
used by subprocesses like GAP or Singular.

(Avoid killing a Sage process with ``kill -9`` from a terminal,
since Sage might not kill child processes, e.g.,
Maple processes, or cleanup temporary files from
``$HOME/.sage/tmp``.)

Your Sage Session
=================

The session is the sequence of input and output
from when you start Sage until you quit. Sage logs all Sage input,
via IPython. In fact, if you're using the interactive shell (not the
notebook interface), then at any point you may type ``%history`` (or ``%hist``) to
get a listing of all input lines typed so far. You can type ``?`` at
the Sage prompt to find out more about IPython, e.g.,
"IPython offers numbered prompts ... with input and output
caching. All input is saved and can be retrieved as variables (besides
the usual arrow key recall). The following GLOBAL variables always
exist (so don't overwrite them!)":

.. CODE-BLOCK:: text

      _:  previous input (interactive shell and notebook)
      __: next previous input (interactive shell only)
      _oh : list of all inputs (interactive shell only)

Here is an example:

.. skip

::

    sage: factor(100)
     _1 = 2^2 * 5^2
    sage: kronecker_symbol(3,5)
     _2 = -1
    sage: %hist   #This only works from the interactive shell, not the notebook.
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    sage: _oh
     _4 = {1: 2^2 * 5^2, 2: -1}
    sage: _i1
     _5 = 'factor(ZZ(100))\n'
    sage: eval(_i1)
     _6 = 2^2 * 5^2
    sage: %hist
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    4: _oh
    5: _i1
    6: eval(_i1)
    7: %hist

We omit the output numbering in the rest of this tutorial and the
other Sage documentation.

You can also store a list of input from session in a macro for that
session.

.. skip

::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: M = ModularSymbols(37)
    sage: %hist
    1: E = EllipticCurve([1,2,3,4,5])
    2: M = ModularSymbols(37)
    3: %hist
    sage: %macro em 1-2
    Macro `em` created. To execute, type its name (without quotes).


.. skip

::

    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field
    sage: E = 5
    sage: M = None
    sage: em
    Executing Macro...
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field

When using the interactive shell, any UNIX shell command can be
executed from Sage by prefacing it by an exclamation point ``!``. For
example,

.. skip

::

    sage: !ls
    auto  example.sage glossary.tex  t  tmp  tut.log  tut.tex

returns the listing of the current directory.

The ``PATH`` has the Sage bin directory at the front, so if you run ``gp``,
``gap``, ``singular``, ``maxima``, etc., you get the versions included
with Sage.

.. skip

::

    sage: !gp
    Reading GPRC: /etc/gprc ...Done.

                               GP/PARI CALCULATOR Version 2.2.11 (alpha)
                      i686 running linux (ix86/GMP-4.1.4 kernel) 32-bit version
    ...
    sage: !singular
                         SINGULAR                             /  Development
     A Computer Algebra System for Polynomial Computations   /   version 3-0-1
                                                           0<
         by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   October 2005
    FB Mathematik der Universitaet, D-67653 Kaiserslautern    \

Logging Input and Output
========================

Logging your Sage session is not the same as saving it (see
:ref:`section-save` for that). To log input (and optionally output) use the
``logstart`` command. Type ``logstart?`` for more details. You can use
this command to log all input you type, all output, and even play
back that input in a future session (by simply reloading the log
file).

.. skip

.. CODE-BLOCK:: shell-session

    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: logstart setup
    Activating auto-logging. Current session state plus future input saved.
    Filename       : setup
    Mode           : backup
    Output logging : False
    Timestamping   : False
    State          : active
    sage: E = EllipticCurve([1,2,3,4,5]).minimal_model()
    sage: F = QQ^3
    sage: x,y = QQ['x,y'].gens()
    sage: G = E.gens()
    sage:
    Exiting Sage (CPU time 0m0.61s, Wall time 0m50.39s).
    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: load("setup")
    Loading log file <setup> one line at a time...
    Finished replaying log file <setup>
    sage: E
    Elliptic Curve defined by y^2 + x*y  = x^3 - x^2 + 4*x + 3 over Rational
    Field
    sage: x*y
    x*y
    sage: G
    [(2 : 3 : 1)]

If you use Sage in the Linux KDE
terminal ``konsole`` then you can save your session as follows: after
starting Sage in ``konsole``, select "settings", then "history...",
then "set unlimited". When you are ready to save your session,
select "edit" then "save history as..." and type in a name to save
the text of your session to your computer. After saving this file,
you could then load it into an editor, such as xemacs, and print
it.

Paste Ignores Prompts
=====================

Suppose you are reading a session of Sage or Python computations
and want to copy them into Sage. But there are annoying ``>>>`` or
``sage:`` prompts to worry about. In fact, you can copy and paste an
example, including the prompts if you want, into Sage. In other
words, by default the Sage parser strips any leading ``>>>`` or
``sage:`` prompt before passing it to Python. For example,

.. skip

::

    sage: 2^10
    1024
    sage: sage: sage: 2^10
    1024
    sage: >>> 2^10
    1024

Timing Commands
===============

If you place the ``%time`` command at the beginning of an input line,
the time the command takes to run will be displayed after the
output. For example, we can compare the running time for a certain
exponentiation operation in several ways. The timings below will
probably be much different on your computer, or even between
different versions of Sage. First, native Python:

.. skip

::

    sage: %time a = int(1938)^int(99484)
    CPU times: user 0.66 s, sys: 0.00 s, total: 0.66 s
    Wall time: 0.66

This means that 0.66 seconds total were taken, and the "Wall time",
i.e., the amount of time that elapsed on your wall clock, is also
0.66 seconds. If your computer is heavily loaded with other
programs, the wall time may be much larger than the CPU time.

It's also possible to use the ``timeit`` function to try to get
timing over a large number of iterations of a command.  This gives
slightly different information, and requires the input of a string
with the command you want to time.

.. skip

::

    sage: timeit("int(1938)^int(99484)")
    5 loops, best of 3: 44.8 ms per loop

Next we time exponentiation using the native Sage Integer type,
which is implemented (in Cython) using the GMP library:

.. skip

::

    sage: %time a = 1938^99484
    CPU times: user 0.04 s, sys: 0.00 s, total: 0.04 s
    Wall time: 0.04

Using the PARI C-library interface:

.. skip

::

    sage: %time a = pari(1938)^pari(99484)
    CPU times: user 0.05 s, sys: 0.00 s, total: 0.05 s
    Wall time: 0.05

GMP is better, but only slightly (as expected, since the version of
PARI built for Sage uses GMP for integer arithmetic).

You can also time a block of commands using
the ``cputime`` command, as illustrated below:

::

    sage: t = cputime()
    sage: a = int(1938)^int(99484)
    sage: b = 1938^99484
    sage: c = pari(1938)^pari(99484)
    sage: cputime(t)                       # somewhat random output
    0.64

.. skip

::

    sage: cputime?
    ...
        Return the time in CPU second since Sage started, or with optional
        argument t, return the time since time t.
        INPUT:
            t -- (optional) float, time in CPU seconds
        OUTPUT:
            float -- time in CPU seconds

The ``walltime`` command behaves just like the ``cputime`` command,
except that it measures wall time.

We can also compute the above power in some of the computer algebra
systems that Sage includes. In each case we execute a trivial command in
the system, in order to start up the server for that program. The
most relevant time is the wall time. However, if there is a
significant difference between the wall time and the CPU time then
this may indicate a performance issue worth looking into.

.. skip

::

    sage: time 1938^99484;
    CPU times: user 0.01 s, sys: 0.00 s, total: 0.01 s
    Wall time: 0.01
    sage: gp(0)
    0
    sage: time g = gp('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: maxima(0)
    0
    sage: time g = maxima('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.30
    sage: kash(0)
    0
    sage: time g = kash('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: mathematica(0)
            0
    sage: time g = mathematica('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.03
    sage: maple(0)
    0
    sage: time g = maple('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.11
    sage: libgap(0)
    0
    sage: time g = libgap.eval('1938^99484;')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 1.02

Note that GAP and Maxima are the slowest in this test (this was run
on the machine ``sage.math.washington.edu``). Because of the pexpect
interface overhead, it is perhaps unfair to compare these to Sage,
which is the fastest.

Other IPython tricks
====================

As noted above, Sage uses IPython as its front end, and so you can use
any of IPython's commands and features.  You can read the `full
IPython documentation <http://ipython.scipy.org/moin/Documentation>`_.
Meanwhile, here are some fun tricks -- these are called "Magic
commands" in IPython:

- You can use ``%edit`` (or ``%ed`` or ``ed``) to open an editor, if
  you want to type in some complex code.  Before you start Sage, make
  sure that the :envvar:`EDITOR` environment variable is set to your
  favorite editor (by putting ``export EDITOR=/usr/bin/emacs`` or
  ``export EDITOR=/usr/bin/vim`` or something similar in the
  appropriate place, like a ``.profile`` file).  From the Sage prompt,
  executing ``%edit`` will open up the named editor.  Then within the
  editor you can define a function:

  .. CODE-BLOCK:: python

    def some_function(n):
        return n**2 + 3*n + 2

  Save and quit from the editor.  For the rest of your Sage session,
  you can then use ``some_function``.  If you want to modify it, type
  ``%edit some_function`` from the Sage prompt.

- If you have a computation and you want to modify its output for
  another use, perform the computation and type ``%rep``: this will
  place the output from the previous command at the Sage prompt, ready
  for you to edit it. ::

    sage: f(x) = cos(x)
    sage: f(x).derivative(x)
    -sin(x)

  At this point, if you type ``%rep`` at the Sage prompt, you will get
  a new Sage prompt, followed by ``-sin(x)``, with the cursor at the
  end of the line.

For more, type ``%quickref`` to get a quick reference guide to
IPython.  As of this writing (April 2011), Sage uses version 0.9.1 of
IPython, and the `documentation for its magic commands
<http://ipython.org/ipython-doc/dev/interactive/tutorial.html#magic-functions>`_
is available online. Various slightly advanced aspects of magic command system are documented `here <http://ipython.org/ipython-doc/stable/interactive/reference.html#magic-command-system>`_ in IPython.

Errors and Exceptions
=====================

When something goes wrong, you will usually see a Python
"exception". Python even tries to suggest what raised the
exception. Often you see the name of the exception, e.g.,
:class:`NameError` or :class:`ValueError` (see the Python Library Reference [PyLR]_
for a complete list of exceptions). For example,

.. skip

::

    sage: 3_2
    ------------------------------------------------------------
       File "<console>", line 1
         ZZ(3)_2
               ^
    SyntaxError: invalid ...

    sage: EllipticCurve([0,infinity])
    ------------------------------------------------------------
    Traceback (most recent call last):
    ...
    TypeError: Unable to coerce Infinity (<class 'sage...Infinity'>) to Rational

The interactive debugger is sometimes useful for understanding what
went wrong. You can toggle it on or off using ``%pdb`` (the
default is off). The prompt ``ipdb>`` appears if an exception is
raised and the debugger is on. From within the debugger, you can
print the state of any local variable, and move up and down the
execution stack. For example,

.. skip

::

    sage: %pdb
    Automatic pdb calling has been turned ON
    sage: EllipticCurve([1,infinity])
    ---------------------------------------------------------------------------
    <class 'exceptions.TypeError'>             Traceback (most recent call last)
    ...

    ipdb>

For a list of commands in the debugger, type ``?`` at the ``ipdb>``
prompt:

.. CODE-BLOCK:: text

    ipdb> ?

    Documented commands (type help <topic>):
    ========================================
    EOF    break  commands   debug    h       l     pdef   quit    tbreak
    a      bt     condition  disable  help    list  pdoc   r       u
    alias  c      cont       down     ignore  n     pinfo  return  unalias
    args   cl     continue   enable   j       next  pp     s       up
    b      clear  d          exit     jump    p     q      step    w
    whatis where

    Miscellaneous help topics:
    ==========================
    exec  pdb

    Undocumented commands:
    ======================
    retval  rv

Type Ctrl-D or ``quit`` to return to Sage.

.. _section-tabcompletion:

Reverse Search and Tab Completion
=================================

Reverse search:
Type the beginning of a command, then ``Ctrl-p`` (or just hit the up
arrow key) to go back to each line you have entered that begins in
that way. This works even if you completely exit Sage and restart
later. You can also do a reverse search through the history using
``Ctrl-r``. All these features use the ``readline`` package, which is
available on most flavors of Linux.

To illustrate tab completion,
first create the three dimensional vector space
:math:`V=\QQ^3` as follows:

::

    sage: V = VectorSpace(QQ,3)
    sage: V
    Vector space of dimension 3 over Rational Field

You can also use the following more concise notation:

::

    sage: V = QQ^3

Then it is easy to list all member functions for :math:`V` using tab
completion. Just type ``V.``, then type the :kbd:`Tab` key on your
keyboard:

.. skip

::

    sage: V.[tab key]
    V._VectorSpace_generic__base_field
    ...
    V.ambient_space
    V.base_field
    V.base_ring
    V.basis
    V.coordinates
    ...
    V.zero_vector

If you type the first few letters of a function, then the :kbd:`Tab` key,
you get only functions that begin as indicated.

.. skip

::

    sage: V.i[tab key]
    V.is_ambient  V.is_dense    V.is_full     V.is_sparse

If you wonder what a particular function does, e.g., the
coordinates function, type ``V.coordinates?`` for help or
``V.coordinates??`` for the source code, as explained in the next
section.


Integrated Help System
======================

Sage features an integrated help facility. Type a function name
followed by ? for the documentation for that function.

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates?
    Type:           instancemethod
    Base Class:     <class 'instancemethod'>
    String Form:    <bound method FreeModule_ambient_field.coordinates of Vector
    space of dimension 3 over Rational Field>
    Namespace:      Interactive
    File:           /home/was/s/local/lib/python2.4/site-packages/sage/modules/f
    ree_module.py
    Definition:     V.coordinates(self, v)
    Docstring:
        Write v in terms of the basis for self.

        Returns a list c such that if B is the basis for self, then

                sum c_i B_i = v.

        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: M = FreeModule(IntegerRing(), 2); M0,M1=M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinates(2*M0-M1)
            [2, -1]

As shown above, the output tells you the type of the object, the
file in which it is defined, and a useful description of the
function with examples that you can paste into your current
session. Almost all of these examples are regularly automatically
tested to make sure they work and behave exactly as claimed.

Another feature that is very much in the spirit of the open source
nature of Sage is that if ``f`` is a Python function, then typing ``f??``
displays the source code that defines ``f``. For example,

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates??
    Type:           instancemethod
    ...
    Source:
    def coordinates(self, v):
            """
            Write $v$ in terms of the basis for self.
            ...
            """
            return self.coordinate_vector(v).list()

This tells us that all the ``coordinates`` function does is call the
``coordinate_vector`` function and change the result into a list.
What does the ``coordinate_vector`` function do?

.. skip

::

    sage: V = QQ^3
    sage: V.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            ...
            return self.ambient_vector_space()(v)

The ``coordinate_vector`` function coerces its input into the
ambient space, which has the effect of computing the vector of
coefficients of :math:`v` in terms of :math:`V`. The space
:math:`V` is already ambient since it's just :math:`\QQ^3`.
There is also a ``coordinate_vector`` function for subspaces, and
it's different. We create a subspace and see:

.. skip

::

    sage: V = QQ^3; W = V.span_of_basis([V.0, V.1])
    sage: W.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            """
             ...
            """
            # First find the coordinates of v wrt echelon basis.
            w = self.echelon_coordinate_vector(v)
            # Next use transformation matrix from echelon basis to
            # user basis.
            T = self.echelon_to_user_matrix()
            return T.linear_combination_of_rows(w)

(If you think the implementation is inefficient, please sign up to
help optimize linear algebra.)

You may also type ``help(command_name)`` or ``help(class)`` for a
manpage-like help file about a given class.

.. skip

::

    sage: help(VectorSpace)
    Help on function VectorSpace in module sage.modules.free_module:

    VectorSpace(K, dimension_or_basis_keys=None, sparse=False, inner_product_matrix=None, *,
                with_basis='standard', dimension=None, basis_keys=None, **args)
    EXAMPLES:

    The base can be complicated, as long as it is a field.

    ::

        sage: V = VectorSpace(FractionField(PolynomialRing(ZZ,'x')),3)
        sage: V
        Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in x
         over Integer Ring
        sage: V.basis()
        [
        (1, 0, 0),
        (0, 1, 0),
    --More--

When you type ``q`` to exit the help system, your session appears
just as it was. The help listing does not clutter up your session,
unlike the output of ``function_name?`` sometimes does. It's
particularly helpful to type ``help(module_name)``. For example,
vector spaces are defined in ``sage.modules.free_module``, so type
``help(sage.modules.free_module)`` for documentation about that
whole module. When viewing documentation using help, you can search
by typing ``/`` and in reverse by typing ``?``.

Saving and Loading Individual Objects
=====================================

Suppose you compute a matrix or worse, a complicated space of
modular symbols, and would like to save it for later use. What can
you do? There are several approaches that computer algebra systems
take to saving individual objects.


#. **Save your Game:** Only support saving and loading of complete
   sessions (e.g., GAP, Magma).

#. **Unified Input/Output:** Make every object print in a way that
   can be read back in (GP/PARI).

#. **Eval**: Make it easy to evaluate arbitrary code in the
   interpreter (e.g., Singular, PARI).


Because Sage uses Python, it takes a different approach, which is that
every object can be serialized, i.e., turned into a string from
which that object can be recovered. This is in spirit similar to
the unified I/O approach of PARI, except it doesn't have the
drawback that objects print to screen in too complicated of a way.
Also, support for saving and loading is (in most cases) completely
automatic, requiring no extra programming; it's simply a feature of
Python that was designed into the language from the ground up.

Almost all Sage objects x can be saved in compressed form to disk using
``save(x, filename)`` (or in many cases ``x.save(filename)``). To load
the object back in, use ``load(filename)``.

.. skip

::

    sage: A = MatrixSpace(QQ,3)(range(9))^2
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]
    sage: save(A, 'A')

You should now quit Sage and restart. Then you can get ``A`` back:

.. skip

::

    sage: A = load('A')
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]

You can do the same with more complicated objects, e.g., elliptic
curves. All data about the object that is cached is stored with the
object. For example,

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: v = E.anlist(100000)              # takes a while
    sage: save(E, 'E')
    sage: quit

The saved version of ``E`` takes 153 kilobytes, since it stores the
first 100000 :math:`a_n` with it.

.. skip

::

    ~/tmp$ ls -l E.sobj
    -rw-r--r--  1 was was 153500 2006-01-28 19:23 E.sobj
    ~/tmp$ sage [...]
    sage: E = load('E')
    sage: v = E.anlist(100000)              # instant!

(In Python, saving and loading is accomplished using
the ``cPickle`` module.   In particular, a Sage object ``x``
can be saved via ``cPickle.dumps(x, 2)``.  Note the ``2``!)

Sage cannot save and load individual objects created in some other
computer algebra systems, e.g., GAP, Singular, Maxima, etc. They
reload in a state marked "invalid". In GAP, though many objects
print in a form from which they can be reconstructed, many don't,
so reconstructing from their print representation is purposely not
allowed.

.. skip

::

    sage: a = libgap(2)
    sage: a.save('a')
    sage: load('a')
    Traceback (most recent call last):
    ...
    ValueError: The session in which this object was defined is no longer
    running.

GP/PARI objects can be saved and loaded since their print
representation is enough to reconstruct them.

.. skip

::

    sage: a = gp(2)
    sage: a.save('a')
    sage: load('a')
    2

Saved objects can be re-loaded later on computers with different
architectures or operating systems, e.g., you could save a huge
matrix on 32-bit OS X and reload it on 64-bit Linux, find the
echelon form, then move it back. Also, in many cases you can even
load objects into versions of Sage that are different than the versions
they were saved in, as long as the code for that object isn't too
different. All the attributes of the objects are saved, along with
the class (but not source code) that defines the object. If that
class no longer exists in a new version of Sage, then the object can't be
reloaded in that newer version. But you could load it in an old
version, get the objects dictionary (with ``x.__dict__``), and
save the dictionary, and load that into the newer version.

Saving as Text
--------------

You can also save the ASCII text representation of objects to a
plain text file by simply opening a file in write mode and writing
the string representation of the object (you can write many objects
this way as well). When you're done writing objects, close the
file.

.. skip

::

    sage: R.<x,y> = PolynomialRing(QQ,2)
    sage: f = (x+y)^7
    sage: o = open('file.txt','w')
    sage: o.write(str(f))
    sage: o.close()

.. _section-save:

Saving and Loading Complete Sessions
====================================

Sage has very flexible support for saving and loading complete
sessions.

The command ``save_session(sessionname)`` saves all the variables
you've defined in the current session as a dictionary in the given
``sessionname``. (In the rare case when a variable does not support
saving, it is simply not saved to the dictionary.) The resulting
file is an ``.sobj`` file and can be loaded just like any other
object that was saved. When you load the objects saved in a
session, you get a dictionary whose keys are the variables names
and whose values are the objects.

You can use the ``load_session(sessionname)`` command to load the
variables defined in ``sessionname`` into the current session. Note
that this does not wipe out variables you've already defined in
your current session; instead, the two sessions are merged.

First we start Sage and define some variables.

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: M = ModularSymbols(37)
    sage: a = 389
    sage: t = M.T(2003).matrix(); t.charpoly().factor()
     _4 = (x - 2004) * (x - 12)^2 * (x + 54)^2

Next we save our session, which saves each of the above variables
into a file. Then we view the file, which is about 3K in size.

.. skip

::

    sage: save_session('misc')
    Saving a
    Saving M
    Saving t
    Saving E
    sage: quit
    was@form:~/tmp$ ls -l misc.sobj
    -rw-r--r--  1 was was 2979 2006-01-28 19:47 misc.sobj

Finally we restart Sage, define an extra variable, and load our saved
session.

.. skip

::

    sage: b = 19
    sage: load_session('misc')
    Loading a
    Loading M
    Loading E
    Loading t

Each saved variable is again available. Moreover, the variable
``b`` was not overwritten.

.. skip

::

    sage: M
    Full Modular Symbols space for Gamma_0(37) of weight 2 with sign 0
    and dimension 5 over Rational Field
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational
    Field
    sage: b
    19
    sage: a
    389

