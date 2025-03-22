.. nodoctest

.. _profiling:

Profiling in Sage
=================

This page lists several methods available in Sage to measure and analyze the
performances of a piece of code. For more general information on profiling, see
:wikipedia:`Profiling_(computer_programming)`.

.. contents:: Table of contents
   :depth: 2
   :class: this-will-duplicate-information-and-it-is-still-useful-here

How long does it take? %time and %timeit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The two IPython magics ``%time`` and ``%timeit`` measure the time it takes to
run a command::

  sage: %time p=random_prime(2**300)
  CPU times: user 152 ms, sys: 0 ns, total: 152 ms
  Wall time: 150 ms

  sage: %timeit p=random_prime(2**300)
  10 loops, best of 3: 62.2 ms per loop


Note that while ``%time`` only runs the command once, ``%timeit`` tries to
return a more meaningful value over several runs.

For more information see ``%timeit?`` or `this page
<https://ipython.org/ipython-doc/dev/interactive/magics.html#magic-timeit>`__.

Note that Sage provides a :class:`timeit
<sage.misc.sage_timeit_class.SageTimeit>` function which also runs in the Sage
notebook.


Python-level function calls: %prun
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``%prun``, you can obtain the list of all Python functions involved in a
computation, as well as the time spent on each of them::

  sage: %prun  _=random_prime(2**500)
        468 function calls in 0.439 seconds

  Ordered by: internal time

  ncalls  tottime  percall  cumtime  percall filename:lineno(function)
      32    0.438    0.014    0.438    0.014 {method 'is_prime' of 'sage.rings.integer.Integer' objects}
      32    0.001    0.000    0.439    0.014 arith.py:407(is_prime)
      32    0.000    0.000    0.001    0.000 random.py:175(randrange)
      32    0.000    0.000    0.000    0.000 random.py:244(_randbelow)
   ...

The most time-consuming functions should appear on the top. A description of the
different columns is `available here
<https://docs.python.org/3/library/profile.html#instant-user-s-manual>`_.

.. NOTE::

   You may want to sort this list differently, e.g: use ``%prun -s cumulative``
   for decreasing cumulative time.

Alternatively, you can "save" this data to a :class:`~pstats.Stats` object for
further inspection::

  sage: %prun -r random_prime(2**500)
  sage: stats_object = _
  sage: stats_object.total_calls
  2547

For more information see ``%prun?`` or `this page
<https://ipython.org/ipython-doc/dev/interactive/magics.html#magic-prun>`__.

**Visualize the statistics:** you can obtain a more graphical output
with `RunSnakeRun <https://pypi.org/project/RunSnakeRun/>`_. Passing
``-D <file>`` to ``%prun`` will dump its output, in binary, to
``<file>``::

  sage: %prun -D random_prime.prof random_prime(2**500)

Now, in a terminal, you can open that file in ``runsnake``::

  $ runsnake random_prime.prof

RunSnakeRun does not need access to the Sage library, so the simplest
way to install it is via pip or your distribution's package manager.

Python-level line-by-line profiling: %lprun
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With `line_profiler <https://pypi.org/project/line-profiler>`_ and its
``%lprun`` magic, you can find out which lines of one (or many) functions are
the most time-consuming. The syntax is the following::

  %lprun -f function1 -f function2 code_to_run

This will display the line-by-line analysis of ``function1`` and ``function2``
when ``code_to_run`` is executed::

  sage: %lprun -f random_prime random_prime(2**500)
  Line #      Hits         Time  Per Hit   % Time  Line Contents
  ==============================================================
  1193                                           def random_prime(n, proof=None, lbound=2):
  ...                                                ...
  1251                                               # since we don't want current_randstate to get
  1252                                               # pulled when you say "from sage.arith.all import *".
  1253         1           11     11.0      0.0      from sage.misc.randstate import current_randstate
  1254         1            7      7.0      0.0      from sage.structure.proof.proof import get_flag
  1255         1            6      6.0      0.0      proof = get_flag(proof, "arithmetic")
  1256         1           17     17.0      0.0      n = ZZ(n)
  ...

In order to install ``line_profiler`` you must first run the following command::

  [user@localhost ~] sage -pip install "line_profiler"

C-level function calls: %crun
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``%crun``, you can obtain the list of all C functions involved in a
computation, as well as the time spent on each of them. You will need to have
`the Google performance analysis tools <https://github.com/gperftools/gperftools>`_
installed on your system::

  sage: %crun p=random_prime(2**500)
  PROFILE: interrupts/evictions/bytes = 45/0/18344
  Total: 45 samples
         0   0.0%   0.0%       35  77.8% PyEval_EvalCode
         0   0.0%   0.0%       35  77.8% PyEval_EvalCodeEx
         0   0.0%   0.0%       35  77.8% PyEval_EvalFrameEx
         0   0.0%   0.0%       35  77.8% PyObject_Call
         0   0.0%   0.0%       35  77.8% PyRun_StringFlags
         0   0.0%   0.0%       35  77.8% __Pyx_PyObject_Call.constprop.73
  ...

For more information on ``%crun``, see :mod:`sage.misc.gperftools`.

C-level line-by-line profiling: perf (Linux only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your code is written in C or in Cython, you can find out line-by-line which
are the most costly using `perf <https://perf.wiki.kernel.org>`_
(included in the Ubuntu package ``linux-tools``).

The easiest way to use it is to run some (very long) computation in Sage, and to
type in a console

.. CODE-BLOCK:: shell-session

  [user@localhost ~] sudo perf top

Select the entry that interests you, and press ``Enter``. The ``annotate``
command will show you:

* the CPU instructions
* the source code
* the associated time

.. CODE-BLOCK:: text

        │     *         cdef unsigned long word = (<unsigned long>1) << (v & self.radix_mod_mask)
        │     *         return (self.edges[place] & word) >> (v & self.radix_mod_mask)             # <<<<<<<<<<<<<<
        │     *
        │     *     cpdef bint has_arc(self, int u, int v) except -1:
        │     */
        │      __pyx_r = (((__pyx_v_self->edges[__pyx_v_place]) & __pyx_v_word) >> (__pyx_v_v & __pyx_v_self->radix_mod_mask));
  10.88 │      movslq %esi,%rsi
   6.52 │      and    (%rdi,%rsi,8),%rax
  12.84 │      shr    %cl,%rax


.. NOTE::

  * press ``s`` to toggle source code view
  * press ``H`` to cycle through hottest instructions
  * press ``h`` for help

Alternatively, or if you have no ``sudo`` privileges, you can record the statistics
of a specific process into a file ``perf.data`` from its PID. Then, visualize
the result using ``perf report``:

.. CODE-BLOCK:: shell-session

  [user@localhost ~] perf record -p PID
  [user@localhost ~] perf report --vmlinux vmlinux
