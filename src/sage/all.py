"""
all.py -- much of sage is imported into this module, so you don't
          have to import everything individually.

TESTS:

This is to test :trac:`10570`. If the number of stackframes at startup
changes due to a patch you made, please check that this was an
intended effect of your patch.

::

    sage: import gc
    sage: import inspect
    sage: from sage import *
    sage: frames = [x for x in gc.get_objects() if inspect.isframe(x)]

We exclude the known files and check to see that there are no others::

    sage: import os
    sage: allowed = [os.path.join("lib","python","threading.py")]
    sage: allowed.append(os.path.join("lib","python","multiprocessing"))
    sage: allowed.append(os.path.join("sage","doctest"))
    sage: allowed.append(os.path.join("bin","sage-runtests"))
    sage: allowed.append(os.path.join("site-packages","IPython"))
    sage: allowed.append(os.path.join("bin","sage-ipython"))
    sage: allowed.append("<ipython console>")
    sage: allowed.append("<doctest sage.all[3]>")
    sage: allowed.append(os.path.join("sage","combinat","species","generating_series.py"))
    sage: for i in frames:
    ....:     filename, lineno, funcname, linelist, indx = inspect.getframeinfo(i)
    ....:     for nm in allowed:
    ....:         if nm in filename:
    ....:             break
    ....:     else:
    ....:         print(filename)
    ....:

Check that the Sage Notebook is not imported at startup (see
:trac:`15335`)::

    sage: sagenb
    Traceback (most recent call last):
    ...
    NameError: name 'sagenb' is not defined

Check lazy import of ``interacts``::

    sage: type(interacts)
    <type 'sage.misc.lazy_import.LazyImport'>
    sage: interacts
    <module 'sage.interacts.all' from '...'>
"""

#*****************************************************************************
#       Copyright (C) 2005-2012 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Future statements which apply to this module. We delete the
# future globals because we do not want these to appear in the sage.all
# namespace. This deleting does not affect the parsing of this module.
from __future__ import absolute_import, division, print_function
del absolute_import, division, print_function

import os
import sys
import operator
import math

from sage.env import SAGE_ROOT, SAGE_SRC, SAGE_DOC_SRC, SAGE_LOCAL, DOT_SAGE, SAGE_ENV

# Add SAGE_SRC at the end of sys.path to enable Cython tracebacks
# (which use paths relative to SAGE_SRC)
sys.path.append(SAGE_SRC)


###################################################################

# This import also setups the interrupt handler
from cysignals.signals import (AlarmInterrupt, SignalError,
        sig_on_reset as sig_on_count)

from time                import sleep

import sage.misc.lazy_import
from sage.misc.all       import *         # takes a while
from sage.typeset.all    import *
from sage.repl.all       import *

from sage.misc.sh import sh

from sage.libs.all       import *
from sage.data_structures.all import *
from sage.doctest.all    import *

from sage.structure.all  import *
from sage.rings.all      import *
from sage.arith.all      import *
from sage.matrix.all     import *

# This must come before Calculus -- it initializes the Pynac library.
import sage.symbolic.pynac

from sage.modules.all    import *
from sage.monoids.all    import *
from sage.algebras.all   import *
from sage.modular.all    import *
from sage.sat.all        import *
from sage.schemes.all    import *
from sage.graphs.all     import *
from sage.groups.all     import *
from sage.databases.all  import *
from sage.categories.all import *
from sage.sets.all       import *
from sage.probability.all import *
from sage.interfaces.all import *

from sage.symbolic.all   import *

from sage.functions.all  import *
from sage.calculus.all   import *

from sage.server.all     import *
import sage.tests.all as tests

from sage.crypto.all     import *
import sage.crypto.mq as mq

from sage.plot.all       import *
from sage.plot.plot3d.all     import *

from sage.coding.all     import *
from sage.combinat.all   import *

from sage.lfunctions.all import *

from sage.geometry.all   import *
from sage.geometry.triangulation.all   import *
from sage.geometry.riemannian_manifolds.all   import *

from sage.dynamics.all   import *

from sage.homology.all   import *

from sage.quadratic_forms.all import *

from sage.gsl.all        import *

from sage.games.all      import *

from sage.media.all      import *

from sage.logic.all      import *

from sage.numerical.all  import *

from sage.stats.all      import *
import sage.stats.all as stats

import sage.finance.all  as finance

from sage.parallel.all   import *

from sage.ext.fast_callable  import fast_callable
from sage.ext.fast_eval      import fast_float

sage.misc.lazy_import.lazy_import('sage.sandpiles.all', '*', globals())

from sage.tensor.all     import *

from sage.matroids.all   import *

from sage.game_theory.all import *

from sage.knots.all import *

from sage.manifolds.all import *

from cysignals.alarm import alarm, cancel_alarm

# Lazily import notebook functions and interacts (#15335)
lazy_import('sagenb.notebook.notebook_object', 'notebook')
lazy_import('sagenb.notebook.notebook_object', 'inotebook')
lazy_import('sagenb.notebook.sage_email', 'email')
lazy_import('sage.interacts', 'all', 'interacts')
lazy_import('sage.interacts.decorator', 'interact')
from sage.interacts.debugger import debug

from copy import copy, deepcopy

# The code executed here uses a large amount of Sage components
from sage.rings.qqbar import _init_qqbar
_init_qqbar()


###########################################################
#### WARNING:
# DO *not* import numpy / matplotlib / networkx here!!
# Each takes a surprisingly long time to initialize,
# and that initialization should be done more on-the-fly
# when they are first needed.
###########################################################

CC = ComplexField()
QQ = RationalField()
RR = RealField()  # default real field
ZZ = IntegerRing()

true = True
false = False
oo = infinity

from sage.misc.copying import license
copying = license
copyright = license

_cpu_time_ = cputime()
_wall_time_ = walltime()

def quit_sage(verbose=True):
    """
    If you use Sage in library mode, you should call this function
    when your application quits.

    It makes sure any child processes are also killed, etc.
    """
    if verbose:
        t1 = cputime(_cpu_time_)
        t1m = int(t1) // 60
        t1s = t1 - t1m * 60
        t2 = walltime(_wall_time_)
        t2m = int(t2) // 60
        t2s = t2 - t2m * 60
        print("Exiting Sage (CPU time %sm%.2fs, Wall time %sm%.2fs)." %
              (t1m, t1s, t2m, t2s))

    import gc
    gc.collect()

    from sage.interfaces.quit import expect_quitall
    expect_quitall(verbose=verbose)

    import sage.matrix.matrix_mod2_dense
    sage.matrix.matrix_mod2_dense.free_m4ri()

    import sage.libs.flint.flint
    sage.libs.flint.flint.free_flint_stack()

    # stop the twisted reactor
    try:
       from twisted.internet import reactor
       if reactor.running:
          reactor.callFromThread(reactor.stop)
    except ImportError:
       pass

    # Free globally allocated mpir integers.
    import sage.rings.integer
    sage.rings.integer.free_integer_pool()
    import sage.algebras.quatalg.quaternion_algebra_element
    sage.algebras.quatalg.quaternion_algebra_element._clear_globals()

    from sage.libs.all import symmetrica
    symmetrica.end()

# A deprecation(20442) warning will be given when this module is
# imported, in particular when these functions are used.
lazy_import("sage.ext.interactive_constructors_c", ["inject_on", "inject_off"])

sage.structure.sage_object.register_unpickle_override('sage.categories.category', 'Sets', Sets)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'HeckeModules', HeckeModules)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'Objects', Objects)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'Rings', Rings)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'Fields', Fields)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'VectorSpaces', VectorSpaces)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'Schemes_over_base', sage.categories.schemes.Schemes_over_base)
sage.structure.sage_object.register_unpickle_override('sage.categories.category_types', 'ModularAbelianVarieties', ModularAbelianVarieties)

# Cache the contents of star imports.
sage.misc.lazy_import.save_cache_file()


### Debugging for Singular, see trac #10903
# from sage.libs.singular.ring import poison_currRing
# sys.settrace(poison_currRing)


# Write a file indicating that Sage was started up successfully.
# This is called by the sage-starts script.
def _write_started_file():
    """
    Write a ``sage-started.txt`` file if it does not exist.  The
    contents of this file do not matter, only its existence.

    The non-existence of this file will be used as a trigger to run
    ``sage-starts`` during the Sage build.

    TESTS:

    Check that the file exists when Sage is running::

        sage: started_file = os.path.join(SAGE_LOCAL, 'etc', 'sage-started.txt')
        sage: os.path.isfile(started_file)
        True
    """
    started_file = os.path.join(SAGE_LOCAL, 'etc', 'sage-started.txt')

    # Current time with a resolution of 1 second
    import datetime
    t = datetime.datetime.now().replace(microsecond=0)

    O = open(started_file, 'w')
    O.write("Sage {} was started at {}\n".format(sage.version.version, t))
    O.close()


# Set a new random number seed as the very last thing
# (so that printing initial_seed() and using that seed
# in set_random_seed() will result in the same sequence you got at
# Sage startup).
set_random_seed()

import warnings
warnings.filters.remove(('ignore', None, DeprecationWarning, None, 0))
# Ignore all deprecations from IPython etc.
warnings.filterwarnings('ignore',
    module='.*(IPython|ipykernel|jupyter_client|jupyter_core|nbformat|notebook|storemagic)')
# but not those that have OUR deprecation warnings
warnings.filterwarnings('default',
    '[\s\S]*See http://trac.sagemath.org/[0-9]* for details.')

# From now on it is ok to resolve lazy imports
sage.misc.lazy_import.finish_startup()
