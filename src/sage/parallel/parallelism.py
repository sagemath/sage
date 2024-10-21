# sage_setup: distribution = sagemath-categories
r"""
Parallelization control

This module defines the singleton class :class:`Parallelism` to govern the
parallelization of computations in some specific topics. It allows the user to
set the number of processes to be used for parallelization.

Some examples of use are provided in the documentation of
:meth:`sage.tensor.modules.comp.Components.contract`.

AUTHORS:

- Marco Mancini, Eric Gourgoulhon, Michal Bejger (2015): initial version
"""

# *****************************************************************************
#       Copyright (C) 2015 Marco Mancini <marco.mancini@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.fast_methods import Singleton
from sage.parallel.ncpus import ncpus
from sage.rings.integer import Integer


class Parallelism(Singleton, SageObject):
    r"""
    Singleton class for managing the number of processes used in parallel
    computations involved in various fields.

    EXAMPLES:

    The number of processes is initialized to 1 (no parallelization) for
    each field::

        sage: Parallelism()
        Number of processes for parallelization:
         - linbox computations: 1
         - tensor computations: 1

    Using 4 processes to parallelize tensor computations::

        sage: Parallelism().set('tensor', nproc=4)
        sage: Parallelism()
        Number of processes for parallelization:
         - linbox computations: 1
         - tensor computations: 4
        sage: Parallelism().get('tensor')
        4

    Using 6 processes to parallelize all types of computations::

        sage: Parallelism().set(nproc=6)
        sage: Parallelism()
        Number of processes for parallelization:
         - linbox computations: 6
         - tensor computations: 6

    Using all the cores available on the computer to parallelize tensor
    computations::

        sage: Parallelism().set('tensor')
        sage: Parallelism()  # random (depends on the computer)
        Number of processes for parallelization:
         - linbox computations: 1
         - tensor computations: 8

    Using all the cores available on the computer to parallelize all types
    of computations::

        sage: Parallelism().set()
        sage: Parallelism()  # random (depends on the computer)
        Number of processes for parallelization:
         - linbox computations: 8
         - tensor computations: 8

    Switching off all parallelizations::

        sage: Parallelism().set(nproc=1)
    """
    def __init__(self):
        r"""
        Construct the single instance of class Parallelism (singleton model).

        TESTS::

            sage: par = Parallelism()
            sage: par
            Number of processes for parallelization:
             - linbox computations: 1
             - tensor computations: 1

        Test of the singleton character::

            sage: Parallelism() is par
            True

        The test suite is passed::

            sage: TestSuite(par).run()
        """
        self._default = ncpus()
        # default number of proc. used in parallelizations

        self._nproc = {'tensor': 1, 'linbox': 1}
        # dict. of number of processes to be used
        # (keys: computational field)

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: Parallelism()._repr_()
            'Number of processes for parallelization:\n - linbox computations: 1\n - tensor computations: 1'
        """
        resu = "Number of processes for parallelization:\n"
        for field in sorted(self._nproc):
            resu += " - {} computations: {}\n".format(field, self._nproc[field])
        return resu[:-1]

    def reset(self):
        r"""
        Put the singleton object ``Parallelism()`` in the same state as
        immediately after its creation.

        EXAMPLES:

        State of ``Parallelism()`` just after its creation::

            sage: Parallelism()
            Number of processes for parallelization:
             - linbox computations: 1
             - tensor computations: 1
            sage: Parallelism().get_default()  # random (depends on the computer)
            8

        Changing some values::

            sage: Parallelism().set_default(6)
            sage: Parallelism().set()
            sage: Parallelism()
            Number of processes for parallelization:
             - linbox computations: 6
             - tensor computations: 6
            sage: Parallelism().get_default()
            6

        Back to the initial state::

            sage: Parallelism().reset()
            sage: Parallelism()
            Number of processes for parallelization:
             - linbox computations: 1
             - tensor computations: 1
            sage: Parallelism().get_default()  # random (depends on the computer)
            8
        """
        self._default = ncpus()
        for field in self._nproc:
            self._nproc[field] = 1

    def set(self, field=None, nproc=None):
        r"""
        Set the number of processes to be launched for parallel computations
        regarding some specific field.

        INPUT:

        - ``field`` -- (default: ``None``) string specifying the computational
          field for which the number of parallel processes is to be set; if
          ``None``, all fields are considered
        - ``nproc`` -- (default: ``None``) number of processes to be used for
          parallelization; if ``None``, the number of processes will be set to
          the default value, which, unless redefined by :meth:`set_default`,
          is the total number of cores found on the computer.

        EXAMPLES:

        The default is a single processor (no parallelization)::

            sage: Parallelism()
            Number of processes for parallelization:
             - linbox computations: 1
             - tensor computations: 1

        Asking for parallelization on 4 cores in tensor algebra::

            sage: Parallelism().set('tensor', nproc=4)
            sage: Parallelism()
            Number of processes for parallelization:
             - linbox computations: 1
             - tensor computations: 4

        Using all the cores available on the computer::

            sage: Parallelism().set('tensor')
            sage: Parallelism()  # random (depends on the computer)
            Number of processes for parallelization:
             - linbox computations: 1
             - tensor computations: 8

        Using 6 cores in all parallelizations::

            sage: Parallelism().set(nproc=6)
            sage: Parallelism()
            Number of processes for parallelization:
             - linbox computations: 6
             - tensor computations: 6

        Using all the cores available on the computer in all parallelizations::

            sage: Parallelism().set()
            sage: Parallelism()  # random (depends on the computer)
            Number of processes for parallelization:
             - linbox computations: 8
             - tensor computations: 8

        Switching off the parallelization::

            sage: Parallelism().set(nproc=1)
            sage: Parallelism()
            Number of processes for parallelization:
             - linbox computations: 1
             - tensor computations: 1
        """
        if field is None:
            for fi in self._nproc:
                self.set(field=fi, nproc=nproc)
        else:
            if field not in self._nproc:
                raise KeyError("entry for field {} is not ".format(field) +
                               "implemented in Parallelism")
            if nproc is None:
                self._nproc[field] = self._default
            else:
                if not isinstance(nproc, (int, Integer)):
                    raise TypeError("nproc must be integer")
                self._nproc[field] = nproc

    def get(self, field):
        r"""
        Return the number of processes which will be used in parallel
        computations regarding some specific field.

        INPUT:

        - ``field`` -- string specifying the part of Sage involved in
          parallel computations

        OUTPUT:

        - number of processes used in parallelization of computations
          pertaining to ``field``

        EXAMPLES:

        The default is a single process (no parallelization)::

            sage: Parallelism().reset()
            sage: Parallelism().get('tensor')
            1

        Asking for parallelization on 4 cores::

            sage: Parallelism().set('tensor', nproc=4)
            sage: Parallelism().get('tensor')
            4
        """
        if field not in self._nproc:
            raise KeyError("entry for field {} is not ".format(field) +
                           "implemented in Parallelism()")
        return self._nproc[field]

    def get_all(self):
        r"""
        Return the number of processes which will be used in parallel
        computations in all fields

        OUTPUT:

        - dictionary of the number of processes, with the computational fields
          as keys

        EXAMPLES::

            sage: Parallelism().reset()
            sage: Parallelism().get_all()
            {'linbox': 1, 'tensor': 1}

        Asking for parallelization on 4 cores::

            sage: Parallelism().set(nproc=4)
            sage: Parallelism().get_all()
            {'linbox': 4, 'tensor': 4}
        """
        return self._nproc

    def set_default(self, nproc=None):
        r"""
        Set the default number of processes to be launched in parallel
        computations.

        INPUT:

        - ``nproc`` -- (default: ``None``) default number of processes;
          if ``None``, the number of processes will be set to the total number
          of cores found on the computer.

        EXAMPLES:

        A priori the default number of process for parallelization is the
        total number of cores found on the computer::

            sage: Parallelism().get_default()  # random (depends on the computer)
            8

        Changing it thanks to ``set_default``::

            sage: Parallelism().set_default(nproc=4)
            sage: Parallelism().get_default()
            4

        Setting it back to the total number of cores available on the computer::

            sage: Parallelism().set_default()
            sage: Parallelism().get_default()  # random (depends on the computer)
            8
        """
        if nproc is None:
            self._default = ncpus()
        else:
            if not isinstance(nproc, (int, Integer)):
                raise TypeError("nproc must be integer")
            self._default = nproc

    def get_default(self):
        r"""
        Return the default number of processes to be launched in parallel
        computations.

        EXAMPLES:

        A priori, the default number of process for parallelization is the
        total number of cores found on the computer::

            sage: Parallelism().reset()
            sage: Parallelism().get_default()  # random (depends on the computer)
            8

        It can be changed via :meth:`set_default`::

            sage: Parallelism().set_default(nproc=4)
            sage: Parallelism().get_default()
            4
        """
        return self._default
