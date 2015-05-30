"""
Configuration for CryptoMiniSat.

EXAMPLE:

We construct a new configuration object::

    sage: from sage.sat.solvers.cryptominisat import SolverConf    # optional - cryptominisat
    sage: s = SolverConf()                                         # optional - cryptominisat
    sage: s.doxorsubsumption                                       # optional - cryptominisat
    True

and modify it such that we disable xor subsumption::

    sage: s.doxorsubsumption = False                               # optional - cryptominisat

Finally, we pass it on to CryptoMiniSat::

    sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat # optional - cryptominisat
    sage: cms = CryptoMiniSat(s)                                   # optional - cryptominisat

Note that we can achieve the same effect by::

    sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat # optional - cryptominisat
    sage: cms = CryptoMiniSat(doxorsubsumption=False)              # optional - cryptominisat

AUTHORS:

- Martin Albrecht (2012): first version
"""
##############################################################################
#  Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

###
# NOTE: To add new options edit solverconf_helper.cpp which maintains
# a map from names to C++ attributes.
###

from libc.stdint cimport uint32_t, uint64_t


cdef class SolverConf(object):
    """
    Configuration for cls:`CryptoMiniSat`

    This class implements an interface to the C++ SolverConf class
    which allows to configure the behaviour of the cls:`CryptoMiniSat`
    solver.
    """
    def __cinit__(self, **kwds):
        """
        Construct a new configuration object.

        INPUT:

        - ``kwds`` - see string representation of any instance of this
                     class for a list of options.

        EXAMPLE::

            sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
            sage: s = SolverConf()                                      # optional - cryptominisat
            sage: s                                                     # optional - cryptominisat
            random_var_freq:   0.001000 # the frequency with which the decision heuristic tries to choose a random variable.        (default 0.02)
            ...
            greedyunbound:       True # if set, then variables will be greedily unbounded (set to l_undef). this is experimental
                 origseed:          0 #
        """
        self._conf = new SolverConfC()
        self._nopts = setup_map(self._map, self._conf[0], 100)

        for k,v in kwds.iteritems():
            self[k] = v

    def __dealloc__(self):
        """
        TESTS::

        sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
        sage: for i in range(100):                                  # optional - cryptominisat
        ...      s = SolverConf()
        ...      del s

        """
        del self._conf

    def __setitem__(self, name, value):
        """
        Set options using dictionary notation.

        EXAMPLE::

            sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
            sage: s = SolverConf()                                      # optional - cryptominisat
            sage: s['polarity_mode'] = 3                                # optional - cryptominisat
            sage: s['polarity_mode']                                    # optional - cryptominisat
            3

        TESTS::

            sage: s['foobar'] = 1.2                                     # optional - cryptominisat
            Traceback (most recent call last):
            ...
            AttributeError: SolverConf has no option 'foobar'

        """
        for i in range(self._nopts):
            name_i = self._map[i].name.lower()
            if name != name_i:
                continue
            if self._map[i].type == t_int:
                (<int*>self._map[i].target)[0] = value
            elif self._map[i].type == t_float:
                (<float*>self._map[i].target)[0] = value
            elif self._map[i].type == t_double:
                (<double*>self._map[i].target)[0] = value
            elif self._map[i].type == t_Var:
                (<uint32_t*>self._map[i].target)[0] = value
            elif self._map[i].type == t_bool:
                (<bint*>self._map[i].target)[0] = value
            elif self._map[i].type == t_uint32_t:
                (<uint32_t*>self._map[i].target)[0] = value
            elif self._map[i].type == t_uint64_t:
                (<uint64_t*>self._map[i].target)[0] = value
            else:
                raise NotImplementedError("Type %d of CryptoMiniSat is not supported yet."%self._map[i].type)
            return
        raise AttributeError("SolverConf has no option '%s'"%name)

    def __setattr__(self, name, value):
        """
        Set options using attributes.

        EXAMPLE::

            sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
            sage: s = SolverConf()                                      # optional - cryptominisat
            sage: s.dovarelim = False                                   # optional - cryptominisat
            sage: s.dovarelim                                           # optional - cryptominisat
            False

        TESTS::

            sage: s.foobar = 1.2                                        # optional - cryptominisat
            Traceback (most recent call last):
            ...
            AttributeError: SolverConf has no option 'foobar'
        """
        self[name] = value

    def __getitem__(self, name):
        """
        Read options using dictionary notation.

        EXAMPLE::

            sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
            sage: s = SolverConf()                                      # optional - cryptominisat
            sage: s['simpstartmult']                                    # optional - cryptominisat
            300.0

        TESTS::

            sage: s['foobar']                                           # optional - cryptominisat
            Traceback (most recent call last):
            ...
            AttributeError: SolverConf has no option 'foobar'

        """
        for i in range(self._nopts):
            name_i = self._map[i].name.lower()
            if name_i != name:
                continue
            if self._map[i].type == t_int:
                return int((<int*>self._map[i].target)[0])
            elif self._map[i].type == t_float:
                return float((<float*>self._map[i].target)[0])
            elif self._map[i].type == t_double:
                return float((<double*>self._map[i].target)[0])
            elif self._map[i].type == t_Var:
                return int((<uint32_t*>self._map[i].target)[0])
            elif self._map[i].type == t_bool:
                return bool((<bint*>self._map[i].target)[0])
            elif self._map[i].type == t_uint32_t:
                return int((<uint32_t*>self._map[i].target)[0])
            elif self._map[i].type == t_uint64_t:
                return int((<uint64_t*>self._map[i].target)[0])
            else:
                raise NotImplementedError("Type %d of CryptoMiniSat is not supported."%self._map[i].type)
        raise AttributeError("SolverConf has no option '%s'"%name)

    def __getattr__(self, name):
        """
        Read options using dictionary notation.

        EXAMPLE::

            sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
            sage: s = SolverConf()                                      # optional - cryptominisat
            sage: s.restrictpickbranch                                  # optional - cryptominisat
            0

        TESTS::

            sage: s.foobar                                              # optional - cryptominisat
            Traceback (most recent call last):
            ...
            AttributeError: SolverConf has no option 'foobar'
       """

        return self[name]

    def __repr__(self):
        """
        Print the current configuration.

        EXAMPLE::

            sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
            sage: s = SolverConf(); s                                   # optional - cryptominisat
            random_var_freq:   0.001000 # the frequency with which the decision heuristic tries to choose a random variable.        (default 0.02)
            ...
            greedyunbound:       True # if set, then variables will be greedily unbounded (set to l_undef). this is experimental
                 origseed:          0 #
        """
        rep = []
        for i in range(self._nopts):
            name = self._map[i].name.lower()
            doc = self._map[i].doc.lower()
            if self._map[i].type == t_int:
                val = "%10d"%(<int*>self._map[i].target)[0]
            elif self._map[i].type == t_float:
                val = "%10f"%(<float*>self._map[i].target)[0]
            elif self._map[i].type == t_double:
                val = "%10f"%(<double*>self._map[i].target)[0]
            elif self._map[i].type == t_Var:
                val = "%10d"%(<uint32_t*>self._map[i].target)[0]
            elif self._map[i].type == t_bool:
                val = "%10s"%bool((<bint*>self._map[i].target)[0])
            elif self._map[i].type == t_uint32_t:
                val = "%10d"%(<uint32_t*>self._map[i].target)[0]
            elif self._map[i].type == t_uint64_t:
                val = "%10d"%(<uint64_t*>self._map[i].target)[0]
            else:
                val = "UNKNOWN"
            rep.append("%20s: %10s # %50s"%(name,val,doc))
        return "\n".join(rep)

    def trait_names(self):
        """
        Return list of all option names.

        EXAMPLE::

            sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
            sage: s = SolverConf()                                      # optional - cryptominisat
            sage: s.trait_names()                                       # optional - cryptominisat
            ['random_var_freq', 'clause_decay', 'restart_first', 'restart_inc', 'learntsize_factor', 'learntsize_inc',
             'expensive_ccmin', 'polarity_mode', 'verbosity', 'restrictpickbranch', 'simpburstsconf', 'simpstartmult',
             'simpstartmmult', 'doperformpresimp', 'failedlitmultiplier', 'dofindxors', 'dofindeqlits',
             'doregfindeqlits', 'doreplace', 'doconglxors', 'doheuleprocess', 'doschedsimp', 'dosatelite',
             'doxorsubsumption', 'dohyperbinres', 'doblockedclause', 'dovarelim', 'dosubsume1', 'doclausvivif',
             'dosortwatched', 'dominimlearntmore', 'dominimlmorerecur', 'dofailedlit', 'doremuselessbins',
             'dosubswbins', 'dosubswnonexistbins', 'doremuselesslbins', 'doprintavgbranch', 'docacheotfssr',
             'docacheotfssrset', 'doextendedscc', 'docalcreach', 'dobxor', 'dootfsubsume', 'maxconfl', 'isplain',
             'maxrestarts', 'needtodumplearnts', 'needtodumporig', 'maxdumplearntssize', 'libraryusage', 'greedyunbound', 'origseed']
        """
        return [self._map[i].name.lower() for i in range(self._nopts)]

    def __richcmp__(self, other, op):
        """
        TESTS::

            sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
            sage: s = SolverConf()                                      # optional - cryptominisat
            sage: t = copy(s)                                           # optional - cryptominisat
            sage: t == s                                                # optional - cryptominisat
            True
            sage: t.dobxor = False                                      # optional - cryptominisat
            sage: t == s                                                # optional - cryptominisat
            False
            sage: t < s                                                 # optional - cryptominisat
            Traceback (most recent call last):
            ...
            TypeError: Configurations are not ordered.
        """
        if op not in (2,3):
            raise TypeError("Configurations are not ordered.")
        res = all(self[name] == other[name] for name in self.trait_names())
        if op == 2: # ==
            return res
        if op == 3: # !=
            return not res

    def __copy__(self):
        """
        Return a copy.

        EXAMPLE::

            sage: from sage.sat.solvers.cryptominisat import SolverConf # optional - cryptominisat
            sage: s = SolverConf()                                      # optional - cryptominisat
            sage: t = copy(s)                                           # optional - cryptominisat
            sage: t.verbosity = 1                                       # optional - cryptominisat
            sage: t['verbosity']                                        # optional - cryptominisat
            1
            sage: s.verbosity                                           # optional - cryptominisat
            0
        """
        cdef SolverConf other = SolverConf.__new__(SolverConf)
        other._conf = new SolverConfC()
        other._nopts = setup_map(other._map, other._conf[0], 100)
        for name in self.trait_names():
            other[name] = self[name]
        return other
