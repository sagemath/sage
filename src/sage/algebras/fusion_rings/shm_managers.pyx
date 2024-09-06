r"""
Shared Memory Managers for F-Symbol Attributes

This module provides an implementation for shared dictionary like
state attributes required by the orthogonal F-matrix solver.

Currently, the attributes only work when the base field of the
:class:`FMatrix` factory is a cyclotomic field.
"""

# ****************************************************************************
#  Copyright (C) 2021 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython
cimport numpy as np
from cysignals.memory cimport sig_malloc
from multiprocessing import shared_memory
from sage.algebras.fusion_rings.poly_tup_engine cimport poly_to_tup, tup_fixes_sq, _flatten_coeffs
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.polynomial.multi_polynomial cimport MPolynomial_libsingular
from sage.rings.polynomial.polydict cimport ETuple

import numpy as np
from os import getpid

cdef class KSHandler:
    r"""
    A shared memory backed dict-like structure to manage the
    ``_ks`` attribute of an F-matrix.

    This structure implements a representation of the known squares dictionary
    using a structured NumPy array backed by a contiguous shared memory
    object.

    The structure mimics a dictionary of ``(idx, known_sq)`` pairs. Each
    integer index corresponds to a variable and each ``known_sq`` is an
    element of the F-matrix factory's base cyclotomic field.

    Each cyclotomic coefficient is stored as a list of numerators and a
    list of denominators representing the rational coefficients. The
    structured array also maintains ``known`` attribute that indicates
    whether the structure contains an entry corresponding to the given index.

    The parent process should construct this object without a
    ``name`` attribute. Children processes use the ``name`` attribute,
    accessed via ``self.shm.name`` to attach to the shared memory block.

    INPUT:

    - ``n_slots`` -- the total number of F-symbols
    - ``field`` -- F-matrix's base cyclotomic field
    - ``use_mp`` -- boolean indicating whether to construct a shared
      memory block to back ``self``
    - ``init_data`` -- dictionary or :class:`KSHandler` object containing
      known squares for initialization, e.g., from a solver checkpoint
    - ``name`` -- the name of a shared memory object (used by child processes
        for attaching)

    .. NOTE::

        To properly dispose of shared memory resources,
        ``self.shm.unlink()`` must be called before exiting.

    .. WARNING::

        This structure *cannot* modify an entry that
        has already been set.

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.shm_managers import KSHandler
        sage: # Create shared data structure
        sage: f = FusionRing("A1", 2).get_fmatrix(inject_variables=True, new=True)
        creating variables fx1..fx14
        Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
        sage: n = f._poly_ring.ngens()
        sage: f.start_worker_pool()
        sage: ks = KSHandler(n, f._field, use_mp=True)
        sage: # In the same shell or in a different shell, attach to fvars
        sage: name = ks.shm.name
        sage: ks2 = KSHandler(n, f._field, name=name, use_mp=True)
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: eqns = [fx1**2 - 4, fx3**2 + f._field.gen()**4 - 1/19*f._field.gen()**2]
        sage: ks.update([poly_to_tup(p) for p in eqns])
        sage: for idx, sq in ks.items():
        ....:     print("Index: {}, square: {}".format(idx, sq))
        ....:
        Index: 1, square: 4
        Index: 3, square: -zeta32^4 + 1/19*zeta32^2
        sage: ks.shm.unlink()
        sage: f.shutdown_worker_pool()
    """
    def __init__(self, n_slots, field, use_mp=False, init_data=None, name=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.fusion_rings.shm_managers import KSHandler
            sage: # Create shared data structure
            sage: f = FusionRing("A1", 2).get_fmatrix(inject_variables=True, new=True)
            creating variables fx1..fx14
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
            sage: n = f._poly_ring.ngens()
            sage: f.start_worker_pool()
            sage: ks = KSHandler(n, f._field, use_mp=True)
            sage: TestSuite(ks).run()
            sage: ks.shm.unlink()
            sage: f.shutdown_worker_pool()
        """
        cdef int n, d
        if init_data is None:
            init_data = {}
        self.field = field
        n = n_slots
        d = self.field.degree()
        ks_t = np.dtype([
            ('known', 'bool', (1, )),
            ('nums', 'i8', (d, )),
            ('denoms', 'u8', (d, ))
        ])
        self.obj_cache = [None]*n
        if use_mp:
            if name is None:
                self.shm = shared_memory.SharedMemory(create=True, size=n*ks_t.itemsize)
            else:
                self.shm = shared_memory.SharedMemory(name=name)
            self.ks_dat = np.ndarray((n, ), dtype=ks_t, buffer=self.shm.buf)
        else:
            self.ks_dat = np.ndarray((n, ), dtype=ks_t)
        if name is None:
            self.ks_dat['known'] = np.zeros((n, 1), dtype='bool')
            self.ks_dat['nums'] = np.zeros((n, d), dtype='i8')
            self.ks_dat['denoms'] = np.ones((n, d), dtype='u8')
        # Populate initializer data
        for idx, sq in init_data.items():
            self.setitem(idx, sq)

    @cython.nonecheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef NumberFieldElement_absolute get(self, int idx):
        r"""
        Retrieve the known square corresponding to the given index,
        if it exists.
        """
        if not self.ks_dat['known'][idx]:
            raise KeyError('index {} does not correspond to a known square'.format(idx))
        if self.obj_cache[idx] is not None:
            return self.obj_cache[idx]
        cdef int d
        cdef list rat
        cdef Py_ssize_t i
        cdef np.ndarray[np.int64_t, ndim=1] nums = self.ks_dat['nums'][idx]
        # cdef np.int64_t[::1] num_view = nums
        cdef np.ndarray[np.uint64_t, ndim=1] denoms = self.ks_dat['denoms'][idx]
        # cdef np.uint64_t[::1] denom_view = denoms
        cdef np.int64_t num
        cdef np.uint64_t denom
        cdef NumberFieldElement_absolute cyc_coeff
        cdef Rational quo
        d = self.field.degree()
        rat = list()
        for i in range(d):
            num = nums[i]
            denom = denoms[i]
            quo = Integer(num) / Integer(denom)
            rat.append(quo)
        cyc_coeff = self.field(rat)
        self.obj_cache[idx] = cyc_coeff
        return cyc_coeff

    cpdef update(self, list eqns):
        r"""
        Update ``self``'s ``shared_memory``-backed dictionary of known
        squares. Keys are variable indices and corresponding values
        are the squares.

        EXAMPLES::

            sage: f = FusionRing("B5", 1).get_fmatrix()
            sage: f._reset_solver_state()
            sage: for idx, sq in f._ks.items():
            ....:     k
            ....:
            sage: f.get_orthogonality_constraints()
            [fx0^2 - 1,
             fx1^2 - 1,
             fx2^2 - 1,
             fx3^2 - 1,
             fx4^2 - 1,
             fx5^2 - 1,
             fx6^2 - 1,
             fx7^2 - 1,
             fx8^2 - 1,
             fx9^2 - 1,
             fx10^2 + fx12^2 - 1,
             fx10*fx11 + fx12*fx13,
             fx10*fx11 + fx12*fx13,
             fx11^2 + fx13^2 - 1]
             sage: f.get_orthogonality_constraints(output=False)
             sage: f._ks.update(f.ideal_basis)
             sage: for idx, sq in f._ks.items():
             ....:     print(idx, "-->", sq)
             ....:
             0 --> 1
             1 --> 1
             2 --> 1
             3 --> 1
             4 --> 1
             5 --> 1
             6 --> 1
             7 --> 1
             8 --> 1
             9 --> 1

        .. WARNING::

            This method assumes every polynomial in ``eqns`` is *monic*.
        """
        cdef ETuple lm
        cdef list rhs
        cdef Py_ssize_t i, idx
        cdef tuple eq_tup
        for i in range(len(eqns)):
            eq_tup = eqns[i]
            if tup_fixes_sq(eq_tup):
                rhs = [-v for v in eq_tup[-1][1]]
                # eq_tup is guaranteed univariate, so we extract variable idx from lm
                lm = eq_tup[0][0]
                idx = lm._data[0]
                try:
                    self.setitem(idx, rhs)
                except OverflowError:
                    print("KS overflowed on index {} with value {}".format(idx, self.field(rhs)))

    @cython.nonecheck(False)
    @cython.wraparound(False)
    @cython.infer_types(False)
    cdef setitem(self, int idx, rhs):
        """
        Create an entry corresponding to the given index.

        The ``rhs`` parameter may be a cyclotomic coefficient or its
        list/tuple representation.
        """
        cdef Py_ssize_t i
        cdef np.ndarray[np.int64_t, ndim=1] nums = self.ks_dat['nums'][idx]
        cdef np.ndarray[np.uint64_t, ndim=1] denoms = self.ks_dat['denoms'][idx]
        cdef np.int64_t num
        cdef np.uint64_t denom
        cdef Rational quo
        self.ks_dat['known'][idx] = True
        if not isinstance(rhs, list):
            rhs = rhs._coefficients()
        for i in range(len(rhs)):
            quo = rhs[i]
            num = quo.numerator()
            denom = quo.denominator()
            if num > 2**32:
                print("WARNING: Large num encountered in KS", num)
            if denom > 2**32:
                print("WARNING: Large denom encountered in KS", denom)
            nums[i] = num
            denoms[i] = denom

    cdef bint contains(self, int idx) noexcept:
        r"""
        Determine whether ``self`` contains entry corresponding to given
        ``idx``.
        """
        return self.ks_dat[idx]['known']

    def __eq__(self, KSHandler other):
        r"""
        Test for equality.

        TESTS::

            sage: f = FusionRing("C2", 2).get_fmatrix()
            sage: f._reset_solver_state()
            sage: f.get_orthogonality_constraints(output=False)
            sage: from sage.algebras.fusion_rings.shm_managers import KSHandler
            sage: n = f._poly_ring.ngens()
            sage: f.start_worker_pool()
            sage: ks = KSHandler(n, f._field, use_mp=True, init_data=f._ks)
            sage: # In the same shell or in a different one, attach to shared memory handler
            sage: name = ks.shm.name
            sage: k2 = KSHandler(n, f._field, name=name, use_mp=True)
            sage: ks == k2
            True
            sage: ks.shm.unlink()
            sage: f.shutdown_worker_pool()
        """
        return all(other.get(idx) == sq for idx, sq in self.items())

    def __reduce__(self):
        r"""
        Provide pickling / unpickling support for ``self``.

        TESTS::

            sage: f = FusionRing("A3", 1).get_fmatrix()
            sage: f._reset_solver_state()
            sage: loads(dumps(f._ks)) == f._ks
            True
            sage: f.find_orthogonal_solution(verbose=False)    # long time
            sage: loads(dumps(f._ks)) == f._ks
            True
        """
        d = {i: sq for i, sq in self.items()}
        return make_KSHandler, (self.ks_dat.size, self.field, d)

    def items(self):
        r"""
        Iterate through existing entries using Python dict-style syntax.

        EXAMPLES::

            sage: f = FusionRing("A3", 1).get_fmatrix()
            sage: f._reset_solver_state()
            sage: f.get_orthogonality_constraints(output=False)
            sage: f._ks.update(f.ideal_basis)
            sage: for idx, sq in f._ks.items():
            ....:     print("Index: {}, sq: {}".format(idx, sq))
            ....:
            Index: 0, sq: 1
            Index: 1, sq: 1
            Index: 2, sq: 1
            Index: 3, sq: 1
            Index: 4, sq: 1
            ...
            Index: 25, sq: 1
            Index: 26, sq: 1
        """
        cdef Py_ssize_t i
        for i in range(self.ks_dat.size):
            if self.ks_dat['known'][i]:
                yield i, self.get(i)


def make_KSHandler(n_slots, field, init_data):
    r"""
    Provide pickling / unpickling support for :class:`KSHandler`.

    TESTS::

        sage: f = FusionRing("B4", 1).get_fmatrix()
        sage: f._reset_solver_state()
        sage: loads(dumps(f._ks)) == f._ks                 # indirect doctest
        True
        sage: f.find_orthogonal_solution(verbose=False)    # long time
        sage: loads(dumps(f._ks)) == f._ks                 # indirect doctest
        True
    """
    return KSHandler(n_slots, field, init_data=init_data)


cdef class FvarsHandler:
    r"""
    A shared memory backed dict-like structure to manage the
    ``_fvars`` attribute of an F-matrix.

    This structure implements a representation of the F-symbols dictionary
    using a structured NumPy array backed by a contiguous shared memory
    object.

    The monomial data is stored in the ``exp_data`` structure. Monomial
    exponent data is stored contiguously and ``ticks`` are used to
    indicate different monomials.

    Coefficient data is stored in the ``coeff_nums`` and ``coeff_denom``
    arrays. The ``coeff_denom`` array stores the value
    ``d = coeff.denominator()`` for each cyclotomic coefficient. The
    ``coeff_nums`` array stores the values
    ``c.numerator() * d for c in coeff._coefficients()``, the abridged
    list representation of the cyclotomic coefficient ``coeff``.

    Each entry also has a boolean ``modified`` attribute, indicating
    whether it has been modified by the parent process. Entry retrieval
    is cached in each process, so each process must check whether
    entries have been modified before attempting retrieval.

    The parent process should construct this object without a
    ``name`` attribute. Children processes use the ``name`` attribute,
    accessed via ``self.shm.name`` to attach to the shared memory block.

    INPUT:

    - ``n_slots`` -- number of generators of the underlying polynomial ring
    - ``field`` -- base field for polynomial ring
    - ``idx_to_sextuple`` -- map relating a single integer index to a sextuple
      of ``FusionRing`` elements
    - ``init_data`` -- dictionary or :class:`FvarsHandler` object containing
      known squares for initialization, e.g., from a solver checkpoint
    - ``use_mp`` -- integer indicating the number of child processes
      used for multiprocessing; if running serially, use 0
    - ``pids_name`` -- the name of a ``ShareableList`` containing the
      process ``pid``'s for every process in the pool (including the
      parent process)
    - ``name`` -- the name of a shared memory object
      (used by child processes for attaching)
    - ``max_terms`` -- maximum number of terms in each entry; since
      we use contiguous C-style memory blocks, the size of the block
      must be known in advance
    - ``n_bytes`` -- the number of bytes that should be allocated for
      each numerator and each denominator stored by the structure

    .. NOTE::

        To properly dispose of shared memory resources,
        ``self.shm.unlink()`` must be called before exiting.

    .. NOTE::

        If you ever encounter an :exc:`OverflowError` when running the
        :meth:`FMatrix.find_orthogonal_solution` solver, consider
        increasing the parameter ``n_bytes``.

    .. WARNING::

        The current data structure supports up to `2^16` entries,
        with each monomial in each entry having at most 254
        nonzero terms. On average, each of the ``max_terms`` monomials
        can have at most 30 terms.

    EXAMPLES::

        sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
        sage: # Create shared data structure
        sage: f = FusionRing("A2", 1).get_fmatrix(inject_variables=True, new=True)
        creating variables fx1..fx8
        Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7
        sage: f.start_worker_pool()
        sage: n_proc = f.pool._processes
        sage: pids_name = f._pid_list.shm.name
        sage: fvars = FvarsHandler(8, f._field, f._idx_to_sextuple, use_mp=n_proc, pids_name=pids_name)
        sage: # In the same shell or in a different shell, attach to fvars
        sage: name = fvars.shm.name
        sage: fvars2 = FvarsHandler(8, f._field, f._idx_to_sextuple, name=name , use_mp=n_proc, pids_name=pids_name)
        sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
        sage: rhs = tuple((exp, tuple(c._coefficients())) for exp, c in poly_to_tup(fx5**5))
        sage: fvars[f2, f1, f2, f2, f0, f0] = rhs
        sage: f._tup_to_fpoly(fvars2[f2, f1, f2, f2, f0, f0])
        fx5^5
        sage: fvars.shm.unlink()
        sage: f.shutdown_worker_pool()
    """
    def __init__(self, n_slots, field, idx_to_sextuple, init_data=None,
                 use_mp=0,
                 pids_name=None, name=None, max_terms=20, n_bytes=32):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
            sage: # Create shared data structure
            sage: f = FusionRing("A2", 1).get_fmatrix(inject_variables=True, new=True)
            creating variables fx1..fx8
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7
            sage: f.start_worker_pool()
            sage: n_proc = f.pool._processes
            sage: pids_name = f._pid_list.shm.name
            sage: fvars = FvarsHandler(8, f._field, f._idx_to_sextuple, use_mp=n_proc, pids_name=pids_name)
            sage: TestSuite(fvars).run(skip='_test_pickling')
            sage: fvars.shm.unlink()
            sage: f.shutdown_worker_pool()
        """
        self.field = field
        self.obj_cache = dict()
        cdef int d = self.field.degree()
        self.bytes = n_bytes
        cdef int slots = self.bytes // 8
        cdef int n_proc = use_mp + 1
        if init_data is None:
            init_data = {}
        self.fvars_t = np.dtype([
            ('modified', np.int8, (n_proc, )),
            ('ticks', 'u1', (max_terms, )),
            ('exp_data', 'u2', (max_terms*30, )),
            ('coeff_nums', np.int64, (max_terms, d, slots)),
            ('coeff_denom', np.uint64, (max_terms, d, slots))
        ])
        self.sext_to_idx = {s: i for i, s in idx_to_sextuple.items()}
        self.ngens = n_slots
        if use_mp:
            if name is None:
                self.shm = shared_memory.SharedMemory(create=True, size=self.ngens*self.fvars_t.itemsize)
            else:
                self.shm = shared_memory.SharedMemory(name=name)
            self.fvars = np.ndarray((self.ngens, ), dtype=self.fvars_t, buffer=self.shm.buf)
            self.pid_list = shared_memory.ShareableList(name=pids_name)
            self.child_id = -1
        else:
            self.fvars = np.ndarray((self.ngens, ), dtype=self.fvars_t)
            self.child_id = 0
        # Populate with initialziation data
        for sextuple, fvar in init_data.items():
            if isinstance(fvar, MPolynomial_libsingular):
                fvar = _flatten_coeffs(poly_to_tup(fvar))
            if isinstance(fvar, NumberFieldElement_absolute):
                fvar = ((ETuple({}, self.ngens), tuple(fvar._coefficients())), )
            if isinstance(fvar, tuple):
                transformed = list()
                for exp, c in fvar:
                    if isinstance(c, NumberFieldElement_absolute):
                        transformed.append((exp, tuple(c._coefficients())))
                if transformed:
                    fvar = tuple(transformed)
            self[sextuple] = fvar

    @cython.nonecheck(False)
    @cython.wraparound(False)
    @cython.boundscheck(False)
    def __getitem__(self, sextuple):
        r"""
        Retrieve a record from the shared memory data structure by
        unflattening its representation and constructing relevant Python
        objects.

        This method returns a tuple of ``(ETuple, coeff)`` pairs,
        where ``coeff`` is an element of ``self.field``.

        EXAMPLES::

            sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
            sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
            sage: f = FusionRing("B7", 1).get_fmatrix(inject_variables=True, new=True)
            creating variables fx1..fx14
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
            sage: f.start_worker_pool()
            sage: n_proc = f.pool._processes
            sage: pids_name = f._pid_list.shm.name
            sage: fvars = FvarsHandler(14, f._field, f._idx_to_sextuple, use_mp=n_proc, pids_name=pids_name)
            sage: rhs = tuple((exp, tuple(c._coefficients()))
            ....:             for exp, c in poly_to_tup(1/8*fx0**15 - 23/79*fx2*fx13**3 - 799/2881*fx1*fx2**5*fx10))
            sage: fvars[(f1, f2, f1, f2, f2, f2)] = rhs
            sage: rhs = tuple((exp, tuple(c._coefficients())) for exp, c in poly_to_tup(f._poly_ring.zero()))
            sage: fvars[f2, f2, f2, f2, f0, f0] = rhs
            sage: rhs = tuple((exp, tuple(c._coefficients())) for exp, c in poly_to_tup(-1/19*f._poly_ring.one()))
            sage: fvars[f2, f1, f2, f1, f2, f2] = rhs
            sage: s, t, r = (f1, f2, f1, f2, f2, f2), (f2, f2, f2, f2, f0, f0), (f2, f1, f2, f1, f2, f2)
            sage: f._tup_to_fpoly(fvars[s]) == 1/8*fx0**15 - 23/79*fx2*fx13**3 - 799/2881*fx1*fx2**5*fx10
            True
            sage: f._tup_to_fpoly(fvars[t]) == 0
            True
            sage: f._tup_to_fpoly(fvars[r]) == -1/19
            True
            sage: fvars.shm.unlink()
            sage: f.shutdown_worker_pool()

        .. NOTE::

            This method implements caching. Only the parent process is allowed
            to modify the shared fvars structure. Each process builds its own
            cache, so each process must update its cache before retrieving a
            modified entry, tagged via its ``modified`` property.
        """
        if sextuple not in self.sext_to_idx:
            raise KeyError('invalid sextuple {}'.format(sextuple))
        cdef Py_ssize_t idx = self.sext_to_idx[sextuple]
        # Each process builds its own cache, so each process must know
        # whether the entry it wants to retrieve has been modified.
        # Each process needs to know where to look, so we use an index
        # every process agrees on. The pid_list[0] belongs to the main process.
        if self.child_id < 0:
            self.child_id = self.pid_list.index(getpid())
        if idx in self.obj_cache:
            if self.fvars['modified'][idx, self.child_id]:
                del self.obj_cache[idx]
            else:
                return self.obj_cache[idx]
        cdef ETuple e, exp
        cdef int count, nnz
        cdef Integer num
        cdef list poly_tup, rats
        cdef NumberFieldElement_absolute cyc_coeff
        cdef Py_ssize_t cum, i, j, k
        cdef Rational quo
        cdef tuple ret
        # Define memory views to reduce Python overhead and ensure correct typing
        cdef np.ndarray[np.uint8_t, ndim=1] ticks = self.fvars['ticks'][idx]
        cdef np.ndarray[np.uint16_t, ndim=1] exp_data = self.fvars['exp_data'][idx]
        cdef np.ndarray[np.int64_t, ndim=3] nums = self.fvars['coeff_nums'][idx]
        cdef np.ndarray[np.uint64_t, ndim=3] denoms = self.fvars['coeff_denom'][idx]
        cdef np.ndarray[np.int8_t, ndim=1] modified = self.fvars['modified'][idx]
        e = ETuple({}, self.ngens)
        poly_tup = list()
        cum = 0
        count = np.count_nonzero(ticks)
        for i in range(count):
            # Construct new ETuple for each monomial
            exp = e._new()
            # Handle constant coeff
            nnz = ticks[i] if ticks[i] < 255 else 0
            exp._nonzero = nnz
            if nnz:
                exp._data = <int*>sig_malloc(sizeof(int)*nnz*2)
                for j in range(2*nnz):
                    exp._data[j] = <int>exp_data[cum]
                    cum += 1

            # Construct cyclotomic field coefficient
            rats = list()
            for k in range(self.field.degree()):
                num = Integer(list(nums[i, k]), 2**63)
                denom = Integer(list(denoms[i, k]), 2**64)
                quo = num / denom
                rats.append(quo)
            cyc_coeff = self.field(rats)
            poly_tup.append((exp, cyc_coeff))
        ret = tuple(poly_tup)
        # Cache object and reset modified
        self.obj_cache[idx] = ret
        modified[self.child_id] = 0
        return ret

    @cython.nonecheck(False)
    @cython.wraparound(False)
    def __setitem__(self, sextuple, fvar):
        r"""
        Given a sextuple of labels and a tuple of ``(ETuple, cyc_coeff)`` pairs,
        create or overwrite an entry in the shared data structure
        corresponding to the given sextuple.

        EXAMPLES::

            sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
            sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
            sage: f = FusionRing("A3", 1).get_fmatrix(inject_variables=True, new=True)
            creating variables fx1..fx27
            Defining fx0, ..., fx26
            sage: f.start_worker_pool()
            sage: n_proc = f.pool._processes
            sage: pids_name = f._pid_list.shm.name
            sage: fvars = FvarsHandler(27, f._field, f._idx_to_sextuple, use_mp=n_proc, pids_name=pids_name)
            sage: rhs = tuple((exp, tuple(c._coefficients()))
            ....:             for exp, c in poly_to_tup(1/8*fx0**15 - 23/79*fx2*fx21**3 - 799/2881*fx1*fx2**5*fx10))
            sage: fvars[(f3, f2, f1, f2, f1, f3)] = rhs
            sage: rhs = tuple((exp, tuple(c._coefficients())) for exp, c in poly_to_tup(f._poly_ring.zero()))
            sage: fvars[f3, f2, f3, f0, f1, f1] = rhs
            sage: rhs = tuple((exp, tuple(c._coefficients())) for exp, c in poly_to_tup(-1/19*f._poly_ring.one()))
            sage: fvars[f3, f3, f3, f1, f2, f2] = rhs
            sage: s, t, r = (f3, f2, f1, f2, f1, f3), (f3, f2, f3, f0, f1, f1), (f3, f3, f3, f1, f2, f2)
            sage: f._tup_to_fpoly(fvars[s]) == 1/8*fx0**15 - 23/79*fx2*fx21**3 - 799/2881*fx1*fx2**5*fx10
            True
            sage: f._tup_to_fpoly(fvars[t]) == 0
            True
            sage: f._tup_to_fpoly(fvars[r]) == -1/19
            True
            sage: fvars.shm.unlink()
            sage: f.shutdown_worker_pool()
        """
        cdef ETuple exp
        cdef Integer num, denom
        cdef tuple coeff_tup
        cdef Py_ssize_t cum, i, idx, j, k, t
        cdef Rational r
        idx = self.sext_to_idx[sextuple]
        # Clear entry before inserting
        self.fvars[idx] = np.zeros((1, ), dtype=self.fvars_t)
        # Define memory views to reduce Python overhead and ensure correct typing
        cdef np.ndarray[np.uint8_t, ndim=1] ticks = self.fvars['ticks'][idx]
        cdef np.ndarray[np.uint16_t, ndim=1] exp_data = self.fvars['exp_data'][idx]
        cdef np.ndarray[np.int64_t, ndim=3] nums = self.fvars['coeff_nums'][idx]
        cdef np.ndarray[np.uint64_t, ndim=3] denoms = self.fvars['coeff_denom'][idx]
        cdef np.ndarray[np.int8_t, ndim=1] modified = self.fvars['modified'][idx]
        cdef list digits
        # Initialize denominators to 1
        denoms[:, :, 0] = 1
        cum = 0
        i = 0
        for exp, coeff_tup in fvar:
            # Handle constant coefficient
            if exp._nonzero > 0:
                ticks[i] = exp._nonzero
            else:
                ticks[i] = -1
            for j in range(2*exp._nonzero):
                exp_data[cum] = exp._data[j]
                cum += 1
            k = 0
            for r in coeff_tup:
                num, denom = r.as_integer_ratio()
                if abs(num) > 2**63 or denom > 2**63:
                    print("Large integers encountered in FvarsHandler", num, denom)
                if abs(num) < 2**63:
                    nums[i, k, 0] = num
                else:
                    digits = num.digits(2**63)
                    # assert len(digits) <= self.bytes // 8, "Numerator {} is too large for shared FvarsHandler. Use at least {} bytes...".format(num, num.nbits()//8+1)
                    for t in range(len(digits)):
                        nums[i, k, t] = <np.int64_t>digits[t]
                if denom < 2**64:
                    denoms[i, k, 0] = denom
                else:
                    digits = denom.digits(2**64)
                    # assert len(digits) <= self.bytes // 8, "Denominator {} is too large for shared FvarsHandler. Use at least {} bytes...".format(denom, denom.nbits()//8+1)
                    for t in range(len(digits)):
                        denoms[i, k, t] = <np.uint64_t>digits[t]
                k += 1
            i += 1
        modified[:] = 1

    def __reduce__(self):
        r"""
        Provide pickling / unpickling support for ``self``.

        TESTS::

            sage: f = FusionRing("F4", 1).get_fmatrix()
            sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
            sage: n = f._poly_ring.ngens()
            sage: f.start_worker_pool()
            sage: n_proc = f.pool._processes
            sage: pids_name = f._pid_list.shm.name
            sage: fvars = FvarsHandler(n, f._field, f._idx_to_sextuple, init_data=f._fvars, use_mp=n_proc, pids_name=pids_name)
            sage: for s, fvar in loads(dumps(fvars)).items():
            ....:     assert f._fvars[s] == f._tup_to_fpoly(fvar)
            ....:
            sage: fvars.shm.unlink()
            sage: f.shutdown_worker_pool()
        """
        n = self.fvars.size
        idx_map = {i: s for s, i in self.sext_to_idx.items()}
        d = {s: fvar for s, fvar in self.items()}
        return make_FvarsHandler, (n, self.field, idx_map, d)

    def items(self):
        r"""
        Iterate through key-value pairs in the data structure as if it
        were a Python dict.

        As in a Python dict, the key-value pairs are yielded in no particular
        order.

        EXAMPLES::

            sage: f = FusionRing("G2", 1).get_fmatrix(inject_variables=True, new=True)
            creating variables fx1..fx5
            Defining fx0, fx1, fx2, fx3, fx4
            sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
            sage: shared_fvars = FvarsHandler(5, f._field, f._idx_to_sextuple, init_data=f._fvars)
            sage: for sextuple, fvar in shared_fvars.items():
            ....:     if sextuple == (f1, f1, f1, f1, f1, f1):
            ....:         f._tup_to_fpoly(fvar)
            ....:
            fx4
        """
        for sextuple in self.sext_to_idx:
            yield sextuple, self[sextuple]


def make_FvarsHandler(n, field, idx_map, init_data):
    r"""
    Provide pickling / unpickling support for :class:`FvarsHandler`.

    TESTS::

        sage: f = FusionRing("G2", 1).get_fmatrix()
        sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
        sage: n = f._poly_ring.ngens()
        sage: f.start_worker_pool()
        sage: n_proc = f.pool._processes
        sage: pids_name = f._pid_list.shm.name
        sage: fvars = FvarsHandler(n, f._field, f._idx_to_sextuple, init_data=f._fvars, use_mp=n_proc, pids_name=pids_name)
        sage: for s, fvar in loads(dumps(fvars)).items():            # indirect doctest
        ....:     assert f._fvars[s] == f._tup_to_fpoly(fvar)
        ....:
        sage: fvars.shm.unlink()
        sage: f.shutdown_worker_pool()
    """
    return FvarsHandler(n, field, idx_map, init_data=init_data)
