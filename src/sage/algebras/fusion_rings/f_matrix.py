r"""
The F-Matrix of a Fusion Ring
"""
# ****************************************************************************
#  Copyright (C) 2019 Daniel Bump <bump at match.stanford.edu>
#                     Guillermo Aboumrad <gh_willieab>
#                     Travis Scrimshaw <tcscrims at gmail.com>
#                     Galit Anikeeva <physicstravels@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import deepcopy
from ctypes import cast, py_object
from itertools import product, zip_longest
from multiprocessing import Pool, cpu_count, set_start_method, shared_memory
import numpy as np
from os import getpid, remove
import pickle

from sage.algebras.fusion_rings.fast_parallel_fmats_methods import (
    _backward_subs, _solve_for_linear_terms,
    executor
)
from sage.algebras.fusion_rings.poly_tup_engine import (
    apply_coeff_map, constant_coeff,
    compute_known_powers,
    get_variables_degrees, variables,
    poly_to_tup, _tup_to_poly, tup_to_univ_poly,
    _unflatten_coeffs,
    poly_tup_sortkey,
    resize
)
from sage.algebras.fusion_rings.shm_managers import KSHandler, FvarsHandler
from sage.graphs.graph import Graph
from sage.matrix.constructor import matrix
from sage.misc.misc import get_main_globals
from sage.rings.ideal import Ideal
from sage.structure.sage_object import SageObject
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polydict import ETuple
from sage.rings.qqbar import AA, QQbar, number_field_elements_from_algebraics


class FMatrix(SageObject):
    r"""
    An F-matrix for a :class:`FusionRing`.

    INPUT:

    - ``FR`` -- a :class:`FusionRing`
    - ``fusion_label`` -- (optional) a string used to label basis elements
      of the :class:`FusionRing` associated to ``self``
      (see :meth:`FusionRing.fusion_labels`)
    - ``var_prefix`` -- (optional) a string indicating the desired prefix
      for variables denoting F-symbols to be solved
    - ``inject_variables`` -- boolean (default: ``False``); whether to inject
      variables (:class:`FusionRing` basis element labels and F-symbols) into
      the global namespace

    The :class:`FusionRing` or Verlinde algebra is the
    Grothendieck ring of a modular tensor category [BaKi2001]_.
    Such categories arise in conformal field theory or in the
    representation theories of affine Lie algebras, or
    quantum groups at roots of unity. They have applications
    to low dimensional topology and knot theory, to conformal
    field theory and to topological quantum computing. The
    :class:`FusionRing` captures much information about a fusion
    category, but to complete the picture, the F-matrices or
    6j-symbols are needed. For example these are required in
    order to construct braid group representations. This
    can be done using the :class:`FusionRing` method
    :meth:`FusionRing.get_braid_generators`, which uses
    the F-matrix.

    We only undertake to compute the F-matrix if the
    :class:`FusionRing` is *multiplicity free* meaning that
    the Fusion coefficients `N^{ij}_k` are bounded
    by 1. For Cartan Types `X_r` and level `k`,
    the multiplicity-free cases are given by the
    following table.

    +------------------------+----------+
    | Cartan Type            | `k`      |
    +========================+==========+
    | `A_1`                  | any      |
    +------------------------+----------+
    | `A_r, r\geq 2`         | `\leq 2` |
    +------------------------+----------+
    | `B_r, r\geq 2`         | `\leq 2` |
    +------------------------+----------+
    | `C_2`                  | `\leq 2` |
    +------------------------+----------+
    | `C_r, r\geq 3`         | `\leq 1` |
    +------------------------+----------+
    | `D_r, r\geq 4`         | `\leq 2` |
    +------------------------+----------+
    | `G_2, F_4, E_6, E_7`   | `\leq 2` |
    +------------------------+----------+
    | `E_8`                  | `\leq 3` |
    +------------------------+----------+

    Beyond this limitation, computation of the F-matrix
    can involve very large systems of equations. A
    rule of thumb is that this code can compute the
    F-matrix for systems with `\leq 14` simple objects
    (primary fields) on a machine with 16 GB of memory.
    (Larger examples can be quite time consuming.)

    The :class:`FusionRing` and its methods capture much
    of the structure of the underlying tensor category.
    But an important aspect that is not encoded in the
    fusion ring is the associator, which is a homomorphism
    `(A\otimes B)\otimes C\to A\otimes(B\otimes C)` that
    requires an additional tool, the F-matrix or 6j-symbol.
    To specify this, we fix a simple object `D`
    and represent the transformation

    .. MATH::

        \text{Hom}(D, (A\otimes B)\otimes C)
        \to \text{Hom}(D, A\otimes(B\otimes C))

    by a matrix `F^{ABC}_D`. This depends on a pair of
    additional simple objects `X` and `Y`. Indeed, we can
    get a basis for `\text{Hom}(D, (A\otimes B)\otimes C)`
    indexed by simple objects `X` in which the corresponding
    homomorphism factors through `X\otimes C`, and similarly
    `\text{Hom}(D, A\otimes(B\otimes C))` has a basis indexed
    by `Y`, in which the basis vector factors through `A\otimes Y`.

    See [TTWL2009]_ for an introduction to this topic,
    [EGNO2015]_ Section 4.9 for a precise mathematical
    definition, and [Bond2007]_ Section 2.5 and [Ab2022]_ for discussions
    of how to compute the F-matrix. In addition to
    [Bond2007]_, worked out F-matrices may be found in
    [RoStWa2009]_ and [CHW2015]_.

    The F-matrix is only determined up to a *gauge*. This
    is a family of embeddings `C \to A\otimes B` for
    simple objects `A, B, C` such that `\text{Hom}(C, A\otimes B)`
    is nonzero. Changing the gauge changes the F-matrix though
    not in a very essential way. By varying the gauge it is
    possible to make the F-matrices unitary, or it is possible
    to make them cyclotomic.

    Due to the large number of equations we may fail to find a
    Groebner basis if there are too many variables.

    EXAMPLES::

        sage: I = FusionRing("E8", 2, conjugate=True)
        sage: I.fusion_labels(["i0", "p", "s"], inject_variables=True)
        sage: f = I.get_fmatrix(inject_variables=True); f
        creating variables fx1..fx14
        Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
        F-Matrix factory for The Fusion Ring of Type E8 and level 2 with Integer Ring coefficients

    We have injected two sets of variables to the global namespace.
    We created three variables ``i0, p, s`` to represent the
    primary fields (simple elements) of the :class:`FusionRing`. Creating
    the :class:`FMatrix` factory also created variables
    ``fx1, fx2, ..., fx14`` in order to solve the hexagon and pentagon
    equations describing the F-matrix. Since we called :class:`FMatrix`
    with the parameter ``inject_variables=True``, these have been injected
    into the global namespace. This is not necessary for the code to work
    but if you want to run the code experimentally you may want access
    to these variables.

    EXAMPLES::

        sage: f.fmatrix(s, s, s, s)
        [fx10 fx11]
        [fx12 fx13]

    The F-matrix has not been computed at this stage, so
    the F-matrix `F^{sss}_s` is filled with variables
    ``fx10``, ``fx11``, ``fx12``, ``fx13``. The task is
    to solve for these.

    As explained above The F-matrix `(F^{ABC}_D)_{X, Y}`
    two other variables `X` and `Y`. We have methods to
    tell us (depending on `A, B, C, D`) what the possibilities
    for these are. In this example with `A=B=C=D=s`
    both `X` and `Y` are allowed to be `i_0` or `s`.

    ::

        sage: f.f_from(s, s, s, s), f.f_to(s, s, s, s)
        ([i0, p], [i0, p])

    The last two statements show that the possible values of
    `X` and `Y` when `A = B = C = D = s` are `i_0` and `p`.

    The F-matrix is computed by solving the so-called
    pentagon and hexagon equations. The *pentagon equations*
    reflect the Mac Lane pentagon axiom in the definition
    of a monoidal category. The hexagon relations
    reflect the axioms of a *braided monoidal category*,
    which are constraints on both the F-matrix and on
    the R-matrix. Optionally, orthogonality constraints
    may be imposed to obtain an orthogonal F-matrix.

    ::

        sage: sorted(f.get_defining_equations("pentagons"))[1:3]
        [fx9*fx12 - fx2*fx13, fx4*fx11 - fx2*fx13]
        sage: sorted(f.get_defining_equations("hexagons"))[1:3]
        [fx6 - 1, fx2 + 1]
        sage: sorted(f.get_orthogonality_constraints())[1:3]
        [fx10*fx11 + fx12*fx13, fx10*fx11 + fx12*fx13]

    There are two methods available to compute an F-matrix.
    The first, :meth:`find_cyclotomic_solution` uses only
    the pentagon and hexagon relations. The second,
    :meth:`find_orthogonal_solution` uses additionally
    the orthogonality relations. There are some differences
    that should be kept in mind.

    :meth:`find_cyclotomic_solution` currently works only with
    smaller examples. For example the :class:`FusionRing` for `G_2`
    at level 2 is too large. When it is available, this method
    produces an F-matrix whose entries are in the same
    cyclotomic field as the underlying :class:`FusionRing`. ::

        sage: f.find_cyclotomic_solution()
        Setting up hexagons and pentagons...
        Finding a Groebner basis...
        Solving...
        Fixing the gauge...
        adding equation... fx1 - 1
        adding equation... fx11 - 1
        Done!

    We now have access to the values of the F-matrix using
    the methods :meth:`fmatrix` and :meth:`fmat`::

        sage: f.fmatrix(s, s, s, s)
        [(-1/2*zeta128^48 + 1/2*zeta128^16)                                  1]
        [                               1/2  (1/2*zeta128^48 - 1/2*zeta128^16)]
        sage: f.fmat(s, s, s, s, p, p)
        (1/2*zeta128^48 - 1/2*zeta128^16)

    :meth:`find_orthogonal_solution` is much more powerful
    and is capable of handling large cases, sometimes
    quickly but sometimes (in larger cases) after hours of
    computation. Its F-matrices are not always in the
    cyclotomic field that is the base ring of the underlying
    :class:`FusionRing`, but sometimes in an extension field adjoining
    some square roots. When this happens, the :class:`FusionRing` is
    modified, adding an attribute ``_basecoer`` that is
    a coercion from the cyclotomic field to the field
    containing the F-matrix. The field containing the F-matrix
    is available through :meth:`field`. ::

        sage: f = FusionRing("B3", 2).get_fmatrix()
        sage: f.find_orthogonal_solution(verbose=False, checkpoint=True)     # not tested (~100 s)
        sage: all(v in CyclotomicField(56) for v in f.get_fvars().values())  # not tested
        True

        sage: f = FusionRing("G2", 2).get_fmatrix()
        sage: f.find_orthogonal_solution(verbose=False) # long time (~11 s)
        sage: f.field()                                 # long time
        Algebraic Field
    """
    def __init__(self, fusion_ring, fusion_label='f', var_prefix='fx', inject_variables=False):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: f = FusionRing("B3", 2).get_fmatrix()
            sage: TestSuite(f).run(skip='_test_pickling')
        """
        self._FR = fusion_ring
        if inject_variables and (self._FR._fusion_labels is None):
            self._FR.fusion_labels(fusion_label, inject_variables=True)
        if not self._FR.is_multiplicity_free():
            raise NotImplementedError("FMatrix is only available for multiplicity free FusionRings")
        # Set up F-symbols entry by entry
        n_vars = self.findcases()
        self._poly_ring = PolynomialRing(self._FR.field(), n_vars, var_prefix)
        if inject_variables:
            print("creating variables %s%s..%s%s" % (var_prefix, 1, var_prefix, n_vars))
            self._poly_ring.inject_variables(get_main_globals())
        self._idx_to_sextuple, self._fvars = self.findcases(output=True)

        # Base field attributes
        self._field = self._FR.field()
        r = self._field.defining_polynomial().roots(ring=QQbar, multiplicities=False)[0]
        self._qqbar_embedding = self._field.hom([r], QQbar)

        # Warm starting
        self._chkpt_status = -1

        # Multiprocessing attributes
        self.mp_thresh = 10000
        self.pool = None

    #######################
    #   Class utilities   #
    #######################

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: FusionRing("B2", 1).get_fmatrix()
            F-Matrix factory for The Fusion Ring of Type B2 and level 1 with Integer Ring coefficients
        """
        return "F-Matrix factory for %s" % self._FR

    def clear_equations(self):
        r"""
        Clear the list of equations to be solved.

        EXAMPLES::

            sage: f = FusionRing("E6", 1).get_fmatrix()
            sage: f.get_defining_equations('hexagons', output=False)
            sage: len(f.ideal_basis)
            6
            sage: f.clear_equations()
            sage: len(f.ideal_basis) == 0
            True
        """
        self.ideal_basis = []

    def clear_vars(self):
        r"""
        Reset the F-symbols.

        EXAMPLES::

            sage: f = FusionRing("C4", 1).get_fmatrix()
            sage: fvars = f.get_fvars()
            sage: some_key = sorted(fvars)[0]
            sage: fvars[some_key]
            fx0
            sage: fvars[some_key] = 1
            sage: f.get_fvars()[some_key]
            1
            sage: f.clear_vars()
            sage: f.get_fvars()[some_key]
            fx0
        """
        self._fvars = {t: self._poly_ring.gen(idx) for idx, t in self._idx_to_sextuple.items()}
        self._solved = [False] * self._poly_ring.ngens()

    def _reset_solver_state(self):
        r"""
        Reset solver state and clear relevant cache.

        Used to ensure state variables are the same for each
        orthogonal solver run.

        EXAMPLES::

            sage: f = FusionRing("G2", 1).get_fmatrix()
            sage: f._reset_solver_state()
            sage: K = f.field()
            sage: len(f._nnz.nonzero_positions())
            1
            sage: f.find_orthogonal_solution(verbose=False)
            sage: K == f.field()
            False
            sage: f._reset_solver_state()
            sage: K == f.field()
            True
            sage: f.FR()._basecoer is None
            True
            sage: f._poly_ring.base_ring() == K
            True
            sage: sum(f._solved) == 0
            True
            sage: len(f.ideal_basis) == 0
            True
            sage: for k, v in f._ks.items():
            ....:     k
            sage: len(f._nnz.nonzero_positions()) == 1
            True
            sage: all(len(x.q_dimension.cache) == 0 for x in f.FR().basis())
            True
            sage: len(f.FR().r_matrix.cache) == 0
            True
            sage: len(f.FR().s_ij.cache) == 0
            True
        """
        self._FR._basecoer = None
        self._field = self._FR.field()
        self._non_cyc_roots = []
        self._poly_ring = self._poly_ring.change_ring(self._field)
        self._chkpt_status = -1
        self.clear_vars()
        self.clear_equations()
        n = self._poly_ring.ngens()
        self._var_degs = [0] * n
        self._kp = {}
        self._ks = KSHandler(n, self._field)
        self._singles = self.get_fvars_by_size(1, indices=True)
        self._nnz = self._get_known_nonz()

        # Clear relevant caches
        [x.q_dimension.clear_cache() for x in self._FR.basis()]
        self._FR.r_matrix.clear_cache()
        self._FR.s_ij.clear_cache()

    def fmat(self, a, b, c, d, x, y, data=True):
        r"""
        Return the F-Matrix coefficient `(F^{a, b, c}_d)_{x, y}`.

        EXAMPLES::

            sage: fr = FusionRing("G2", 1, fusion_labels=("i0", "t"), inject_variables=True)
            sage: f = fr.get_fmatrix()
            sage: [f.fmat(t, t, t, t, x, y) for x in fr.basis() for y in fr.basis()]
            [fx1, fx2, fx3, fx4]
            sage: f.find_cyclotomic_solution(output=True)
            Setting up hexagons and pentagons...
            Finding a Groebner basis...
            Solving...
            Fixing the gauge...
            adding equation... fx2 - 1
            Done!
            {(t, t, t, i0, t, t): 1,
             (t, t, t, t, i0, i0): (-zeta60^14 + zeta60^6 + zeta60^4 - 1),
             (t, t, t, t, i0, t): 1,
             (t, t, t, t, t, i0): (-zeta60^14 + zeta60^6 + zeta60^4 - 1),
             (t, t, t, t, t, t): (zeta60^14 - zeta60^6 - zeta60^4 + 1)}
            sage: [f.fmat(t, t, t, t, x, y) for x in f._FR.basis() for y in f._FR.basis()]
            [(-zeta60^14 + zeta60^6 + zeta60^4 - 1),
             1,
             (-zeta60^14 + zeta60^6 + zeta60^4 - 1),
             (zeta60^14 - zeta60^6 - zeta60^4 + 1)]
        """
        if (self._FR.Nk_ij(a, b, x) == 0
                or self._FR.Nk_ij(x, c, d) == 0
                or self._FR.Nk_ij(b, c, y) == 0
                or self._FR.Nk_ij(a, y, d) == 0):
            return 0

        # Some known zero F-symbols
        if a == self._FR.one():
            return 1 if x == b and y == d else 0

        if b == self._FR.one():
            return 1 if x == a and y == c else 0

        if c == self._FR.one():
            return 1 if x == d and y == b else 0

        if data:
            # Better to use try/except for speed. Somewhat trivial, but worth
            # hours when method is called ~10^11 times
            try:
                return self._fvars[a, b, c, d, x, y]
            except KeyError:
                return 0

        return (a, b, c, d, x, y)

    def fmatrix(self, a, b, c, d):
        r"""
        Return the F-Matrix `F^{a, b, c}_d`.

        INPUT:

        - ``a``, ``b``, ``c``, ``d`` -- basis elements of the associated
          :class:`FusionRing`

        EXAMPLES::

            sage: fr = FusionRing("A1", 2, fusion_labels='c', inject_variables=True)
            sage: f = fr.get_fmatrix(new=True)
            sage: f.fmatrix(c1, c1, c1, c1)
            [fx0 fx1]
            [fx2 fx3]
            sage: f.find_cyclotomic_solution(verbose=False);
            adding equation... fx4 - 1
            adding equation... fx10 - 1
            sage: f.f_from(c1, c1, c1, c1)
            [c0, c2]
            sage: f.f_to(c1, c1, c1, c1)
            [c0, c2]
            sage: f.fmatrix(c1, c1, c1, c1)
            [ (1/2*zeta32^12 - 1/2*zeta32^4) (-1/2*zeta32^12 + 1/2*zeta32^4)]
            [ (1/2*zeta32^12 - 1/2*zeta32^4)  (1/2*zeta32^12 - 1/2*zeta32^4)]
        """
        X = self.f_from(a, b, c, d)
        Y = self.f_to(a, b, c, d)
        return matrix([[self.fmat(a, b, c, d, x, y) for y in Y] for x in X])

    def field(self):
        r"""
        Return the base field containing the F-symbols.

        When ``self`` is initialized, the field is set to be the
        cyclotomic field of the :class:`FusionRing` associated
        to ``self``.

        The field may change after running :meth:`find_orthogonal_solution`.
        At that point, this method could return the
        associated :class:`FusionRing`'s cyclotomic field, an
        appropriate :func:`NumberField` that was computed on the fly
        by the F-matrix solver, or the :class:`QQbar<AlgebraicField>`.

        Depending on the ``CartanType`` of ``self``, the solver may need
        to compute an extension field containing certain square roots that
        do not belong to the associated :class:`FusionRing`'s cyclotomic field.

        In certain cases we revert to :class:`QQbar<AlgebraicField>` because
        the extension field computation does not seem to terminate. See
        :meth:`attempt_number_field_computation` for more details.

        The method :meth:`get_non_cyclotomic_roots` returns a list of
        roots defining the extension of the :class:`FusionRing`'s
        cyclotomic field needed to contain all F-symbols.

        EXAMPLES::

            sage: f = FusionRing("G2", 1).get_fmatrix()
            sage: f.field()
            Cyclotomic Field of order 60 and degree 16
            sage: f.find_orthogonal_solution(verbose=False)
            sage: f.field()
            Number Field in a with defining polynomial y^32 - ... - 22*y^2 + 1
            sage: phi = f.get_qqbar_embedding()
            sage: [phi(r).n() for r in f.get_non_cyclotomic_roots()]
            [-0.786151377757423 - 8.92806368517581e-31*I]

        .. NOTE::

            Consider using ``self.field().optimized_representation()`` to
            obtain an equivalent :func:`NumberField` with a defining
            polynomial with smaller coefficients, for a more efficient
            element representation.
        """
        return self._field

    def FR(self):
        r"""
        Return the :class:`FusionRing` associated to ``self``.

        EXAMPLES::

            sage: f = FusionRing("D3", 1).get_fmatrix()
            sage: f.FR()
            The Fusion Ring of Type D3 and level 1 with Integer Ring coefficients
        """
        return self._FR

    def findcases(self, output=False):
        r"""
        Return unknown F-matrix entries.

        If run with ``output=True``,
        this returns two dictionaries; otherwise it just returns the
        number of unknown values.

        EXAMPLES::

            sage: f = FusionRing("G2", 1, fusion_labels=("i0", "t")).get_fmatrix()
            sage: f.findcases()
            5
            sage: f.findcases(output=True)
             ({0: (t, t, t, i0, t, t),
              1: (t, t, t, t, i0, i0),
              2: (t, t, t, t, i0, t),
              3: (t, t, t, t, t, i0),
              4: (t, t, t, t, t, t)},
             {(t, t, t, i0, t, t): fx0,
              (t, t, t, t, i0, i0): fx1,
              (t, t, t, t, i0, t): fx2,
              (t, t, t, t, t, i0): fx3,
              (t, t, t, t, t, t): fx4})
        """
        i = 0
        if output:
            idx_map = {}
            ret = {}
        id_anyon = self._FR.one()
        for a, b, c, d in product(self._FR.basis(), repeat=4):
            if a == id_anyon or b == id_anyon or c == id_anyon:
                continue
            for x in self.f_from(a, b, c, d):
                for y in self.f_to(a, b, c, d):
                    if output:
                        v = self._poly_ring.gen(i)
                        ret[(a, b, c, d, x, y)] = v
                        idx_map[i] = (a, b, c, d, x, y)
                    i += 1

        return (idx_map, ret) if output else i

    def f_from(self, a, b, c, d):
        r"""
        Return the possible `x` such that there are morphisms
        `d \to x \otimes c \to (a \otimes b) \otimes c`.

        INPUT:

        - ``a``, ``b``, ``c``, ``d`` -- basis elements of the associated
          :class:`FusionRing`

        EXAMPLES::

            sage: fr = FusionRing("A1", 3, fusion_labels='a', inject_variables=True)
            sage: f = fr.get_fmatrix()
            sage: f.fmatrix(a1, a1, a2, a2)
            [fx6 fx7]
            [fx8 fx9]
            sage: f.f_from(a1, a1, a2, a2)
            [a0, a2]
            sage: f.f_to(a1, a1, a2, a2)
            [a1, a3]
        """
        return [x for x in self._FR.basis()
                if self._FR.Nk_ij(a, b, x) != 0 and self._FR.Nk_ij(x, c, d) != 0]

    def f_to(self, a, b, c, d):
        r"""
        Return the possible `y` such that there are morphisms
        `d \to a \otimes y \to a \otimes (b \otimes c)`.

        INPUT:

        - ``a``, ``b``, ``c``, ``d`` -- basis elements of the associated
          :class:`FusionRing`

        EXAMPLES::

            sage: b22 = FusionRing("B2", 2)
            sage: b22.fusion_labels("b", inject_variables=True)
            sage: B = b22.get_fmatrix()
            sage: B.fmatrix(b2, b4, b2, b4)
            [fx266 fx267 fx268]
            [fx269 fx270 fx271]
            [fx272 fx273 fx274]
            sage: B.f_from(b2, b4, b2, b4)
            [b1, b3, b5]
            sage: B.f_to(b2, b4, b2, b4)
            [b1, b3, b5]
        """
        return [y for y in self._FR.basis()
                if self._FR.Nk_ij(b, c, y) != 0 and self._FR.Nk_ij(a, y, d) != 0]

    ####################
    #   Data getters   #
    ####################

    def get_fvars(self):
        r"""
        Return a dictionary of F-symbols.

        The keys are sextuples `(a, b, c, d, x, y)` of basis elements of
        ``self.FR()`` and the values are the corresponding F-symbols
        `(F^{a, b, c}_d)_{xy}`.

        These values reflect the current state of a solver's computation.

        EXAMPLES::

            sage: f = FusionRing("A2", 1).get_fmatrix(inject_variables=True)
            creating variables fx1..fx8
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7
            sage: f.get_fvars()[(f1, f1, f1, f0, f2, f2)]
            fx0
            sage: f.find_orthogonal_solution(verbose=False)
            sage: f.get_fvars()[(f1, f1, f1, f0, f2, f2)]
            1
        """
        return self._fvars

    def get_poly_ring(self):
        r"""
        Return the polynomial ring whose generators denote the desired F-symbols.

        EXAMPLES::

            sage: f = FusionRing("B6", 1).get_fmatrix()
            sage: f.get_poly_ring()
            Multivariate Polynomial Ring in fx0, ..., fx13 over
             Cyclotomic Field of order 96 and degree 32
        """
        return self._poly_ring

    # TODO: this method is incredibly slow... improve by keeping track of the cyclotomic polynomials, NOT their roots in QQbar
    def get_non_cyclotomic_roots(self):
        r"""
        Return a list of roots that define the extension of the associated
        :class:`FusionRing`'s base
        :class:`Cyclotomic field<sage.rings.number_field.number_field.CyclotomicFieldFactory>`,
        containing all the F-symbols.

        OUTPUT:

        The list of non-cyclotomic roots is given as a list of elements of the
        field returned by :meth:`field()`.

        If ``self.field() == self.FR().field()`` then this method
        returns an empty list.

        EXAMPLES::

            sage: f = FusionRing("E6", 1).get_fmatrix()
            sage: f.find_orthogonal_solution(verbose=False)
            sage: f.field() == f.FR().field()
            True
            sage: f.get_non_cyclotomic_roots()
            []
            sage: f = FusionRing("G2", 1).get_fmatrix()
            sage: f.find_orthogonal_solution(verbose=False)
            sage: f.field() == f.FR().field()
            False
            sage: phi = f.get_qqbar_embedding()
            sage: [phi(r).n() for r in f.get_non_cyclotomic_roots()]
            [-0.786151377757423 - 8.92806368517581e-31*I]

        When ``self.field()`` is a ``NumberField``, one may use
        :meth:`get_qqbar_embedding` to embed the resulting values into
        :class:`QQbar<AlgebraicField>`.
        """
        return sorted(set(self._non_cyc_roots))

    def get_qqbar_embedding(self):
        r"""
        Return an embedding from the base field containing F-symbols (the
        associated :class:`FusionRing`'s
        :class:`Cyclotomic field<sage.rings.number_field.number_field.CyclotomicFieldFactory>`,
        a :func:`NumberField`, or :class:`QQbar<AlgebraicField>`) into
        :class:`QQbar<AlgebraicField>`.

        This embedding is useful for getting a better sense for the
        F-symbols, particularly when they are computed as elements of a
        :func:`NumberField`. See also :meth:`get_non_cyclotomic_roots`.

        EXAMPLES::

            sage: fr = FusionRing("G2", 1)
            sage: f = fr.get_fmatrix(fusion_label='g', inject_variables=True, new=True)
            creating variables fx1..fx5
            Defining fx0, fx1, fx2, fx3, fx4
            sage: f.find_orthogonal_solution()
            Computing F-symbols for The Fusion Ring of Type G2 and level 1 with Integer Ring coefficients with 5 variables...
            Set up 10 hex and orthogonality constraints...
            Partitioned 10 equations into 2 components of size:
            [4, 1]
            Elimination epoch completed... 0 eqns remain in ideal basis
            Hex elim step solved for 4 / 5 variables
            Set up 0 reduced pentagons...
            Pent elim step solved for 4 / 5 variables
            Partitioned 0 equations into 0 components of size:
            []
            Partitioned 1 equations into 1 components of size:
            [1]
            Computing appropriate NumberField...
            sage: phi = f.get_qqbar_embedding()
            sage: phi(f.fmat(g1, g1, g1, g1, g1, g1)).n()
            -0.618033988749895 + 1.46674215951686e-29*I
        """
        return self._qqbar_embedding

    def get_coerce_map_from_fr_cyclotomic_field(self):
        r"""
        Return a coercion map from the associated :class:`FusionRing`'s
        cyclotomic field into the base field containing all F-symbols
        (this could be the :class:`FusionRing`'s
        :class:`Cyclotomic field<sage.rings.number_field.number_field.CyclotomicFieldFactory>`,
        a :func:`NumberField`, or :class:`QQbar<AlgebraicField>`).

        EXAMPLES::

            sage: f = FusionRing("G2", 1).get_fmatrix()
            sage: f.find_orthogonal_solution(verbose=False)
            sage: f.FR().field()
            Cyclotomic Field of order 60 and degree 16
            sage: f.field()
            Number Field in a with defining polynomial y^32 - ... - 22*y^2 + 1
            sage: phi = f.get_coerce_map_from_fr_cyclotomic_field()
            sage: phi.domain() == f.FR().field()
            True
            sage: phi.codomain() == f.field()
            True

        When F-symbols are computed as elements of the associated
        :class:`FusionRing`'s base
        :class:`Cyclotomic field<sage.rings.number_field.number_field.CyclotomicFieldFactory>`,
        we have ``self.field() == self.FR().field()`` and this
        returns the identity map on ``self.field()``. ::

            sage: f = FusionRing("A2", 1).get_fmatrix()
            sage: f.find_orthogonal_solution(verbose=False)
            sage: phi = f.get_coerce_map_from_fr_cyclotomic_field()
            sage: f.field()
            Cyclotomic Field of order 48 and degree 16
            sage: f.field() == f.FR().field()
            True
            sage: phi.domain() == f.field()
            True
            sage: phi.is_identity()
            True
        """
        # If base field is different from associated FusionRing's CyclotomicField,
        # return coercion map
        try:
            return self._coerce_map_from_cyc_field
        # Otherwise, return identity map CyclotomicField <-> CyclotomicField
        except AttributeError:
            F = self._FR.field()
            return F.hom([F.gen()], F)

    def get_fvars_in_alg_field(self):
        r"""
        Return F-symbols as elements of the :class:`QQbar<AlgebraicField>`.

        This method uses the embedding defined by
        :meth:`get_qqbar_embedding` to coerce
        F-symbols into :class:`QQbar<AlgebraicField>`.

        EXAMPLES::

            sage: fr = FusionRing("G2", 1)
            sage: f = fr.get_fmatrix(fusion_label='g', inject_variables=True, new=True)
            creating variables fx1..fx5
            Defining fx0, fx1, fx2, fx3, fx4
            sage: f.find_orthogonal_solution(verbose=False)
            sage: f.field()
            Number Field in a with defining polynomial y^32 - ... - 22*y^2 + 1
            sage: f.get_fvars_in_alg_field()
            {(g1, g1, g1, g0, g1, g1): 1,
             (g1, g1, g1, g1, g0, g0): 0.61803399? + 0.?e-8*I,
             (g1, g1, g1, g1, g0, g1): -0.7861514? + 0.?e-8*I,
             (g1, g1, g1, g1, g1, g0): -0.7861514? + 0.?e-8*I,
             (g1, g1, g1, g1, g1, g1): -0.61803399? + 0.?e-8*I}
        """
        return {sextuple: self._qqbar_embedding(fvar) for sextuple, fvar in self._fvars.items()}

    def get_radical_expression(self):
        """
        Return a radical expression of F-symbols.

        EXAMPLES::

            sage: f = FusionRing("G2", 1).get_fmatrix()
            sage: f.FR().fusion_labels("g", inject_variables=True)
            sage: f.find_orthogonal_solution(verbose=False)
            sage: radical_fvars = f.get_radical_expression()       # long time (~1.5s)
            sage: radical_fvars[g1, g1, g1, g1, g1, g0]            # long time
            -sqrt(1/2*sqrt(5) - 1/2)
        """
        return {sextuple: val.radical_expression() for sextuple, val in self.get_fvars_in_alg_field().items()}

    #######################
    #   Private helpers   #
    #######################

    def _get_known_vals(self):
        r"""
        Construct a dictionary of ``idx``, ``known_val`` pairs used for
        substituting into remaining equations.

        EXAMPLES::

            sage: f = FusionRing("D4", 1).get_fmatrix()
            sage: f._reset_solver_state()
            sage: len(f._get_known_vals()) == 0
            True
            sage: f.find_orthogonal_solution(verbose=False)
            sage: len(f._get_known_vals()) == f._poly_ring.ngens()
            True
        """
        return {i: self._fvars[s] for i, s in self._idx_to_sextuple.items() if self._solved[i]}

    def _get_known_nonz(self):
        r"""
        Construct an :class:`ETuple` indicating positions of
        known nonzero variables.

        .. NOTE::

            MUST be called after ``self._ks = _get_known_sq()``.
            This method is called by the constructor of ``self``.

        EXAMPLES::

            sage: f = FusionRing("D5", 1).get_fmatrix()  # indirect doctest
            sage: f._reset_solver_state()
            sage: f._nnz
            (100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
        """
        nonz = {idx: 100 for idx in self._singles}
        for idx, _ in self._ks.items():
            nonz[idx] = 100
        return ETuple(nonz, self._poly_ring.ngens())

    ##############################
    #   Variables partitioning   #
    ##############################

    def largest_fmat_size(self):
        r"""
        Get the size of the largest F-matrix `F^{abc}_d`.

        EXAMPLES::

            sage: f = FusionRing("B3", 2).get_fmatrix()
            sage: f.largest_fmat_size()
            4
        """
        return max(self.fmatrix(*tup).nrows() for tup in product(self._FR.basis(), repeat=4))

    def get_fvars_by_size(self, n, indices=False):
        r"""
        Return the set of F-symbols that are entries of an `n \times n` matrix
        `F^{a, b, c}_d`.

        INPUT:

        - ``n`` -- positive integer
        - ``indices`` -- boolean (default: ``False``)

        If ``indices`` is ``False`` (default),
        this method returns a set of sextuples `(a, b, c, d, x, y)` identifying
        the corresponding F-symbol. Each sextuple is a key in the
        dictionary returned by :meth:`get_fvars`.

        Otherwise the method returns a list of integer indices that
        internally identify the F-symbols. The ``indices=True`` option is
        meant for internal use.

        EXAMPLES::

            sage: f = FusionRing("A2", 2).get_fmatrix(inject_variables=True)
            creating variables fx1..fx287
            Defining fx0, ..., fx286
            sage: f.largest_fmat_size()
            2
            sage: f.get_fvars_by_size(2)
            {(f2, f2, f2, f4, f1, f1),
             (f2, f2, f2, f4, f1, f5),
             ...
             (f4, f4, f4, f4, f4, f0),
             (f4, f4, f4, f4, f4, f4)}
        """
        var_set = set()
        one = self._FR.one()
        for a, b, c, d in product(self._FR.basis(), repeat=4):
            X = self.f_from(a, b, c, d)
            Y = self.f_to(a, b, c, d)
            if len(X) == n and len(Y) == n:
                for x in X:
                    for y in Y:
                        # Discard trivial 1x1 F-matrix
                        trivial = a == one and x == b and y == d
                        trivial |= b == one and x == a and y == c
                        trivial |= c == one and x == d and y == b
                        if not trivial:
                            var_set.add((a, b, c, d, x, y))
        if indices:
            sext_to_idx = {v: k for k, v in self._idx_to_sextuple.items()}
            return {sext_to_idx[fx] for fx in var_set}
        return var_set

    ############################
    #   Checkpoint utilities   #
    ############################

    def save_fvars(self, filename):
        r"""
        Save computed F-symbols for later use.

        INPUT:

        - ``filename`` -- string specifying the name of the pickle file
          to be used

        The current directory is used unless an absolute path to a file in
        a different directory is provided.

        .. NOTE::

            This method should only be used *after* successfully running one
            of the solvers, e.g. :meth:`find_cyclotomic_solution` or
            :meth:`find_orthogonal_solution`.

        When used in conjunction with :meth:`load_fvars`, this method may
        be used to restore state of an :class:`FMatrix` object at the end
        of a successful F-matrix solver run.

        EXAMPLES::

            sage: f = FusionRing("A2", 1).get_fmatrix(new=True)
            sage: f.find_orthogonal_solution(verbose=False)
            sage: fvars = f.get_fvars()
            sage: K = f.field()
            sage: filename = f.get_fr_str() + "_solver_results.pickle"
            sage: f.save_fvars(filename)
            sage: del f
            sage: f2 = FusionRing("A2", 1).get_fmatrix(new=True)
            sage: f2.load_fvars(filename)
            sage: fvars == f2.get_fvars()
            True
            sage: K == f2.field()
            True
            sage: os.remove(filename)
        """
        final_state = [self._fvars,
                       self._non_cyc_roots,
                       self.get_coerce_map_from_fr_cyclotomic_field(),
                       self._qqbar_embedding]
        with open(filename, 'wb') as f:
            pickle.dump(final_state, f)

    def load_fvars(self, filename):
        r"""
        Load previously computed F-symbols from a pickle file.

        See :meth:`save_fvars` for more information.

        EXAMPLES::

            sage: f = FusionRing("A2", 1).get_fmatrix(new=True)
            sage: f.find_orthogonal_solution(verbose=False)
            sage: fvars = f.get_fvars()
            sage: K = f.field()
            sage: filename = f.get_fr_str() + "_solver_results.pickle"
            sage: f.save_fvars(filename)
            sage: del f
            sage: f2 = FusionRing("A2", 1).get_fmatrix(new=True)
            sage: f2.load_fvars(filename)
            sage: fvars == f2.get_fvars()
            True
            sage: K == f2.field()
            True
            sage: os.remove(filename)

        .. NOTE::

            :meth:`save_fvars`. This method does not work with intermediate
            checkpoint pickles; it only works with pickles containing *all*
            F-symbols, i.e. those created by :meth:`save_fvars` and by
            specifying an optional ``save_results`` parameter for
            :meth:`find_orthogonal_solution`.
        """
        with open(filename, 'rb') as f:
            self._fvars, self._non_cyc_roots, self._coerce_map_from_cyc_field, self._qqbar_embedding = pickle.load(f)
        # Update state attributes
        self._chkpt_status = 7
        self._solved = [True for v in self._fvars]
        self._field = self._qqbar_embedding.domain()

    def get_fr_str(self):
        r"""
        Auto-generate an identifying key for saving results.

        EXAMPLES::

            sage: f = FusionRing("B3", 1).get_fmatrix()
            sage: f.get_fr_str()
            'B31'
        """
        ct = self._FR.cartan_type()
        return ct.letter + str(ct.n) + str(self._FR.fusion_level())

    def _checkpoint(self, do_chkpt, status, verbose=True):
        r"""
        Pickle current solver state.

        EXAMPLES::

            sage: f = FusionRing("A1", 3).get_fmatrix(new=True)
            sage: f._reset_solver_state()
            sage: f.get_orthogonality_constraints(output=False)
            sage: f.get_defining_equations('hexagons', output=False)
            sage: f.ideal_basis = f._par_graph_gb(verbose=False)
            sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_tup_sortkey, poly_to_tup
            sage: f.ideal_basis.sort(key=poly_tup_sortkey)
            sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
            sage: n = f._poly_ring.ngens()
            sage: f._fvars = FvarsHandler(n, f._field, f._idx_to_sextuple, init_data=f._fvars)
            sage: f._triangular_elim(verbose=False)
            sage: f._update_reduction_params()
            sage: f._checkpoint(do_chkpt=True, status=2)
            Checkpoint 2 reached!
            sage: del f
            sage: f = FusionRing("A1", 3).get_fmatrix(new=True)
            sage: f.find_orthogonal_solution(warm_start='fmatrix_solver_checkpoint_A13.pickle')
            Computing F-symbols for The Fusion Ring of Type A1 and level 3 with Integer Ring coefficients with 71 variables...
            Set up 121 reduced pentagons...
            Elimination epoch completed... 18 eqns remain in ideal basis
            Elimination epoch completed... 5 eqns remain in ideal basis
            Pent elim step solved for 64 / 71 variables
            Partitioned 5 equations into 1 components of size:
            [4]
            Elimination epoch completed... 0 eqns remain in ideal basis
            Partitioned 6 equations into 6 components of size:
            [1, 1, 1, 1, 1, 1]
            Computing appropriate NumberField...
            sage: f._chkpt_status == 7
            True
            sage: sum(f._solved) == f._poly_ring.ngens()
            True
            sage: os.remove("fmatrix_solver_checkpoint_A13.pickle")
            sage: f = FusionRing("A1", 2).get_fmatrix(new=True)
            sage: f._reset_solver_state()
            sage: f.get_orthogonality_constraints(output=False)
            sage: f.get_defining_equations('hexagons', output=False)
            sage: f.ideal_basis = f._par_graph_gb(verbose=False)
            sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_tup_sortkey
            sage: f.ideal_basis.sort(key=poly_tup_sortkey)
            sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
            sage: n = f._poly_ring.ngens()
            sage: f._fvars = FvarsHandler(n, f._field, f._idx_to_sextuple, init_data=f._fvars)
            sage: f._triangular_elim(verbose=False)
            sage: f._update_reduction_params()
            sage: f.get_defining_equations('pentagons', output=False)
            sage: f.ideal_basis.sort(key=poly_tup_sortkey)
            sage: f._triangular_elim(verbose=False)
            sage: f._checkpoint(do_chkpt=True, status=4)
            Checkpoint 4 reached!
            sage: del f
            sage: f = FusionRing("A1", 2).get_fmatrix(new=True)
            sage: f.find_orthogonal_solution(warm_start='fmatrix_solver_checkpoint_A12.pickle')
            Computing F-symbols for The Fusion Ring of Type A1 and level 2 with Integer Ring coefficients with 14 variables...
            Partitioned 0 equations into 0 components of size:
            []
            Partitioned 2 equations into 2 components of size:
            [1, 1]
            sage: f._chkpt_status == 7
            True
            sage: sum(f._solved) == f._poly_ring.ngens()
            True
            sage: os.remove("fmatrix_solver_checkpoint_A12.pickle")
        """
        if not do_chkpt:
            return
        filename = "fmatrix_solver_checkpoint_" + self.get_fr_str() + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump([self._fvars, list(self._solved), self._ks, self.ideal_basis, status], f)
        if verbose:
            print(f"Checkpoint {status} reached!")

    def _restore_state(self, filename):
        r"""
        Load solver state from file. Use this method both for warm-starting
        :meth:`find_orthogonal_solution` and to load pickled results.

        EXAMPLES::

            sage: f = FusionRing("A1", 2).get_fmatrix(new=True)
            sage: f._reset_solver_state()
            sage: f.get_orthogonality_constraints(output=False)
            sage: f.get_defining_equations('hexagons', output=False)
            sage: f.ideal_basis = f._par_graph_gb(verbose=False)
            sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_tup_sortkey, poly_to_tup
            sage: f.ideal_basis.sort(key=poly_tup_sortkey)
            sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
            sage: n = f._poly_ring.ngens()
            sage: f._fvars = FvarsHandler(n, f._field, f._idx_to_sextuple, init_data=f._fvars)
            sage: f._triangular_elim(verbose=False)
            sage: f._update_reduction_params()
            sage: fvars = f._fvars
            sage: ib = f.ideal_basis
            sage: solved = f._solved
            sage: ks = f._ks
            sage: status = f._chkpt_status
            sage: f._checkpoint(do_chkpt=True, status=2)
            Checkpoint 2 reached!
            sage: del f
            sage: f = FusionRing("A1", 2).get_fmatrix(new=True)
            sage: f._reset_solver_state()
            sage: f._restore_state("fmatrix_solver_checkpoint_A12.pickle")
            sage: for sextuple, fvar in fvars.items():
            ....:     assert fvar == f._fvars[sextuple]
            ....:
            sage: ib == f.ideal_basis
            True
            sage: ks == f._ks
            True
            sage: solved == f._solved
            True
            sage: 2 == f._chkpt_status
            True
            sage: os.remove("fmatrix_solver_checkpoint_A12.pickle")

        TESTS::

            sage: f = FusionRing("A1", 3).get_fmatrix(new=True)
            sage: f.find_orthogonal_solution(save_results='test.pickle', verbose=False)   # long time
            sage: del f
            sage: f = FusionRing("A1", 3).get_fmatrix(new=True)
            sage: f.find_orthogonal_solution(warm_start='test.pickle')                   # long time
            sage: f._chkpt_status == 7                                                   # long time
            True
            sage: os.remove("test.pickle")                                               # long time
        """
        with open(filename, 'rb') as f:
            state = pickle.load(f)
        # Loading saved results pickle
        if len(state) == 4:
            self.load_fvars(filename)
            self._chkpt_status = 7
            return
        self._fvars, self._solved, self._ks, self.ideal_basis, self._chkpt_status = state
        self._update_reduction_params()

    #################
    #   MapReduce   #
    #################

    def start_worker_pool(self, processes=None):
        """
        Initialize a ``multiprocessing`` worker pool for parallel processing,
        which may be used e.g. to set up defining equations using
        :meth:`get_defining_equations`.

        This method sets ``self``'s ``pool`` attribute. The worker
        pool may be used time and again. Upon initialization, each process
        in the pool attaches to the necessary shared memory resources.

        When you are done using the worker pool, use
        :meth:`shutdown_worker_pool` to close the pool and properly dispose
        of shared memory resources.

        INPUT:

        - ``processes`` -- integer indicating the number of workers
          in the pool; if left unspecified, the number of workers
          equals the number of processors available

        OUTPUT: boolean; whether a worker pool was successfully initialized

        EXAMPLES::

            sage: f = FusionRing("G2", 1).get_fmatrix(new=True)
            sage: f.start_worker_pool()
            sage: he = f.get_defining_equations('hexagons')
            sage: sorted(he)
            [fx0 - 1,
             fx2*fx3 + (zeta60^14 + zeta60^12 - zeta60^6 - zeta60^4 + 1)*fx4^2 + (zeta60^6)*fx4,
             fx1*fx3 + (zeta60^14 + zeta60^12 - zeta60^6 - zeta60^4 + 1)*fx3*fx4 + (zeta60^14 - zeta60^4)*fx3,
             fx1*fx2 + (zeta60^14 + zeta60^12 - zeta60^6 - zeta60^4 + 1)*fx2*fx4 + (zeta60^14 - zeta60^4)*fx2,
             fx1^2 + (zeta60^14 + zeta60^12 - zeta60^6 - zeta60^4 + 1)*fx2*fx3 + (-zeta60^12)*fx1]
            sage: pe = f.get_defining_equations('pentagons')
            sage: f.shutdown_worker_pool()

        .. WARNING::

            This method is needed to initialize the worker pool using the
            necessary shared memory resources. Simply using the
            ``multiprocessing.Pool`` constructor will not work with our
            class methods.

        .. WARNING::

            Failure to call :meth:`shutdown_worker_pool` may result in a memory
            leak, since shared memory resources outlive the process that created
            them.
        """
        try:
            set_start_method('fork')
        except RuntimeError:
            pass
        if not hasattr(self, '_nnz'):
            self._reset_solver_state()
        # Set up shared memory resource handlers
        n_proc = cpu_count() if processes is None else processes
        self._pid_list = shared_memory.ShareableList([0] * (n_proc+1))
        pids_name = self._pid_list.shm.name
        self._solved = shared_memory.ShareableList(self._solved)
        s_name = self._solved.shm.name
        self._var_degs = shared_memory.ShareableList(self._var_degs)
        vd_name = self._var_degs.shm.name
        n = self._poly_ring.ngens()
        self._ks = KSHandler(n, self._field, use_mp=True, init_data=self._ks)
        ks_names = self._ks.shm.name
        self._shared_fvars = FvarsHandler(n, self._field, self._idx_to_sextuple, use_mp=n_proc, pids_name=pids_name, init_data=self._fvars)
        fvar_names = self._shared_fvars.shm.name
        # Initialize worker pool processes
        args = (id(self), s_name, vd_name, ks_names, fvar_names, n_proc, pids_name)

        def init(fmats_id, solved_name, vd_name, ks_names, fvar_names, n_proc, pids_name):
            """
            Connect worker process to shared memory resources
            """
            fmats_obj = cast(fmats_id, py_object).value
            fmats_obj._solved = shared_memory.ShareableList(name=solved_name)
            fmats_obj._var_degs = shared_memory.ShareableList(name=vd_name)
            n = fmats_obj._poly_ring.ngens()
            K = fmats_obj._field
            fmats_obj._fvars = FvarsHandler(n, K, fmats_obj._idx_to_sextuple, name=fvar_names, use_mp=n_proc, pids_name=pids_name)
            fmats_obj._ks = KSHandler(n, K, name=ks_names, use_mp=True)

        self.pool = Pool(processes=n_proc, initializer=init, initargs=args)
        self._pid_list[0] = getpid()
        for i, p in enumerate(self.pool._pool):
            self._pid_list[i + 1] = p.pid
        # return True

    def shutdown_worker_pool(self):
        r"""
        Shutdown the given worker pool and dispose of shared memory resources
        created when the pool was set up using :meth:`start_worker_pool`.

        .. WARNING::

            Failure to call this method after using :meth:`start_worker_pool`
            to create a process pool may result in a memory
            leak, since shared memory resources outlive the process that
            created them.

        EXAMPLES::

            sage: f = FusionRing("A1", 3).get_fmatrix(new=True)
            sage: f.start_worker_pool()
            sage: he = f.get_defining_equations('hexagons')
            sage: f.shutdown_worker_pool()
        """
        if self.pool is not None:
            self.pool.close()
            self.pool = None
            self._solved.shm.unlink()
            self._var_degs.shm.unlink()
            self._ks.shm.unlink()
            self._shared_fvars.shm.unlink()
            self._pid_list.shm.unlink()
            del self.__dict__['_shared_fvars']

    def _map_triv_reduce(self, mapper, input_iter, worker_pool=None, chunksize=None, mp_thresh=None):
        r"""
        Apply the given mapper to each element of the given input iterable and
        return the results (with no duplicates) in a list.

        INPUT:

        - ``mapper`` -- string specifying the name of a function defined in
          the ``fast_parallel_fmats_methods`` module

        .. NOTE::

            If ``worker_pool`` is not provided, function maps and reduces on a
            single process.
            If ``worker_pool`` is provided, the function attempts to determine
            whether it should use multiprocessing based on the length of the
            input iterable. If it can't determine the length of the input
            iterable then it uses multiprocessing with the default chunksize of
            `1` unless a chunksize is provided.

        EXAMPLES::

            sage: f = FusionRing("A1", 2).get_fmatrix()
            sage: f._reset_solver_state()
            sage: len(f._map_triv_reduce('get_reduced_hexagons', [(0, 1, False)]))
            11
            sage: f.start_worker_pool()
            sage: mp_params = [(i, f.pool._processes, True) for i in range(f.pool._processes)]
            sage: len(f._map_triv_reduce('get_reduced_pentagons', mp_params, worker_pool=f.pool, chunksize=1, mp_thresh=0))
            33
            sage: f.shutdown_worker_pool()
        """
        if mp_thresh is None:
            mp_thresh = self.mp_thresh
        # Compute multiprocessing parameters
        if worker_pool is not None:
            try:
                n = len(input_iter)
            except (TypeError, ValueError, AttributeError):
                n = mp_thresh + 1
            if chunksize is None:
                chunksize = n // (worker_pool._processes**2) + 1
        no_mp = worker_pool is None or n < mp_thresh
        # Map phase
        input_iter = zip_longest([], input_iter, fillvalue=(mapper, id(self)))
        if no_mp:
            mapped = map(executor, input_iter)
        else:
            mapped = worker_pool.imap_unordered(executor, input_iter, chunksize=chunksize)
        # Reduce phase
        results = set()
        for child_eqns in mapped:
            if child_eqns is not None:
                results.update(child_eqns)
        results = list(results)
        return results

    ########################
    #   Equations set up   #
    ########################

    def get_orthogonality_constraints(self, output=True):
        r"""
        Get equations imposed on the F-matrix by orthogonality.

        INPUT:

        - ``output`` -- boolean

        OUTPUT:

        If ``output=True``, orthogonality constraints are returned as
        polynomial objects.

        Otherwise, the constraints are appended to ``self.ideal_basis``.
        They are stored in the internal tuple representation. The
        ``output=False`` option is meant mostly for internal use by the
        F-matrix solver.

        EXAMPLES::

            sage: f = FusionRing("B4", 1).get_fmatrix()
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
        """
        eqns = []
        for tup in product(self._FR.basis(), repeat=4):
            mat = self.fmatrix(*tup)
            eqns.extend((mat.T * mat - matrix.identity(mat.nrows())).coefficients())
        if output:
            return eqns
        self.ideal_basis.extend([poly_to_tup(eq) for eq in eqns])

    def get_defining_equations(self, option, output=True):
        r"""
        Get the equations defining the ideal generated by the hexagon or
        pentagon relations.

        INPUT:

        - ``option`` -- string determining equations to be set up:

          * ``'hexagons'`` -- get equations imposed on the F-matrix by
            the hexagon relations in the definition of a braided category

          * ``'pentagons'`` -- get equations imposed on the F-matrix by
            the pentagon relations in the definition of a monoidal category

        - ``output`` -- boolean (default: ``True``); whether
          results should be returned, where the equations will be polynomials.
          Otherwise, the constraints are appended to ``self.ideal_basis``.
          Constraints are stored in the internal tuple representation. The
          ``output=False`` option is meant only for internal use by the
          F-matrix solver. When computing the hexagon equations with the
          ``output=False`` option, the initial state of the F-symbols is used.

        .. NOTE::

            To set up the defining equations using parallel processing,
            use :meth:`start_worker_pool` to initialize multiple processes
            *before* calling this method.

        EXAMPLES::

            sage: f = FusionRing("B2", 1).get_fmatrix()
            sage: sorted(f.get_defining_equations('hexagons'))
            [fx7 + 1,
             fx6 - 1,
             fx2 + 1,
             fx0 - 1,
             fx11*fx12 + (-zeta32^8)*fx13^2 + (zeta32^12)*fx13,
             fx10*fx12 + (-zeta32^8)*fx12*fx13 + (zeta32^4)*fx12,
             fx10*fx11 + (-zeta32^8)*fx11*fx13 + (zeta32^4)*fx11,
             fx10^2 + (-zeta32^8)*fx11*fx12 + (-zeta32^12)*fx10,
             fx4*fx9 + fx7,
             fx3*fx8 - fx6,
             fx1*fx5 + fx2]
            sage: pe = f.get_defining_equations('pentagons')
            sage: len(pe)
            33
        """
        if not hasattr(self, '_nnz'):
            self._reset_solver_state()
        n_proc = self.pool._processes if self.pool is not None else 1
        params = [(child_id, n_proc, output) for child_id in range(n_proc)]
        eqns = self._map_triv_reduce('get_reduced_' + option, params,
                                     worker_pool=self.pool, chunksize=1,
                                     mp_thresh=0)
        if output:
            F = self._field
            for i, eq_tup in enumerate(eqns):
                eqns[i] = _unflatten_coeffs(F, eq_tup)
            return [self._tup_to_fpoly(p) for p in eqns]
        self.ideal_basis.extend(eqns)

    ############################
    #   Equations processing   #
    ############################

    def _tup_to_fpoly(self, eq_tup):
        r"""
        Assemble a polynomial object from its tuple representation.

        .. WARNING::

            This method avoids implicit casting when constructing a
            polynomial object, and may therefore lead to SEGFAULTs.
            It is meant for internal use by the F-matrix solver.

        This method is a left inverse of
        :meth:`sage.algebras.fusion_rings.poly_tup_engine.poly_to_tup`.

        EXAMPLES::

            sage: f = FusionRing("C3", 1).get_fmatrix()
            sage: f.start_worker_pool()
            sage: he = f.get_defining_equations('hexagons')
            sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
            sage: all(f._tup_to_fpoly(poly_to_tup(h)) for h in he)
            True
            sage: f.shutdown_worker_pool()
        """
        return _tup_to_poly(eq_tup, parent=self._poly_ring)

    def _update_reduction_params(self, eqns=None):
        r"""
        Update reduction parameters that are solver state attributes.

        EXAMPLES::

            sage: f = FusionRing("A1", 3).get_fmatrix()
            sage: f._reset_solver_state()
            sage: f.get_orthogonality_constraints(output=False)
            sage: f.start_worker_pool()
            sage: f.get_defining_equations('hexagons', output=False)
            sage: f.ideal_basis = f._par_graph_gb(verbose=False)
            sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_tup_sortkey, poly_to_tup
            sage: f.ideal_basis.sort(key=poly_tup_sortkey)
            sage: f.mp_thresh = 0
            sage: f._fvars = f._shared_fvars
            sage: f._triangular_elim(verbose=False)  # indirect doctest
            sage: f.ideal_basis
            []
            sage: f.shutdown_worker_pool()
        """
        if eqns is None:
            eqns = self.ideal_basis
        self._ks.update(eqns)
        for i, d in enumerate(get_variables_degrees(eqns, self._poly_ring.ngens())):
            self._var_degs[i] = d
        self._nnz = self._get_known_nonz()
        self._kp = compute_known_powers(self._var_degs, self._get_known_vals(), self._field.one())

    def _triangular_elim(self, eqns=None, verbose=True):
        r"""
        Perform triangular elimination of linear terms in two-term equations
        until no such terms exist.

        .. NOTE::

            For optimal usage of triangular elimination, pass in a
            *sorted* list of equations.

        EXAMPLES::

            sage: f = FusionRing("D3", 1).get_fmatrix()
            sage: f.get_defining_equations('hexagons', output=False)
            sage: f.get_orthogonality_constraints(output=False)
            sage: gb = f._par_graph_gb(verbose=False)
            sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_tup_sortkey, poly_to_tup
            sage: f.ideal_basis = sorted(gb, key=poly_tup_sortkey)
            sage: from sage.algebras.fusion_rings.shm_managers import FvarsHandler
            sage: n = f._poly_ring.ngens()
            sage: f._fvars = FvarsHandler(n, f._field, f._idx_to_sextuple, init_data=f._fvars)
            sage: f._triangular_elim()
            Elimination epoch completed... 0 eqns remain in ideal basis
            sage: f.ideal_basis
            []
        """
        if eqns is None:
            eqns = self.ideal_basis

        while True:
            linear_terms_exist = _solve_for_linear_terms(self, eqns)
            if not linear_terms_exist:
                break
            _backward_subs(self)
            # Compute new reduction params and update eqns
            self._update_reduction_params(eqns=eqns)
            if self.pool is not None and len(eqns) > self.mp_thresh:
                n = self.pool._processes
                chunks = [[] for i in range(n)]
                for i, eq_tup in enumerate(eqns):
                    chunks[i % n].append(eq_tup)
                eqns = chunks
            else:
                eqns = [eqns]
            eqns = self._map_triv_reduce('update_reduce', eqns, worker_pool=self.pool, mp_thresh=0)
            eqns.sort(key=poly_tup_sortkey)
            if verbose:
                print("Elimination epoch completed... {} eqns remain in ideal basis".format(len(eqns)))
        self.ideal_basis = eqns

    #####################
    #   Graph methods   #
    #####################

    def equations_graph(self, eqns=None):
        r"""
        Construct a graph corresponding to the given equations.

        Every node corresponds to a variable and nodes are connected when
        the corresponding variables appear together in an equation.

        INPUT:

        - ``eqns`` -- list of polynomials

        Each polynomial is either an object in the ring returned by
        :meth:`get_poly_ring` or it is a tuple of pairs representing
        a polynomial using the internal representation.

        If no list of equations is passed, the graph is built from the
        polynomials in ``self.ideal_basis``. In this case the method assumes
        the internal representation of a polynomial as a tuple of pairs is
        used.

        This method is crucial to :meth:`find_orthogonal_solution`. The
        hexagon equations, obtained using :meth:`get_defining_equations`,
        define a disconnected graph that breaks up into many small components.
        The :meth:`find_orthogonal_solution` solver exploits this when
        undertaking a Groebner basis computation.

        OUTPUT:

        A ``Graph`` object. If a list of polynomial objects was given,
        the set of nodes in the output graph is the subset polynomial
        ring generators appearing in the equations.

        If the internal representation was used, the set of nodes is
        the subset of indices corresponding to polynomial ring generators.
        This option is meant for internal use by the F-matrix solver.

        EXAMPLES::

            sage: f = FusionRing("A3", 1).get_fmatrix()
            sage: f.get_poly_ring().ngens()
            27
            sage: he = f.get_defining_equations('hexagons')
            sage: graph = f.equations_graph(he)
            sage: graph.connected_components_sizes()
            [6, 3, 3, 3, 3, 3, 3, 1, 1, 1]
        """
        if eqns is None:
            eqns = self.ideal_basis

        G = Graph()
        if not eqns:
            return G

        # Eqns could be a list of poly objects or poly tuples stored in internal repn
        if isinstance(eqns[0], tuple):
            G.add_vertices([x for eq_tup in eqns for x in variables(eq_tup)])
        else:
            G.add_vertices([x for eq in eqns for x in eq.variables()])
        for eq in eqns:
            # Eqns could be a list of poly objects or poly tuples stored in internal repn
            if isinstance(eq, tuple):
                s = list(variables(eq))
            else:
                s = list(eq.variables())
            for x in s:
                for y in s:
                    if y != x:
                        G.add_edge(x, y)
        return G

    def _partition_eqns(self, eqns=None, verbose=True):
        r"""
        Partition equations corresponding to edges in a disconnected graph.

        OUTPUT:

        This method returns a dictionary of (c, e) pairs, where
        c is a tuple denoting a connected component in the graph produced
        by calling :meth:`equations_graph` with the given ``eqns`` and
        e is a list of all equations with variables in c.

        EXAMPLES::

            sage: f = FusionRing("C2", 1).get_fmatrix()
            sage: f.get_defining_equations('hexagons', output=False)
            sage: partition = f._partition_eqns()
            Partitioned 11 equations into 5 components of size:
            [4, 3, 3, 3, 1]
            sage: from sage.algebras.fusion_rings.poly_tup_engine import variables
            sage: for c, e in partition.items():
            ....:     assert set(v for eq_tup in e for v in variables(eq_tup)) == set(c)
            sage: vars_in_partition = set()
            sage: eqns_in_partition = set()
            sage: for c, e in partition.items():
            ....:     vars_in_partition.update(c)
            ....:     eqns_in_partition.update(e)
            sage: vars_in_partition == set(v for eq_tup in f.ideal_basis for v in variables(eq_tup))
            True
            sage: eqns_in_partition == set(f.ideal_basis)
            True
            sage: from itertools import product
            sage: for e1, e2 in product(partition.values(), repeat=2):
            ....:     assert e1 == e2 or set(e1).isdisjoint(set(e2))
        """
        if eqns is None:
            eqns = self.ideal_basis
        graph = self.equations_graph(eqns)
        partition = {tuple(c): [] for c in graph.connected_components(sort=True)}
        for eq_tup in eqns:
            partition[tuple(graph.connected_component_containing_vertex(variables(eq_tup)[0], sort=True))].append(eq_tup)
        if verbose:
            print("Partitioned {} equations into {} components of size:".format(len(eqns), graph.connected_components_number()))
            print(graph.connected_components_sizes())
        return partition

    def _par_graph_gb(self, eqns=None, term_order='degrevlex', largest_comp=45, verbose=True):
        r"""
        Compute a Groebner basis for a list of equations partitioned
        according to their corresponding graph.

        .. NOTE::

            If the graph has more than 50 components, this method computes the
            Groebner basis in parallel when a ``worker_pool`` is provided.

            This method will refuse to find a Groebner basis for a component
            of size larger than 60, since such a calculation does not seem to
            terminate.

        EXAMPLES::

            sage: f = FusionRing("F4", 1).get_fmatrix()
            sage: f._reset_solver_state()
            sage: f.get_orthogonality_constraints(output=False)
            sage: f.start_worker_pool()
            sage: f.get_defining_equations('hexagons', output=False)
            sage: gb = f._par_graph_gb()
            Partitioned 10 equations into 2 components of size:
            [4, 1]
            sage: from sage.algebras.fusion_rings.poly_tup_engine import _unflatten_coeffs
            sage: ret = [f._tup_to_fpoly(_unflatten_coeffs(f.field(), t)) for t in gb]
            sage: ret.sort(); ret
            [fx4 + (-zeta80^24 + zeta80^16),
             fx2 - fx3,
             fx1 + (zeta80^24 - zeta80^16),
             fx0 - 1,
             fx3^2 + (zeta80^24 - zeta80^16)]
            sage: f.shutdown_worker_pool()
        """
        if eqns is None:
            eqns = self.ideal_basis
        small_comps = []
        temp_eqns = []
        for comp, comp_eqns in self._partition_eqns(eqns=eqns, verbose=verbose).items():
            # Check if component is too large to process
            if len(comp) > largest_comp:
                temp_eqns.extend(comp_eqns)
            else:
                small_comps.append(comp_eqns)
        input_iter = zip_longest(small_comps, [], fillvalue=term_order)
        small_comp_gb = self._map_triv_reduce('compute_gb', input_iter, worker_pool=self.pool, chunksize=1, mp_thresh=50)
        ret = small_comp_gb + temp_eqns
        return ret

    def _get_component_variety(self, var, eqns):
        r"""
        Translate equations in each connected component to smaller polynomial
        rings so we can call built-in variety method.

        INPUT:

        - ``var`` -- list of variable indices
        - ``eqns`` -- list of polynomial equations in the internal
          tuple of pairs representation

        EXAMPLES::

            sage: f = FusionRing("G2", 2).get_fmatrix(new=True)
            sage: f.start_worker_pool()
            sage: f.get_defining_equations('hexagons', output=False)                     # long time
            sage: f.shutdown_worker_pool()
            sage: partition = f._partition_eqns()                                        # long time
            Partitioned 327 equations into 35 components of size:
            [27, 27, 27, 24, 24, 16, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
             9, 9, 6, 6, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1]
            sage: c = (216, 292, 319)
            sage: from sage.algebras.fusion_rings.poly_tup_engine import poly_to_tup
            sage: eqns = partition[c] + [poly_to_tup(f._poly_ring.gen(216)-1)]           # long time
            sage: f._get_component_variety(c, eqns)                                      # long time
            [{216: -1, 292: -1, 319: 1}]
        """
        # Define smaller poly ring in component vars
        R = PolynomialRing(self._FR.field(), len(var), 'a', order='lex')

        # Zip tuples into R and compute Groebner basis
        idx_map = {old: new for new, old in enumerate(sorted(var))}
        nvars = len(var)
        eqns = [_unflatten_coeffs(self._field, eq_tup) for eq_tup in eqns]
        polys = [_tup_to_poly(resize(eq_tup, idx_map, nvars), parent=R) for eq_tup in eqns]
        var_in_R = Ideal(sorted(polys)).variety(ring=AA)

        # Change back to fmats poly ring and append to temp_eqns
        inv_idx_map = {v: k for k, v in idx_map.items()}
        return [{inv_idx_map[i]: value for i, (key, value) in enumerate(sorted(soln.items()))} for soln in var_in_R]

    #######################
    #   Solution method   #
    #######################

    # TODO: this can probably be improved by constructing a set of defining polynomials
    # and checking, one by one, if it's irreducible over the current field.
    # If it is, we construct an extension. Perhaps it's best to go one by one here...
    def attempt_number_field_computation(self):
        r"""
        Based on the ``CartanType`` of ``self`` and data
        known on March 17, 2021, determine whether to attempt
        to find a :func:`NumberField` containing all the F-symbols.

        This method is used by :meth:`find_orthogonal_solution`
        to determine a field containing all F-symbols.
        See :meth:`field` and :meth:`get_non_cyclotomic_roots`.

        For certain :class:`fusion rings <FusionRing>`, the number field
        computation does not terminate in reasonable time.
        In these cases, we report F-symbols as elements
        of the :class:`QQbar<AlgebraicField>`.

        EXAMPLES::

            sage: f = FusionRing("F4", 2).get_fmatrix()
            sage: f.attempt_number_field_computation()
            False
            sage: f = FusionRing("G2", 1).get_fmatrix()
            sage: f.attempt_number_field_computation()
            True

        .. NOTE::

            In certain cases, F-symbols are found in the associated
            :class:`FusionRing`'s cyclotomic field and a
            :func:`NumberField` computation is not needed. In these
            cases this method returns ``True`` but the
            :meth:`find_orthogonal_solution` solver does *not*
            undertake a :func:`NumberField` computation.
        """
        ct = self._FR.cartan_type()
        k = self._FR._k
        # Don't try when k is large and odd for SU(2)_k
        if ct.letter == 'A':
            if ct.n == 1 and k >= 9 and k % 2:
                return False
        if ct.letter == 'C':
            if ct.n >= 9 and ct.n % 2 and k == 1:
                return False
        if ct.letter == 'E':
            if ct.n < 8 and k == 2:
                return False
            if ct.n == 8 and k == 3:
                return False
        if ct.letter == 'F' and k == 2:
            return False
        if ct.letter == 'G' and k == 2:
            return False
        return True

    def _get_explicit_solution(self, eqns=None, verbose=True):
        r"""
        Construct an explicit solution of ``self``.

        When this method is called, the solution is already found in
        terms of Groeber basis. A few degrees of freedom remain.
        By specializing the free variables and back substituting, a
        solution in the base field is now obtained.

        EXAMPLES::

            sage: f = FusionRing("A1", 3).get_fmatrix()  # indirect doctest
            sage: f.find_orthogonal_solution()           # long time
            Computing F-symbols for The Fusion Ring of Type A1 and level 3 with Integer Ring coefficients with 71 variables...
            Set up 134 hex and orthogonality constraints...
            Partitioned 134 equations into 17 components of size:
            [12, 12, 6, 6, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1]
            Elimination epoch completed... 10 eqns remain in ideal basis
            Elimination epoch completed... 0 eqns remain in ideal basis
            Hex elim step solved for 51 / 71 variables
            Set up 121 reduced pentagons...
            Elimination epoch completed... 18 eqns remain in ideal basis
            Elimination epoch completed... 5 eqns remain in ideal basis
            Pent elim step solved for 64 / 71 variables
            Partitioned 5 equations into 1 components of size:
            [4]
            Elimination epoch completed... 0 eqns remain in ideal basis
            Partitioned 6 equations into 6 components of size:
            [1, 1, 1, 1, 1, 1]
            Computing appropriate NumberField...
        """
        if eqns is None:
            eqns = self.ideal_basis
        # Don't add square fixers when warm starting from a late-stage checkpoint
        if self._chkpt_status < 5:
            n = self._poly_ring.ngens()
            one = self._field.one()
            for fx, rhs in self._ks.items():
                if not self._solved[fx]:
                    lt = (ETuple({fx: 2}, n), one)
                    eqns.append((lt, (ETuple({}, n), -rhs)))
        eqns_partition = self._partition_eqns(verbose=verbose)

        F = self._field
        R = F['x']
        numeric_fvars = {}
        non_cyclotomic_roots = []
        must_change_base_field = False
        phi = F.hom([F.gen()], F)
        for comp, part in eqns_partition.items():
            # If component has only one equation in a single variable, get a root
            if len(comp) == 1 and len(part) == 1:
                # Attempt to find cyclotomic root
                univ_poly = tup_to_univ_poly(part[0], R)
                roots = univ_poly.roots(multiplicities=False)
                if roots:
                    numeric_fvars[comp[0]] = roots[0]
                else:
                    # A real solution is preferred
                    roots = univ_poly.roots(ring=AA, multiplicities=False)
                    if not roots:
                        roots = univ_poly.roots(ring=QQbar, multiplicities=False)
                    non_cyclotomic_roots.append((comp[0], roots[0]))
                    must_change_base_field = True
            # Otherwise, compute the component variety and select a point to obtain a numerical solution
            else:
                sols = self._get_component_variety(comp, part)
                for fx, rhs in sols[0].items():
                    non_cyclotomic_roots.append((fx, rhs))
                must_change_base_field = True

        if must_change_base_field:
            # Attempt to compute smallest number field containing all the F-symbols
            # If calculation takes too long, we use QQbar as the base field
            if self.attempt_number_field_computation():
                if verbose:
                    print("Computing appropriate NumberField...")
                roots = [self._FR.field().gen()] + [r[1] for r in non_cyclotomic_roots]
                self._field, bf_elts, self._qqbar_embedding = number_field_elements_from_algebraics(roots, minimal=True)
            else:
                self._field = QQbar
                bf_elts = [self._qqbar_embedding(F.gen())]
                bf_elts += [rhs for fx, rhs in non_cyclotomic_roots]
                self._qqbar_embedding = lambda x: x
            self._non_cyc_roots = bf_elts[1:]

            # Embed cyclotomic field into newly constructed base field
            cyc_gen_as_bf_elt = bf_elts.pop(0)
            phi = self._FR.field().hom([cyc_gen_as_bf_elt], self._field)
            self._coerce_map_from_cyc_field = phi
            numeric_fvars = {k: phi(v) for k, v in numeric_fvars.items()}
            for i, elt in enumerate(bf_elts):
                numeric_fvars[non_cyclotomic_roots[i][0]] = elt
            # Update polynomial ring
            self._poly_ring = self._poly_ring.change_ring(self._field)

        # Ensure all F-symbols are known
        for fx in numeric_fvars:
            self._solved[fx] = True
        nvars = self._poly_ring.ngens()
        assert sum(self._solved) == nvars, "Some F-symbols are still missing...{}".format([self._poly_ring.gen(fx) for fx in range(nvars) if not self._solved[fx]])

        # Backward substitution step. Traverse variables in reverse lexicographical order. (System is in triangular form)
        self._fvars = {sextuple: apply_coeff_map(rhs, phi) for sextuple, rhs in self._fvars.items()}
        for fx, rhs in numeric_fvars.items():
            self._fvars[self._idx_to_sextuple[fx]] = ((ETuple({}, nvars), rhs), )
        _backward_subs(self, flatten=False)
        self._fvars = {sextuple: constant_coeff(rhs, self._field) for sextuple, rhs in self._fvars.items()}

        # Update base field attributes
        self._FR._field = self.field()
        self._FR._basecoer = self.get_coerce_map_from_fr_cyclotomic_field()
        if self._FR._basecoer:
            self._FR.r_matrix.clear_cache()

    def find_orthogonal_solution(self, checkpoint=False, save_results='', warm_start='', use_mp=True, verbose=True):
        r"""
        Solve the the hexagon and pentagon relations, along with
        orthogonality constraints, to evaluate an orthogonal F-matrix.

        INPUT:

        - ``checkpoint`` -- boolean (default: ``False``); whether
          the computation should be checkpointed. Depending on the associated
          ``CartanType``, the computation may take hours to complete. For
          large examples, checkpoints are recommended. This method supports
          "warm" starting, so the calculation may be resumed from a checkpoint,
          using the ``warm_start`` option.

          Checkpoints store necessary state in the pickle file
          ``"fmatrix_solver_checkpoint_" + key + ".pickle"``, where ``key``
          is the result of :meth:`get_fr_str`.

          Checkpoint pickles are automatically deleted when the solver exits
          a successful run.

        - ``save_results`` -- (optional) a string indicating the name of a
          pickle file in which to store calculated F-symbols for later use.

          If ``save_results`` is not provided (default), F-matrix results
          are not stored to file.

          The F-symbols may be saved to file after running the solver using
          :meth:`save_fvars`.

        - ``warm_start`` -- (optional) a string indicating the name of a pickle
          file containing checkpointed solver state. This file must have been
          produced by a previous call to the solver using the ``checkpoint``
          option.

          If no file name is provided, the calculation begins from scratch.

        - ``use_mp`` -- boolean (default: ``True``); whether to use
          multiprocessing to speed up calculation. The default value
          ``True`` is highly recommended, since parallel processing yields
          results much more quickly.

        - ``verbose`` -- boolean (default: ``True``); whether the
          solver should print out intermediate progress reports

        OUTPUT:

        This method returns ``None``. If the solver runs successfully, the
        results may be accessed through various methods, such as
        :meth:`get_fvars`, :meth:`fmatrix`, :meth:`fmat`, etc.

        EXAMPLES::

            sage: f = FusionRing("B5", 1).get_fmatrix(fusion_label='b', inject_variables=True)
            creating variables fx1..fx14
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
            sage: f.find_orthogonal_solution()
            Computing F-symbols for The Fusion Ring of Type B5 and level 1 with Integer Ring coefficients with 14 variables...
            Set up 25 hex and orthogonality constraints...
            Partitioned 25 equations into 5 components of size:
            [4, 3, 3, 3, 1]
            Elimination epoch completed... 0 eqns remain in ideal basis
            Hex elim step solved for 10 / 14 variables
            Set up 7 reduced pentagons...
            Elimination epoch completed... 0 eqns remain in ideal basis
            Pent elim step solved for 12 / 14 variables
            Partitioned 0 equations into 0 components of size:
            []
            Partitioned 2 equations into 2 components of size:
            [1, 1]
            sage: f.fmatrix(b2, b2, b2, b2)
            [ 1/2*zeta80^30 - 1/2*zeta80^10 -1/2*zeta80^30 + 1/2*zeta80^10]
            [ 1/2*zeta80^30 - 1/2*zeta80^10  1/2*zeta80^30 - 1/2*zeta80^10]
            sage: f.fmat(b2, b2, b2, b2, b0, b1)
            -1/2*zeta80^30 + 1/2*zeta80^10

        Every F-matrix `F^{a, b, c}_d` is orthogonal and in many cases real.
        We may use :meth:`fmats_are_orthogonal` and :meth:`fvars_are_real`
        to obtain correctness certificates.

        EXAMPLES::

            sage: f.fmats_are_orthogonal()
            True

        In any case, the F-symbols are obtained as elements of the associated
        :class:`FusionRing`'s
        :class:`Cyclotomic field<sage.rings.number_field.number_field.CyclotomicFieldFactory>`,
        a computed :func:`NumberField`, or :class:`QQbar<AlgebraicField>`.
        Currently, the field containing the F-symbols is determined based
        on the ``CartanType`` associated to ``self``.

        .. SEEALSO::

            :meth:`attempt_number_field_computation`
        """
        if self._poly_ring.ngens() == 0:
            return
        self._reset_solver_state()

        # Resume computation from checkpoint
        if warm_start:
            self._restore_state(warm_start)
            # Loading from a pickle with solved F-symbols
            if self._chkpt_status > 5:
                return
        if use_mp:
            self.start_worker_pool()
        if verbose:
            print("Computing F-symbols for {} with {} variables...".format(self._FR, self._poly_ring.ngens()))

        if self._chkpt_status < 1:
            # Set up hexagon equations and orthogonality constraints
            self.get_orthogonality_constraints(output=False)
            self.get_defining_equations('hexagons', output=False)
            # Report progress
            if verbose:
                print("Set up {} hex and orthogonality constraints...".format(len(self.ideal_basis)))

        # Unzip _fvars and link to shared_memory structure if using multiprocessing
        if use_mp:  # and loads_shared_memory:
            self._fvars = self._shared_fvars
        else:
            n = self._poly_ring.ngens()
            self._fvars = FvarsHandler(n, self._field, self._idx_to_sextuple, init_data=self._fvars)
        self._checkpoint(checkpoint, 1, verbose=verbose)

        if self._chkpt_status < 2:
            # Set up equations graph. Find GB for each component in parallel. Eliminate variables
            self.ideal_basis = self._par_graph_gb(verbose=verbose)
            self.ideal_basis.sort(key=poly_tup_sortkey)
            self._triangular_elim(verbose=verbose)
            # Report progress
            if verbose:
                print("Hex elim step solved for {} / {} variables".format(sum(self._solved), len(self._poly_ring.gens())))
        self._checkpoint(checkpoint, 2, verbose=verbose)

        if self._chkpt_status < 3:
            # Set up pentagon equations in parallel
            self.get_defining_equations('pentagons', output=False)
            # Report progress
            if verbose:
                print("Set up {} reduced pentagons...".format(len(self.ideal_basis)))
        self._checkpoint(checkpoint, 3, verbose=verbose)

        if self._chkpt_status < 4:
            # Simplify and eliminate variables
            self.ideal_basis.sort(key=poly_tup_sortkey)
            self._triangular_elim(verbose=verbose)
            # Report progress
            if verbose:
                print("Pent elim step solved for {} / {} variables".format(sum(self._solved), len(self._poly_ring.gens())))
        self._checkpoint(checkpoint, 4, verbose=verbose)

        # Try adding degrevlex gb -> elim loop until len(ideal_basis) does not change

        # Set up new equations graph and compute variety for each component
        if self._chkpt_status < 5:
            self.ideal_basis = self._par_graph_gb(term_order='lex', verbose=verbose)
            self.ideal_basis.sort(key=poly_tup_sortkey)
            self._triangular_elim(verbose=verbose)
        self._checkpoint(checkpoint, 5, verbose=verbose)
        self.shutdown_worker_pool()

        # Find numeric values for each F-symbol
        self._get_explicit_solution(verbose=verbose)
        # The calculation was successful, so we may delete checkpoints
        self._chkpt_status = 7
        self.clear_equations()
        if checkpoint:
            remove("fmatrix_solver_checkpoint_" + self.get_fr_str() + ".pickle")
        if save_results:
            self.save_fvars(save_results)

    #########################
    #   Cyclotomic method   #
    #########################

    def _fix_gauge(self, algorithm=""):
        r"""
        Fix the gauge by forcing F-symbols not already fixed to equal `1`.

        .. NOTE::

            This method should be used *after* adding hexagon and pentagon
            equations to ``self.ideal_basis``.

        EXAMPLES::

            sage: f = FusionRing("A3", 1).get_fmatrix()
            sage: f._reset_solver_state()                             # long time
            sage: f._var_to_sextuple = {f._poly_ring.gen(i): s for i, s in f._idx_to_sextuple.items()}  # long time
            sage: eqns = f.get_defining_equations("hexagons")+f.get_defining_equations("pentagons")     # long time
            sage: f.ideal_basis = set(Ideal(eqns).groebner_basis())   # long time
            sage: _, _ = f._substitute_degree_one()                   # long time
            sage: f._fix_gauge()                                      # long time
            adding equation... fx1 - 1
            adding equation... fx18 - 1
            adding equation... fx21 - 1
        """
        while not all(self._solved):
            # Get a variable that has not been fixed
            # In ascending index order, for consistent results
            for i, var in enumerate(self._poly_ring.gens()):
                if not self._solved[i]:
                    break

            # Fix var = 1, substitute, and solve equations
            self.ideal_basis.add(var - 1)
            print("adding equation...", var - 1)
            self.ideal_basis = set(Ideal(list(self.ideal_basis)).groebner_basis(algorithm=algorithm))
            self._substitute_degree_one()
            self._update_equations()

    def _substitute_degree_one(self, eqns=None):
        r"""
        Substitute known value from linear univariate polynomial and
        solve, following [Bond2007]_ p.37, for two-term linear equation
        for one of the variables. See also [Ab2022]_.

        EXAMPLES::

            sage: fr = FusionRing("D3", 1)
            sage: f = fr.get_fmatrix(inject_variables=True, new=True)
            creating variables fx1..fx27
            Defining fx0, ..., fx26
            sage: f._reset_solver_state()
            sage: f._var_to_sextuple = {f._poly_ring.gen(i): s for i, s in f._idx_to_sextuple.items()}
            sage: f.ideal_basis = [fx0 - 8, fx4**2 - 3, fx4 + fx10 + 3, fx4 + fx9]
            sage: _, _ = f._substitute_degree_one()
            sage: f._fvars[f._var_to_sextuple[fx0]]
            8
            sage: f._fvars[f._var_to_sextuple[fx4]]
            -fx9
        """
        if eqns is None:
            eqns = self.ideal_basis

        new_knowns = set()
        useless = set()
        for eq in eqns:
            if eq.degree() == 1 and sum(eq.degrees()) <= 2 and eq.lm() not in self._solved:
                self._fvars[self._var_to_sextuple[eq.lm()]] = -sum(c * m for c, m in zip(eq.coefficients()[1:], eq.monomials()[1:])) / eq.lc()
                # Add variable to set of known values and remove this equation
                new_knowns.add(eq.lm())
                useless.add(eq)

        # Update fvars depending on other variables
        for idx, fx in enumerate(self._poly_ring.gens()):
            if fx in new_knowns:
                self._solved[idx] = fx
        for sextuple, rhs in self._fvars.items():
            d = {var: self._fvars[self._var_to_sextuple[var]] for var in rhs.variables() if var in self._solved}
            if d:
                self._fvars[sextuple] = rhs.subs(d)
        return new_knowns, useless

    def _update_equations(self):
        r"""
        Perform backward substitution on equations in ``self.ideal_basis``.

        EXAMPLES::

            sage: fr = FusionRing("D3", 1)
            sage: f = fr.get_fmatrix(inject_variables=True, new=True)
            creating variables fx1..fx27
            Defining fx0, ..., fx26
            sage: f._reset_solver_state()
            sage: f._var_to_sextuple = {f._poly_ring.gen(i): s for i, s in f._idx_to_sextuple.items()}
            sage: f.ideal_basis = [fx0 - 8, fx4 + fx9, fx4**2 + fx3 - fx9**2]
            sage: _, _ = f._substitute_degree_one()
            sage: f._update_equations()
            sage: f.ideal_basis
            {fx3}
        """
        special_values = {known: self._fvars[self._var_to_sextuple[known]] for known in self._solved if known}
        self.ideal_basis = {eq.subs(special_values) for eq in self.ideal_basis}
        self.ideal_basis.discard(0)

    def find_cyclotomic_solution(self, equations=None, algorithm='', verbose=True, output=False):
        r"""
        Solve the hexagon and pentagon relations to evaluate the F-matrix.

        This method (omitting the orthogonality constraints) produces
        output in the cyclotomic field, but it is very limited in the size
        of examples it can handle: for example, `G_2` at level 2 is
        too large for this method. You may use :meth:`find_orthogonal_solution`
        to solve much larger examples.

        INPUT:

        - ``equations`` -- (optional) a set of equations to be
          solved; defaults to the hexagon and pentagon equations
        - ``algorithm`` -- (optional) algorithm to compute Groebner Basis
        - ``output`` -- boolean (default: ``False``); output a dictionary of
          F-matrix values; this may be useful to see but may be omitted
          since this information will be available afterwards via the
          :meth:`fmatrix` and :meth:`fmat` methods.

        EXAMPLES::

            sage: fr = FusionRing("A2", 1, fusion_labels='a', inject_variables=True)
            sage: f = fr.get_fmatrix(inject_variables=True)
            creating variables fx1..fx8
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7
            sage: f.find_cyclotomic_solution(output=True)
            Setting up hexagons and pentagons...
            Finding a Groebner basis...
            Solving...
            Fixing the gauge...
            adding equation... fx4 - 1
            Done!
            {(a2, a2, a2, a0, a1, a1): 1,
             (a2, a2, a1, a2, a1, a0): 1,
             (a2, a1, a2, a2, a0, a0): 1,
             (a2, a1, a1, a1, a0, a2): 1,
             (a1, a2, a2, a2, a0, a1): 1,
             (a1, a2, a1, a1, a0, a0): 1,
             (a1, a1, a2, a1, a2, a0): 1,
             (a1, a1, a1, a0, a2, a2): 1}

        After you successfully run :meth:`find_cyclotomic_solution` you may
        check the correctness of the F-matrix by running
        :meth:`get_defining_equations` with ``option='hexagons'`` and
        ``option='pentagons'``. These should return empty lists
        of equations.

        EXAMPLES::

            sage: f.get_defining_equations("hexagons")
            []
            sage: f.get_defining_equations("pentagons")
            []
        """
        if self._poly_ring.ngens() == 0:
            return
        self._reset_solver_state()
        self._var_to_sextuple = {self._poly_ring.gen(i): s for i, s in self._idx_to_sextuple.items()}

        if equations is None:
            if verbose:
                print("Setting up hexagons and pentagons...")
            equations = self.get_defining_equations("hexagons") + self.get_defining_equations("pentagons")
        if verbose:
            print("Finding a Groebner basis...")
        self.ideal_basis = set(Ideal(equations).groebner_basis(algorithm=algorithm))
        if verbose:
            print("Solving...")
        self._substitute_degree_one()
        if verbose:
            print("Fixing the gauge...")
        self._fix_gauge(algorithm=algorithm)
        if verbose:
            print("Done!")
        if output:
            return self._fvars

    #####################
    #   Verifications   #
    #####################

    def fmats_are_orthogonal(self):
        r"""
        Verify that all F-matrices are orthogonal.

        This method should always return ``True`` when called after running
        :meth:`find_orthogonal_solution`.

        EXAMPLES::

            sage: f = FusionRing("D4", 1).get_fmatrix()
            sage: f.find_orthogonal_solution(verbose=False)
            sage: f.fmats_are_orthogonal()
            True
        """
        is_orthog = []
        for a, b, c, d in product(self._FR.basis(), repeat=4):
            mat = self.fmatrix(a, b, c, d)
            is_orthog.append(mat.T * mat == matrix.identity(mat.nrows()))
        return all(is_orthog)

    def fvars_are_real(self):
        r"""
        Test whether all F-symbols are real.

        EXAMPLES::

            sage: f = FusionRing("A1", 3).get_fmatrix()
            sage: f.find_orthogonal_solution(verbose=False) # long time
            sage: f.fvars_are_real()                        # not tested (cypari issue in doctesting framework)
            True
        """
        try:
            for k, v in self._fvars.items():
                AA(self._qqbar_embedding(v))
        except ValueError:
            print("the F-symbol {} (key {}) has a nonzero imaginary part".format(v, k))
            return False
        return True

    def certify_pentagons(self, use_mp=True, verbose=False):
        r"""
        Obtain a certificate of satisfaction for the pentagon equations,
        up to floating-point error.

        This method converts the computed F-symbols (available through
        :meth:`get_fvars`) to native Python floats and then checks whether
        the pentagon equations are satisfied using floating point arithmetic.

        When ``self.FR().basis()`` has many elements, verifying satisfaction
        of the pentagon relations exactly using :meth:`get_defining_equations`
        with ``option="pentagons"`` may take a long time. This method is
        faster, but it cannot provide mathematical guarantees.

        EXAMPLES::

            sage: f = FusionRing("C3", 1).get_fmatrix()
            sage: f.find_orthogonal_solution()        # long time
            Computing F-symbols for The Fusion Ring of Type C3 and level 1 with Integer Ring coefficients with 71 variables...
            Set up 134 hex and orthogonality constraints...
            Partitioned 134 equations into 17 components of size:
            [12, 12, 6, 6, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1]
            Elimination epoch completed... 10 eqns remain in ideal basis
            Elimination epoch completed... 0 eqns remain in ideal basis
            Hex elim step solved for 51 / 71 variables
            Set up 121 reduced pentagons...
            Elimination epoch completed... 18 eqns remain in ideal basis
            Elimination epoch completed... 5 eqns remain in ideal basis
            Pent elim step solved for 64 / 71 variables
            Partitioned 5 equations into 1 components of size:
            [4]
            Elimination epoch completed... 0 eqns remain in ideal basis
            Partitioned 6 equations into 6 components of size:
            [1, 1, 1, 1, 1, 1]
            Computing appropriate NumberField...
            sage: f.certify_pentagons()  is None      # not tested (cypari issue in doctesting framework), long time (~1.5s)
            True
        """
        fvars_copy = deepcopy(self._fvars)
        self._fvars = {sextuple: float(rhs) for sextuple, rhs in self.get_fvars_in_alg_field().items()}
        if use_mp:
            pool = Pool()
        else:
            pool = None
        n_proc = pool._processes if pool is not None else 1
        params = [(child_id, n_proc, verbose) for child_id in range(n_proc)]
        pe = self._map_triv_reduce('pent_verify', params, worker_pool=pool, chunksize=1, mp_thresh=0)
        if np.all(np.isclose(np.array(pe), 0, atol=1e-7)):
            if verbose:
                print("Found valid F-symbols for {}".format(self._FR))
            pe = None
        else:
            if verbose:
                print("Something went wrong. Pentagons remain.")
        self._fvars = fvars_copy
        return pe
