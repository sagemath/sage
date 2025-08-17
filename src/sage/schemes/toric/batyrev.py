r"""
Smooth Fano toric varieties.

This module provides support for smooth Fano toric varieties, which have been
classified by Batyrev in dimensions up to 4.
"""
from collections import defaultdict, Counter
from sage.geometry.cone import Cone
from sage.geometry.fan import Fan, FaceFan
from sage.geometry.lattice_polytope import LatticePolytope
from sage.geometry.toric_lattice import ToricLattice
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.arith.misc import GCD as gcd
from sage.schemes.toric.variety import (
    DEFAULT_PREFIX,
    ToricVariety_field,
    normalize_names,
)
from sage.schemes.toric.fano_variety import FanoToricVariety_field
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.categories.fields import Fields

_Fields = Fields()


from .poly_db_3d import BATYREV_3FOLD_LOOKUP, CLASS_BY_LABEL_3FOLD    
from .poly_db_4d import BATYREV_4FOLD_LOOKUP, CLASS_BY_LABEL_4FOLD


# TODO: add the method to fix a basis of the Picard group, and express divisors in this basis
class SmoothFanoToricVariety_field(FanoToricVariety_field, metaclass=ClasscallMetaclass):
    r"""
    Construct a smooth Fano toric variety from a smooth reflexive polytope
    **or directly from a smooth complete fan**.

    INPUT:

    - ``Delta`` -- reflexive :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The fan of the
      constructed CPR-Fano toric variety will be a crepant subdivision of the
      *normal fan* of ``Delta``. Either ``Delta`` or ``Delta_polar`` must be
      given, but not both at the same time, since one is completely determined
      by another via :meth:`polar
      <sage.geometry.lattice_polytope.LatticePolytopeClass.polar>` method.

    - ``Delta_polar`` -- reflexive :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The fan of the
      constructed CPR-Fano toric variety will be a crepant subdivision of the
      *face fan* of ``Delta_polar``. Either ``Delta`` or ``Delta_polar`` must
      be given, but not both at the same time, since one is completely
      determined by another via :meth:`polar
      <sage.geometry.lattice_polytope.LatticePolytopeClass.polar>` method.

    - ``fan`` -- a **smooth complete** :class:`~sage.geometry.fan.Fan`. When
      provided, the constructor first computes the **anticanonical polytope**
      :math:`\Delta = \{ m \in M_\mathbb{R} : \langle m, v_\rho\rangle \ge -1\}`
      in :math:`M=N^\vee`, checks boundedness (Fano) and reflexivity (Gorenstein),
      and then sets ``Delta_polar := \Delta^\circ \subset N_\mathbb{R}``.

    - ``coordinate_names`` -- names of variables for the coordinate ring, see
      :func:`~sage.schemes.toric.variety.normalize_names`
      for acceptable formats. If not given, indexed variable names will be
      created automatically.

    - ``names`` -- an alias of ``coordinate_names`` for internal
      use. You may specify either ``names`` or ``coordinate_names``,
      but not both.

    - ``coordinate_name_indices`` -- list of integers, indices for indexed
      variables. If not given, the index of each variable will coincide with
      the index of the corresponding point of ``Delta_polar``.

    - ``base_ring`` / ``base_field`` -- base field of the Fano toric variety
      (default: `\QQ`). If both are given, ``base_field`` takes precedence.

    - ``check`` -- if True (default), validate inputs (e.g. smoothness of the
      face fan of ``Delta_polar``).

    OUTPUT: :class:`smooth Fano toric variety <SmoothFanoToricVariety_field>`

    Examples
    --------
    We start with the product of two projective lines, which is smooth::

        sage: from sage.geometry import lattice_polytope
        sage: P1xP1 = SmoothFanoToricVariety(Delta_polar=lattice_polytope.cross_polytope(2))
        sage: P1xP1
        2-d smooth Fano toric variety covered by 4 affine patches
        sage: P1xP1.is_smooth()
        True

    Construct :math:`\mathbb{P}^2` **from its fan** and let the constructor
    compute the anticanonical polytope automatically::

        sage: N = ToricLattice(2)
        sage: e1, e2 = N.basis()
        sage: rays = [e1, e2, -(e1+e2)]
        sage: cones = [(0,1), (1,2), (2,0)]
        sage: fan_P2 = Fan(cones, rays=rays, check=True)
        sage: X = SmoothFanoToricVariety(fan=fan_P2)
        sage: X.dimension_relative(), X.picard_rank(), X.fan().nrays()
        (2, 1, 3)
        sage: X.is_smooth()
        True
    """

    @staticmethod
    def __classcall_private__(
        cls,
        fan=None,
        Delta=None,
        Delta_polar=None,
        coordinate_names=None,
        names=None,
        coordinate_name_indices=None,
        base_ring=None,
        base_field=None,
        check=True,
    ):
        if names is not None:
            if coordinate_names is not None:
                raise ValueError("You must not specify both coordinate_names and names!")
            coordinate_names = names

        # AL: allow the user to specify a smooth complete fan directly
        if fan is not None:
            if not fan.is_complete():
                raise ValueError("fan must be complete!")
            if not fan.is_smooth():
                raise ValueError("fan must be smooth!")

            Delta, Delta_polar = cls.anticanonical_polytope_from_fan(fan, require_reflexive=True)

            # Normalize / choose base field
            if base_field is not None:
                base_ring = base_field
            if base_ring is None:
                base_ring = QQ
            elif base_ring not in _Fields:
                raise TypeError("need a field to construct a Fano toric variety!\n Got %s" % base_ring)

            return typecall(
                cls,
                Delta_polar=Delta_polar,   # in N_R
                fan=fan,
                coordinate_names=coordinate_names,
                coordinate_name_indices=coordinate_name_indices,
                base_field=base_ring,
            )

        if (Delta is None) == (Delta_polar is None):
            raise ValueError("exactly one of Delta and Delta_polar must be given")

        if Delta_polar is None:
            Delta_polar = Delta.polar()

        if check and not Delta_polar.is_reflexive():
            raise ValueError("Delta_polar must be reflexive")

        fan = FaceFan(Delta_polar)

        if check and not fan.is_smooth():
            raise ValueError("The face fan of Delta_polar is not smooth")

        # Check/normalize base_field
        if base_field is not None:
            base_ring = base_field
        if base_ring is None:
            base_ring = QQ
        elif base_ring not in _Fields:
            raise TypeError("need a field to construct a Fano toric variety!\n Got %s" % base_ring)

        return typecall(
            cls,
            Delta_polar=Delta_polar,
            fan=fan,
            coordinate_names=coordinate_names,
            coordinate_name_indices=coordinate_name_indices,
            base_field=base_ring,
        )

    def __init__(
        self, Delta_polar, fan, coordinate_names, coordinate_name_indices, base_field
    ):
        r"""
        See :class:`SmoothFanoToricVariety_field` for documentation.

        Use ``SmoothFanoToricVariety`` to construct a smooth Fano toric variety.
        """
        FanoToricVariety_field.__init__(self,
            Delta_polar, fan, coordinate_names, coordinate_name_indices, base_field
        )
        self.Sigma = self.Delta_polar.face_fan()
        self.n = self.dimension()
        self.V, self.max_cones = self._matrix_from_rays_and_cones(self.Sigma)
        self.Q = self._gale_dual(self.V)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: from sage.geometry import lattice_polytope
            sage: P1xP1 = SmoothFanoToricVariety(Delta_polar=lattice_polytope.cross_polytope(2))
            sage: print(P1xP1._repr_())
            2-d smooth Fano toric variety covered by 4 affine patches
        """
        return "%d-d smooth Fano toric variety covered by %d affine patches" % (
            self.dimension_relative(),
            self.fan().ngenerating_cones(),
        )
    
    def _cohomology_dim(self, divisor, degree):
        r"""
        Return ``dim H^{degree}(X, O_X(divisor))`` as an ``int``.
        
        INPUT:
        - ``divisor`` -- a divisor on ``X``, given as a list of integers, the coefficients of the torus-invariant divisors in terms of the rays of the fan. Here we use the divisor to represent the line bundle O_X(divisor).
        """
        L = self.divisor(divisor)
        H = L.cohomology(dim=True)
        return H[degree]

    def _h0_dim(self, divisor):
        r"""
        Return ``dim H^0(X, O_X(divisor))`` without computing higher cohomology.

        INPUT:
        - ``divisor`` -- a divisor on ``X``, given as a list of integers, the coefficients
          of the torus-invariant divisors in the order of the fan’s rays.

        METHOD:
        Build the polyhedron in ``M_R`` given by inequalities
        ``<m, v_rho> + a_rho >= 0`` and count its lattice points.
        Falls back to ``cohomology(dim=True)[0]`` if needed.
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        from sage.rings.rational_field import QQ

        # fast path: build P_D = { m : <m, v_rho> >= -a_rho } in M_R
        try:
            a = list(divisor)                               # coefficients a_rho
            N = self.fan().lattice()
            M = N.dual()
            rays = [N(v) for v in self.fan().rays()]        # primitive generators v_rho

            # polyhedron expects rows [const, coeffs...], representing const + <coeffs, m> >= 0
            ieqs = [[QQ(ai)] + [QQ(c) for c in v] for ai, v in zip(a, rays)]
            P = Polyhedron(ieqs=ieqs)                        # lives in M_R

            # robust boundedness check (avoid backend .is_bounded quirks):
            if P.n_rays() != 0 or P.n_lines() != 0:
                raise RuntimeError("Unbounded section polyhedron; falling back to cohomology().")

            # count lattice points. Prefer the fast counter if available.
            try:
                return int(P.integral_points_count())
            except Exception:
                # fallback: explicit enumeration (can be slower on large polytopes).
                return len(list(P.integral_points()))

        except Exception:
            # safe fallback: use the existing general routine.
            return self._cohomology_dim(divisor, 0)
        
    def _hom_dim(self, L1, L2):
      r"""Return ``dim Hom(L1, L2) = dim H^0(X, O_X(L2-L1))`` via a cheap H^0 counter."""
      return self._h0_dim(L2 - L1)
    
    def picard_rank(self) -> int:
        r"""
        Return the Picard rank of ``self``.
        TODO: 
        - we can implement a more general method for any toric variety. The current
        implementation only works for smooth complete toric varieties.
        - add testcases
        """
        ambient_dim = self.fan().lattice().ngens()
        from sage.modules.free_module import VectorSpace
        ambient_space = VectorSpace(QQ, ambient_dim)
        dim_cone = ambient_space.span(self.fan().rays()).dimension()
        return self.fan().nrays() - dim_cone
    
    def table_invariants(self):
      n = self.n
      X = ToricVariety(self.Sigma)

      # c1^n and (if n>=4) c1^2 c2 via Chow ring
      try:
          c1n, c1_2_c2 = self._chern_numbers_toric(X)
      except Exception:
          # fallback for c1^n only using polar volume (requires Fano)
          Delta_dual = self.Delta.polar()
          c1n = self._degree_c1_power_n(Delta_dual)
          c1_2_c2 = None

      # b2 (=\rho) and b4 from Chow ring when available
      try:
          b2, b4 = self._betti_numbers_toric(X)
      except Exception:
          b2 = self.V.ncols() - n
          b4 = None

      h0 = self._h0_of_anticanonical(self.Delta)
      try:
          aV = self._automorphism_dimension(self.Sigma)
      except Exception:
          aV = n  # conservative lower bound

      if n == 3:
          return (3, int(c1n), int(b2), int(h0), int(aV))
      elif n == 4:
          if c1_2_c2 is None or b4 is None:
              raise RuntimeError("Need Chow ring support to compute c1^2 c2 and b4 in dimension 4.")
          return (4, int(c1n), int(c1_2_c2), int(b2), int(b4), int(h0), int(aV))
      else:
          raise ValueError("This method standardizes invariants for 3- and 4-folds.")

    def identify_batyrev_label(self, lookup=BATYREV_3FOLD_LOOKUP):
        inv = self.table_invariants()
        return [lab for lab, tpl in lookup.items() if tuple(tpl) == tuple(inv)]

    def classify(self, lookup=BATYREV_3FOLD_LOOKUP):
        """
        Combine invariants with structure tests to pick a unique label
        among candidates sharing the same invariants.
        """
        inv = self.table_invariants()
        cands = self.identify_batyrev_label(lookup)
        prod, factors = self.is_cartesian_product_of_projective_spaces()
        pb, k = self.is_projective_bundle_over_projective_space()
        blow, center_dim = self.is_blow_up_of_projective_space()

        # heuristic disambiguation if multiple candidates share invariants
        label = None
        if len(cands) == 1:
            label = cands[0]
        elif cands:
            # filter by coarse class
            def label_class(L):
                return CLASS_BY_LABEL.get(L, None)
            if prod:
                pool = [L for L in cands if label_class(L) == "product"]
                label = pool[0] if pool else None
            elif pb:
                pool = [L for L in cands if label_class(L) == "projective_bundle"]
                label = pool[0] if pool else None
            elif blow:
                pool = [L for L in cands if label_class(L) == "blow_up"]
                label = pool[0] if pool else None

        return {
            "invariants": inv,
            "candidates": cands,
            "label": label,
            "structure": {
                "is_product": bool(prod), "factors": factors,
                "is_projective_bundle_over_P": bool(pb), "fiber_dimension": k,
                "is_blow_up_of_P": bool(blow), "center_dim": center_dim,
            },
        }

    def is_cartesian_product_of_projective_spaces(self):
        """
        Recognize X ≅ Π P^{a_i} from Gale-dual blocks and cone combinatorics.
        Returns (bool, [a1,a2,...]) where ai are the factor dimensions if True.
        """
        n, m = self.V.nrows(), self.V.ncols()
        Q = self.Q

        if Q.nrows() == 0:
            if self.is_fan_of_projective_space(self.V, self.max_cones):
                return True, [n]
            return False, None

        reps, col_class = [], [-1]*m
        for j in range(m):
            qj = vector(QQ, Q.column(j))
            assigned = False
            for k, rk in enumerate(reps):
                if qj in Span([rk]):
                    col_class[j] = k; assigned = True; break
            if not assigned:
                reps.append(qj); col_class[j] = len(reps)-1

        blocks = defaultdict(list)
        for j, cls in enumerate(col_class):
            blocks[cls].append(j)
        blocks = list(blocks.values())
        if any(len(B) < 2 for B in blocks):
            return False, None

        for c in self.max_cones:
            cnt = Counter(col_class[j] for j in c)
            for b, B in enumerate(blocks):
                if cnt[b] != len(B) - 1:
                    return False, None

        dims = [len(B)-1 for B in blocks]
        if sum(dims) != n:
            return False, None
        return True, dims

    def is_projective_bundle_over_projective_space(self):
        """
        Detect X ≅ P^k-bundle over P^{n-k}.
        Returns (bool, k) with k ≥ 1 if True.
        """
        prod, _ = self.is_cartesian_product_of_projective_spaces()
        if prod:
            return False, None

        n, m = self.V.nrows(), self.V.ncols()
        Q = self.Q

        reps, col_class = [], [-1]*m
        for j in range(m):
            qj = vector(QQ, Q.column(j))
            assigned = False
            for k, rk in enumerate(reps):
                if qj in Span([rk]):
                    col_class[j] = k; assigned = True; break
            if not assigned:
                reps.append(qj); col_class[j] = len(reps)-1

        blocks = defaultdict(list)
        for j, cls in enumerate(col_class):
            blocks[cls].append(j)
        blocks = list(blocks.values())

        for S in blocks:
            k = len(S) - 1
            if k < 1 or k >= n:
                continue
            B = [j for j in range(m) if j not in S]
            if len(B) != (n-k)+1:
                continue
            Sset, Bset = set(S), set(B)
            ok = True
            for c in self.max_cones:
                s_in = len([j for j in c if j in Sset])
                b_in = len([j for j in c if j in Bset])
                if not (s_in == len(S)-1 and b_in == len(B)-1):
                    ok = False; break
            if not ok:
                continue
            # Base must be P^{n-k}
            Vb = self.V[:, B]
            base_max = list({tuple(sorted([j for j in c if j in Bset])) for c in self.max_cones})
            if self._is_fan_of_projective_space(Vb, base_max):
                return True, k
        return False, None

    def is_blow_up_of_projective_space(self):
        """
        Detect a single toric blow-up of P^n.
        Returns (bool, center_dim) where center_dim in {0,1} for n=3 (point/line).
        """
        n, m = self.V.nrows(), self.V.ncols()
        for r in range(m):
            keep = [j for j in range(m) if j != r]
            if len(keep) != n+1:
                continue
            base_cones = sorted({tuple(sorted([j for j in c if j in keep]))
                                for c in self.max_cones if r not in c})
            if len(base_cones) != n+1:
                continue
            Vb = self.V[:, keep]
            if self._is_fan_of_projective_space(Vb, base_cones):
                center_dim = None
                if n == 3:
                    used = [c for c in self.max_cones if r in c]
                    partners = set()
                    for c in used:
                        for j in c:
                            if j != r: partners.add(j)
                    # 3 partners -> point, 2 partners -> line
                    center_dim = 0 if len(partners) == 3 else 1
                return True, center_dim
        return False, None

    @staticmethod
    def _matrix_from_rays_and_cones(fan):
        rays = [vector(ZZ, r.primitive_generator()) for r in fan.rays()]
        V = matrix(ZZ, list(zip(*rays)))  # n x m
        r2i = {tuple(r): i for i, r in enumerate(rays)}
        max_cones = []
        for s in fan.maximal_cones():
            cols = [r2i[tuple(vector(ZZ, rr))] for rr in s.rays()]
            max_cones.append(tuple(cols))
        return V, max_cones

    @staticmethod
    def _gale_dual(V):
        U, D, W = matrix(ZZ, V.transpose()).smith_form()
        r = V.rank()
        Q = W[r:, :]
        if Q.nrows() == 0:
            return matrix(ZZ, 0, V.ncols())
        for i in range(Q.nrows()):
            g = gcd(list(Q.row(i)))
            if g > 1:
                Q.set_row(i, Q.row(i) // g)
        return Q

    @staticmethod
    def _is_fan_of_projective_space(V, max_cones):
        n, m = V.nrows(), V.ncols()
        if m != n+1:
            return False
        if not (all(len(c) == n for c in max_cones) and len(max_cones) == m):
            return False
        dets = {abs(V[:, list(J)].determinant()) for J in combinations(range(m), n)}
        if dets != {1}:
            return False
        s = sum(V.columns())
        return all(x == 0 for x in s)

    @staticmethod
    def _h0_of_anticanonical(Delta):
        return Delta.integral_points_count()

    @staticmethod
    def _degree_c1_power_n(Delta_dual):
        n = Delta_dual.ambient_dim()
        vol = Delta_dual.volume(measure='induced_lattice')
        return ZZ(factorial(n) * vol)

    @staticmethod
    def _betti_numbers_toric(X):
        A, _ = X.chow_ring()
        n = X.dimension()
        def rank_deg(k): return len(A.homogeneous_component(k).basis())
        b2 = rank_deg(1)
        b4 = rank_deg(2) if n >= 4 else None
        return b2, b4

    @staticmethod
    def _chern_numbers_toric(X):
        A, divs = X.chow_ring()
        one = A(1)
        c_tot = one
        for D in divs.values():
            c_tot *= (one + D)
        n = X.dimension()
        c = [A(0)]*(n+1)
        for k in range(n+1):
            c[k] = c_tot.homogeneous_component(k)
        c1, c2 = c[1], (c[2] if n >= 2 else None)
        deg = A.degree()
        c1n = deg(c1**n)
        c1_2_c2 = deg((c1**2)*c2) if n >= 4 else None
        return ZZ(c1n), (ZZ(c1_2_c2) if c1_2_c2 is not None else None)

    @staticmethod
    def _automorphism_dimension(fan):
        # a(V) = n + #Demazure roots
        V, _ = SmoothToricFanoVariety._matrix_from_rays_and_cones(fan)
        n, m = V.nrows(), V.ncols()
        rays = [V.column(i) for i in range(m)]
        roots = set()
        for j, vj in enumerate(rays):
            A, b = [], []
            for i, vi in enumerate(rays):
                if i == j:
                    A.append(list(vi)); b.append(-1)     # <m, vj> = -1
                else:
                    A.append([-x for x in vi]); b.append(0)  # <m, vi> ≥ 0
            P = Polyhedron(ieqs=[[bi]+Ai for Ai,bi in zip(A,b)])
            for w in P.integral_points():
                roots.add(tuple(w))
        return n + len(roots)
    
    def _is_cartesian_product_of_projective_spaces(self):
        r"""
        Check if the toric variety is a cartesian product of projective spaces. This will determine the method to compute the cohomologically zero line bundles.
        """
        # TODO: if no neat algorithm related to the reflexive polytope, we can just enrich the library with its original name and geometric properties. 
        if self.picard_rank != self.toric_variety.dimension():
            return False
        # Quick ray count check: P^m has m+1 rays.  
        # product of r copies ⇒ total rays = Σ(m_i+1). We only test consistency here, not sufficiency.
        rays = len(self.toric_variety.fan().rays())
        return rays == self.toric_variety.dimension() + self.picard_rank

        
    def _is_projective_bundle_over_projective_space(self):
        # AL: do we also want to return the splitting type of the bundle?
        if self.picard_rank == 2:
            # all smooth toric varieties of Picard rank 2 are projective bundles over projective spaces
            return True
        # TODO: investigate so called Cayley test

    
    def _is_blow_up_of_projective_space(self):
        # Rough proxy: smooth Fano with Picard rank > 1 but small.
        return 1 < self.picard_rank() <= self.dimension() + 1


    @staticmethod
    # TODO: something wrong with this method for some testcase
    def anticanonical_polytope_from_fan(fan, require_reflexive=True):
        r"""
        Construct the anticanonical polytope from a complete fan.

        Given a complete fan in a lattice :math:`N`, define
        \[
            \Delta \;=\; \{\, m \in M_\mathbb{R} \mid \langle m, v_\rho\rangle \ge -1
            \ \text{for each primitive ray generator } v_\rho \in N \,\}
            \subset M_\mathbb{R},\quad M = N^\vee.
        \]
        Then:
          * :math:`\Delta` is bounded iff :math:`-K_X` is ample (i.e. the variety is Fano).
          * If the variety is Gorenstein Fano (e.g. smooth Fano), then :math:`\Delta` is reflexive
            and its polar :math:`\Delta^\circ \subset N_\mathbb{R}` equals the convex hull of the
            primitive ray generators.

        Parameters
        ----------
        fan : :class:`~sage.geometry.fan.Fan`
            A complete fan in a lattice :math:`N`.
        require_reflexive : bool, default True
            If True, raise an error when :math:`\Delta` is not reflexive (i.e. the variety is not
            Gorenstein Fano).

        Returns
        -------
        (Delta, Delta_polar) : tuple
            ``Delta`` is a :class:`~sage.geometry.lattice_polytope.LatticePolytope` in :math:`M`.
            ``Delta_polar`` is its polar, a lattice polytope in :math:`N`.

        Examples
        --------
        A Fano–Gorenstein example: :math:`\mathbb{P}^2`. Its fan has rays
        :math:`(1,0)`, :math:`(0,1)`, :math:`(-1,-1)` and three 2D cones::

        
            sage: N = ToricLattice(2)
            sage: e1, e2 = N.basis()
            sage: rays = [e1, e2, -(e1+e2)]
            sage: cones = [(0,1), (1,2), (2,0)]
            sage: fan_P2 = Fan(cones, rays=rays, check=True)
            sage: Delta, Delta_polar = SmoothFanoToricVariety_field.anticanonical_polytope_from_fan(fan_P2, require_reflexive=True)
            sage: Delta.is_reflexive()
            True
            sage: len(list(Delta_polar.vertices()))
            3
            sage: set(map(tuple, Delta_polar.vertices())) == set(map(tuple, fan_P2.rays()))
            True

        A smooth but **non-Fano** example: the Hirzebruch surface :math:`\mathbb{F}_2`.
        The standard smooth complete fan uses rays :math:`(0,1)`, :math:`(1,0)`, :math:`(0,-1)`,
        :math:`(-1,2)` in cyclic order. Here :math:`-K` is not ample, so :math:`\Delta` is unbounded
        and the construction fails with a clear error::

            sage: N = ToricLattice(2)
            sage: e1, e2 = N.basis()                   # e1=(1,0), e2=(0,1)
            sage: rays = [e2, e1, -e2, -e1 + 2*e2]     # (0,1), (1,0), (0,-1), (-1,2)
            sage: cones = [(0,1), (1,2), (2,3), (3,0)]
            sage: fan_F2 = Fan(cones, rays=rays, check=True)
        """
        if not fan.is_complete():
            raise ValueError("fan must be complete.")

        N = fan.lattice()
        M = N.dual()

        # primitive ray generators in N
        rays_N = [N(v) for v in fan.rays()]

        # inequalities:  <m, v> + 1 >= 0   encoded as [const, coeffs...]
        ieqs = [[QQ(1)] + [QQ(c) for c in v] for v in rays_N]

        P = Polyhedron(ieqs=ieqs)  # in M_R
        # TODO: to check if the polyhedron is bounded
        vertices = []
        for v in P.vertices():
            vec = v.vector()
            if require_reflexive and any(c.denominator() != 1 for c in vec):
                print(vec)
                raise ValueError("The anticanonical polyhedron is not a lattice polytope: the variety is not Gorenstein Fano.")
            vertices.append(M(vec))
        Delta = LatticePolytope(vertices)
        if require_reflexive and not Delta.is_reflexive():
            raise ValueError("Delta is not reflexive: the fan does not define a Gorenstein Fano variety.")

        Delta_polar = Delta.polar()  # in N_R
        return Delta, Delta_polar


def SmoothFanoToricVariety(*arg, **kwds):
    """
    Convenience constructor.

    Examples
    --------
    Using a reflexive polar polytope::

        sage: from sage.geometry import lattice_polytope
        sage: X = SmoothFanoToricVariety(Delta_polar=lattice_polytope.cross_polytope(2))
        sage: X.dimension_relative(), X.fan().nrays()
        (2, 4)

    Using a smooth complete fan (the constructor checks Fano/Gorenstein via the anticanonical polytope)::

        sage: N = ToricLattice(2)
        sage: e1, e2 = N.basis()
        sage: rays = [e1, e2, -(e1+e2)]
        sage: cones = [(0,1), (1,2), (2,0)]
        sage: fan_P2 = Fan(cones, rays=rays, check=True)
        sage: X = SmoothFanoToricVariety(fan=fan_P2)
        sage: X.is_smooth(), X.picard_rank()
        (True, 1)
    """
    return SmoothFanoToricVariety_field(*arg, **kwds)
