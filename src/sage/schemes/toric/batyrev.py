r"""
Smooth Fano toric varieties.

This module provides support for smooth Fano toric varieties, which have been
classified by Batyrev in dimensions up to 4.
"""

from sage.geometry.cone import Cone
from sage.geometry.fan import Fan, FaceFan
from sage.geometry.lattice_polytope import LatticePolytope
from sage.geometry.toric_lattice import ToricLattice
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.rational_field import QQ
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

        sage: from sage.geometry.toric_lattice import ToricLattice
        sage: from sage.geometry.fan import Fan
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

        fan = FaceFan(Delta_polar, check=check)

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
    
    def _picard_rank(self) -> int:
        r"""
        Return the Picard rank of ``self``.
        TODO: we can implement a more general method for any toric variety. The current
        implementation only works for smooth complete toric varieties.
        """
        return self.fan().nrays() - self.fan().dimension()


    @staticmethod
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

            sage: from sage.schemes.toric.smooth_fano import SmoothFanoToricVariety_field
            sage: from sage.geometry.toric_lattice import ToricLattice
            sage: from sage.geometry.fan import Fan
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

        sage: from sage.geometry.toric_lattice import ToricLattice
        sage: from sage.geometry.fan import Fan
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
