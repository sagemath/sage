"""
Neighbors
"""

from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from copy import deepcopy
from sage.matrix.constructor import matrix


# ############################################################################
# Routines used for understanding p-neighbors and computing classes in a genus
# ############################################################################


def find_primitive_p_divisible_vector__random(self, p):
    """
    Find a random `p`-primitive vector in `L/pL` whose value is `p`-divisible.

    .. NOTE::

        Since there are about `p^{(n-2)}` of these lines, we have a `1/p`
        chance of randomly finding an appropriate vector.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [10,1,4])
        sage: v = Q.find_primitive_p_divisible_vector__random(5)
        sage: tuple(v) in ((1, 0), (1, 1), (2, 0), (2, 2), (3, 0), (3, 3), (4, 0), (4, 4))
        True
        sage: 5.divides(Q(v))
        True
        sage: Q = QuadraticForm(QQ, matrix.diagonal([1,1,1,1]))
        sage: v = Q.find_primitive_p_divisible_vector__random(2)
        sage: Q(v)
        2
    """
    n = self.dim()
    v = vector([ZZ.random_element(p) for _ in range(n)])

    # Repeatedly choose random vectors, and evaluate until the value
    # is p-divisible.
    k = 0
    while k < 1000:
        k = k + 1
        a = self(v)
        if a in ZZ and (a % p == 0) and (v != 0):
            return v
        else:
            v[ZZ.random_element(n)] = ZZ.random_element(p)
            # Replace a random entry and try again.
    raise RuntimeError("unable to find a p divisible vector")


def find_primitive_p_divisible_vector__next(self, p, v=None):
    """
    Find the next `p`-primitive vector (up to scaling) in `L/pL` whose
    value is `p`-divisible, where the last vector returned was `v`.  For
    an initial call, no `v` needs to be passed.

    Return vectors whose last nonzero entry is normalized to 0 or 1 (so no
    lines are counted repeatedly).  The ordering is by increasing the
    first non-normalized entry.  If we have tested all (lines of)
    vectors, then return None.

    OUTPUT: vector or None

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [10,1,4])
        sage: v = Q.find_primitive_p_divisible_vector__next(5); v
        (1, 1)
        sage: v = Q.find_primitive_p_divisible_vector__next(5, v); v
        (1, 0)
        sage: v = Q.find_primitive_p_divisible_vector__next(5, v); v
        sage: v = Q.find_primitive_p_divisible_vector__next(2) ; v
        (0, 1)
        sage: v = Q.find_primitive_p_divisible_vector__next(2, v) ; v
        (1, 0)
        sage: Q = QuadraticForm(QQ, matrix.diagonal([1,1,1,1]))
        sage: v = Q.find_primitive_p_divisible_vector__next(2)
        sage: Q(v)
        2
    """
    # Initialize
    n = self.dim()
    if v is None:
        w = vector(ZZ, [0] * (n - 1) + [1])
        a = self(w)
        if a in ZZ and (a % p == 0):
            return w
    else:
        w = deepcopy(v)

    # Handle n = 1 separately.
    if n <= 1:
        raise NotImplementedError("Sorry -- Not implemented yet!")

    # Look for the last nonzero entry (which must be 1)
    nz = n - 1
    while w[nz] == 0:
        nz += -1

    # Test that the last nonzero entry is 1 (to detect tampering).
    if w[nz] != 1:
        print("Warning: The input vector to QuadraticForm.find_primitive_p_divisible_vector__next() is not normalized properly.")

    # Look for the next vector, until w == 0
    while True:

        # Look for the first non-maximal (non-normalized) entry
        ind = 0
        while (ind < nz) and (w[ind] == p - 1):
            ind += 1

        # Increment
        if ind < nz:
            w[ind] += 1
            for j in range(ind):
                w[j] = 0
        else:
            for j in range(ind + 1):    # Clear all entries
                w[j] = 0

            if nz != 0:
                # Move the nonzero normalized index over by one, or
                # return the zero vector
                w[nz - 1] = 1
                nz += -1

        # Test for zero vector
        if w == 0:
            return None

        # Test for p-divisibility
        a = self(w)
        if a in ZZ and (a % p == 0):
            return w


def find_p_neighbor_from_vec(self, p, y, return_matrix=False):
    r"""
    Return the `p`-neighbor of ``self`` defined by ``y``.

    Let `(L,q)` be a lattice with `b(L,L) \subseteq \ZZ` which is maximal at `p`.
    Let `y \in L` with `b(y,y) \in p^2\ZZ` then the `p`-neighbor of
    `L` at `y` is given by
    `\ZZ y/p + L_y` where `L_y = \{x \in L | b(x,y) \in p \ZZ \}`
    and `b(x,y) = q(x+y)-q(x)-q(y)` is the bilinear form associated to `q`.

    INPUT:

    - ``p`` -- a prime number
    - ``y`` -- a vector with `q(y) \in p \ZZ`
    - ``odd`` -- boolean (default: ``False``); if `p=2`, return also odd neighbors
    - ``return_matrix`` -- boolean (default: ``False``); return
      the transformation matrix instead of the quadratic form

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: v = vector([0,2,1,1])
        sage: X = Q.find_p_neighbor_from_vec(3, v); X
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 0 0 0 ]
        [ * 1 4 4 ]
        [ * * 5 12 ]
        [ * * * 9 ]
        sage: B = Q.find_p_neighbor_from_vec(3, v, return_matrix=True)
        sage: Q(B) == X
        True

    Since the base ring and the domain are not yet separate,
    for rational, half integral forms we just pretend
    the base ring is `\ZZ`::

        sage: # needs sage.libs.pari
        sage: Q = QuadraticForm(QQ, matrix.diagonal([1,1,1,1]))
        sage: v = vector([1,1,1,1])
        sage: Q.find_p_neighbor_from_vec(2, v)
        Quadratic form in 4 variables over Rational Field with coefficients:
        [ 1/2 1 1 1 ]
        [ * 1 1 2 ]
        [ * * 1 2 ]
        [ * * * 2 ]
    """
    p = ZZ(p)
    if not p.divides(self(y)):
        raise ValueError(f"y={y} must be of square divisible by p={p}")
    if self.base_ring() not in [ZZ, QQ]:
        raise NotImplementedError("the base ring of this form must be the integers or the rationals")

    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

    n = self.dim()
    G = self.Hessian_matrix()
    R = self.base_ring()
    odd = False
    if R is QQ:
        odd = True
        if G.denominator() != 1:
            raise ValueError("the associated bilinear form q(x+y)-q(x)-q(y) must be integral.")
    b = y * G * y
    if not b % p == 0:
        raise ValueError("y^2 must be divisible by p=%s" % p)
    y_dual = y * G
    if p != 2 and b % p**2 != 0:
        for k in range(n):
            if y_dual[k] % p != 0:
                z = (ZZ**n).gen(k)
                break
        else:
            raise ValueError("either y is not primitive or self is not maximal at %s" % p)
        z *= (2 * y * G * z).inverse_mod(p)
        y = y - b * z
        # assert y*G*y % p^2 == 0
    if p == 2:
        val = b.valuation(p)
        if val <= 1:
            raise ValueError("y=%s must be of square divisible by 2" % y)
        if val == 2 and not odd:
            # modify it to have square 4
            for k in range(n):
                if y_dual[k] % p != 0:
                    z = (ZZ**n).gen(k)
                    break
            else:
                raise ValueError("either y is not primitive or self is not even, maximal at 2")
            y += 2 * z
            # assert y*G*y % 8 == 0

    y_dual = G * y
    Ly = y_dual.change_ring(GF(p)).column().kernel().matrix().lift()
    B = Ly.stack(p * matrix.identity(n))
    # the rows of B now generate L_y = { x in L | (x,y)=0 mod p}
    B = y.row().stack(p * B)
    B = B.hermite_form()[:n, :] / p
    # the rows of B generate ZZ * y/p + L_y
    # by definition this is the p-neighbor of L at y
    # assert B.det().abs() == 1

    if return_matrix:
        return B.T

    QF = self.parent()
    Gnew = (B * G * B.T).change_ring(R)
    return QF(Gnew)


def neighbor_iteration(seeds, p, mass=None, max_classes=None,
                       algorithm=None, max_neighbors=1000, verbose=False):
    r"""
    Return all classes in the `p`-neighbor graph of ``self``.

    Starting from the given seeds, this function successively
    finds `p`-neighbors until no new quadratic form (class) is obtained.

    INPUT:

    - ``seeds`` -- list of quadratic forms in the same genus

    - ``p`` -- a prime number

    - ``mass`` -- (optional) a rational number; the mass of this genus

    - ``max_classes`` -- (default: ``1000``) break the computation when ``max_classes`` are found

    - ``algorithm`` -- (optional) one of ``'orbits'``, ``'random'``, ``'exhaustion'``

    - ``max_random_trys`` -- (default: ``1000``) the maximum number of neighbors
      computed for a single lattice

    OUTPUT: list of quadratic forms

    EXAMPLES::

        sage: from sage.quadratic_forms.quadratic_form__neighbors import neighbor_iteration
        sage: Q = QuadraticForm(ZZ, 3, [1, 0, 0, 2, 1, 3])
        sage: Q.det()
        46

        sage: # needs sage.symbolic
        sage: mass = Q.conway_mass()
        sage: g1 = neighbor_iteration([Q], 3,   # long time
        ....:                         mass=mass, algorithm='random')
        sage: g2 = neighbor_iteration([Q], 3, algorithm='exhaustion')   # long time
        sage: g3 = neighbor_iteration([Q], 3, algorithm='orbits')                       # needs sage.libs.gap
        sage: mass == sum(1/q.number_of_automorphisms() for q in g1)    # long time
        True
        sage: mass == sum(1/q.number_of_automorphisms() for q in g2)    # long time
        True
        sage: mass == sum(1/q.number_of_automorphisms() for q in g3)                    # needs sage.libs.gap
        True

    TESTS::

        sage: from sage.quadratic_forms.quadratic_form__neighbors import neighbor_iteration
        sage: Q = QuadraticForm(ZZ, 3, [1, 0, 0, 2, 1, 3])
        sage: g = neighbor_iteration([Q], 3, mass=Q.conway_mass(), max_classes=2)       # needs sage.symbolic
        ...
        UserWarning: reached the maximum number of isometry classes=2.
        Increase the optional argument max_classes to obtain more.
        Warning: not all classes in the genus were found
        sage: neighbor_iteration([Q], 3,                                                # needs sage.symbolic
        ....:                    mass=Q.conway_mass(), max_neighbors=0, algorithm='random')
        Warning: not all classes in the genus were found
        []
    """
    from sage.quadratic_forms.quadratic_form import QuadraticForm
    from warnings import warn
    p = ZZ(p)
    if max_classes is None:
        max_classes = 1000
    if not all(isinstance(s, QuadraticForm) for s in seeds):
        raise ValueError("seeds must be a list of quadratic forms")
    if algorithm is None:
        n = seeds[0].dim()
        if p**n > ZZ(2)**18:
            # too many lines to compute the orbits fast
            algorithm = 'random'
        else:
            algorithm = 'orbits'

    if algorithm == 'orbits':
        def p_divisible_vectors(Q, max_neighbors):
            yield from iter(v.lift() for v in Q.orbits_lines_mod_p(p)
                            if v != 0 and Q(v.lift()).valuation(p) > 0)

    elif algorithm == 'exhaustion':
        def p_divisible_vectors(Q, max_neighbors):
            k = 0
            v = Q.find_primitive_p_divisible_vector__next(p)
            while k < max_neighbors:
                k = k + 1
                v = Q.find_primitive_p_divisible_vector__next(p, v)
                if v is not None:
                    yield v
    elif algorithm == 'random':
        def p_divisible_vectors(Q, max_neighbors):
            k = 0
            while k < max_neighbors:
                k += 1
                v = Q.find_primitive_p_divisible_vector__random(p)
                yield v
    else:
        raise ValueError("unknown algorithm")
    waiting_list = list(seeds)
    isom_classes = []
    mass_count = QQ.zero()
    n_isom_classes = ZZ.zero()
    while waiting_list and mass != mass_count and n_isom_classes < max_classes:
        # find all p-neighbors of Q
        Q = waiting_list.pop()
        for v in p_divisible_vectors(Q, max_neighbors):
            Q_neighbor = Q.find_p_neighbor_from_vec(p, v)
            if not any(Q_neighbor.is_globally_equivalent_to(S) for S in isom_classes):
                Q_neighbor = Q_neighbor.lll()
                isom_classes.append(Q_neighbor)
                waiting_list.append(Q_neighbor)
                n_isom_classes += 1
                mass_count += Q_neighbor.number_of_automorphisms()**(-1)
                if verbose:
                    print(max_neighbors)
                    print(len(waiting_list))
                if mass_count == mass or n_isom_classes >= max_classes:
                    break

    if len(isom_classes) >= max_classes:
        warn("reached the maximum number of isometry classes=%s. Increase the optional argument max_classes to obtain more." % max_classes)

    if mass is not None:
        assert mass_count <= mass
        if mass_count < mass:
            print("Warning: not all classes in the genus were found")
    return isom_classes


def orbits_lines_mod_p(self, p):
    r"""
    Let `(L, q)` be a lattice. This returns representatives of the
    orbits of lines in `L/pL` under the orthogonal group of `q`.

    INPUT:

    - ``p`` -- a prime number

    OUTPUT: list of vectors over ``GF(p)``

    EXAMPLES::

        sage: from sage.quadratic_forms.quadratic_form__neighbors import orbits_lines_mod_p
        sage: Q = QuadraticForm(ZZ, 3, [1, 0, 0, 2, 1, 3])
        sage: Q.orbits_lines_mod_p(2)                                                   # needs sage.libs.gap sage.libs.pari
        [(0, 0, 1),
         (0, 1, 0),
         (0, 1, 1),
         (1, 0, 0),
         (1, 0, 1),
         (1, 1, 0),
         (1, 1, 1)]
    """
    from sage.libs.gap.libgap import libgap
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

    # careful the self.automorphism_group() acts from the left
    # but in gap we act from the right!! --> transpose
    gens = self.automorphism_group().gens()
    gens = [g.matrix().transpose().change_ring(GF(p)) for g in gens]
    orbs = libgap.function_factory("""function(gens, p)
        local one, G, reps, V, n, orb;
        one:= One(GF(p));
        G:=Group(List(gens, g -> g*one));
        n:= Size(gens[1]);
        V:= GF(p)^n;
        orb:= OrbitsDomain(G, V, OnLines);
        reps:= List(orb, g->g[1]);
        return reps;
        end;""")
    orbs_reps = orbs(gens, p)
    M = GF(p)**self.dim()
    return [M(m.sage()) for m in orbs_reps if not m.IsZero()]
