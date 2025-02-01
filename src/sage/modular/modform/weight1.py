# sage.doctest: needs sage.libs.pari
r"""
Weight 1 modular forms

This module contains routines for computing weight 1 modular forms, using
George Schaeffer's "Hecke stability" algorithm (detailed in [Sch2015]_). These
functions are mostly for internal use; a more convenient interface is offered
by the usual ModularForms and CuspForms constructors.

AUTHORS:

- David Loeffler (2017-11): first version
"""

from sage.misc.cachefunc import cached_function
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.misc.verbose import verbose
from sage.structure.sequence import Sequence
from sage.modular.arithgroup.congroup_gamma0 import Gamma0_constructor as Gamma0
from sage.modular.arithgroup.congroup_gammaH import GammaH_constructor as GammaH
from sage.modular.arithgroup.arithgroup_generic import ArithmeticSubgroup


@cached_function
def modular_ratio_space(chi):
    r"""
    Compute the space of 'modular ratios', i.e. meromorphic modular forms f
    level N and character chi such that f * E is a holomorphic cusp form for
    every Eisenstein series E of weight 1 and character 1/chi.

    Elements are returned as `q`-expansions up to precision R, where R is one
    greater than the weight 3 Sturm bound.

    EXAMPLES::

        sage: chi = DirichletGroup(31,QQ).0
        sage: sage.modular.modform.weight1.modular_ratio_space(chi)
        [q - 8/3*q^3 + 13/9*q^4 + 43/27*q^5 - 620/81*q^6 + 1615/243*q^7 + 3481/729*q^8 + O(q^9),
         q^2 - 8/3*q^3 + 13/9*q^4 + 70/27*q^5 - 620/81*q^6 + 1858/243*q^7 + 2752/729*q^8 + O(q^9)]
    """
    from sage.modular.modform.constructor import EisensteinForms, CuspForms

    if chi(-1) == 1:
        return []
    N = chi.modulus()
    chi = chi.minimize_base_ring()
    K = chi.base_ring()
    R = Gamma0(N).sturm_bound(3) + 1
    verbose("Working precision is %s" % R, level=1)
    verbose("Coeff field is %s" % K, level=1)

    V = K**R
    I = V
    d = I.rank()

    t = verbose("Calculating Eisenstein forms in weight 1...",level=1)
    B0 = EisensteinForms(~chi, 1).q_echelon_basis(prec=R)
    B = [b + B0[0] for b in B0]
    verbose("Done (dimension %s)" % len(B),level=1,t=t)

    t = verbose("Calculating in weight 2...", level=1)
    C = CuspForms(Gamma0(N), 2).q_echelon_basis(prec=R)
    verbose("Done (dimension %s)" % len(C), t=t, level=1)

    t = verbose("Computing candidate space", level=1)
    for b in B:
        quots = (c/b for c in C)
        W = V.span(V(x.padded_list(R)) for x in quots)
        I = I.intersection(W)
        if I.rank() < d:
            verbose(" Cut down to dimension %s" % I.rank(), level=1)
            d = I.rank()
        if I.rank() == 0:
            break

    verbose("Done: intersection is %s-dimensional" % I.dimension(), t=t, level=1)

    A = PowerSeriesRing(K, 'q')
    return [A(x.list()).add_bigoh(R) for x in I.gens()]


def modular_ratio_to_prec(chi, qexp, prec):
    r"""
    Given a `q`-expansion of a modular ratio up to sufficient precision to
    determine it uniquely, compute it to greater precision.

    EXAMPLES::

        sage: from sage.modular.modform.weight1 import modular_ratio_to_prec
        sage: R.<q> = QQ[[]]
        sage: modular_ratio_to_prec(DirichletGroup(31,QQ).0, q-q^2-q^5-q^7+q^8+O(q^9), 20)
        q - q^2 - q^5 - q^7 + q^8 + q^9 + q^10 + q^14 - q^16 - q^18 - q^19 + O(q^20)
    """
    if prec <= qexp.prec():
        return qexp.add_bigoh(prec)
    from sage.modular.modform.constructor import EisensteinForms, CuspForms
    C = CuspForms(chi.level(), 2, base_ring=qexp.base_ring())
    B = EisensteinForms(~chi, 1).gen(0).qexp(prec)
    qexp = qexp.add_bigoh(C.sturm_bound())
    fB = qexp * B
    fB_elt = C(fB, check=False)
    return fB_elt.qexp(prec) / B


@cached_function
def hecke_stable_subspace(chi, aux_prime=ZZ(2)):
    r"""
    Compute a `q`-expansion basis for `S_1(\chi)`.

    Results are returned as `q`-expansions to a certain fixed (and fairly high)
    precision. If more precision is required this can be obtained with
    :func:`modular_ratio_to_prec`.

    EXAMPLES::

        sage: from sage.modular.modform.weight1 import hecke_stable_subspace
        sage: hecke_stable_subspace(DirichletGroup(59, QQ).0)
        [q - q^3 + q^4 - q^5 - q^7 - q^12 + q^15 + q^16 + 2*q^17 - q^19 - q^20 + q^21 + q^27 - q^28 - q^29 + q^35 + O(q^40)]
    """
    # Deal quickly with the easy cases.
    if chi(-1) == 1:
        return []
    N = chi.modulus()
    H = chi.kernel()
    G = GammaH(N, H)
    try:
        if ArithmeticSubgroup.dimension_cusp_forms(G, 1) == 0:
            verbose("no wt 1 cusp forms for N=%s, chi=%s by Riemann-Roch" % (N, chi._repr_short_()), level=1)
            return []
    except NotImplementedError:
        pass

    from sage.modular.modform.constructor import EisensteinForms
    chi = chi.minimize_base_ring()
    K = chi.base_ring()

    # Auxiliary prime for Hecke stability method
    l = aux_prime
    while l.divides(N):
        l = l.next_prime()
    verbose("Auxiliary prime: %s" % l, level=1)

    # Compute working precision
    R = l*Gamma0(N).sturm_bound(l + 2)

    t = verbose("Computing modular ratio space", level=1)
    mrs = modular_ratio_space(chi)

    t = verbose("Computing modular ratios to precision %s" % R, level=1)
    qexps = [modular_ratio_to_prec(chi, f, R) for f in mrs]
    verbose("Done", t=t, level=1)

    # We want to compute the largest subspace of I stable under T_l. To do
    # this, we compute I intersect T_l(I) modulo q^(R/l), and take its preimage
    # under T_l, which is then well-defined modulo q^R.

    from sage.modular.modform.hecke_operator_on_qexp import hecke_operator_on_qexp

    t = verbose("Computing Hecke-stable subspace", level=1)
    A = PowerSeriesRing(K, 'q')
    r = R // l
    V = K**R
    W = K**r
    Tl_images = [hecke_operator_on_qexp(f, l, 1, chi) for f in qexps]
    qvecs = [V(x.padded_list(R)) for x in qexps]
    qvecs_trunc = [W(x.padded_list(r)) for x in qexps]
    Tvecs = [W(x.padded_list(r)) for x in Tl_images]

    I = V.submodule(qvecs)
    Iimage = W.span(qvecs_trunc)
    TlI = W.span(Tvecs)
    Jimage = Iimage.intersection(TlI)
    J = I.Hom(W)(Tvecs).inverse_image(Jimage)

    verbose("Hecke-stable subspace is %s-dimensional" % J.dimension(), t=t, level=1)

    if J.rank() == 0:
        return []

    # The theory does not guarantee that J is exactly S_1(chi), just that it is
    # intermediate between S_1(chi) and M_1(chi). In every example I know of,
    # it is equal to S_1(chi), but just for honesty, we check this anyway.
    t = verbose("Checking cuspidality", level=1)
    JEis = V.span(V(x.padded_list(R)) for x in EisensteinForms(chi, 1).q_echelon_basis(prec=R))
    D = JEis.intersection(J)
    if D.dimension() != 0:
        raise ArithmeticError("Got non-cuspidal form!")
    verbose("Done", t=t, level=1)
    qexps = Sequence(A(x.list()).add_bigoh(R) for x in J.gens())
    return qexps


@cached_function
def dimension_wt1_cusp_forms(chi):
    r"""
    Return the dimension of the space of cusp forms of weight 1 and character chi.

    EXAMPLES::

        sage: chi = DirichletGroup(59, QQ).0
        sage: sage.modular.modform.weight1.dimension_wt1_cusp_forms(chi)
        1
    """
    return len(hecke_stable_subspace(chi))


@cached_function
def dimension_wt1_cusp_forms_gH(group):
    r"""
    Return the dimension of the space of cusp forms of weight 1 for the given
    group (which should be of GammaH type). Computed by summing over Galois
    orbits of characters modulo H.

    EXAMPLES::

        sage: sage.modular.modform.weight1.dimension_wt1_cusp_forms_gH(GammaH(31, [7]))
        1
    """
    chis = [g.minimize_base_ring() for g in group.characters_mod_H(galois_orbits=True)]
    return sum(dimension_wt1_cusp_forms(chi) * chi.base_ring().degree() for chi in chis)
