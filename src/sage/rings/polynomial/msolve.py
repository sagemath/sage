r"""
Solution of polynomial systems using msolve

`msolve <https://msolve.lip6.fr/>`_ is a multivariate polynomial system solver
based on Gröbner bases.

This module provide implementations of some operations on polynomial ideals
based on msolve.

Note that the :ref:`optional package msolve <spkg_msolve>` must be installed.

.. SEEALSO::

    - :mod:`sage.features.msolve`
    - :mod:`sage.rings.polynomial.multi_polynomial_ideal`
"""

import os
import tempfile
import subprocess

import sage.structure.proof.proof

from sage.features.msolve import msolve
from sage.misc.converting_dict import KeyConvertingDict
from sage.misc.sage_eval import sage_eval
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.rational_field import QQ
from sage.rings.real_arb import RealBallField
from sage.rings.real_double import RealDoubleField_class
from sage.rings.real_mpfr import RealField_class
from sage.rings.real_mpfi import RealIntervalField_class, RealIntervalField
from sage.structure.sequence import Sequence


def _run_msolve(ideal, options):
    r"""
    Internal utility function
    """

    base = ideal.base_ring()
    if not (base is QQ or isinstance(base, FiniteField) and
            base.is_prime_field() and base.characteristic() < 2**31):
        raise NotImplementedError(f"unsupported base field: {base}")

    # Run msolve

    drlpolring = ideal.ring().change_ring(order='degrevlex')
    polys = ideal.change_ring(drlpolring).gens()
    msolve_in = tempfile.NamedTemporaryFile(mode='w',
                                            encoding='ascii', delete=False)
    command = [msolve().absolute_filename(), "-f", msolve_in.name] + options
    try:
        print(",".join(drlpolring.variable_names()), file=msolve_in)
        print(base.characteristic(), file=msolve_in)
        print(*(pol._repr_().replace(" ", "") for pol in polys),
                sep=',\n', file=msolve_in)
        msolve_in.close()
        msolve_out = subprocess.run(command, capture_output=True, text=True)
    finally:
        os.unlink(msolve_in.name)
    msolve_out.check_returncode()

    return msolve_out.stdout


def groebner_basis_degrevlex(ideal, proof=True):
    r"""
    Compute a degrevlex Gröbner basis using msolve

    EXAMPLES::

        sage: from sage.rings.polynomial.msolve import groebner_basis_degrevlex

        sage: R.<a,b,c> = PolynomialRing(GF(101), 3, order='lex')
        sage: I = sage.rings.ideal.Katsura(R,3)
        sage: gb = groebner_basis_degrevlex(I); gb # optional - msolve
        [a + 2*b + 2*c - 1, b*c - 19*c^2 + 10*b + 40*c,
        b^2 - 41*c^2 + 20*b - 20*c, c^3 + 28*c^2 - 37*b + 13*c]
        sage: gb.universe() is R # optional - msolve
        False
        sage: gb.universe().term_order() # optional - msolve
        Degree reverse lexicographic term order
        sage: ideal(gb).transformed_basis(other_ring=R) # optional - msolve
        [c^4 + 38*c^3 - 6*c^2 - 6*c, 30*c^3 + 32*c^2 + b - 14*c,
        a + 2*b + 2*c - 1]

    Gröbner bases over the rationals require `proof=False`::

        sage: R.<x, y> = PolynomialRing(QQ, 2)
        sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])
        sage: I.groebner_basis(algorithm='msolve') # optional - msolve
        Traceback (most recent call last):
        ...
        ValueError: msolve relies on heuristics; please use proof=False
        sage: I.groebner_basis(algorithm='msolve', proof=False) # optional - msolve
        [x*y - 1, x^2 + y^2 - 4*x - 2*y + 4, y^3 - 2*y^2 + x + 4*y - 4]

    TESTS::

        sage: R.<foo, bar> = PolynomialRing(GF(536870909), 2)
        sage: I = Ideal([ foo^2 - 1, bar^2 - 1 ])
        sage: I.groebner_basis(algorithm='msolve') # optional - msolve
        [bar^2 - 1, foo^2 - 1]
    """

    if ideal.base_ring() is QQ and sage.structure.proof.proof.get_flag(proof, "polynomial"):
            raise ValueError("msolve relies on heuristics; please use proof=False")

    drlpolring = ideal.ring().change_ring(order='degrevlex')
    msolve_out = _run_msolve(ideal, ["-g", "2"])
    gbasis = sage_eval(msolve_out[:-2], locals=drlpolring.gens_dict())
    return Sequence(gbasis)


def variety(ideal, ring, *, proof=True):
    r"""
    Compute the variety of a zero-dimensional ideal using msolve.

    Part of the initial implementation was loosely based on the example
    interfaces available as part of msolve, with the authors' permission.

    EXAMPLES::

        sage: from sage.rings.polynomial.msolve import variety
        sage: p = 536870909
        sage: R.<x, y> = PolynomialRing(GF(p), 2, order='lex')
        sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])
        sage: sorted(variety(I, GF(p^2), proof=False), key=str) # optional - msolve
        [{x: 1, y: 1},
         {x: 254228855*z2 + 114981228, y: 232449571*z2 + 402714189},
         {x: 267525699, y: 473946006},
         {x: 282642054*z2 + 154363985, y: 304421338*z2 + 197081624}]

    TESTS::

        sage: p = 536870909
        sage: R.<x, y> = PolynomialRing(GF(p), 2, order='lex')
        sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])

        sage: sorted(I.variety(algorithm='msolve', proof=False), key=str) # optional - msolve
        [{x: 1, y: 1}, {x: 267525699, y: 473946006}]

        sage: K.<a> = GF(p^2)
        sage: sorted(I.variety(K, algorithm='msolve', proof=False), key=str) # optional - msolve
        [{x: 1, y: 1},
         {x: 118750849*a + 194048031, y: 510295713*a + 18174854},
         {x: 267525699, y: 473946006},
         {x: 418120060*a + 75297182, y: 26575196*a + 44750050}]

        sage: R.<x, y> = PolynomialRing(GF(2147483659), 2, order='lex')
        sage: ideal([x, y]).variety(algorithm='msolve', proof=False)
        Traceback (most recent call last):
        ...
        NotImplementedError: unsupported base field: Finite Field of size 2147483659

        sage: R.<x, y> = PolynomialRing(QQ, 2, order='lex')
        sage: I = Ideal([ x*y - 1, (x-2)^2 + (y-1)^2 - 1])

        sage: I.variety(algorithm='msolve', proof=False) # optional - msolve
        [{x: 1, y: 1}]
        sage: I.variety(RealField(100), algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.7692923542386314152404094643, y: 0.36110308052864737763464656216},
         {x: 1.0000000000000000000000000000, y: 1.0000000000000000000000000000}]
        sage: I.variety(RealIntervalField(100), algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.76929235423863141524040946434?, y: 0.361103080528647377634646562159?},
         {x: 1, y: 1}]
        sage: I.variety(RBF, algorithm='msolve', proof=False) # optional - msolve
        [{x: [2.76929235423863 +/- 2.08e-15], y: [0.361103080528647 +/- 4.53e-16]},
         {x: 1.000000000000000, y: 1.000000000000000}]
        sage: I.variety(RDF, algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.7692923542386314, y: 0.36110308052864737}, {x: 1.0, y: 1.0}]
        sage: I.variety(AA, algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.769292354238632?, y: 0.3611030805286474?},
         {x: 1.000000000000000?, y: 1.000000000000000?}]
        sage: I.variety(QQbar, algorithm='msolve', proof=False) # optional - msolve
        [{x: 2.769292354238632?, y: 0.3611030805286474?},
         {x: 1, y: 1},
         {x: 0.11535382288068429? + 0.5897428050222055?*I, y: 0.3194484597356763? - 1.633170240915238?*I},
         {x: 0.11535382288068429? - 0.5897428050222055?*I, y: 0.3194484597356763? + 1.633170240915238?*I}]
        sage: I.variety(ComplexField(100))
        [{y: 1.0000000000000000000000000000, x: 1.0000000000000000000000000000},
         {y: 0.36110308052864737763464656216, x: 2.7692923542386314152404094643},
         {y: 0.31944845973567631118267671892 - 1.6331702409152376561188467320*I, x: 0.11535382288068429237979526783 + 0.58974280502220550164728074602*I},
         {y: 0.31944845973567631118267671892 + 1.6331702409152376561188467320*I, x: 0.11535382288068429237979526783 - 0.58974280502220550164728074602*I}]

        sage: Ideal(x^2 + y^2 - 1, x - y).variety(RBF, algorithm='msolve', proof=False) # optional - msolve
        [{x: [-0.707106781186547 +/- 6.29e-16], y: [-0.707106781186547 +/- 6.29e-16]},
         {x: [0.707106781186547 +/- 6.29e-16], y: [0.707106781186547 +/- 6.29e-16]}]
        sage: sorted(Ideal(x^2 - 1, y^2 - 1).variety(QQ, algorithm='msolve', proof=False), key=str) # optional - msolve
        [{x: -1, y: -1}, {x: -1, y: 1}, {x: 1, y: -1}, {x: 1, y: 1}]
        sage: Ideal(x^2-1, y^2-2).variety(CC, algorithm='msolve', proof=False) # optional - msolve
        [{x: 1.00000000000000, y: 1.41421356237310},
         {x: -1.00000000000000, y: 1.41421356237309},
         {x: 1.00000000000000, y: -1.41421356237309},
         {x: -1.00000000000000, y: -1.41421356237310}]

        sage: Ideal([x, y, x + y]).variety(algorithm='msolve', proof=False) # optional - msolve
        [{x: 0, y: 0}]

        sage: Ideal([x, y, x + y - 1]).variety(algorithm='msolve', proof=False) # optional - msolve
        []
        sage: Ideal([x, y, x + y - 1]).variety(RR, algorithm='msolve', proof=False) # optional - msolve
        []

        sage: Ideal([x*y - 1]).variety(QQbar, algorithm='msolve', proof=False) # optional - msolve
        Traceback (most recent call last):
        ...
        ValueError: positive-dimensional ideal

        sage: R.<x, y> = PolynomialRing(RR, 2, order='lex')
        sage: Ideal(x, y).variety(algorithm='msolve', proof=False)
        Traceback (most recent call last):
        ...
        NotImplementedError: unsupported base field: Real Field with 53 bits of precision

        sage: R.<x, y> = PolynomialRing(QQ, 2, order='lex')
        sage: Ideal(x, y).variety(ZZ, algorithm='msolve', proof=False)
        Traceback (most recent call last):
        ...
        ValueError: no coercion from base field Rational Field to output ring Integer Ring
    """

    proof = sage.structure.proof.proof.get_flag(proof, "polynomial")
    if proof:
        raise ValueError("msolve relies on heuristics; please use proof=False")

    base = ideal.base_ring()
    if ring is None:
        ring = base
    if not ring.has_coerce_map_from(base):
        raise ValueError(
            f"no coercion from base field {base} to output ring {ring}")

    if isinstance(ring, (RealIntervalField_class, RealBallField,
                         RealField_class, RealDoubleField_class)):
        parameterization = False
        options = ["-p", str(ring.precision())]
    else:
        parameterization = True
        options = ["-P", "1"]

    msolve_out = _run_msolve(ideal, options)

    # Interpret output

    try:
        data = sage_eval(msolve_out[:-2])
    except SyntaxError:
        raise NotImplementedError(f"unsupported msolve output format: {data}")

    dim = data[0]
    if dim == -1:
        return []
    elif dim > 0:
        raise ValueError("positive-dimensional ideal")
    else:
        assert dim.is_zero()

    out_ring = ideal.ring().change_ring(ring)

    if parameterization:

        def to_poly(p, d=1, *, upol=PolynomialRing(base, 't')):
            assert len(p[1]) == p[0] + 1 or p == [-1, [0]]
            return upol(p[1])/d

        try:
            [char, nvars, deg, vars, _, [one, [elim, den, param]]] = data[1]
        except (IndexError, ValueError):
            raise NotImplementedError(
                f"unsupported msolve output format: {data}")
        assert char == ideal.base_ring().characteristic()
        assert one.is_one()
        assert len(vars) == nvars
        ringvars = out_ring.variable_names()
        assert sorted(vars[:len(ringvars)]) == sorted(ringvars)
        vars = [out_ring(name) for name in vars[:len(ringvars)]]
        elim = to_poly(elim)
        # Criterion suggested by Mohab Safey El Din to avoid cases where there
        # is no rational parameterization or where the one returned by msolve
        # has a significant probability of being incorrect.
        if deg >= char > 0 or 0 < char <= 2**17 and deg != elim.degree():
            raise NotImplementedError(f"characteristic {char} too small")
        den = to_poly(den)
        # As of msolve 0.4.4, param is of the form [pol, denom] in char 0, but
        # [pol] in char p > 0. My understanding is that both cases will
        # eventually use the same format, so let's not be too picky.
        param = [to_poly(*f) for f in param]
        elim_roots = elim.roots(ring, multiplicities=False)
        variety = []
        for rt in elim_roots:
            den_of_rt = den(rt)
            point = [-p(rt) / den_of_rt for p in param]
            if len(param) != len(vars):
                point.append(rt)
            assert len(point) == len(vars)
            variety.append(point)

    else:

        if len(data[1]) < 2 or len(data[1]) != data[1][0] + 1:
            raise NotImplementedError(
                f"unsupported msolve output format: {data}")
        if isinstance(ring, (RealIntervalField_class, RealBallField)):
            to_out_ring = ring
        else:
            assert isinstance(ring, (RealField_class, RealDoubleField_class))
            myRIF = RealIntervalField(ring.precision())
            to_out_ring = lambda iv: ring.coerce(myRIF(iv).center())
        vars = out_ring.gens()
        variety = [[to_out_ring(iv) for iv in point]
                   for l in data[1][1:]
                   for point in l]

    return [KeyConvertingDict(out_ring, zip(vars, point)) for point in variety]
