# coding: utf-8
r"""
Solution of polynomial systems using msolve

`msolve <https://msolve.lip6.fr/>`_ is a multivariate polynomial system solver
developed mainly by Jérémy Berthomieu (Sorbonne University), Christian Eder
(TU Kaiserslautern), and Mohab Safey El Din (Sorbonne University).

This module provide implementations of some operations on polynomial ideals
based on msolve. Currently the only supported operation is the computation of
the variety of zero-dimensional ideal over the rationals.

Note that msolve must be installed separately.

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
from sage.rings.rational_field import QQ
from sage.rings.real_arb import RealBallField
from sage.rings.real_double import RealDoubleField_class
from sage.rings.real_mpfr import RealField_class
from sage.rings.real_mpfi import RealIntervalField_class, RealIntervalField


def _variety(ideal, ring, proof):
    r"""
    Compute the variety of a zero-dimensional ideal using msolve.

    Part of the initial implementation was loosely based on the example
    interfaces available as part of msolve, with the authors' permission.

    TESTS::

        sage: K.<x, y> = PolynomialRing(QQ, 2, order='lex')
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

        sage: K.<x, y> = PolynomialRing(RR, 2, order='lex')
        sage: Ideal(x, y).variety(algorithm='msolve', proof=False)
        Traceback (most recent call last):
        ...
        NotImplementedError: unsupported base field: Real Field with 53 bits of precision

        sage: K.<x, y> = PolynomialRing(QQ, 2, order='lex')
        sage: Ideal(x, y).variety(ZZ, algorithm='msolve', proof=False)
        Traceback (most recent call last):
        ...
        ValueError: no coercion from base field Rational Field to output ring Integer Ring
    """

    # Normalize and check input

    base = ideal.base_ring()
    if ring is None:
        ring = base
    proof = sage.structure.proof.proof.get_flag(proof, "polynomial")
    if proof:
        raise ValueError("msolve relies on heuristics; please use proof=False")
    # As of msolve 0.2.4, prime fields seem to be supported, by I cannot
    # make sense of msolve's output in the positive characteristic case.
    # if not (base is QQ or isinstance(base, FiniteField) and
    #         base.is_prime_field() and base.characteristic() < 2**31):
    if base is not QQ:
        raise NotImplementedError(f"unsupported base field: {base}")
    if not ring.has_coerce_map_from(base):
        raise ValueError(
            f"no coercion from base field {base} to output ring {ring}")

    # Run msolve

    msolve().require()

    drlpolring = ideal.ring().change_ring(order='degrevlex')
    polys = ideal.change_ring(drlpolring).gens()
    msolve_in = tempfile.NamedTemporaryFile(mode='w',
                                            encoding='ascii', delete=False)
    command = ["msolve", "-f", msolve_in.name]
    if isinstance(ring, (RealIntervalField_class, RealBallField,
                         RealField_class, RealDoubleField_class)):
        parameterization = False
        command += ["-p", str(ring.precision())]
    else:
        parameterization = True
        command += ["-P", "1"]
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

    # Interpret output

    data = sage_eval(msolve_out.stdout[:-2])

    dim = data[0]
    if dim == -1:
        return []
    elif dim > 0:
        raise ValueError("positive-dimensional ideal")
    else:
        assert dim.is_zero()

    out_ring = ideal.ring().change_ring(ring)

    if parameterization:

        def to_poly(p, upol=PolynomialRing(base, 't')):
            assert len(p[1]) == p[0] + 1
            return upol(p[1])

        if len(data) != 3:
            raise NotImplementedError(
                f"unsupported msolve output format: {data}")
        [dim1, nvars, _, vars, _, [one, elim, den, param]] = data[1]
        assert dim1.is_zero()
        assert one.is_one()
        assert len(vars) == nvars
        ringvars = out_ring.variable_names()
        assert sorted(vars[:len(ringvars)]) == sorted(ringvars)
        vars = [out_ring(name) for name in vars[:len(ringvars)]]
        elim = to_poly(elim)
        den = to_poly(den)
        param = [to_poly(f)/d for [f, d] in param]
        elim_roots = elim.roots(ring, multiplicities=False)
        variety = []
        for rt in elim_roots:
            den_of_rt = den(rt)
            point = [-p(rt)/den_of_rt for p in param]
            if len(param) != len(vars):
                point.append(rt)
            assert len(point) == len(vars)
            variety.append(point)

    else:

        if len(data) != 2 or data[1][0] != 1:
            raise NotImplementedError(
                f"unsupported msolve output format: {data}")
        _, [_, variety] = data
        if isinstance(ring, (RealIntervalField_class, RealBallField)):
            to_out_ring = ring
        else:
            assert isinstance(ring, (RealField_class, RealDoubleField_class))
            myRIF = RealIntervalField(ring.precision())
            to_out_ring = lambda iv: ring.coerce(myRIF(iv).center())
        vars = out_ring.gens()
        variety = [[to_out_ring(iv) for iv in point]
                   for point in variety]

    return [KeyConvertingDict(out_ring, zip(vars, point)) for point in variety]

