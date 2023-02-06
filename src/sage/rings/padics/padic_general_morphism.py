r"""
TODO
"""
# ****************************************************************************
#       Copyright (C)      2019 David Roe <roed.math@gmail.com>
#                     2019-2022 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.rings.morphism import RingHomomorphism

class pAdicGeneralMorphism(RingHomomorphism):
    """
    A homomorphism from a relative extension L/K of p-adic rings or fields to an arbitrary ring A.

    One can specify such a homomorphism by giving either
    * a homomorphism from K into A together with the image in A of the generator of L/K.
    * a homomorphism from the backend representation to A.

    TODO: When A is non-exact, describe precision.

    INPUT:

        - ``parent`` -- the Homset containing this morphism
        - ``im_gen`` -- the image of the generator of L/K (a length one list is also accepted)
        - ``base_hom`` -- the homomorphism from the base ring (defaults to the coercion map from K to A)
        - ``backend_hom`` -- the homomorphism on the two-step extension representating L/K
        - ``check`` -- whether to check that ``im_gen`` defines a valid homomorphism
    """
    def __init__(self, parent, im_gen=None, base_hom=None, backend_hom=None, check=True):
        # We specify either im_gen or backend_hom but not both
        L = parent.domain()
        K = L.base_ring()
        A = parent.codomain()
        Ax = A['x']
        B = L._backend()
        if backend_hom is None:
            if base_hom is None:
                base_hom = A.coerce_map_from(K)
                if base_hom is None:
                    raise ValueError("Must specify homomorphism on the base ring")
            elif base_hom.domain() is not K or base_hom.codomain() is not A:
                raise ValueError("Base homomorphism does not have correct domain/codomain")
            if im_gen is None:
                if self.degree() > 1:
                    raise ValueError("Must specify the image of the generator")
            elif check:
                # Should check be done here or in the hom constructions down below?
                if im_gen.parent() is not A:
                    raise ValueError
                pol = Ax([base_hom(c) for c in L.defining_polynomial()])
                if pol(im_gen) != 0:
                    raise ValueError("relations do not all (canonically) map to 0 under map determined by images of generators")

            def get_image(pol):
                if im_gen is None:
                    assert pol.degree() < 1
                    return pol[0]
                else:
                    return pol(im_gen)
            if B.absolute_e() == 1:
                # B is unramified
                if B.absolute_f() == 1:
                    # B is a trivial extension of Qp or Zp
                    # Every hom from Qp or Zp is a natural one
                    backend_hom = B.hom(A)
                else:
                    u = L(B.gen())
                    u_image = get_image(Ax([base_hom(c) for c in u.polynomial()]))
                    backend_hom = B.hom([u_image], check=False)
            else:
                if B.absolute_f() == 1:
                    # B is a totally ramified extension of Qp or Zp
                    pi = L(B.uniformizer())
                    pi_image = get_image(Ax([base_hom(c) for c in pi.polynomial()]))
                    backend_hom = B.hom([pi_image], check=False)
                else:
                    # Two step extension
                    U = B.base_ring()
                    u = L(B(U.gen()))
                    pi = L(B.uniformizer())
                    u_image = get_image(Ax([base_hom(c) for c in u.polynomial()]))
                    pi_image = get_image(Ax([base_hom(c) for c in pi.polynomial()]))
                    backend_hom = B.hom([pi_image], base_hom=U.hom([u_image]), check=False)
        else:
            if im_gen is not None or base_hom is not None:
                raise NotImplementedError("Cannot specify both im_gen/base_hom and backend_hom")
