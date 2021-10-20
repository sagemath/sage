r"""
Morphisms between finitely generated modules over a PID

AUTHOR:

- William Stein, 2009
"""
# *************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
# *************************************************************************

from sage.categories.morphism import Morphism, is_Morphism
from .fgp_module import DEBUG
from sage.structure.richcmp import richcmp, op_NE
from sage.misc.cachefunc import cached_method

class FGP_Morphism(Morphism):
    """
    A morphism between finitely generated modules over a PID.

    EXAMPLES:

    An endomorphism::

        sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
        sage: Q = V/W; Q
        Finitely generated module V/W over Integer Ring with invariants (4, 12)
        sage: phi = Q.hom([Q.0+3*Q.1, -Q.1]); phi
        Morphism from module over Integer Ring with invariants (4, 12) to module with invariants (4, 12) that sends the generators to [(1, 3), (0, 11)]
        sage: phi(Q.0) == Q.0 + 3*Q.1
        True
        sage: phi(Q.1) == -Q.1
        True

    A morphism between different modules V1/W1 ---> V2/W2 in
    different ambient spaces::

        sage: V1 = ZZ^2; W1 = V1.span([[1,2],[3,4]]); A1 = V1/W1
        sage: V2 = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W2 = V2.span([2*V2.0+4*V2.1, 9*V2.0+12*V2.1, 4*V2.2]); A2=V2/W2
        sage: A1
        Finitely generated module V/W over Integer Ring with invariants (2)
        sage: A2
        Finitely generated module V/W over Integer Ring with invariants (4, 12)
        sage: phi = A1.hom([2*A2.0])
        sage: phi(A1.0)
        (2, 0)
        sage: 2*A2.0
        (2, 0)
        sage: phi(2*A1.0)
        (0, 0)

    TESTS::

        sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W
        sage: phi = Q.hom([Q.0,Q.0 + 2*Q.1])
        sage: loads(dumps(phi)) == phi
        True
    """
    def __init__(self, parent, phi, check=True):
        """
        A morphism between finitely generated modules over a PID.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: phi = Q.hom([Q.0+3*Q.1, -Q.1]); phi
            Morphism from module over Integer Ring with invariants (4, 12) to module with invariants (4, 12) that sends the generators to [(1, 3), (0, 11)]
            sage: phi(Q.0) == Q.0 + 3*Q.1
            True
            sage: phi(Q.1) == -Q.1
            True

        For full documentation, see :class:`FGP_Morphism`.
        """
        Morphism.__init__(self, parent)
        M = parent.domain()
        N = parent.codomain()
        if isinstance(phi, FGP_Morphism):
            if check:
                if phi.parent() != parent:
                    raise TypeError
            phi = phi._phi
            check = False # no need

        # input: phi is a morphism from MO = M.optimized().V() to N.V()
        # that sends MO.W() to N.W()
        if check:
            if not is_Morphism(phi) and M == N:
                A = M.optimized()[0].V()
                B = N.V()
                s = M.base_ring()(phi) * B.coordinate_module(A).basis_matrix()
                phi = A.Hom(B)(s)

            MO, _ = M.optimized()
            if phi.domain() != MO.V():
                raise ValueError("domain of phi must be the covering module for the optimized covering module of the domain")
            if phi.codomain() != N.V():
                raise ValueError("codomain of phi must be the covering module the codomain.")
            # check that MO.W() gets sent into N.W()
            # todo (optimize): this is slow:
            for x in MO.W().basis():
                if phi(x) not in N.W():
                    raise ValueError("phi must send optimized submodule of M.W() into N.W()")
        self._phi = phi

    def _repr_(self):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: phi = Q.hom([Q.0+3*Q.1, -Q.1])
            sage: phi._repr_()
            'Morphism from module over Integer Ring with invariants (4, 12) to module with invariants (4, 12) that sends the generators to [(1, 3), (0, 11)]'
        """
        return "Morphism from module over %s with invariants %s to module with invariants %s that sends the generators to %s"%(
            self.domain().base_ring(), self.domain().invariants(), self.codomain().invariants(),
            list(self.im_gens()))

    @cached_method
    def im_gens(self):
        """
        Return tuple of the images of the generators of the domain
        under this morphism.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W
            sage: phi = Q.hom([Q.0,Q.0 + 2*Q.1])
            sage: phi.im_gens()
            ((1, 0), (1, 2))
            sage: phi.im_gens() is phi.im_gens()
            True
        """
        return tuple([self(x) for x in self.domain().gens()])

    def _richcmp_(self, right, op):
        """
        Comparison of ``self`` and ``right``.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ)
            sage: W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: phi = Q.hom([Q.0,Q.0 + 2*Q.1])
            sage: phi.im_gens()
            ((1, 0), (1, 2))
            sage: phi.im_gens() is phi.im_gens()
            True
            sage: phi == phi
            True
            sage: psi = Q.hom([Q.0,Q.0 - 2*Q.1])
            sage: phi == psi
            False
            sage: psi = Q.hom([Q.0,Q.0 - 2*Q.1])
            sage: phi < psi
            True
            sage: psi >= phi
            True
            sage: psi = Q.hom([Q.0,Q.0 + 2*Q.1])
            sage: phi == psi
            True
        """
        a = (self.domain(), self.codomain())
        b = (right.domain(), right.codomain())
        if a != b:
            return (op == op_NE)
        return richcmp(self.im_gens(), right.im_gens(), op)

    def __add__(self, right):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q=V/W; phi = Q.hom([2*Q.0, Q.1]); phi
            Morphism from module over Integer Ring with invariants (4, 12) to module with invariants (4, 12) that sends the generators to [(2, 0), (0, 1)]
            sage: phi + phi
            Morphism from module over Integer Ring with invariants (4, 12) to module with invariants (4, 12) that sends the generators to [(0, 0), (0, 2)]
        """
        if not isinstance(right, FGP_Morphism):  # todo: implement using coercion model
            right = self.parent()(right)
        return FGP_Morphism(self.parent(), self._phi + right._phi, check=DEBUG)

    def __sub__(self, right):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q=V/W; phi = Q.hom([2*Q.0, Q.1])
            sage: phi - phi
            Morphism from module over Integer Ring with invariants (4, 12) to module with invariants (4, 12) that sends the generators to [(0, 0), (0, 0)]
        """
        if not isinstance(right, FGP_Morphism):  # todo: implement using coercion model
            right = self.parent()(right)
        return FGP_Morphism(self.parent(), self._phi - right._phi, check=DEBUG)

    def __neg__(self):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q=V/W; phi = Q.hom([2*Q.0, Q.1])
            sage: -phi
            Morphism from module over Integer Ring with invariants (4, 12) to module with invariants (4, 12) that sends the generators to [(2, 0), (0, 11)]
        """
        return FGP_Morphism(self.parent(), -self._phi, check=DEBUG)

    def __call__(self, x):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W
            sage: phi = Q.hom([Q.0+3*Q.1, -Q.1])
            sage: phi(Q.0) == Q.0 + 3*Q.1
            True

        We compute the image of some submodules of the domain::

            sage: phi(Q)
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: phi(Q.submodule([Q.0]))
            Finitely generated module V/W over Integer Ring with invariants (4)
            sage: phi(Q.submodule([Q.1]))
            Finitely generated module V/W over Integer Ring with invariants (12)
            sage: phi(W/W)
            Finitely generated module V/W over Integer Ring with invariants ()

        We try to evaluate on a module that is not a submodule of the domain, which raises a ValueError::

            sage: phi(V/W.scale(2))
            Traceback (most recent call last):
            ...
            ValueError: x must be a submodule or element of the domain

        We evaluate on an element of the domain that is not in the V
        for the optimized representation of the domain::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: O, X = Q.optimized()
            sage: O.V()
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [ 0  6  1]
            [ 0 -2  0]
            sage: phi = Q.hom([Q.0, 4*Q.1])
            sage: x = Q(V.0); x
            (0, 8)
            sage: x == 8*Q.1
            True
            sage: x in O.V()
            False
            sage: phi(x)
            (0, 8)
            sage: phi(8*Q.1)
            (0, 8)
            sage: phi(8*Q.1) == phi(x)
            True
        """
        from .fgp_module import is_FGP_Module
        if is_FGP_Module(x):
            if not x.is_submodule(self.domain()):
                raise ValueError("x must be a submodule or element of the domain")
            # perhaps can be optimized with a matrix multiply; but note
            # the subtlety of optimized representations.
            return self.codomain().submodule([self(y) for y in x.smith_form_gens()])
        else:
            C = self.codomain()
            D = self.domain()
            O, X = D.optimized()
            x = D(x)
            if O is D:
                x = x.lift()
            else:
                # Now we have to transform x so that it is in the optimized representation.
                x = D.V().coordinate_vector(x.lift()) * X
            return C(self._phi(x))

    def kernel(self):
        """
        Compute the kernel of this homomorphism.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.hom([0, Q.1]).kernel()
            Finitely generated module V/W over Integer Ring with invariants (4)
            sage: A = Q.hom([Q.0, 0]).kernel(); A
            Finitely generated module V/W over Integer Ring with invariants (12)
            sage: Q.1 in A
            True
            sage: phi = Q.hom([Q.0-3*Q.1, Q.0+Q.1])
            sage: A = phi.kernel(); A
            Finitely generated module V/W over Integer Ring with invariants (4)
            sage: phi(A)
            Finitely generated module V/W over Integer Ring with invariants ()
        """
        # The kernel is just got by taking the inverse image of the submodule W
        # of the codomain quotient object.
        V = self._phi.inverse_image(self.codomain().W())
        D = self.domain()
        V = D.W() + V
        return D._module_constructor(V, D.W(), check=DEBUG)

    def inverse_image(self, A):
        """
        Given a submodule A of the codomain of this morphism, return
        the inverse image of A under this morphism.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: phi = Q.hom([0, Q.1])
            sage: phi.inverse_image(Q.submodule([]))
            Finitely generated module V/W over Integer Ring with invariants (4)
            sage: phi.kernel()
            Finitely generated module V/W over Integer Ring with invariants (4)
            sage: phi.inverse_image(phi.codomain())
            Finitely generated module V/W over Integer Ring with invariants (4, 12)

            sage: phi.inverse_image(Q.submodule([Q.0]))
            Finitely generated module V/W over Integer Ring with invariants (4)
            sage: phi.inverse_image(Q.submodule([Q.1]))
            Finitely generated module V/W over Integer Ring with invariants (4, 12)

            sage: phi.inverse_image(ZZ^3)
            Traceback (most recent call last):
            ...
            TypeError: A must be a finitely generated quotient module
            sage: phi.inverse_image(ZZ^3 / W.scale(2))
            Traceback (most recent call last):
            ...
            ValueError: A must be a submodule of the codomain
        """
        from .fgp_module import is_FGP_Module
        if not is_FGP_Module(A):
            raise TypeError("A must be a finitely generated quotient module")
        if not A.is_submodule(self.codomain()):
            raise ValueError("A must be a submodule of the codomain")
        V = self._phi.inverse_image(A.V())
        D = self.domain()
        V = D.W() + V
        return D._module_constructor(V, D.W(), check=DEBUG)

    def image(self):
        """
        Compute the image of this homomorphism.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = V/W; Q
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.hom([Q.0+3*Q.1, -Q.1]).image()
            Finitely generated module V/W over Integer Ring with invariants (4, 12)
            sage: Q.hom([3*Q.1, Q.1]).image()
            Finitely generated module V/W over Integer Ring with invariants (12)
        """
        V = self._phi.image() + self.codomain().W()
        W = V.intersection(self.codomain().W())
        return self.codomain()._module_constructor(V, W, check=DEBUG)

    def lift(self, x):
        """
        Given an element x in the codomain of self, if possible find an
        element y in the domain such that self(y) == x.  Raise a ValueError
        if no such y exists.

        INPUT:

        - ``x`` -- element of the codomain of self.

        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q=V/W; phi = Q.hom([2*Q.0, Q.1])
            sage: phi.lift(Q.1)
            (0, 1)
            sage: phi.lift(Q.0)
            Traceback (most recent call last):
            ...
            ValueError: no lift of element to domain
            sage: phi.lift(2*Q.0)
            (1, 0)
            sage: phi.lift(2*Q.0+Q.1)
            (1, 1)
            sage: V = span([[5, -1/2]],ZZ); W = span([[20,-2]],ZZ); Q = V/W; phi=Q.hom([2*Q.0])
            sage: x = phi.image().0; phi(phi.lift(x)) == x
            True

        """
        x = self.codomain()(x)

        # We view self as a map V/W --> V'/W', where V/W is the
        # optimized representation (which is fine to work with since
        # there is a lift to the optimized representation if and only
        # if there is a lift to the non-optimized representation).
        CD = self.codomain()
        A = self._phi.matrix()
        try:
            H, U = self.__lift_data
        except AttributeError:
            # Get the matrix of self: V --> V' wrt the basis for V and V'.

            # Stack it on top of the basis for W'.
            Wp = CD.V().coordinate_module(CD.W()).basis_matrix()
            B = A.stack(Wp)

            # Compute Hermite form of C with transformation
            H, U = B.hermite_form(transformation=True)
            self.__lift_data = H, U

        # write x in terms of the basis for V.
        w = CD.V().coordinate_vector(x.lift())

        # Solve z*H = w.
        try:
            z = H.solve_left(w)
            if z.denominator() != 1:
                raise ValueError
        except ValueError:
            raise ValueError("no lift of element to domain")

        # Write back in terms of rows of B, and delete rows not corresponding to A,
        # since those corresponding to relations
        v = (z*U)[:A.nrows()]

        # Take the linear combination that v defines.
        y = v*self.domain().optimized()[0].V().basis_matrix()

        # Return the finitely generated module element defined by y.
        y = self.domain()(y)
        assert self(y) == x, "bug in phi.lift()"
        return y

from sage.categories.homset import Homset

import sage.misc.weak_dict
_fgp_homset = sage.misc.weak_dict.WeakValueDictionary()


def FGP_Homset(X, Y):
    """
    EXAMPLES::

        sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W
        sage: Q.Hom(Q)           # indirect doctest
        Set of Morphisms from Finitely generated module V/W over Integer Ring with invariants (4, 12) to Finitely generated module V/W over Integer Ring with invariants (4, 12) in Category of modules over Integer Ring
        sage: True # Q.Hom(Q) is Q.Hom(Q)
        True
        sage: type(Q.Hom(Q))
        <class 'sage.modules.fg_pid.fgp_morphism.FGP_Homset_class_with_category'>
    """
    key = (X, Y)
    try:
        return _fgp_homset[key]
    except KeyError:
        pass
    H = FGP_Homset_class(X, Y)
    # Caching breaks tests in fgp_module.
    # _fgp_homset[key] = H
    return H


class FGP_Homset_class(Homset):
    """
    Homsets of :class:`~sage.modules.fg_pid.fgp_module.FGP_Module`

    TESTS::

        sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W
        sage: H = Hom(Q,Q); H    # indirect doctest
        Set of Morphisms from Finitely generated module V/W over Integer Ring with invariants (4, 12) to Finitely generated module V/W over Integer Ring with invariants (4, 12) in Category of modules over Integer Ring
        sage: type(H)
        <class 'sage.modules.fg_pid.fgp_morphism.FGP_Homset_class_with_category'>
    """
    Element = FGP_Morphism

    def __init__(self, X, Y, category=None):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W
            sage: type(Q.Hom(Q))
            <class 'sage.modules.fg_pid.fgp_morphism.FGP_Homset_class_with_category'>
        """
        if category is None:
            from sage.modules.free_module import is_FreeModule
            if is_FreeModule(X) and is_FreeModule(Y):
                from sage.categories.all import FreeModules
                category = FreeModules(X.base_ring())
            else:
                from sage.categories.all import Modules
                category = Modules(X.base_ring())
        Homset.__init__(self, X, Y, category)

    def _coerce_map_from_(self, S):
        """
        EXAMPLES::

            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W
            sage: phi = Q.hom([Q.0,Q.0 + 2*Q.1]);  psi = loads(dumps(phi))
            sage: phi.parent()._coerce_map_from_(psi.parent())
            True
            sage: phi.parent()._coerce_map_from_(Q.Hom(ZZ^3))
            False
        """
        # We define this so that morphisms in equal parents canonically coerce,
        # since otherwise, e.g., the dumps(loads(...)) doctest above would fail.
        if isinstance(S, FGP_Homset_class) and S == self:
            return True
        if self.is_endomorphism_set():
            R = self.domain().base_ring()
            return R == S or bool(R._coerce_map_from_(S))
        return False

    def __call__(self, x):
        """
        Convert x into an fgp morphism.

        EXAMPLES::

            sage: V = span([[1/2,0,0],[3/2,2,1],[0,0,1]],ZZ); W = V.span([V.0+2*V.1, 9*V.0+2*V.1, 4*V.2])
            sage: Q = V/W; H = Q.Hom(Q)
            sage: H(3)
            Morphism from module over Integer Ring with invariants (4, 16) to module with invariants (4, 16) that sends the generators to [(3, 0), (0, 3)]
        """
        return self.element_class(self, x)
