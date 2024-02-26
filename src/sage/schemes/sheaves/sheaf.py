r"""
Coherent sheaves

EXAMPLES:

We define the Fermat cubic surface in \PP^2 and construct its structure sheaf::

    sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
    sage: X = P2.subscheme(x^4 + y^4 + z^4)
    sage: sh = X.structure_sheaf()

AUTHORS:

- Kwankyu Lee (2024-01-22): initial version

"""

from functools import cached_property
from sage.structure.sage_object import SageObject
from sage.modules.free_module import FreeModule


class Sheaf(SageObject):
    r"""
    Coherent sheaf on a projective scheme.

    INPUT:

    - ``scheme`` -- the base scheme on which the sheaf is defined

    - ``module`` -- a free module or its quotient by a submodule

    - ``twist`` -- (default: 0) an integer

    This class constructs the coherent sheaf `\tilde M(n)` if `M` is the
    ``module`` and `n` is the ``twist``.
    """
    def __init__(self, scheme, module, twist=0):
        """
        Initialize ``self``.
        """
        try:
            if module.is_ambient():
                module = module.quotient(module.zero_submodule())
        except AttributeError:
            pass

        assert module.cover() == module.free_cover()

        self._base_scheme = scheme
        self._module = module
        self._twist = twist

    @cached_property
    def _cohomology(self):
        """
        Return an object that computes the cohomology.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: sheaf = X.structure_sheaf()
            sage: c = sheaf._cohomology
            sage: c.H(1).dimension()
            3
            sage: c.h(1)
            3
        """
        raise NotImplementedError('_cohomology is not implemented')

    def _repr_(self):
        sheaf = 'Twisted Sheaf' if self._twist else 'Sheaf'
        return f'{sheaf} on {self._base_scheme}'

    def base_scheme(self):
        """
        Return the base scheme on which this sheaf is defined.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: s = X.structure_sheaf()
            sage: s.base_scheme()
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^4 + y^4 + z^4
        """
        return self._base_scheme

    def defining_twist(self):
        """
        Return the integer by which the module defining this coherent sheaf is
        twisted.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: s = X.structure_sheaf(3)
            sage: s.defining_twist()
            3
        """
        return self._twist

    def defining_module(self):
        """
        Return the module defining this coherent sheaf.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: X.structure_sheaf()
            Sheaf on Closed subscheme of Projective Space of dimension 2 over
            Rational Field defined by: x^4 + y^4 + z^4
        """
        return self._module

    def cohomology_group(self, r=0):
        """
        Return the `r`-th cohomology as a vector space.

        INPUT:

        - ``r`` -- (default: 0) a non-negative integer

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: s = X.structure_sheaf()
            sage: s.cohomology_group(1).dimension()
            3
        """
        return self._cohomology.H(r)

    def cohomology(self, r=0):
        """
        Return the dimension of the `r`-th cohomology as a vector space.

        INPUT:

        - ``r`` -- (default: 0) a non-negative integer

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: sheaf = X.structure_sheaf()
            sage: sheaf.cohomology(0)
            1
            sage: sheaf.cohomology(1)
            3
            sage: sheaf.cohomology(2)
            0
        """
        return self._cohomology.h(r)

    def twist(self, t=0):
        r"""
        Return the twisted sheaf `F(n)` of this sheaf `F`.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: sh = X.structure_sheaf()
            sage: sh.twist(1).euler_characteristic()
            2
            sage: sh.twist(2).euler_characteristic()
            6
        """
        return type(self)(self._base_scheme, self._module, self._twist + t)

    def euler_characteristic(self):
        """
        Return the Euler characteristic of this coherent sheaf.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: sh = X.structure_sheaf()
            sage: sh.euler_characteristic()
            -2
        """
        d = self._base_scheme.dimension()
        chi = 0
        for r in range(d + 1):  # for Grothendieck's vanishing theorem
            d = self.cohomology(r)
            if r % 2:
                chi = chi - d
            else:
                chi = chi + d
        return chi


class Sheaf_on_projective_space(Sheaf):

    @cached_property
    def _cohomology(self):
        """
        This property keeps the cohomology object for this sheaf.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P2.structure_sheaf()._cohomology
            Maruyama Complex induced from S(0) <-- 0
        """
        from sage.schemes.sheaves.cohomology import MaruyamaComplex
        return MaruyamaComplex(self._module, twist=self._twist)


class Sheaf_on_projective_subscheme(Sheaf):

    @cached_property
    def _cohomology(self):
        """
        This property keeps the cohomology object for this sheaf.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: X.structure_sheaf()._cohomology
            Maruyama Complex induced from S(0) <-- S(-4) <-- 0
        """
        return self.image_to_ambient_space()._cohomology

    def image_to_ambient_space(self):
        """
        Return the direct image of this sheaf to the ambient space.

        The image is with respect to the inclusion morphism from the base
        scheme into the projective Space.

        INPUT:

        - ``twist`` -- (default: `0`) an integer

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(x^4 + y^4 + z^4)
            sage: X.structure_sheaf().image_to_ambient_space()
            Sheaf on Projective Space of dimension 2 over Rational Field
        """
        X = self._base_scheme
        A = X.ambient_space()
        S = A.coordinate_ring()

        d = self._module.degree()
        M = FreeModule(S, d)
        I = X.defining_polynomials()
        J = self._module.relations().gens()
        G = [f * M.gen(i) for i in range(d) for f in I] + [v.change_ring(S) for v in J]
        N = M.submodule(G)
        return A.coherent_sheaf(M.quotient(N), twist=self._twist)
