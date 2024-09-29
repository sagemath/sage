from sage.misc.lazy_attribute import lazy_attribute

from sage.categories.modules import Modules
from sage.categories.category_types import Category_over_base_ring
from sage.categories.homsets import Homsets

from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

class OreModules(Category_over_base_ring):
    r"""
    Category of Ore modules.
    """
    @staticmethod
    def __classcall_private__(cls, ring, twist):
        r"""
        Normalize the input and call the init function.

        INPUT:

        - ``ring`` -- a commutative ring, the base ring of
          the Ore modules

        - ``twist`` -- a twisting morphism/derivation or a
          Ore polynomial ring

        TESTS::

            sage: from sage.categories.ore_modules import OreModules
            sage: K.<a> = GF(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: cat = OreModules(K, Frob)
            sage: cat
            Category of Ore modules over Finite Field in a of size 5^3 twisted by a |--> a^5

            sage: S = cat.ore_ring('y')
            sage: cat is OreModules(K, S)
            True
        """
        if isinstance(twist, OrePolynomialRing):
            ore = twist.change_var('x')
            if ore.base_ring() is not ring:
                raise ValueError("base rings do not match")
        else:
            ore = OrePolynomialRing(ring, twist, names='x', polcast=False)
        return cls.__classcall__(cls, ore)

    def __init__(self, ore):
        r"""
        Initialize this category.

        TESTS::

            sage: from sage.categories.ore_modules import OreModules
            sage: K.<a> = GF(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: cat = OreModules(K, Frob)

            sage: TestSuite(cat).run()
        """
        base = ore.base_ring()
        Category_over_base_ring.__init__(self, base)
        self._ore = ore

    def __reduce__(self):
        r"""
        Return the arguments which were used to create this instance.

        This method is needed for pickling.

        TESTS::

            sage: from sage.categories.ore_modules import OreModules
            sage: K.<a> = GF(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: cat = OreModules(K, Frob)
            sage: cat2 = loads(dumps(cat))  # indirect doctest
            sage: cat is cat2
            True
        """
        return OreModules, (self.base_ring(), self._ore)

    def super_categories(self):
        r"""
        Return the immediate super categories of this category.

        EXAMPLES::

            sage: from sage.categories.ore_modules import OreModules
            sage: K.<a> = GF(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: cat = OreModules(K, Frob)
            sage: cat.super_categories()
            [Category of vector spaces over Finite Field in a of size 5^3]
        """
        return [Modules(self.base())]

    def _repr_object_names(self):
        r"""
        Return a string representation naming the objects
        in this category.

        EXAMPLES::

            sage: from sage.categories.ore_modules import OreModules
            sage: K.<a> = GF(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: cat = OreModules(K, Frob)
            sage: cat._repr_object_names()
            'Ore modules over Finite Field in a of size 5^3 twisted by a |--> a^5'
        """
        return "Ore modules over %s %s" % (self.base_ring(), self._ore._repr_twist())

    def ore_ring(self, var='x'):
        r"""
        Return the underlying Ore polynomial ring.

        INPUT:

        - ``var`` (default; ``x``) -- the variable name

        EXAMPLES::

            sage: from sage.categories.ore_modules import OreModules
            sage: K.<a> = GF(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: cat = OreModules(K, Frob)
            sage: cat.ore_ring()
            Ore Polynomial Ring in x over Finite Field in a of size 5^3 twisted by a |--> a^5

            sage: cat.ore_ring('y')
            Ore Polynomial Ring in y over Finite Field in a of size 5^3 twisted by a |--> a^5
        """
        return self._ore.change_var(var)

    def twisting_morphism(self):
        r"""
        Return the underlying twisting morphism.

        EXAMPLES::

            sage: from sage.categories.ore_modules import OreModules
            sage: K.<a> = GF(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: cat = OreModules(K, Frob)
            sage: cat.twisting_morphism()
            Frobenius endomorphism a |--> a^5 on Finite Field in a of size 5^3

        If the twising morphism is the identity, nothing is returned::

            sage: R.<t> = QQ[]
            sage: d = R.derivation()
            sage: cat = OreModules(R, d)
            sage: cat.twisting_morphism()
        """
        return self._ore.twisting_morphism()

    def twisting_derivation(self):
        r"""
        Return the underlying twisting derivation.

        EXAMPLES::

            sage: from sage.categories.ore_modules import OreModules
            sage: R.<t> = QQ[]
            sage: d = R.derivation()
            sage: cat = OreModules(R, d)
            sage: cat.twisting_derivation()
            d/dt

        If the twising derivation is zero, nothing is returned::

            sage: K.<a> = GF(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: cat = OreModules(K, Frob)
            sage: cat.twisting_derivation()
        """
        return self._ore.twisting_derivation()
