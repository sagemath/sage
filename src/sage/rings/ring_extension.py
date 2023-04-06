r"""
Extension of rings

Sage offers the possibility to work with ring extensions `L/K` as
actual parents and perform meaningful operations on them and their
elements.

The simplest way to build an extension is to use the method
:meth:`sage.categories.commutative_rings.CommutativeRings.ParentMethods.over` on the top ring,
that is `L`.
For example, the following line constructs the extension of
finite fields `\mathbf{F}_{5^4}/\mathbf{F}_{5^2}`::

    sage: GF(5^4).over(GF(5^2))
    Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base

By default, Sage reuses the canonical generator of the top ring
(here `z_4 \in \mathbf{F}_{5^4}`), together with its name. However,
the user can customize them by passing in appropriate arguments::

    sage: F = GF(5^2)
    sage: k = GF(5^4)
    sage: z4 = k.gen()
    sage: K.<a> = k.over(F, gen = 1-z4)
    sage: K
    Field in a with defining polynomial x^2 + z2*x + 4 over its base

The base of the extension is available via the method :meth:`base` (or
equivalently :meth:`base_ring`)::

    sage: K.base()
    Finite Field in z2 of size 5^2

It is also possible to build an extension on top of another extension,
obtaining this way a tower of extensions::

    sage: L.<b> = GF(5^8).over(K)
    sage: L
    Field in b with defining polynomial x^2 + (4*z2 + 3*a)*x + 1 - a over its base
    sage: L.base()
    Field in a with defining polynomial x^2 + z2*x + 4 over its base
    sage: L.base().base()
    Finite Field in z2 of size 5^2

The method :meth:`bases` gives access to the complete list of rings in
a tower::

    sage: L.bases()
    [Field in b with defining polynomial x^2 + (4*z2 + 3*a)*x + 1 - a over its base,
     Field in a with defining polynomial x^2 + z2*x + 4 over its base,
     Finite Field in z2 of size 5^2]

Once we have constructed an extension (or a tower of extensions), we
have interesting methods attached to it. As a basic example, one can
compute a basis of the top ring over any base in the tower::

    sage: L.basis_over(K)
    [1, b]
    sage: L.basis_over(F)
    [1, a, b, a*b]

When the base is omitted, the default is the natural base of the extension::

    sage: L.basis_over()
    [1, b]

The method :meth:`sage.rings.ring_extension_element.RingExtensionWithBasis.vector`
computes the coordinates of an element according to the above basis::

    sage: u = a + 2*b + 3*a*b
    sage: u.vector()   # over K
    (a, 2 + 3*a)
    sage: u.vector(F)
    (0, 1, 2, 3)

One can also compute traces and norms with respect to any base of the tower::

    sage: u.trace()           # over K
    (2*z2 + 1) + (2*z2 + 1)*a
    sage: u.trace(F)
    z2 + 1
    sage: u.trace().trace()   # over K, then over F
    z2 + 1

    sage: u.norm()            # over K
    (z2 + 1) + (4*z2 + 2)*a
    sage: u.norm(F)
    2*z2 + 2

And minimal polynomials::

    sage: u.minpoly()
    x^2 + ((3*z2 + 4) + (3*z2 + 4)*a)*x + (z2 + 1) + (4*z2 + 2)*a
    sage: u.minpoly(base=F)
    x^4 + (4*z2 + 4)*x^3 + x^2 + (z2 + 1)*x + 2*z2 + 2


AUTHOR:

- Xavier Caruso (2019)
"""

# **************************************************************************
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#                  2022 Julian Rüth <julian.rueth@fsfe.org>
#
#    This program is free softwGare: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ***************************************************************************


from sage.misc.cachefunc import cached_method
from sage.cpython.getattr import AttributeErrorMessage
from sage.cpython.getattr import dir_with_other_class
from sage.misc.latex import latex, latex_variable_name

from sage.structure.factory import UniqueFactory
from sage.structure.category_object import normalize_names
from sage.categories.map import Map
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields
from sage.rings.ring import CommutativeAlgebra
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity

from sage.rings.ring_extension_element import (
    RingExtensionElement, RingExtensionFractionFieldElement, RingExtensionWithBasisElement)
from sage.rings.ring_extension_morphism import (
    RingExtensionHomomorphism, RingExtensionHomomorphism_baseinclusion, RingExtensionBackendIsomorphism, RingExtensionBackendReverseIsomorphism,
    are_different_morphisms, MapFreeModuleToRelativeRing, MapRelativeRingToFreeModule)
from sage.rings.ring_extension_conversion import (
    backend_parent, to_backend, from_backend)


# Helper functions
##################

def tower_bases(ring, degree):
    r"""
    Return the list of bases of ``ring`` (including itself); if
    degree is ``True``, restrict to finite extensions and return
    in addition the degree of ``ring`` over each base.

    INPUT:

    - ``ring`` -- a commutative ring

    - ``degree`` -- a boolean

    EXAMPLES::

        sage: from sage.rings.ring_extension import tower_bases
        sage: S.<x> = QQ[]
        sage: T.<y> = S[]
        sage: tower_bases(T, False)
        ([Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field,
          Univariate Polynomial Ring in x over Rational Field,
          Rational Field],
         [])
        sage: tower_bases(T, True)
        ([Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field],
         [1])

        sage: K.<a> = Qq(5^2)
        sage: L.<w> = K.extension(x^3 - 5)
        sage: tower_bases(L, True)
        ([5-adic Eisenstein Extension Field in w defined by x^3 - 5 over its base field,
          5-adic Unramified Extension Field in a defined by x^2 + 4*x + 2,
          5-adic Field with capped relative precision 20],
         [1, 3, 6])
    """
    bases = [ ]
    degrees = [ ]
    base = ring
    deg = 1
    while True:
        bases.append(base)
        if degree:
            degrees.append(deg)
            try:
                d = base.relative_degree()
            except AttributeError:
                try:
                    d = base.degree()
                except AttributeError:
                    break
            if d is Infinity:
                break
            deg *= d
        newbase = base._base
        if newbase is base:
            break
        base = newbase
    return bases, degrees


def common_base(K, L, degree):
    """
    Return a common base on which ``K`` and ``L`` are defined.

    INPUT:

    - ``K`` -- a commutative ring

    - ``L`` -- a commutative ring

    - ``degree`` -- a boolean; if true, return the degree of
      ``K`` and ``L`` over their common base

    EXAMPLES::

        sage: from sage.rings.ring_extension import common_base

        sage: common_base(GF(5^3), GF(5^7), False)
        Finite Field of size 5
        sage: common_base(GF(5^3), GF(5^7), True)
        (Finite Field of size 5, 3, 7)

        sage: common_base(GF(5^3), GF(7^5), False)
        Traceback (most recent call last):
        ...
        NotImplementedError: unable to find a common base

    When ``degree`` is set to ``True``, we only look up for bases on
    which both ``K`` and ``L`` are finite::

        sage: S.<x> = QQ[]
        sage: common_base(S, QQ, False)
        Rational Field
        sage: common_base(S, QQ, True)
        Traceback (most recent call last):
        ...
        NotImplementedError: unable to find a common base

    """
    bases_K, degrees_K = tower_bases(K, degree)
    bases_L, degrees_L = tower_bases(L, degree)
    base = None
    for iL in range(len(bases_L)):
        try:
            iK = bases_K.index(bases_L[iL])
            base = bases_L[iL]
            break
        except ValueError:
            pass
    if base is None:
        raise NotImplementedError("unable to find a common base")
    if degree:
        return base, degrees_K[iK], degrees_L[iL]
    else:
        return base


def _generators(ring, base):
    r"""
    Return the generators of ``ring`` over ``base``.

    INPUT:

    - ``ring`` -- a backend ring

    - ``base`` -- a base ring of ``ring`` or ``None`` which is equivalent to
      ``ring.base_ring()``

    EXAMPLES::

        sage: from sage.rings.ring_extension import _generators
        sage: S.<x> = QQ[]
        sage: T.<y> = S[]

        sage: _generators(T, S)
        (y,)
        sage: _generators(T, QQ)
        (y, x)
    """
    gens = tuple()

    while ring is not ring.base_ring() and (base is None or not base.has_coerce_map_from(ring)):
        gens += tuple(ring.gens())
        ring = ring.base_ring()

    if base is not None:
        gens = tuple([x for x in gens if x not in base])

    return gens


def variable_names(ring, base):
    r"""
    Return the variable names of the generators of ``ring``
    over ``base``.

    INPUT:

    - ``ring`` -- a commutative ring

    - ``base`` -- a commutative ring

    EXAMPLES::

        sage: from sage.rings.ring_extension import variable_names
        sage: S.<x> = QQ[]
        sage: T.<y> = S[]

        sage: variable_names(T, S)
        ('y',)
        sage: variable_names(T, QQ)
        ('y', 'x')
    """
    names = tuple()
    while ring is not ring.base_ring() and (base is None or not base.has_coerce_map_from(ring)):
        gens = ring.gens()
        vars = ring.variable_names()
        if len(gens) != len(vars):
            raise NotImplementedError("cannot figure out the variable names")
        if base is None:
            names += tuple(vars)
        else:
            for gen, var in zip(gens, vars):
                if gen not in base:
                    names += (var,)
        ring = ring.base_ring()
    return names


# Factory
#########

class RingExtensionFactory(UniqueFactory):
    """
    Factory for ring extensions.

    TESTS::

        sage: E = QQ.over(ZZ)
        sage: QQ.over(ZZ) is E
        True

        sage: K.<a> = QQ.extension(x^2 - 2)
        sage: E = K.over(QQ)
        sage: E
        Field in a with defining polynomial x^2 - 2 over its base

        sage: E2.<b> = K.over(QQ)
        sage: E2 is E
        False
    """
    def create_key_and_extra_args(self, ring, defining_morphism=None, gens=None, names=None, canonical_backend=None, constructors=None):
        """
        Create a key and return it together with a list of constructors
        of the object.

        INPUT:

        - ``ring`` -- a commutative ring

        - ``defining_morphism`` -- a ring homomorphism or a commutative
          ring or ``None`` (default: ``None``); the defining morphism of
          this extension or its base (if it coerces to ``ring``)

        - ``gens`` -- a list of generators of this extension (over its base)
          or ``None`` (default: ``None``);

        - ``names`` -- a list or a tuple of variable names or ``None``
          (default: ``None``)

        - ``canonical_backend`` -- whether the association between this ring extension is canonical,
          or if there are choices involved (such as choosing a root of some defining polynomial).
          This affects the existence of a coercion from the backend and coercions between ring extensions.
          Defaults to True if the defining morphism is the coercion map.

        - ``constructors`` -- a list of constructors; each constructor
          is a pair ``(class, arguments)`` where ``class`` is the class
          implementing the extension and ``arguments`` is the dictionary
          of arguments to pass in to init function

        TESTS::

            sage: from sage.rings.ring_extension import RingExtension
            sage: RingExtension.create_key_and_extra_args(QQ, ZZ)
            ((Natural morphism:
                From: Integer Ring
                To:   Rational Field,
              (),
              (),
              True),
             {'constructors': [(<class 'sage.rings.ring_extension.RingExtension_generic'>,
                {'print_options': {'print_elements_as': None, 'print_parent_as': None}})]})

            sage: RingExtension.create_key_and_extra_args(GF(5^4), GF(5^2), names=('a',))
            ((Ring morphism:
                From: Finite Field in z2 of size 5^2
                To:   Finite Field in z4 of size 5^4
                Defn: z2 |--> z4^3 + z4^2 + z4 + 3,
              (z4,),
              ('a',),
              True),
             {'constructors': [(<class 'sage.rings.ring_extension.RingExtensionWithGen'>,
                {'gen': z4, 'names': ('a',)}),
               (<class 'sage.rings.ring_extension.RingExtension_generic'>,
                {'print_options': {'print_elements_as': None, 'print_parent_as': None}})]})

        We can fall back on the generic construction even in the presence of generators,
        since RingExtensionWithGen only supports finite extensions::

            sage: R.<x> = GF(5)['x'].over()
            sage: type(R)
            <class 'sage.rings.ring_extension.RingExtension_generic_with_category'>
        """
        use_generic_constructor = True
        print_as = None

        if defining_morphism is None:
            base = ring.base_ring()
        elif isinstance(defining_morphism, Map):
            base = defining_morphism.domain()
        elif defining_morphism in CommutativeRings():
            base = defining_morphism
            defining_morphism = None
        else:
            raise TypeError(f"ring extension cannot be created from defining morphism {defining_morphism}")

        # We compute the defining morphism
        if defining_morphism is None:
            defining_morphism = ring.coerce_map_from(base)
            if defining_morphism is None and isinstance(base, RingExtension_generic):
                backend_base = base._backend
                if ring.has_coerce_map_from(backend_base):
                    bHr = base.Hom(ring)
                    defining_morphism = bHr.__make_element_class__(RingExtensionHomomorphism)(bHr, ring.coerce_map_from(backend_base))
            if defining_morphism is None:
                raise ValueError(f"No coercion map from {base} to {ring}")
            if canonical_backend is None:
                canonical_backend = True
        else:
            if defining_morphism.domain() is not base:
                defining_morphism = defining_morphism.extend_domain(base)
            if defining_morphism.codomain() is not ring:
                defining_morphism = defining_morphism.extend_codomain(ring)
            if canonical_backend is None:
                canonical_backend = not are_different_morphisms(defining_morphism, None)

        # We normalize other attributes
        if gens is not None:
            if not isinstance(gens, (list, tuple)):
                raise TypeError("gens must be a list or a tuple")
            gens = tuple(ring(g) for g in gens)
            if names is None:
                raise TypeError("you must specify the names of the generators")
            names = normalize_names(len(gens), names)
            #use_generic_constructor = False
        else:
            gens = _generators(ring, base)
            if names is None:
                try:
                    names = variable_names(ring, base)
                except NotImplementedError:
                    gens = names = None
            else:
                names = normalize_names(len(gens), names)
                #use_generic_constructor = False

        # We figure out what are the best constructors
        if constructors is None:
            constructors = []
            if gens is not None and len(gens) == 1:
                constructors.append((RingExtensionWithGen,
                                     {'gen': gens[0], 'names': names}))
            if use_generic_constructor:
                constructors.append((RingExtension_generic,
                                     {'print_options': {'print_parent_as': print_as,
                                                        'print_elements_as': print_as}}))

        assert defining_morphism.codomain() is ring

        # We build the key and return it
        return (defining_morphism, gens, names, canonical_backend), {'constructors': constructors}

    def create_object(self, version, key, **extra_args):
        """
        Return the object associated to a given key.

        TESTS::

            sage: from sage.rings.ring_extension import RingExtension
            sage: key, extra_args = RingExtension.create_key_and_extra_args(QQ, ZZ, canonical_backend=True)
            sage: RingExtension.create_object((8,9,0), key, **extra_args)
            Rational Field over its base
        """
        if len(key) == 3:
            defining_morphism, gens, names = key
            canonical_backend = True
        else:
            defining_morphism, gens, names, canonical_backend = key
        constructors = extra_args['constructors']
        if len(constructors) == 0:
            raise NotImplementedError("no constructor available for this extension")
        for (constructor, kwargs) in constructors[:-1]:
            try:
                return constructor(defining_morphism, canonical_backend=canonical_backend, **kwargs)
            except (NotImplementedError, ValueError, TypeError):
                pass
        (constructor, kwargs) = constructors[-1]
        return constructor(defining_morphism, canonical_backend=canonical_backend, **kwargs)


RingExtension = RingExtensionFactory("sage.rings.ring_extension.RingExtension")


# General extensions
####################

class RingExtension_generic(CommutativeAlgebra):
    r"""
    A generic class for all ring extensions.

    TESTS::

        sage: Q = QQ.over(ZZ)  # indirect doctest
        sage: Q
        Rational Field over its base

        sage: type(Q)
        <class 'sage.rings.ring_extension.RingExtension_generic_with_category'>

        sage: TestSuite(Q).run()

    """
    Element = RingExtensionElement

    def __init__(self, defining_morphism, print_options={}, import_methods=True, canonical_backend=False, category=None, check=True):
        r"""
        Initialize this ring extension.

        INPUT:

        - ``defining_morphism`` -- a ring homomorphism

        - ``print_options`` -- a dictionary

        - ``import_methods`` -- a boolean (default: ``True``); whether this
          parent (resp. its elements) import the methods of the backend
          parent class (resp. element class)

        - ``canonical_backend`` -- a boolean (default: ``False``)
          Whether the association between this ring extension and its backend is canonical,
          or if there are choices involved (such as choosing a root of some defining polynomial).
          This affects the existence of a coercion from the backend and coercions between ring extensions.

        - ``category`` -- the category for the resulting parent
          (default: ``CommutativeRings()``)

        OUTPUT:

        The extension defined by ``defining_morphism``

        EXAMPLES::

            sage: QQ.over(ZZ)
            Rational Field over its base

            sage: S.<x> = QQ[]
            sage: S.over()  # over QQ
            Univariate Polynomial Ring in x over Rational Field over its base

        TESTS::

            sage: ZZ.over(NN)
            Traceback (most recent call last):
            ...
            TypeError: ring extension cannot be created from defining morphism Non negative integer semiring

            sage: K = GF(5^3)
            sage: L = K.over(K.frobenius_endomorphism()); L
            Finite Field in z3 of size 5^3 over its base
            sage: L(K.gen())
            2*z3^2 + 4*z3 + 4
            sage: K.gen()^5
            2*z3^2 + 4*z3 + 4

        """
        base = defining_morphism.domain()
        ring = defining_morphism.codomain()
        if category is None:
            # Another option would be to set category = CommutativeAlgebras(base)
            # but CommutativeRings() seems safer, especially when dealing with
            # morphisms which do not need to preserve the base
            category = CommutativeRings()
        CommutativeAlgebra.__init__(self, ZZ, category=category)
        self._base = base
        self._backend = ring
        self._backend_defining_morphism = defining_morphism
        bHs = base.Hom(self)
        self._defining_morphism = bHs.__make_element_class__(RingExtensionHomomorphism_baseinclusion)(bHs)
        self._print_options = print_options.copy()
        if 'over' not in self._print_options:
            self._print_options['over'] = ZZ(0)
        self._import_methods = import_methods
        self._canonical_backend = canonical_backend
        self._type = "Ring"
        if self._backend in Fields():
            self._type = "Field"

        if check:
            if (base not in CommutativeRings()
                or ring not in CommutativeRings()
                or not defining_morphism.category_for().is_subcategory(CommutativeRings())):
                raise TypeError("only commutative rings are accepted")

        self.register_coercion(self._defining_morphism.__copy__())

        rHs = ring.Hom(self)
        sHr = self.Hom(ring)
        self._from_backend_morphism = rHs.__make_element_class__(RingExtensionBackendIsomorphism)(rHs)
        self._to_backend_morphism = sHr.__make_element_class__(RingExtensionBackendReverseIsomorphism)(sHr)

        if canonical_backend:
            self.register_coercion(self._from_backend_morphism)
            ring.register_conversion(self._to_backend_morphism)

    def __getattr__(self, name):
        """
        If this extension was created with ``import_methods = True``,
        return a wrapper to the corresponding method of the backend
        parent (if it exists).

        EXAMPLES::

            sage: K.<a> = QQ.extension(x^2 - 2)
            sage: E = K.over()  # over QQ

            sage: hasattr(E, 'automorphisms')
            True
            sage: E.automorphisms()
            [Ring endomorphism of Field in a with defining polynomial x^2 - 2 over its base
               Defn: a |--> a,
             Ring endomorphism of Field in a with defining polynomial x^2 - 2 over its base
               Defn: a |--> -a]
        """
        try:
            return super().__getattr__(name)
        except AttributeError:
            pass
        method = None
        if self._import_methods and hasattr(self._backend, name):
            method = getattr(self._backend, name)
        if not callable(method):
            raise AttributeError(AttributeErrorMessage(self, name))

        def wrapper(*args, **kwargs):
            output = method(*to_backend(args), **to_backend(kwargs))
            return from_backend(output, self)
        wrapper.__doc__ = method.__doc__
        return wrapper

    def __dir__(self):
        """
        Return the list of all the attributes of this extension;
        if the extension was created with ``import_methods = True``,
        concatenate this list with the list of all the methods of
        the backend parent.

        EXAMPLES::

            sage: A.<a> = QQ.extension(x^2 - 2)
            sage: K.<a> = A.over()

            sage: dir(K)
            ['CartesianProduct',
             'Element',
             'Hom',
             ...
             'zeta',
             'zeta_coefficients',
             'zeta_function',
             'zeta_order']
        """
        d = dir_with_other_class(self, self.category().parent_class)
        if not self._import_methods:
            return d
        for name in dir(self._backend):
            if name[0] == "_":
                continue
            try:
                attribute = getattr(self._backend, name)
                if callable(attribute):
                    d.append(name)
            except Exception:
                pass
        return sorted(set(d))

    def construction(self):
        """
        Return the functorial construction of this extension, if defined.

        EXAMPLES::

             sage: E = GF(5^3).over()
             sage: E.construction()

        """
        # One could define a construction functor K' -> K' otimes_K L, but we leave this to another issue
        pass

    def from_base_ring(self, r):
        r"""
        Return the canonical embedding of ``r`` into this extension.

        INPUT:

        - ``r`` -- an element of the base of the ring of this extension

        EXAMPLES::

            sage: k = GF(5)
            sage: K.<u> = GF(5^2).over(k)
            sage: L.<v> = GF(5^4).over(K)

            sage: x = L.from_base_ring(k(2)); x
            2
            sage: x.parent()
            Field in v with defining polynomial x^2 + (3 - u)*x + u over its base

            sage: x = L.from_base_ring(u); x
            u
            sage: x.parent()
            Field in v with defining polynomial x^2 + (3 - u)*x + u over its base
        """
        r = self._base.coerce(r)
        return self.element_class(self, self._backend_defining_morphism(r))

    def print_options(self, **options):
        """
        Update the printing options of this extension.

        INPUT:

        - ``over`` -- an integer or ``Infinity`` (default: ``0``); the maximum
          number of bases included in the printing of this extension

        - ``base`` -- a base over which this extension is finite free;
          elements in this extension will be printed as a linear
          combinaison of a basis of this extension over the given base

        EXAMPLES::

            sage: A.<a> = GF(5^2).over()   # over GF(5)
            sage: B.<b> = GF(5^4).over(A)
            sage: C.<c> = GF(5^12).over(B)
            sage: D.<d> = GF(5^24).over(C)

        Observe what happens when we modify the option ``over``::

            sage: D
            Field in d with defining polynomial x^2 + ((1 - a) + ((1 + 2*a) - b)*c + ((2 + a) + (1 - a)*b)*c^2)*x + c over its base

            sage: D.print_options(over=2)
            sage: D
            Field in d with defining polynomial x^2 + ((1 - a) + ((1 + 2*a) - b)*c + ((2 + a) + (1 - a)*b)*c^2)*x + c over
            Field in c with defining polynomial x^3 + (1 + (2 - a)*b)*x^2 + (2 + 2*b)*x - b over
            Field in b with defining polynomial x^2 + (3 - a)*x + a over its base

            sage: D.print_options(over=Infinity)
            sage: D
            Field in d with defining polynomial x^2 + ((1 - a) + ((1 + 2*a) - b)*c + ((2 + a) + (1 - a)*b)*c^2)*x + c over
            Field in c with defining polynomial x^3 + (1 + (2 - a)*b)*x^2 + (2 + 2*b)*x - b over
            Field in b with defining polynomial x^2 + (3 - a)*x + a over
            Field in a with defining polynomial x^2 + 4*x + 2 over
            Finite Field of size 5

        Now the option ``base``::

            sage: d^2
            -c + ((-1 + a) + ((-1 + 3*a) + b)*c + ((3 - a) + (-1 + a)*b)*c^2)*d

            sage: D.basis_over(B)
            [1, c, c^2, d, c*d, c^2*d]
            sage: D.print_options(base=B)
            sage: d^2
            -c + (-1 + a)*d + ((-1 + 3*a) + b)*c*d + ((3 - a) + (-1 + a)*b)*c^2*d

            sage: D.basis_over(A)
            [1, b, c, b*c, c^2, b*c^2, d, b*d, c*d, b*c*d, c^2*d, b*c^2*d]
            sage: D.print_options(base=A)
            sage: d^2
            -c + (-1 + a)*d + (-1 + 3*a)*c*d + b*c*d + (3 - a)*c^2*d + (-1 + a)*b*c^2*d
        """
        for (name, value) in options.items():
            method = None
            if hasattr(self, '_print_option_' + name):
                method = getattr(self, '_print_option_' + name)
            if not callable(method):
                raise ValueError("option '%s' does not exist" % name)
            self._print_options[name] = method(value)

    def _print_option_over(self, over):
        """
        Check and normalize the print option ``over``

        INPUT:

        - ``over`` -- an integer or ``Infinity``

        OUTPUT:

        The normalized value of ``over``

        TESTS::

            sage: E = QQ.over(ZZ)
            sage: E.print_options(over=-2) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: 'over' must be nonnegative

            sage: E.print_options(over=ZZ)
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce <class 'sage.rings.integer_ring.IntegerRing_class'> to an integer
        """
        if over is not None and over is not Infinity:
            over = ZZ(over)
            if over < 0:
                raise ValueError("'over' must be nonnegative")
        return over

    def _repr_(self, **options):
        r"""
        Return a string representation of this extension.

        INPUT:

        - ``over`` -- an integer, ``Infinity`` or ``None``; the maximum
          number of bases included in the print representation of
          this extension;
          if ``None``, use the print options of this extension.

        EXAMPLES::

            sage: E = QQ.over(ZZ)
            sage: E
            Rational Field over its base

            sage: E._repr_()
            'Rational Field over its base'
            sage: E._repr_(over=Infinity)
            'Rational Field over Integer Ring'
        """
        if 'print_parent_as' in options:
            print_as = options.pop('print_parent_as')
        else:
            print_as = self._print_options.get('print_parent_as')
        if print_as is not None:
            if isinstance(print_as, RingExtension_generic):
                return print_as._repr_(**options)
            else:
                return str(print_as)
        print_options = self._print_options.copy()
        for (name, value) in options.items():
            method = None
            if hasattr(self, '_print_option_' + name):
                method = getattr(self, '_print_option_' + name)
            if not callable(method):
                raise ValueError("option '%s' does not exist" % name)
            print_options[name] = method(value)
        over = print_options.pop('over', None)

        try:
            trivial = self.relative_degree() == 1
        except NotImplementedError:
            trivial = False

        if trivial and self._canonical_backend:
            s = "Trivial extension of "
        else:
            s = self._repr_topring(**print_options)
            if over > 0:
                s += " over "
                over -= 1
            else:
                s += " over its base"
                return s
        if trivial or over >= 0:
            base = self._base
            if isinstance(base, RingExtension_generic):
                s += base._repr_(over=over)
            else:
                s += str(base)
        return s

    def _repr_topring(self, **options):
        r"""
        Return a string representation of top ring of this extension.

        EXAMPLES::

            sage: E = QQ.over(ZZ)
            sage: E._repr_topring()
            'Rational Field'
        """
        return str(self._backend)

    def _latex_(self, **options):
        r"""
        Return a LaTeX representation of this extension.

        - ``over`` -- an integer, ``Infinity`` or ``None``; the maximum
          number of bases included in the LaTeX representation of
          this extension;
          if ``None``, use the print options of this extension.

        EXAMPLES::

            sage: E = QQ.over(ZZ)
            sage: latex(E)
            \Bold{Q}

            sage: E._latex_()
            \Bold{Q}
            sage: E._latex_(over=Infinity)
            \Bold{Q} / \Bold{Z}
        """
        if 'print_parent_as' in options:
            print_as = options.pop('print_parent_as')
        else:
            print_as = self._print_options.get('print_parent_as')
        if print_as is not None:
            if isinstance(print_as, RingExtension_generic):
                return print_as._latex_(**options)
            else:
                return latex(print_as)
        print_options = self._print_options.copy()
        for (name, value) in options.items():
            method = None
            if hasattr(self, '_print_option_' + name):
                method = getattr(self, '_print_option_' + name)
            if not callable(method):
                raise ValueError("option '%s' does not exist" % name)
            print_options[name] = method(value)
        over = print_options.pop('over', None)

        try:
            trivial = self.relative_degree() == 1
        except NotImplementedError:
            trivial = False

        s = ""
        if not trivial:
            s = self._latex_topring(**print_options)
            if over > 0:
                s += " / "
                over -= 1
            else:
                return s
        if trivial or over >= 0:
            base = self._base
            if isinstance(base, RingExtension_generic):
                s += base._latex_(over=over)
            else:
                s += latex(base)
        return s

    def _latex_topring(self, **options):
        r"""
        Return a string representation of top ring of this extension.

        EXAMPLES::

            sage: E = QQ.over(ZZ)
            sage: E._latex_topring()
            \Bold{Q}
        """
        return latex(self._backend)

    def _coerce_map_from_(self, other):
        r"""
        Return a coerce map from this extension to ``other`` if defined.

        COERCION MODEL

        If `L/K` is an extension, a coercion map `K \to (L/K)`
        (acting through the defining morphism of `L/K`) is set.

        If `L_1/K_1` and `L_2/K_2` are two extensions, a coercion
        map `(L_1/K_1) \to (L_2/K_2)`` is set when `L_1` coerces to
        `L_2` and `K_1` coerces to `K_2` in such a way that the
        appropriate diagram commutes.

        These rules have the following consequence regarding iterated
        extensions.
        Given two iterated extensions `A = (A_n/\cdots/A_2/A_1)` and
        `B = (B_m/\cdots/B_2/B_1)`, there is a coercion map `A \to B`
        if there exists a strictly increasing function
        `\sigma : \{1,\ldots,n\} \to \{1,\ldots,m\}` and coercion maps
        `A_i \to B_{\sigma(i)}` making all the appropriate diagrams
        commutative.

        .. NOTE::

            In order to avoid discrepancies, it is forbidden to create
            an extension with exotic defining morphisms:
            if (A_n/\cdots/A_2/A_1) is an iterated extension and `i
            \leq j` are two indices such that `A_i` coerces to `A_j`,
            then the composition defining morphism `A_i \to A_{i+1}
            \to \cdots \to A_j` must agree with the coercion map.

        TESTS::

            sage: E1 = GF(3^6).over(GF(3^3))
            sage: E1.coerce_map_from(GF(3^3))  # indirect doctest
            Base injection morphism:
              From: Finite Field in z3 of size 3^3
              To:   Field in z6 with defining polynomial x^2 + (2*z3 + 1)*x + z3 over its base

            sage: E1.coerce_map_from(GF(3))    # indirect doctest
            Composite map:
              From: Finite Field of size 3
              To:   Field in z6 with defining polynomial x^2 + (2*z3 + 1)*x + z3 over its base
              Defn:   Ring morphism:
                      From: Finite Field of size 3
                      To:   Finite Field in z3 of size 3^3
                      Defn: 1 |--> 1
                    then
                      Base injection morphism:
                      From: Finite Field in z3 of size 3^3
                      To:   Field in z6 with defining polynomial x^2 + (2*z3 + 1)*x + z3 over its base

            sage: E2 = GF(3^18).over(GF(3^9))
            sage: E2.coerce_map_from(E1)       # indirect doctest
            Ring morphism:
              From: Field in z6 with defining polynomial x^2 + (2*z3 + 1)*x + z3 over its base
              To:   Field in z18 with defining polynomial x^2 + (z9^8 + 2*z9^7 + z9^5 + 2*z9^4 + z9^2 + z9 + 1)*x + z9 over its base
              Defn: z6 |--> (2*z9^7 + z9^6 + 2*z9^2 + 2*z9) + (z9^8 + 2*z9^7 + 2*z9^6 + z9^5 + z9^3 + z9 + 1)*z18

        A test with iterated extensions::

            sage: A = GF(3^18).over(GF(3^3))   #   simple extension GF(3^3) -> GF(3^18)
            sage: B = GF(3^18).over(E1)        # iterated extension GF(3^3) -> GF(3^6) -> GF(3^18)
            sage: A.has_coerce_map_from(B)
            False
            sage: B.has_coerce_map_from(A)
            True

        """
        if isinstance(other, RingExtension_generic) and self._canonical_backend and other._canonical_backend:
            right = other
            if self._backend.has_coerce_map_from(right._backend) and self._base.has_coerce_map_from(right._base):
                back_map = self._backend.coerce_map_from(right._backend)
                base_map = self._base.coerce_map_from(right._base)
                f = back_map * right._backend_defining_morphism
                g = self._backend_defining_morphism * base_map
                if not are_different_morphisms(f, g):
                    back_map = self._from_backend_morphism * back_map
                    rHs = right.Hom(self)
                    return rHs.__make_element_class__(RingExtensionHomomorphism)(rHs, back_map)

    def base(self):
        r"""
        Return the base of this extension.

        EXAMPLES::

            sage: F = GF(5^2)
            sage: K = GF(5^4).over(F)
            sage: K.base()
            Finite Field in z2 of size 5^2

        In case of iterated extensions, the base is itself an extension::

            sage: L = GF(5^8).over(K)
            sage: L.base()
            Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
            sage: L.base() is K
            True

        .. SEEALSO::

            :meth:`bases`, :meth:`absolute_base`, :meth:`is_defined_over`
        """
        return self._base

    def bases(self):
        r"""
        Return the list of successive bases of this extension
        (including itself).

        EXAMPLES::

            sage: F = GF(5^2).over()  # over GF(5)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(K)

            sage: F.bases()
            [Field in z2 with defining polynomial x^2 + 4*x + 2 over its base,
             Finite Field of size 5]

            sage: K.bases()
            [Field in z4 with defining polynomial x^2 + (3 - z2)*x + z2 over its base,
             Field in z2 with defining polynomial x^2 + 4*x + 2 over its base,
             Finite Field of size 5]

            sage: L.bases()
            [Field in z12 with defining polynomial x^3 + (1 + (2 - z2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base,
             Field in z4 with defining polynomial x^2 + (3 - z2)*x + z2 over its base,
             Field in z2 with defining polynomial x^2 + 4*x + 2 over its base,
             Finite Field of size 5]

        .. SEEALSO::

            :meth:`base`, :meth:`absolute_base`, :meth:`is_defined_over`
        """
        L = [ self ]
        base = self
        while isinstance(base, RingExtension_generic):
            base = base.base_ring()
            L.append(base)
        return L

    def absolute_base(self):
        r"""
        Return the absolute base of this extension.

        By definition, the absolute base of an iterated extension
        `K_n/\cdots K_2/K_1` is the ring `K_1`.

        EXAMPLES::

            sage: F = GF(5^2).over()   # over GF(5)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(K)

            sage: F.absolute_base()
            Finite Field of size 5
            sage: K.absolute_base()
            Finite Field of size 5
            sage: L.absolute_base()
            Finite Field of size 5

        .. SEEALSO::

            :meth:`base`, :meth:`bases`, :meth:`is_defined_over`
        """
        return self.bases()[-1]

    def is_defined_over(self, base):
        r"""
        Return whether or not ``base`` is one of the bases of this
        extension.

        INPUT:

        - ``base`` -- a commutative ring, which might be itself an
          extension

        EXAMPLES::

            sage: A = GF(5^4).over(GF(5^2))
            sage: B = GF(5^12).over(A)

            sage: A.is_defined_over(GF(5^2))
            True
            sage: A.is_defined_over(GF(5))
            False

            sage: B.is_defined_over(A)
            True
            sage: B.is_defined_over(GF(5^4))
            True
            sage: B.is_defined_over(GF(5^2))
            True
            sage: B.is_defined_over(GF(5))
            False

        Note that an extension is defined over itself::

            sage: A.is_defined_over(A)
            True
            sage: A.is_defined_over(GF(5^4))
            True

        .. SEEALSO::

            :meth:`base`, :meth:`bases`, :meth:`absolute_base`
        """
        b = self
        while isinstance(b, RingExtension_generic):
            if b is base or b._backend is base:
                return True
            b = b._base
        return b is base

    def _check_base(self, base):
        r"""
        Check if ``base`` is one of the successive bases of this
        extension and, if it is, normalize it.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        OUTPUT:

        The base ``base`` normalized as a parent appearing in the
        list of bases of this extension as returned by :meth:`bases`.

        EXAMPLES::

            sage: F = GF(5^2)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(K)
            sage: L.bases()
            [Field in z12 with defining polynomial x^3 + (1 + (4*z2 + 2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base,
             Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base,
             Finite Field in z2 of size 5^2]

            sage: L._check_base(K)
            Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base

        The backend field is not supported as a base since (in the presence of
        trivial extensions) it might be both a field in the tower and a backend
        field::

            sage: L._check_base(GF(5^4))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z4 of size 5^4

        When ``base`` is ``None``, the base of the extension is returned::

            sage: L._check_base(None)
            Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
            sage: L._check_base(None) is L.base()
            True

        """
        if base is None:
            return self._base
        b = self
        while isinstance(b, RingExtension_generic):
            if b is base:
                return b
            b = b._base
        if b is base:
            return b
        raise ValueError("not (explicitly) defined over %s" % base)

    def defining_morphism(self, base=None):
        r"""
        Return the defining morphism of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: F = GF(5^2)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(K)

            sage: K.defining_morphism()
            Base injection morphism:
              From: Finite Field in z2 of size 5^2
              To:   Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base

            sage: L.defining_morphism()
            Base injection morphism:
              From: Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
              To:   Field in z12 with defining polynomial x^3 + (1 + (4*z2 + 2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base

        One can also pass in a base over which the extension is explicitly
        defined (see also :meth:`is_defined_over`)::

            sage: L.defining_morphism(F)
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Field in z12 with defining polynomial x^3 + (1 + (4*z2 + 2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base
              Defn:   Base injection morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
                    then
                      Base injection morphism:
                      From: Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
                      To:   Field in z12 with defining polynomial x^3 + (1 + (4*z2 + 2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base

            sage: L.defining_morphism(GF(5))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field of size 5
        """
        base = self._check_base(base)
        return self.coerce_map_from(base)

    def _an_element_(self):
        r"""
        Return an element of this extension.

        TESTS::

            sage: E = QQ.over(ZZ)
            sage: x = E.an_element()  # indirect doctest
            sage: x
            1/2
            sage: x.parent()
            Rational Field over its base

        """
        elt = self._backend.an_element()
        return self.element_class(self, elt)

    def gens(self, base=None):
        r"""
        Return the generators of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``); if omitted,
          use the base of this extension

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: K.gens()
            (a,)
            sage: L.<b> = GF(5^4).over(K)
            sage: L.gens()
            (b,)
            sage: L.gens(GF(5))
            (b, a)

            sage: S.<x> = QQ[]
            sage: T.<y> = S[]
            sage: T.over(S).gens()
            (y,)
            sage: T.over(QQ).gens()
            (y, x)
        """
        self._check_base(base)
        # TODO: This implementation is relying on the fact that there are
        # coercions between the backends that are compatible with the defining
        # morphisms on the frontend.
        return tuple([ self._from_backend_morphism(x) for x in _generators(self._backend, backend_parent(self._base)) ])

    def ngens(self, base=None):
        r"""
        Return the number of generators of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: K = GF(5^2).over()   # over GF(5)
            sage: K.gens()
            (z2,)
            sage: K.ngens()
            1

            sage: L = GF(5^4).over(K)
            sage: L.gens(GF(5))
            (z4, z2)
            sage: L.ngens(GF(5))
            2
        """
        return len(self.gens(base))

    def gen(self, i=0):
        r"""
        Return the first generator of this extension.

        EXAMPLES::

            sage: K = GF(5^2).over()   # over GF(5)
            sage: x =K.gen(); x
            z2

        Observe that the generator lives in the extension::

            sage: x.parent()
            Field in z2 with defining polynomial x^2 + 4*x + 2 over its base
            sage: x.parent() is K
            True
        """
        return self.gens()[i]

    def random_element(self):
        r"""
        Return a random element in this extension.

        EXAMPLES::

            sage: K = GF(5^2).over()   # over GF(5)
            sage: x = K.random_element(); x   # random
            3 + z2

            sage: x.parent()
            Field in z2 with defining polynomial x^2 + 4*x + 2 over its base
            sage: x.parent() is K
            True
        """
        elt = self._backend.random_element()
        return self.element_class(self, elt)

    def _factor_univariate_polynomial(self, f, **kwds):
        """
        Factor a univariate polynomial using code for the backend

        EXAMPLES:

        We need to disable randomness since the defining polynomial of finite
        field extensions is determined randomly, and this affects the output
        below::

            sage: set_random_seed(0)

        ::

            sage: K.<a> = GF(625, base=GF(25))
            sage: R.<x> = K[]
            sage: (x^2 - a^14).factor() # indirect doctest
            (x + 2 + a + 3*a^2 + a^3) * (x + 3 - a + 2*a^2 - a^3)
            sage: a^7
            3 - a + 2*a^2 - a^3
        """
        from sage.structure.factorization import Factorization
        F = f.change_ring(self._to_backend_morphism).factor()
        return Factorization([(g.change_ring(self._from_backend_morphism), e) for (g, e) in F], unit=self._from_backend_morphism(F.unit()))

    def degree_over(self, base=None):
        r"""
        Return the degree of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: F = GF(5^2)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(K)

            sage: K.degree_over(F)
            2
            sage: L.degree_over(K)
            3
            sage: L.degree_over(F)
            6

        If ``base`` is omitted, the degree is computed over the base
        of the extension::

            sage: K.degree_over()
            2
            sage: L.degree_over()
            3

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: K.degree_over(GF(5))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field of size 5
        """
        base = self._check_base(base)
        return self._degree_over(base)

    def _degree_over(self, base):
        r"""
        Return the degree of this extension over ``base``.

        Should be implemented in subclasses.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        TESTS::

            sage: A.<a> = QQ.extension(x^2 - 2)
            sage: B.<b> = QQ.extension(x^6 - 2)
            sage: f = A.hom([b^3])
            sage: E = B.over(f)
            sage: E.degree_over()  # indirect doctest
            3
        """
        if base is self:
            return ZZ(1)
        raise NotImplementedError("degree is not implemented (and maybe not defined) for this extension")

    def degree(self, base):
        r"""
        Return the degree of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        EXAMPLES::

            sage: A = GF(5^4).over(GF(5^2))
            sage: B = GF(5^12).over(A)

            sage: A.degree(GF(5^2))
            2
            sage: B.degree(A)
            3
            sage: B.degree(GF(5^2))
            6

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: A.degree(GF(5))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field of size 5

        .. SEEALSO::

            :meth:`relative_degree`, :meth:`absolute_degree`
        """
        return self.degree_over(base)

    def relative_degree(self):
        r"""
        Return the degree of this extension over its base

        EXAMPLES::

            sage: A = GF(5^4).over(GF(5^2))
            sage: A.relative_degree()
            2

        .. SEEALSO::

            :meth:`degree`, :meth:`absolute_degree`
        """
        return self._degree_over(self._base)

    def absolute_degree(self):
        r"""
        Return the degree of this extension over its absolute base

        EXAMPLES::

            sage: A = GF(5^4).over(GF(5^2))
            sage: B = GF(5^12).over(A)

            sage: A.absolute_degree()
            2
            sage: B.absolute_degree()
            6

        .. SEEALSO::

            :meth:`degree`, :meth:`relative_degree`
        """
        return self._degree_over(self.absolute_base())

    def is_finite_over(self, base=None):
        r"""
        Return whether or not this extension is finite over ``base`` (as a module).

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: K = GF(5^2).over()  # over GF(5)
            sage: L = GF(5^4).over(K)

            sage: L.is_finite_over(K)
            True
            sage: L.is_finite_over(GF(5))
            True

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: L.is_finite_over()
            True
        """
        base = self._check_base(base)
        if base is self:
            return True
        try:
            return self._is_finite_over(base)
        except NotImplementedError:
            pass
        b = self._base
        while b is not base:
            try:
                if self._is_finite_over(b) and b.is_finite_over(base):
                    return True
            except NotImplementedError:
                pass
            b = b._base
        raise NotImplementedError

    def _is_finite_over(self, base):
        r"""
        Return whether or not this extension is finite over ``base``.

        Should be implemented in subclasses.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: K = GF(5^2).over()  # over GF(5)
            sage: K.is_finite_over()  # indirect doctest
            True
        """
        raise NotImplementedError

    def is_free_over(self, base=None):
        r"""
        Return ``True`` if this extension is free (as a module)
        over ``base``

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: K = GF(5^2).over()  # over GF(5)
            sage: L = GF(5^4).over(K)

            sage: L.is_free_over(K)
            True
            sage: L.is_free_over(GF(5))
            True

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: L.is_free_over()
            True
        """
        base = self._check_base(base)
        if base is self or base.is_field():
            return True
        try:
            return self._is_free_over(base)
        except NotImplementedError:
            pass
        b = self._base
        while b is not base:
            try:
                if self._is_free_over(b) and b.is_free_over(base):
                    return True
            except NotImplementedError:
                pass
            b = b._base
        raise NotImplementedError

    def _is_free_over(self, base):
        r"""
        Return whether or not this extension is finite over ``base``.

        Should be implemented in subclasses.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: K = GF(5^2).over()  # over GF(5)
            sage: K.is_free_over()  # indirect doctest
            True
        """
        raise NotImplementedError

    def is_field(self, proof=True):
        r"""
        Return whether or not this extension is a field.

        INPUT:

        - ``proof`` -- a boolean (default: ``False``)

        EXAMPLES::

            sage: K = GF(5^5).over()  # over GF(5)
            sage: K.is_field()
            True

            sage: S.<x> = QQ[]
            sage: A = S.over(QQ)
            sage: A.is_field()
            False

            sage: B = A.fraction_field()
            sage: B.is_field()
            True
        """
        return self._backend.is_field(proof=proof)

    @cached_method
    def fraction_field(self, extend_base=False):
        r"""
        Return the fraction field of this extension.

        INPUT:

        - ``extend_base`` -- a boolean (default: ``False``);

        If ``extend_base`` is ``False``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/L/K`, except
        if `L` is already a field in which case the fraction field
        of `L/K` is `L/K` itself.

        If ``extend_base`` is ``True``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/\textrm{Frac}(K)`
        (provided that the defining morphism extends to the fraction
        fields, i.e. is injective).

        EXAMPLES::

            sage: A.<a> = ZZ.extension(x^2 - 5)
            sage: OK = A.over()   # over ZZ
            sage: OK
            Order in Number Field in a with defining polynomial x^2 - 5 over its base

            sage: K1 = OK.fraction_field()
            sage: K1
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K1.bases()
            [Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base,
             Order in Number Field in a with defining polynomial x^2 - 5 over its base,
             Integer Ring]

            sage: K2 = OK.fraction_field(extend_base=True)
            sage: K2
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K2.bases()
            [Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base,
             Rational Field]

        Note that there is no coercion between `K_1` and `K_2`::

            sage: K1.has_coerce_map_from(K2)
            False
            sage: K2.has_coerce_map_from(K1)
            False

        We check that when the extension is a field, its fraction field does not change::

            sage: K1.fraction_field() is K1
            True
            sage: K2.fraction_field() is K2
            True

        TESTS::

            sage: A = GF(5).over(ZZ)
            sage: A.fraction_field(extend_base=True)
            Traceback (most recent call last):
            ...
            ValueError: the morphism is not injective
        """
        defining_morphism = self._defining_morphism_fraction_field(extend_base)
        if defining_morphism is None:
            return self
        ring = defining_morphism.codomain()
        constructor = RingExtensionFractionField, {'ring': self}
        return RingExtension(ring, defining_morphism, canonical_backend=self._canonical_backend, constructors=[constructor])

    def _defining_morphism_fraction_field(self, extend_base):
        r"""
        Return the defining morphism of the fraction field of this extension.

        This is a helper function.

        INPUT:

        - ``extend_base`` -- a boolean (default: ``False``); see
          :meth:`fraction_field` for more information

        TESTS::

            sage: K = GF(5^2).over()
            sage: K.fraction_field()  # indirect doctest
            Field in z2 with defining polynomial x^2 + 4*x + 2 over its base

            sage: K = QQ.over(ZZ)
            sage: K.fraction_field(extend_base=True)  # indirect doctest
            Rational Field over its base
        """
        if extend_base:
            defining_morphism = self._backend_defining_morphism.extend_to_fraction_field()
            if isinstance(self._base, RingExtension_generic):
                base = self._base.fraction_field(extend_base)
                backend = defining_morphism.codomain()
                bHb = base.Hom(backend)
                defining_morphism = bHb.__make_element_class__(RingExtensionHomomorphism)(bHb, defining_morphism)
        else:
            if self.is_field():
                defining_morphism = None
            else:
                backend = self._backend.fraction_field()
                sHb = self.Hom(backend)
                defining_morphism = sHb.__make_element_class__(RingExtensionHomomorphism)(sHb, backend.coerce_map_from(self._backend))
        return defining_morphism

    def _Hom_(self, codomain, category):
        r"""
        Return the homset from this extension of ``codomain`` is the category ``category``.

        INPUT:

        - ``codomain`` -- a parent

        - ``category`` -- a subcategory of the category of rings

        EXAMPLES::

            sage: F = GF(5^2)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(F)

            sage: K.Hom(L)  # indirect doctest
            Set of Homomorphisms from Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
            to Field in z12 with defining polynomial x^6 + (4*z2 + 3)*x^5 + x^4 + (3*z2 + 1)*x^3 + x^2 + (4*z2 + 1)*x + z2 over its base

            sage: K.Hom(L, category=Sets())
            Set of Morphisms from Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
            to Field in z12 with defining polynomial x^6 + (4*z2 + 3)*x^5 + x^4 + (3*z2 + 1)*x^3 + x^2 + (4*z2 + 1)*x + z2 over its base
            in Category of sets

        """
        from sage.rings.ring_extension_homset import RingExtensionHomset
        if category.is_subcategory(CommutativeRings()):
            return RingExtensionHomset(self, codomain, category)
        raise TypeError("category must be a subcategory of rings")

    def hom(self, im_gens, codomain=None, base_map=None, category=None, check=True):
        r"""
        Return the unique homomorphism from this extension to
        ``codomain`` that sends ``self.gens()`` to the entries
        of ``im_gens`` and induces the map ``base_map`` on the
        base ring.

        INPUT:

        - ``im_gens`` -- the images of the generators of this extension

        - ``codomain`` -- the codomain of the homomorphism; if omitted, it
          is set to the smallest parent containing all the entries of ``im_gens``

        - ``base_map`` -- a map from one of the bases of this extension into
          something that coerces into the codomain; if omitted, coercion maps
          are used

        - ``category`` -- the category of the resulting morphism

        - ``check`` -- a boolean (default: ``True``); whether to verify that the
          images of generators extend to define a map (using only canonical coercions)

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()    # over GF(5)
            sage: L.<b> = GF(5^6).over(K)

        We define (by hand) the relative Frobenius endomorphism of the extension `L/K`::

            sage: L.hom([b^25])
            Ring endomorphism of Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: b |--> 2 + 2*a*b + (2 - a)*b^2

        Defining the absolute Frobenius of `L` can be bit more complicated
        because it is not a homomorphism of `K`-algebras.
        In this case it works since `b` generates `GF(5^6)` over its base `GF(5)`::

            sage: L.hom([b^5])
            Ring endomorphism of Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: b |--> (-1 + a) + (1 + 2*a)*b + a*b^2

        However, it can fail in other cases, if we choose a poor generator for example::

            sage: l = GF(5^6); z6 = l.gen()
            sage: z3 = 4*z6^5 + z6^4 + 4*z6^3 + 4*z6^2 + 3*z6 + 3
            sage: z3.minpoly()
            x^3 + 3*x + 3
            sage: M = GF(5^12); img = z3.minpoly().any_root(M)
            sage: L.<b> = GF(5^6).over(K, gens=[z3])
            sage: b.minpoly()
            x^3 + 3*x + 3
            sage: L.hom([img])
            Traceback (most recent call last):
            ...
            ValueError: There is no coercion from base ring to codomain; try specifying base_map

        What we need is to specify a base map which is compatible with our choice of img::

            sage: imga = a.minpoly().any_root(M)
            sage: base_map = K.hom([imga])
            sage: phi = L.hom([img], base_map=base_map); phi
            Ring morphism:
              From: Field in b with defining polynomial x^3 + 3*x + 3 over its base
              To:   Finite Field in z12 of size 5^12
              Defn: b |--> 4*z12^10 + 3*z12^9 + z12^8 + 4*z12^7 + z12^6 + 2*z12^5 + 3*z12^4 + 4*z12^3 + 3*z12^2 + 3
                    with map on base ring:
                    a |--> 4*z12^11 + 4*z12^10 + 2*z12^9 + 3*z12^7 + 4*z12^6 + 4*z12^5 + 3*z12^4 + 3*z12^3 + 4*z12^2 + 4*z12 + 3

        As a shortcut, we may use the following construction::

            sage: psi = L.hom([img, imga])
            sage: phi == psi
            True
        """
        if isinstance(im_gens, Map):
            if base_map is not None:
                raise ValueError("base_map cannot be set when passing in a morphism")
            if codomain is not None and codomain is not im_gens.codomain():
                im_gens = im_gens.extend_codomain(codomain)
            if im_gens.domain() is not self:
                im_gens = im_gens.extend_domain(domain)
            if category is not None and category is not im_gens.category_for():
                parent = self.Hom(im_gens.codomain(), category=category)
                im_gens = parent(im_gens)
            return im_gens
        elif isinstance(im_gens, (list, tuple)):
            if codomain is None:
                from sage.structure.sequence import Sequence
                codomain = Sequence(im_gens).universe()
            n = self.ngens()
            if n > len(im_gens):
                raise ValueError("Not enough images provided")
            elif n < len(im_gens):
                base_map = self.base_ring().hom(im_gens[n:], codomain=codomain, base_map=base_map, category=category, check=check)
                im_gens = im_gens[:n]
            im_gens = [codomain(g) for g in im_gens]
        parent = self.Hom(codomain, category=category)
        return parent(im_gens, base_map=base_map, check=check)

    def characteristic(self):
        r"""
        Return the characteristic of the extension as a ring.

        OUTPUT:

        A prime number or zero.

        EXAMPLES::

            sage: F = GF(5^2).over()   # over GF(5)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(K)
            sage: F.characteristic()
            5
            sage: K.characteristic()
            5
            sage: L.characteristic()
            5

        ::

            sage: F = RR.over(ZZ)
            sage: F.characteristic()
            0

        ::

            sage: F = GF(11)
            sage: A.<x> = F[]
            sage: K = Frac(F).over(F)
            sage: K.characteristic()
            11

        ::

            sage: E = GF(7).over(ZZ)
            sage: E.characteristic()
            7

        TESTS:

        Ensure issue :trac:`34692` is fixed::

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: k = Frac(FqX)
            sage: K = k.over(FqX)
            sage: K.frobenius_endomorphism()
            Frobenius endomorphism x |--> x^11 of Fraction Field of Univariate Polynomial Ring in X over Finite Field of size 11 over its base
        """
        return self._backend.characteristic()

# Fraction fields
#################

class RingExtensionFractionField(RingExtension_generic):
    r"""
    A class for ring extensions of the form `\textrm{Frac}(A)/A`.

    TESTS::

        sage: Z = ZZ.over()   # over ZZ itself
        sage: Q = Z.fraction_field()
        sage: Q
        Fraction Field of Integer Ring over its base

        sage: type(Q)
        <class 'sage.rings.ring_extension.RingExtensionFractionField_with_category'>

        sage: TestSuite(Q).run()

    """
    Element = RingExtensionFractionFieldElement

    def __init__(self, defining_morphism, ring=None, **kwargs):
        r"""
        Initialize this ring extension.

        INPUT:

        - ``defining_morphism`` -- a ring homomorphism

        - ``ring`` -- the commutative ring whose fraction field is this
          extension

        TESTS::

            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: OK = A.over()
            sage: K = OK.fraction_field()
            sage: K
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 2 over its base

            sage: TestSuite(K).run()

        """
        RingExtension_generic.__init__(self, defining_morphism, **kwargs)
        self._ring = ring or self.base_ring()

    def ring(self):
        r"""
        Return the ring whose fraction field is this extension.

        EXAMPLES::

            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: OK = A.over()
            sage: K = OK.fraction_field()
            sage: K
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 2 over its base

            sage: K.ring()
            Order in Number Field in a with defining polynomial x^2 - 2 over its base
            sage: K.ring() is OK
            True
        """
        return self._ring

    def _repr_topring(self, **options):
        r"""
        Return a string representation of top ring of this extension.

        EXAMPLES::

            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: OK = A.over()
            sage: K = OK.fraction_field()

            sage: K._repr_topring()
            'Fraction Field of Order in Number Field in a with defining polynomial x^2 - 2'
        """
        if isinstance(self._ring, RingExtension_generic):
            sr = self._ring._repr_topring(**options)
        else:
            sr = str(self._ring)
        if self._ring in Fields():
            return sr
        else:
            return "Fraction Field of %s" % sr

    def _latex_topring(self, **options):
        r"""
        Return a LaTeX representation of top ring of this extension.

        EXAMPLES::

            sage: Z = ZZ.over()
            sage: Q = Z.fraction_field()

            sage: Q._latex_topring()
            '\\mathrm{Frac}(\\Bold{Z})'
        """
        if self._ring in Fields():
            return self._ring._latex_topring(**options)
        else:
            return "\\mathrm{Frac}(%s)" % latex(self._ring)


# Finite free extensions
########################

class RingExtensionWithBasis(RingExtension_generic):
    """
    A class for finite free ring extensions equipped
    with a basis.

    TESTS::

        sage: E = GF(5^4).over(GF(5^2))
        sage: E
        Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base

        sage: TestSuite(E).run()
    """
    Element = RingExtensionWithBasisElement

    def __init__(self, defining_morphism, basis, names, check=True, **kwargs):
        r"""
        Initialize this ring extension.

        INPUT:

        - ``defining_morphism`` -- a ring homomorphism

        - ``basis`` -- a tuple of elements in this extension

        - ``names`` -- a tuple of strings or ``None`` (default: ``None``);
          the way the elements of the basis are printed

        - ``check`` -- a boolean (default: ``True``); whether to check if
          ``basis`` is indeed a basis

        TESTS::

            sage: K.<a> = QQ.extension(x^3 - 2)
            sage: E = K.over()
            sage: E
            Field in a with defining polynomial x^3 - 2 over its base

            sage: TestSuite(E).run()
        """
        RingExtension_generic.__init__(self, defining_morphism, **kwargs)
        self._basis = [ self._from_backend_morphism(b) for b in basis ]
        if len(names) != len(self._basis):
            raise ValueError(f"unexpected number of names for basis elements (expected {len(basis)}, got {len(names)})")
        self._basis_names = names
        self._basis_latex_names = [ latex_variable_name(name) for name in names ]
        self._names = tuple(names)
        if check:
            try:
                _ = self.free_module(map=True)
            except (ZeroDivisionError, ArithmeticError):
                raise ValueError("the given family is not a basis")
        if 'base' not in self._print_options:
            self._print_options['base'] = self._base

    def _print_option_base(self, base):
        r"""
        Return a normalized form of the print option ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        TESTS::

            sage: F = GF(5)
            sage: K = GF(5^2).over(F)
            sage: L = GF(5^4).over(K)

            sage: L._print_option_base(F) is F
            True
            sage: L._print_option_base(K) is K
            True

            sage: L._print_option_base(None) is K
            True

            sage: L._print_option_base(GF(5^2)) is K
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z2 of size 5^2

            sage: L._print_option_base(L)
            Traceback (most recent call last):
            ...
            ValueError: base must be strict

            sage: K._print_option_base(L)
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Field in z4 with defining polynomial x^2 + (3 - z2)*x + z2 over its base

        """
        if 'print_elements_as' in self._print_options:
            raise NotImplementedError("printing is handled by an external function or another parent")
        base = self._check_base(base)
        if base is self:
            raise ValueError("base must be strict")
        b = self._base
        while b is not base:
            if not isinstance(b, RingExtensionWithBasis):
                raise NotImplementedError
            b = b.base_ring()
        return base

    def _degree_over(self, base):
        r"""
        Return the degree of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        TESTS::

            sage: A.<a> = QQ.extension(x^2 - 2)
            sage: B.<b> = QQ.extension(x^6 - 2)
            sage: f = A.hom([b^3])
            sage: E = B.over(f)
            sage: E.degree_over()  # indirect doctest
            3
        """
        if base is self:
            return ZZ(1)
        elif base is self._base:
            return len(self._basis)
        else:
            return len(self._basis) * self._base._degree_over(base)

    def dimension(self, base=None):
        """
        The dimension of this ring extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: A.<a> = QQ.extension(x^2 - 2)
            sage: B.<b> = QQ.extension(x^6 - 2)
            sage: f = A.hom([b^3])
            sage: E = B.over(f)
            sage: E.dimension()
            3
        """
        if base is None:
            base = self._base
        return self._degree_over(base)

    def rank(self, base=None):
        """
        The dimension of this ring extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: A.<a> = QQ.extension(x^2 - 2)
            sage: B.<b> = QQ.extension(x^6 - 2)
            sage: f = A.hom([b^3])
            sage: E = B.over(f)
            sage: E.rank()
            3
        """
        return self.dimension()

    def _is_finite_over(self, base):
        r"""
        Return whether or not this extension is finite over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: K = GF(5^2).over()  # over GF(5)
            sage: K.is_finite_over()  # indirect doctest
            True
        """
        if base is self or base is self._base:
            return True
        return self._base._is_finite_over(base)

    def _is_free_over(self, base):
        r"""
        Return whether or not this extension is free over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: K = GF(5^2).over()  # over GF(5)
            sage: K.is_free_over()    # indirect doctest
            True
        """
        if base is self or base is self._base:
            return True
        return self._base._is_free_over(base)

    def basis_over(self, base=None):
        r"""
        Return a basis of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        EXAMPLES::

            sage: F.<a> = GF(5^2).over()  # over GF(5)
            sage: K.<b> = GF(5^4).over(F)
            sage: L.<c> = GF(5^12).over(K)

            sage: L.basis_over(K)
            [1, c, c^2]

            sage: L.basis_over(F)
            [1, b, c, b*c, c^2, b*c^2]

            sage: L.basis_over(GF(5))
            [1, a, b, a*b, c, a*c, b*c, a*b*c, c^2, a*c^2, b*c^2, a*b*c^2]

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: L.basis_over()
            [1, c, c^2]

            sage: K.basis_over()
            [1, b]

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: L.degree_over(GF(5^6))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z6 of size 5^6
        """
        base = self._check_base(base)
        return self._basis_over(base)

    def _basis_over(self, base):
        r"""
        Return a basis of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: A.<a> = QQ.extension(x^3 - 2)
            sage: K.<u> = A.over()
            sage: K.basis_over() # indirect doctest
            [1, u, u^2]
        """
        if base is self:
            return [ self.one() ]
        elif base is self._base:
            return self._basis[:]
        else:
            b = self._base._basis_over(base)
            return [ x*y for x in self._basis for y in b ]

    def free_module(self, base=None, map=True):
        r"""
        Return a free module V over ``base`` which is isomorphic to
        this ring

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        - ``map`` -- boolean (default ``True``); whether to return
          isomorphisms between this ring and V

        OUTPUT:

        - A finite-rank free module V over ``base``

        - The isomorphism from V to this ring corresponding to the
          basis output by the method :meth:`basis_over`
          (only included if ``map`` is ``True``)

        - The reverse isomorphism of the isomorphism above
          (only included if ``map`` is ``True``)

        EXAMPLES::

            sage: F = GF(11)
            sage: K.<a> = GF(11^2).over()
            sage: L.<b> = GF(11^6).over(K)

        Forgetting a part of the multiplicative structure, the field L
        can be viewed as a vector space of dimension 3 over K, equipped
        with a distinguished basis, namely `(1, b, b^2)`::

            sage: V, i, j = L.free_module(K)
            sage: V
            Vector space of dimension 3 over Field in a with defining polynomial x^2 + 7*x + 2 over its base
            sage: i
            Generic map:
              From: Vector space of dimension 3 over Field in a with defining polynomial x^2 + 7*x + 2 over its base
              To:   Field in b with defining polynomial x^3 + (7 + 2*a)*x^2 + (2 - a)*x - a over its base
            sage: j
            Generic map:
              From: Field in b with defining polynomial x^3 + (7 + 2*a)*x^2 + (2 - a)*x - a over its base
              To:   Vector space of dimension 3 over Field in a with defining polynomial x^2 + 7*x + 2 over its base

            sage: j(b)
            (0, 1, 0)
            sage: i((1, a, a+1))
            1 + a*b + (1 + a)*b^2

        Similarly, one can view L as a F-vector space of dimension 6::

            sage: V, i, j, = L.free_module(F)
            sage: V
            Vector space of dimension 6 over Finite Field of size 11

        In this case, the isomorphisms between `V` and `L` are given by the
        basis `(1, a, b, ab, b^2, ab^2)`:

            sage: j(a*b)
            (0, 0, 0, 1, 0, 0)
            sage: i((1,2,3,4,5,6))
            (1 + 2*a) + (3 + 4*a)*b + (5 + 6*a)*b^2

        When ``base`` is omitted, the default is the base of this extension::

            sage: L.free_module(map=False)
            Vector space of dimension 3 over Field in a with defining polynomial x^2 + 7*x + 2 over its base

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: L.degree(GF(11^3))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z3 of size 11^3

        """
        base = self._check_base(base)
        return self._free_module(base, map)

    @cached_method
    def _free_module(self, base, map):
        r"""
        Return a free module V over ``base`` which is isomorphic to
        this ring

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        - ``map`` -- boolean (default ``True``); whether to return
          isomorphisms between this ring and V

        OUTPUT:

        - A finite-rank free module V over ``base``

        - The isomorphism from V to this ring corresponding to the
          basis output by the method :meth:`basis_over`
          (only included if ``map`` is ``True``)

        - The reverse isomorphism of the isomorphism above
          (only included if ``map`` is ``True``)

        TESTS::

            sage: K = GF(7^5).over()
            sage: L = GF(7^15).over(K)
            sage: for base in L.bases():
            ....:     V, i, j = L.free_module(base) # indirect doctest
            ....:     assert([ i(v) for v in V.basis() ] == L.basis_over(base))
            ....:     assert([ j(x) for x in L.basis_over(base) ] == V.basis())

        """
        d = self._degree_over(base)
        if map:
            V = base**d
            VHs = V.Hom(self)
            sHV = self.Hom(V)
            return (V,
                    VHs.__make_element_class__(MapFreeModuleToRelativeRing)(VHs),
                    sHV.__make_element_class__(MapRelativeRingToFreeModule)(sHV))
        else:
            return base**d

    def coordinates(self, x, base=None):
        """
        The coordinates of an element in terms of the basis.

        INPUT:

        - ``x`` -- an element of this ring extension

        - ``base`` -- a commutative ring, one of the (iterative) bases of this ring extension

        EXAMPLES::

            sage: K = GF(7^5).over()
            sage: L = GF(7^15).over(K)
            sage: x = L.random_element()
            sage: coords = L.coordinates(x)
            sage: basis = L.basis_over()
            sage: x == sum(c*v for (c, v) in zip(coords, basis))
            True
        """
        V, from_V, to_V = self.free_module(base, map=True)
        return list(to_V(x))

    @cached_method
    def fraction_field(self, extend_base=False):
        r"""
        Return the fraction field of this extension.

        INPUT:

        - ``extend_base`` -- a boolean (default: ``False``);

        If ``extend_base`` is ``False``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/L/K`, except
        is `L` is already a field in which base the fraction field
        of `L/K` is `L/K` itself.

        If ``extend_base`` is ``True``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/\textrm{Frac}(K)`
        (provided that the defining morphism extends to the fraction
        fields, i.e. is injective).

        EXAMPLES::

            sage: A.<a> = ZZ.extension(x^2 - 5)
            sage: OK = A.over()   # over ZZ
            sage: OK
            Order in Number Field in a with defining polynomial x^2 - 5 over its base

            sage: K1 = OK.fraction_field()
            sage: K1
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K1.bases()
            [Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base,
             Order in Number Field in a with defining polynomial x^2 - 5 over its base,
             Integer Ring]

            sage: K2 = OK.fraction_field(extend_base=True)
            sage: K2
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K2.bases()
            [Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base,
             Rational Field]

        Note that there is no coercion map between `K_1` and `K_2`::

            sage: K1.has_coerce_map_from(K2)
            False
            sage: K2.has_coerce_map_from(K1)
            False

        We check that when the extension is a field, its fraction field does not change::

            sage: K1.fraction_field() is K1
            True
            sage: K2.fraction_field() is K2
            True

        TESTS::

            sage: A = GF(5).over(ZZ)
            sage: A.fraction_field(extend_base=True)
            Traceback (most recent call last):
            ...
            ValueError: the morphism is not injective
        """
        defining_morphism = self._defining_morphism_fraction_field(extend_base)
        if defining_morphism is None:
            return self
        if extend_base:
            basis = self._basis
            names = self._basis_names
            constructor = RingExtensionWithBasis
            kwargs = { 'basis': basis, 'names': names, 'check': False }
        else:
            gen = names = None
            constructor = RingExtensionFractionField
            kwargs = { 'print_options': {'print_elements_as': self.fraction_field(extend_base=True)} }
        ring = defining_morphism.codomain()
        return RingExtension(ring, defining_morphism, gen=gen, names=names, canonical_backend=self._canonical_backend, constructors=[(constructor, kwargs)])


class RingExtensionWithGen(RingExtensionWithBasis):
    """
    A class for finite free ring extensions generated by
    a single element

    TESTS::

        sage: A.<a> = QQ.extension(x^3 - 7)
        sage: K = A.over()

        sage: type(K)
        <class 'sage.rings.ring_extension.RingExtensionWithGen_with_category'>

        sage: TestSuite(K).run()

    """
    def __init__(self, defining_morphism, gen, names, check=True, **kwargs):
        r"""
        Initialize this ring extension.

        INPUT:

        - ``defining_morphism`` -- a ring homomorphism

        - ``gen`` -- a generator of this extension

        - ``names`` -- a tuple of strings or ``None`` (default: ``None``);
          the way the elements of the basis are printed

        - ``check`` -- a boolean (default: ``True``); whether to check if
          ``gen`` is indeed a generator

        TESTS::

            sage: K.<a> = QQ.extension(x^3 + 3*x + 1)
            sage: E = K.over()
            sage: E
            Field in a with defining polynomial x^3 + 3*x + 1 over its base

            sage: TestSuite(E).run()
        """
        backend_base = backend_parent(defining_morphism.domain())
        _, deg_domain, deg_codomain = common_base(backend_base, defining_morphism.codomain(), True)
        degree = deg_codomain // deg_domain

        if len(names) != 1:
            raise ValueError(f"expected exactly one name for the generators of this ring extension but got {names}")

        name = names[0]
        latex_name = latex_variable_name(name)

        if degree == 1:
            basis_names = [""]
            basis_latex_names = [""]
        else:
            basis_names = [ "", name ] + [ f"{name}^{i}" for i in range(2, degree) ]
            basis_latex_names = [ "", latex_name] + [ f"{latex_name}^{i}" for i in range(2, degree) ]

        basis = [gen ** i for i in range(degree)]

        RingExtensionWithBasis.__init__(self, defining_morphism, basis, basis_names, check, **kwargs)

        # Override the names that RingExtensionWithBasis.__init__ set.
        self._name = name
        self._names = (name,)
        self._latex_names = (latex_name,)
        self._basis_latex_names = basis_latex_names

        self._gen = self._backend(gen)

    def _repr_topring(self, **options):
        r"""
        Return a string representation of top ring of this extension.

        EXAMPLES::

            sage: K.<a> = GF(5^3).over()
            sage: K._repr_topring()
            'Field in a with defining polynomial x^3 + 3*x + 3'

            sage: L.<b> = GF(5^9).over(K)
            sage: L._repr_topring()
            'Field in b with defining polynomial x^3 + (1 + 3*a^2)*x^2 + (3 + 2*a + 2*a^2)*x - a'
        """
        return "%s in %s with defining polynomial %s" % (self._type, self._name, self.modulus())

    def _latex_topring(self):
        r"""
        Return a LaTeX representation of top ring of this extension.

        EXAMPLES::

            sage: K.<a> = GF(5^3).over()
            sage: K._latex_topring()
            '\\Bold{F}_{5}[a]'

            sage: L.<b> = GF(5^9).over(K)
            sage: L._latex_topring()
            '\\Bold{F}_{5}[a][b]'
        """
        if isinstance(self._base, RingExtension_generic):
            return "%s[%s]" % (self._base._latex_topring(), self.latex_variable_names()[0])
        else:
            return "%s[%s]" % (latex(self._base), self.latex_variable_names()[0])

    def modulus(self, var='x'):
        r"""
        Return the defining polynomial of this extension, that is the
        minimal polynomial of the given generator of this extension.

        INPUT:

        - ``var`` -- a variable name (default: ``x``)

        EXAMPLES::

            sage: K.<u> = GF(7^10).over(GF(7^2))
            sage: K
            Field in u with defining polynomial x^5 + (6*z2 + 4)*x^4 + (3*z2 + 5)*x^3 + (2*z2 + 2)*x^2 + 4*x + 6*z2 over its base

            sage: P = K.modulus(); P
            x^5 + (6*z2 + 4)*x^4 + (3*z2 + 5)*x^3 + (2*z2 + 2)*x^2 + 4*x + 6*z2
            sage: P(u)
            0

        We can use a different variable name::

            sage: K.modulus('y')
            y^5 + (6*z2 + 4)*y^4 + (3*z2 + 5)*y^3 + (2*z2 + 2)*y^2 + 4*y + 6*z2
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        _, _, j = self.free_module(map=True)
        d = self.relative_degree()
        coeffs = [ -c for c in j(self._from_backend_morphism(self._gen**d)) ] + [ 1 ]
        S = PolynomialRing(self._base, name=var)
        return S(coeffs)

    def gens(self, base=None, check=True):
        r"""
        Return the generators of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)
        - ``check`` -- whether to check the validity of ``base``

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: K.gens()
            (a,)

            sage: L.<b> = GF(5^4).over(K)
            sage: L.gens()
            (b,)
            sage: L.gens(GF(5))
            (b, a)

        TESTS:

        Generators are returned as elements of this ring::

            sage: all(g.parent() is L for g in L.gens(GF(5)))
            True

        ::

            sage: k = GF(3)
            sage: F.<a> = k.extension(2, absolute=False)
            sage: F.gens(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Integer Ring

        """
        if check:
            base = self._check_base(base)

        if base is self:
            return tuple()

        gens = (self._from_backend_morphism(self._gen),)

        if base is self.base():
            return gens

        return gens + tuple([self(gen) for gen in self.base().gens(base=base, check=False)])

    def gen(self, i=0):
        r"""
        Return the first generator of this extension.

        EXAMPLES::

            sage: K = GF(5^2).over()   # over GF(5)
            sage: x =K.gen(); x
            z2

        """
        if i != 0:
            raise ValueError("ring has only a single generator, use gens(base) to get other relative generators")
        return self._from_backend_morphism(self._gen)

    @cached_method
    def fraction_field(self, extend_base=False):
        r"""
        Return the fraction field of this extension.

        INPUT:

        - ``extend_base`` -- a boolean (default: ``False``);

        If ``extend_base`` is ``False``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/L/K`, except
        is `L` is already a field in which base the fraction field
        of `L/K` is `L/K` itself.

        If ``extend_base`` is ``True``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/\textrm{Frac}(K)`
        (provided that the defining morphism extends to the fraction
        fields, i.e. is injective).

        EXAMPLES::

            sage: A.<a> = ZZ.extension(x^2 - 5)
            sage: OK = A.over()   # over ZZ
            sage: OK
            Order in Number Field in a with defining polynomial x^2 - 5 over its base

            sage: K1 = OK.fraction_field()
            sage: K1
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K1.bases()
            [Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base,
             Order in Number Field in a with defining polynomial x^2 - 5 over its base,
             Integer Ring]

            sage: K2 = OK.fraction_field(extend_base=True)
            sage: K2
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K2.bases()
            [Fraction Field of Order in Number Field in a with defining polynomial x^2 - 5 over its base,
             Rational Field]

        Note that there is no coercion map between `K_1` and `K_2`::

            sage: K1.has_coerce_map_from(K2)
            False
            sage: K2.has_coerce_map_from(K1)
            False

        We check that when the extension is a field, its fraction field does not change::

            sage: K1.fraction_field() is K1
            True
            sage: K2.fraction_field() is K2
            True

        TESTS::

            sage: A = GF(5).over(ZZ)
            sage: A.fraction_field(extend_base=True)
            Traceback (most recent call last):
            ...
            ValueError: the morphism is not injective
        """
        defining_morphism = self._defining_morphism_fraction_field(extend_base)
        if defining_morphism is None:
            return self
        if extend_base:
            gen = self._gen
            names = self._names
            constructor = RingExtensionWithGen
            kwargs = { 'gen': gen, 'names': names, 'check': False }
        else:
            gen = names = None
            constructor = RingExtensionFractionField
            kwargs = { 'print_options': {'print_elements_as': self.fraction_field(extend_base=True)} }
        ring = defining_morphism.codomain()
        return RingExtension(ring, defining_morphism, gen=gen, names=names, canonical_backend=self._canonical_backend, constructors=[(constructor, kwargs)])
