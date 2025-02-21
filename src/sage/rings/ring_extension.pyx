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

    sage: GF(5^4).over(GF(5^2))                                                         # needs sage.rings.finite_rings
    Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base

By default, Sage reuses the canonical generator of the top ring
(here `z_4 \in \mathbf{F}_{5^4}`), together with its name. However,
the user can customize them by passing in appropriate arguments::

    sage: # needs sage.rings.finite_rings
    sage: F = GF(5^2)
    sage: k = GF(5^4)
    sage: z4 = k.gen()
    sage: K.<a> = k.over(F, gen=1-z4); K
    Field in a with defining polynomial x^2 + z2*x + 4 over its base

The base of the extension is available via the method :meth:`base` (or
equivalently :meth:`base_ring`)::

    sage: K.base()                                                                      # needs sage.rings.finite_rings
    Finite Field in z2 of size 5^2

It is also possible to build an extension on top of another extension,
obtaining this way a tower of extensions::

    sage: L.<b> = GF(5^8).over(K); L                                                    # needs sage.rings.finite_rings
    Field in b with defining polynomial x^2 + (4*z2 + 3*a)*x + 1 - a over its base
    sage: L.base()                                                                      # needs sage.rings.finite_rings
    Field in a with defining polynomial x^2 + z2*x + 4 over its base
    sage: L.base().base()                                                               # needs sage.rings.finite_rings
    Finite Field in z2 of size 5^2

The method :meth:`bases` gives access to the complete list of rings in
a tower::

    sage: L.bases()                                                                     # needs sage.rings.finite_rings
    [Field in b with defining polynomial x^2 + (4*z2 + 3*a)*x + 1 - a over its base,
     Field in a with defining polynomial x^2 + z2*x + 4 over its base,
     Finite Field in z2 of size 5^2]

Once we have constructed an extension (or a tower of extensions), we
have interesting methods attached to it. As a basic example, one can
compute a basis of the top ring over any base in the tower::

    sage: L.basis_over(K)                                                               # needs sage.rings.finite_rings
    [1, b]
    sage: L.basis_over(F)                                                               # needs sage.rings.finite_rings
    [1, a, b, a*b]

When the base is omitted, the default is the natural base of the extension::

    sage: L.basis_over()                                                                # needs sage.rings.finite_rings
    [1, b]

The method :meth:`sage.rings.ring_extension_element.RingExtensionWithBasis.vector`
computes the coordinates of an element according to the above basis::

    sage: u = a + 2*b + 3*a*b                                                           # needs sage.rings.finite_rings
    sage: u.vector()   # over K                                                         # needs sage.rings.finite_rings
    (a, 2 + 3*a)
    sage: u.vector(F)                                                                   # needs sage.rings.finite_rings
    (0, 1, 2, 3)

One can also compute traces and norms with respect to any base of the tower::

    sage: # needs sage.rings.finite_rings
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

    sage: u.minpoly()                                                                   # needs sage.rings.finite_rings
    x^2 + ((3*z2 + 4) + (3*z2 + 4)*a)*x + (z2 + 1) + (4*z2 + 2)*a
    sage: u.minpoly(F)                                                                  # needs sage.rings.finite_rings
    x^4 + (4*z2 + 4)*x^3 + x^2 + (z2 + 1)*x + 2*z2 + 2


AUTHOR:

- Xavier Caruso (2019)
"""

# ###########################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free softwGare: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************


from sage.misc.fast_methods cimport hash_by_id
from sage.misc.cachefunc import cached_method
from sage.cpython.getattr cimport AttributeErrorMessage
from sage.cpython.getattr import dir_with_other_class
from sage.misc.latex import latex, latex_variable_name

from sage.structure.factory import UniqueFactory
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.structure.category_object import normalize_names
from sage.categories.map cimport Map
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity

from sage.rings.ring_extension_element cimport (
    RingExtensionElement, RingExtensionFractionFieldElement, RingExtensionWithBasisElement)
from sage.rings.ring_extension_morphism cimport (
    RingExtensionHomomorphism, RingExtensionBackendIsomorphism, RingExtensionBackendReverseIsomorphism,
    are_equal_morphisms, MapFreeModuleToRelativeRing, MapRelativeRingToFreeModule)
from sage.rings.ring_extension_conversion cimport (
    backend_parent, backend_morphism, to_backend, from_backend)


# Helper functions
##################

def tower_bases(ring, degree):
    r"""
    Return the list of bases of ``ring`` (including itself); if
    degree is ``True``, restrict to finite extensions and return
    in addition the degree of ``ring`` over each base.

    INPUT:

    - ``ring`` -- a commutative ring

    - ``degree`` -- boolean

    EXAMPLES::

        sage: from sage.rings.ring_extension import tower_bases
        sage: S.<x> = QQ[]
        sage: T.<y> = S[]
        sage: tower_bases(T, False)
        ([Univariate Polynomial Ring in y over
           Univariate Polynomial Ring in x over Rational Field,
          Univariate Polynomial Ring in x over Rational Field,
          Rational Field],
         [])
        sage: tower_bases(T, True)
        ([Univariate Polynomial Ring in y over
           Univariate Polynomial Ring in x over Rational Field],
         [1])

        sage: K.<a> = Qq(5^2)                                                           # needs sage.rings.padics
        sage: L.<w> = K.extension(x^3 - 5)                                              # needs sage.rings.padics
        sage: tower_bases(L, True)                                                      # needs sage.rings.padics
        ([5-adic Eisenstein Extension Field in w defined by x^3 - 5 over its base field,
          5-adic Unramified Extension Field in a defined by x^2 + 4*x + 2,
          5-adic Field with capped relative precision 20],
         [1, 3, 6])
    """
    bases = []
    degrees = []
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
            if d is Infinity: break
            deg *= d
        newbase = base._base
        if newbase is base: break
        base = newbase
    return bases, degrees


def common_base(K, L, degree):
    """
    Return a common base on which ``K`` and ``L`` are defined.

    INPUT:

    - ``K`` -- a commutative ring

    - ``L`` -- a commutative ring

    - ``degree`` -- boolean; if ``True``, return the degree of
      ``K`` and ``L`` over their common base

    EXAMPLES::

        sage: from sage.rings.ring_extension import common_base

        sage: common_base(GF(5^3), GF(5^7), False)                                      # needs sage.rings.finite_rings
        Finite Field of size 5
        sage: common_base(GF(5^3), GF(5^7), True)                                       # needs sage.rings.finite_rings
        (Finite Field of size 5, 3, 7)

        sage: common_base(GF(5^3), GF(7^5), False)                                      # needs sage.rings.finite_rings
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


def generators(ring, base):
    r"""
    Return the generators of ``ring`` over ``base``.

    INPUT:

    - ``ring`` -- a commutative ring

    - ``base`` -- a commutative ring

    EXAMPLES::

        sage: from sage.rings.ring_extension import generators
        sage: S.<x> = QQ[]
        sage: T.<y> = S[]

        sage: generators(T, S)
        (y,)
        sage: generators(T, QQ)
        (y, x)
    """
    gens = tuple()
    while ring is not ring.base_ring() and (base is None or not base.has_coerce_map_from(ring)):
        gens += tuple(ring.gens())
        ring = ring.base_ring()
    if base is None:
        return gens
    else:
        return tuple([x for x in gens if x not in base])


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

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = QQ.extension(x^2 - 2)                                             # needs sage.rings.number_field
        sage: E = K.over(QQ); E                                                         # needs sage.rings.number_field
        Field in a with defining polynomial x^2 - 2 over its base

        sage: E2.<b> = K.over(QQ)                                                       # needs sage.rings.number_field
        sage: E2 is E                                                                   # needs sage.rings.number_field
        False
    """
    def create_key_and_extra_args(self, ring, defining_morphism=None, gens=None, names=None, constructors=None):
        """
        Create a key and return it together with a list of constructors
        of the object.

        INPUT:

        - ``ring`` -- a commutative ring

        - ``defining_morphism`` -- a ring homomorphism or a commutative
          ring or ``None`` (default: ``None``); the defining morphism of
          this extension or its base (if it coerces to ``ring``)

        - ``gens`` -- list of generators of this extension (over its base)
          or ``None`` (default: ``None``)

        - ``names`` -- list or a tuple of variable names or ``None``
          (default: ``None``)

        - ``constructors`` -- list of constructors; each constructor
          is a pair `(class, arguments)` where `class` is the class
          implementing the extension and `arguments` is the dictionary
          of arguments to pass in to init function

        TESTS::

            sage: from sage.rings.ring_extension import RingExtension
            sage: RingExtension.create_key_and_extra_args(QQ, ZZ)
            ((Ring morphism:
                From: Integer Ring
                To:   Rational Field
                Defn: 1 |--> 1, (), ()),
             {'constructors': [(<class 'sage.rings.ring_extension.RingExtension_generic'>,
                               {'is_backend_exposed': True,
                                'print_options': {'print_elements_as': None,
                                                  'print_parent_as': None}})]})

            sage: RingExtension.create_key_and_extra_args(GF(5^4), GF(5^2),             # needs sage.rings.finite_rings
            ....:                                         names=('a',))
            ((Ring morphism:
                From: Finite Field in z2 of size 5^2
                To:   Finite Field in z4 of size 5^4
                Defn: z2 |--> z4^3 + z4^2 + z4 + 3, (z4,), ('a',)),
             {'constructors': [(<class 'sage.rings.ring_extension.RingExtensionWithGen'>,
                               {'gen': z4, 'is_backend_exposed': True, 'names': ('a',)})]})
        """
        use_generic_constructor = True
        is_backend_exposed = True
        print_as = None

        if defining_morphism is None:
            base = ring.base_ring()
        elif isinstance(defining_morphism, Map):
            base = defining_morphism.domain()
        elif defining_morphism in CommutativeRings():
            base = defining_morphism
            defining_morphism = None
        else:
            raise TypeError("only commutative rings are accepted")

        # We compute the defining morphism
        if defining_morphism is None:
            if isinstance(base, RingExtension_generic):
                backend_base = (<RingExtension_generic>base)._backend
                if ring.has_coerce_map_from(backend_base):
                    defining_morphism = RingExtensionHomomorphism(base.Hom(ring), ring.coerce_map_from(backend_base))
            else:
                if ring.has_coerce_map_from(base):
                    defining_morphism = ring.coerce_map_from(base)
            if defining_morphism is None:
                raise ValueError("No coercion map from %s to %s" % (base,ring))
        else:
            if defining_morphism.domain() is not base:
                defining_morphism = defining_morphism.extend_domain(base)
            if defining_morphism.codomain() is not ring:
                defining_morphism = defining_morphism.extend_codomain(ring)
        if not isinstance(defining_morphism, RingExtensionHomomorphism):
            defining_morphism = RingExtensionHomomorphism(defining_morphism.parent(), defining_morphism)
        if isinstance(ring, RingExtension_generic):
            defining_morphism = backend_morphism(defining_morphism, forget='codomain')
            if (<RingExtension_generic>ring)._is_backend_exposed:
                print_as = (<RingExtension_generic>ring)._backend
            else:
                use_generic_constructor = False
                is_backend_exposed = False
            ring = (<RingExtension_generic>ring)._backend

        # We normalize other attributes
        if gens is not None:
            if not isinstance(gens, (list, tuple)):
                raise TypeError("gens must be a list or a tuple")
            gens = tuple(ring(g) for g in gens )
            if names is None:
                raise TypeError("you must specify the names of the generators")
            names = normalize_names(len(gens), names)
            use_generic_constructor = False
        else:
            gens = generators(ring, base)
            if names is None:
                try:
                    names = variable_names(ring, base)
                except NotImplementedError:
                    gens = names = None
            else:
                names = normalize_names(len(gens), names)
                use_generic_constructor = False

        # We figure out what are the best constructors
        if constructors is None:
            constructors = []
            if gens is not None and len(gens) == 1:
                constructors.append((RingExtensionWithGen,
                                     {'gen': gens[0], 'names': names,
                                      'is_backend_exposed': is_backend_exposed}))
            if use_generic_constructor:
                constructors.append((RingExtension_generic,
                                     {'print_options': {'print_parent_as': print_as,
                                                        'print_elements_as': print_as},
                                      'is_backend_exposed': is_backend_exposed}))

        # We build the key and return it
        return (defining_morphism, gens, names), {'constructors': constructors}

    def create_object(self, version, key, **extra_args):
        """
        Return the object associated to a given key.

        TESTS::

            sage: from sage.rings.ring_extension import RingExtension
            sage: key, extra_args = RingExtension.create_key_and_extra_args(QQ, ZZ)
            sage: RingExtension.create_object((8,9,0), key, **extra_args)
            Rational Field over its base
        """
        defining_morphism, gens, names = key
        constructors = extra_args['constructors']
        if len(constructors) == 0:
            raise NotImplementedError("no constructor available for this extension")
        for (constructor, kwargs) in constructors[:-1]:
            try:
                return constructor(defining_morphism, **kwargs)
            except (NotImplementedError, ValueError, TypeError):
                pass
        (constructor, kwargs) = constructors[-1]
        return constructor(defining_morphism, **kwargs)


RingExtension = RingExtensionFactory("sage.rings.ring_extension.RingExtension")


# General extensions
####################

cdef class RingExtension_generic(Parent):
    r"""
    A generic class for all ring extensions.

    TESTS::

        sage: Q = QQ.over(ZZ)  # indirect doctest
        sage: Q
        Rational Field over its base

        sage: type(Q)
        <class 'sage.rings.ring_extension.RingExtension_generic'>

        sage: TestSuite(Q).run()
    """
    Element = RingExtensionElement

    def __init__(self, defining_morphism, print_options={}, import_methods=True, is_backend_exposed=False, category=None):
        r"""
        Initialize this ring extension.

        INPUT:

        - ``defining_morphism`` -- a ring homomorphism

        - ``print_options`` -- dictionary

        - ``import_methods`` -- boolean (default: ``True``); whether this
          parent (resp. its elements) import the methods of the backend
          parent class (resp. element class)

        - ``is_backend_exposed`` -- boolean (default: ``False``); whether
          the backend ring can be exposed to the user

        - ``category`` -- the category for the resulting parent
          (default: ``CommutativeRings()``)

        .. NOTE::

            The attribute ``is_backend_exposed`` is only used for printing;
            when it is ``False``, printing an element like its backend is
            disabled (and a :exc:`RuntimeError` is raised when it would occur).

        OUTPUT: the extension defined by ``defining_morphism``

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
            TypeError: only commutative rings are accepted

            sage: K = GF(5^3)                                                           # needs sage.rings.finite_rings
            sage: K.over(K.frobenius_endomorphism())                                    # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: exotic defining morphism between two rings in the tower; consider using another variable name
        """
        cdef Parent base, ring
        cdef Parent b, backend
        cdef Map f

        base = defining_morphism.domain()
        ring = defining_morphism.codomain()
        if category is None:
            # Another option would be to set category = CommutativeAlgebras(base)
            # but CommutativeRings() seems safer, especially when dealing with
            # morphisms which do not need to preserve the base
            category = CommutativeRings()
        Parent.__init__(self, base=ZZ, category=category)
        self._base = base
        self._backend = ring
        self._backend_defining_morphism = defining_morphism
        self._defining_morphism = RingExtensionHomomorphism(self._base.Hom(self), defining_morphism)
        self._print_options = print_options.copy()
        if 'over' not in self._print_options:
            self._print_options['over'] = ZZ(0)
        self._import_methods = import_methods
        self._is_backend_exposed = is_backend_exposed
        self._type = "Ring"
        if self._backend in Fields():
            self._type = "Field"

        # Some checks
        if (base not in CommutativeRings()
         or ring not in CommutativeRings()
         or not defining_morphism.category_for().is_subcategory(CommutativeRings())):
            raise TypeError("only commutative rings are accepted")
        f = ring.Hom(ring).identity()
        b = self
        while isinstance(b, RingExtension_generic):
            f *= backend_morphism((<RingExtension_generic>b)._backend_defining_morphism)
            b = b._base
            if isinstance(b, RingExtension_generic):
                backend = (<RingExtension_generic>b)._backend
            else:
                backend = b
            if ring.has_coerce_map_from(backend) and not are_equal_morphisms(f, None):
                # TODO: find a better message
                raise ValueError("exotic defining morphism between two rings in the tower; consider using another variable name")

        # We register coercion/conversion maps
        self.register_coercion(self._defining_morphism.__copy__())
        self.register_coercion(RingExtensionBackendIsomorphism(ring.Hom(self)))
        ring.register_conversion(RingExtensionBackendReverseIsomorphism(self.Hom(ring)))

    def __getattr__(self, name):
        """
        If this extension was created with ``import_methods = True``,
        return a wrapper to the corresponding method of the backend
        parent (if it exists).

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
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
            return self.getattr_from_category(name)
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

            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = QQ.extension(x^2 - 2)                                         # needs sage.rings.number_field
            sage: K.<a> = A.over()                                                      # needs sage.rings.number_field

            sage: dir(K)                                                                # needs sage.rings.number_field
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
            if name[0] == "_": continue
            try:
                attribute = getattr(self._backend, name)
                if callable(attribute):
                    d.append(name)
            except Exception:
                pass
        return sorted(set(d))

    def __hash__(self):
        """
        Return a hash of this extension.

        EXAMPLES:

            sage: E = GF(5^3).over()                                                    # needs sage.rings.finite_rings
            sage: hash(E)   # random                                                    # needs sage.rings.finite_rings
            140257667982632
        """
        return hash_by_id(<void *>self)

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        extension.

        TESTS::

            sage: K = GF(7^3).over()                                                    # needs sage.rings.finite_rings
            sage: type(K)                                                               # needs sage.rings.finite_rings
            <class 'sage.rings.ring_extension.RingExtensionWithGen'>
            sage: loads(dumps(K)) is K                                                  # needs sage.rings.finite_rings
            True
        """
        (defining_morphism, gens, names) = self._factory_data[2]
        constructors = self._factory_data[3]['constructors']
        return RingExtension, (self._backend, defining_morphism, gens, names, constructors)

    def construction(self):
        """
        Return the functorial construction of this extension, if defined.

        EXAMPLES::

             sage: E = GF(5^3).over()                                                   # needs sage.rings.finite_rings
             sage: E.construction()                                                     # needs sage.rings.finite_rings
        """
        # One could define a construction functor K' -> K' otimes_K L, but we leave this to another issue
        pass

    def backend(self, force=False):
        """
        Return the backend of this extension.

        INPUT:

        - ``force`` -- boolean (default: ``False``); if ``False``,
          raise an error if the backend is not exposed

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K = GF(5^3)
            sage: E = K.over()
            sage: E
            Field in z3 with defining polynomial x^3 + 3*x + 3 over its base
            sage: E.backend()
            Finite Field in z3 of size 5^3
            sage: E.backend() is K
            True
        """
        if force or self._is_backend_exposed:
            return self._backend
        raise ValueError("backend is not exposed; try force=True")

    def from_base_ring(self, r):
        r"""
        Return the canonical embedding of ``r`` into this extension.

        INPUT:

        - ``r`` -- an element of the base of the ring of this extension

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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
        if r not in self._base:
            raise TypeError("%s is not an element of the base of %s (= %s)" % (r, self._backend, self._base))
        return self.element_class(self, r)

    def print_options(self, **options):
        """
        Update the printing options of this extension.

        INPUT:

        - ``over`` -- integer or ``Infinity`` (default: ``0``); the maximum
          number of bases included in the printing of this extension

        - ``base`` -- a base over which this extension is finite free;
          elements in this extension will be printed as a linear
          combinaison of a basis of this extension over the given base

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A.<a> = GF(5^2).over()   # over GF(5)
            sage: B.<b> = GF(5^4).over(A)
            sage: C.<c> = GF(5^12).over(B)
            sage: D.<d> = GF(5^24).over(C)

        Observe what happens when we modify the option ``over``::

            sage: # needs sage.rings.finite_rings
            sage: D
            Field in d with defining polynomial
             x^2 + ((1 - a) + ((1 + 2*a) - b)*c + ((2 + a) + (1 - a)*b)*c^2)*x + c over its base
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

            sage: # needs sage.rings.finite_rings
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
        Check and normalize the print option ``over``.

        INPUT:

        - ``over`` -- integer or ``Infinity``

        OUTPUT: the normalized value of ``over``

        TESTS::

            sage: E = QQ.over(ZZ)
            sage: E.print_options(over=-2)
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

        - ``over`` -- integer, ``Infinity`` or ``None``; the maximum
          number of bases included in the print representation of
          this extension. If ``None``, use the print options of this extension.

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
        s = self._repr_topring(**print_options)
        if over is None or over == 0:
            s += " over its base"
        else:
            s += " over "
            base = self._base
            if isinstance(base, RingExtension_generic):
                s += base._repr_(over=over-1)
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
        if not self._is_backend_exposed:
            raise RuntimeError("backend is not exposed to the user; cannot print")
        return str(self._backend)

    def _latex_(self, **options):
        r"""
        Return a LaTeX representation of this extension.

        - ``over`` -- integer, ``Infinity`` or ``None``; the maximum
          number of bases included in the LaTeX representation of
          this extension. If ``None``, use the print options of this extension.

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
        s = self._latex_topring(**print_options)
        if over > 0:
            s += " / "
            base = self._base
            if isinstance(base, RingExtension_generic):
                s += base._latex_(over=over-1)
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
        if not self._is_backend_exposed:
            raise RuntimeError("backend is not exposed to the user; cannot print")
        return latex(self._backend)

    cpdef _coerce_map_from_(self, other):
        r"""
        Return a coerce map from this extension to ``other`` if defined.

        COERCION MODEL

        If `L/K` is an extension, a coercion map `K \to (L/K)`
        (acting through the defining morphism of `L/K`) is set.

        If `L_1/K_1` and `L_2/K_2` are two extensions, a coercion
        map `(L_1/K_1) \to (L_2/K_2)` is set when `L_1` coerces to
        `L_2` and `K_1` coerces to `K_2` in such a way that the
        appropriate diagram commutes.

        These rules have the following consequence regarding iterated
        extensions.
        Given two iterated extensions `A = (A_n/\cdots/A_2/A_1)` and
        `B = (B_m/\cdots/B_2/B_1)`, there is a coercion map `A \to B`
        if there exists a strictly increasing function
        `sigma : \{1,\ldots,n\} \to \{1,\ldots,m\}` and coercion maps
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

            sage: # needs sage.rings.finite_rings
            sage: E1 = GF(3^6).over(GF(3^3))
            sage: E1.coerce_map_from(GF(3^3))  # indirect doctest
            Ring morphism:
              From: Finite Field in z3 of size 3^3
              To:   Field in z6 with defining polynomial x^2 + (2*z3 + 1)*x + z3 over its base
              Defn: z3 |--> z3
            sage: E1.coerce_map_from(GF(3))    # indirect doctest
            Ring morphism:
              From: Finite Field of size 3
              To:   Field in z6 with defining polynomial x^2 + (2*z3 + 1)*x + z3 over its base
              Defn: 1 |--> 1
            sage: E2 = GF(3^18).over(GF(3^9))
            sage: E2.coerce_map_from(E1)       # indirect doctest
            Ring morphism:
              From: Field in z6 with defining polynomial x^2 + (2*z3 + 1)*x + z3 over its base
              To:   Field in z18 with defining polynomial x^2 + (z9^8 + 2*z9^7 + z9^5 + 2*z9^4 + z9^2 + z9 + 1)*x + z9 over its base
              Defn: z6 |--> (2*z9^7 + z9^6 + 2*z9^2 + 2*z9) + (z9^8 + 2*z9^7 + 2*z9^6 + z9^5 + z9^3 + z9 + 1)*z18

        A test with iterated extensions::

            sage: # needs sage.rings.finite_rings
            sage: A = GF(3^18).over(GF(3^3))   #   simple extension GF(3^3) -> GF(3^18)
            sage: B = GF(3^18).over(E1)        # iterated extension GF(3^3) -> GF(3^6) -> GF(3^18)
            sage: A.has_coerce_map_from(B)
            False
            sage: B.has_coerce_map_from(A)
            True
        """
        cdef RingExtension_generic right
        if isinstance(other, RingExtension_generic):
            right = <RingExtension_generic>other
            if self._backend.has_coerce_map_from(right._backend) and self._base.has_coerce_map_from(right._base):
                backend = self._backend.coerce_map_from(right._backend)
                f = backend * backend_morphism(right._defining_morphism)
                g = backend_morphism(self._defining_morphism * self._base.coerce_map_from(right._base))
                if are_equal_morphisms(f, g):
                    return RingExtensionHomomorphism(right.Hom(self), backend)

    def base(self):
        r"""
        Return the base of this extension.

        EXAMPLES::

            sage: F = GF(5^2)                                                           # needs sage.rings.finite_rings
            sage: K = GF(5^4).over(F)                                                   # needs sage.rings.finite_rings
            sage: K.base()                                                              # needs sage.rings.finite_rings
            Finite Field in z2 of size 5^2

        In case of iterated extensions, the base is itself an extension::

            sage: L = GF(5^8).over(K)                                                   # needs sage.rings.finite_rings
            sage: L.base()                                                              # needs sage.rings.finite_rings
            Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
            sage: L.base() is K                                                         # needs sage.rings.finite_rings
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

            sage: # needs sage.rings.finite_rings
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
            [Field in z12 with defining polynomial
              x^3 + (1 + (2 - z2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base,
             Field in z4 with defining polynomial x^2 + (3 - z2)*x + z2 over its base,
             Field in z2 with defining polynomial x^2 + 4*x + 2 over its base,
             Finite Field of size 5]

        .. SEEALSO::

            :meth:`base`, :meth:`absolute_base`, :meth:`is_defined_over`
        """
        L = [self]
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

            sage: # needs sage.rings.finite_rings
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

    cpdef is_defined_over(self, base):
        r"""
        Return whether or not ``base`` is one of the bases of this
        extension.

        INPUT:

        - ``base`` -- a commutative ring, which might be itself an
          extension

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A = GF(5^4).over(GF(5^2))
            sage: B = GF(5^12).over(A)
            sage: A.is_defined_over(GF(5^2))
            True
            sage: A.is_defined_over(GF(5))
            False

            sage: # needs sage.rings.finite_rings
            sage: B.is_defined_over(A)
            True
            sage: B.is_defined_over(GF(5^4))
            True
            sage: B.is_defined_over(GF(5^2))
            True
            sage: B.is_defined_over(GF(5))
            False

        Note that an extension is defined over itself::

            sage: A.is_defined_over(A)                                                  # needs sage.rings.finite_rings
            True
            sage: A.is_defined_over(GF(5^4))                                            # needs sage.rings.finite_rings
            True

        .. SEEALSO::

            :meth:`base`, :meth:`bases`, :meth:`absolute_base`
        """
        cdef Parent b
        b = self
        while isinstance(b, RingExtension_generic):
            if b is base or (<RingExtension_generic>b)._backend is base: return True
            b = (<RingExtension_generic>b)._base
        return b is base

    cpdef Parent _check_base(self, Parent base):
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

            sage: # needs sage.rings.finite_rings
            sage: F = GF(5^2)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(K)
            sage: L.bases()
            [Field in z12 with defining polynomial x^3 + (1 + (4*z2 + 2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base,
             Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base,
             Finite Field in z2 of size 5^2]
            sage: L._check_base(K)
            Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
            sage: L._check_base(GF(5^4))
            Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
            sage: L._check_base(GF(5^4)) is K
            True

        When ``base`` is ``None``, the base of the extension is returned::

            sage: L._check_base(None)                                                   # needs sage.rings.finite_rings
            Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
            sage: L._check_base(None) is L.base()                                       # needs sage.rings.finite_rings
            True
        """
        cdef Parent b
        if base is None:
            return self._base
        b = self
        while isinstance(b, RingExtension_generic):
            if b is base or (<RingExtension_generic>b)._backend is base: return b
            b = (<RingExtension_generic>b)._base
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

            sage: # needs sage.rings.finite_rings
            sage: F = GF(5^2)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(K)
            sage: K.defining_morphism()
            Ring morphism:
              From: Finite Field in z2 of size 5^2
              To:   Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
              Defn: z2 |--> z2
            sage: L.defining_morphism()
            Ring morphism:
              From: Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
              To:   Field in z12 with defining polynomial
                    x^3 + (1 + (4*z2 + 2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base
              Defn: z4 |--> z4

        One can also pass in a base over which the extension is explicitly
        defined (see also :meth:`is_defined_over`)::

            sage: L.defining_morphism(F)                                                # needs sage.rings.finite_rings
            Ring morphism:
              From: Finite Field in z2 of size 5^2
              To:   Field in z12 with defining polynomial
                    x^3 + (1 + (4*z2 + 2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base
              Defn: z2 |--> z2
            sage: L.defining_morphism(GF(5))                                            # needs sage.rings.finite_rings
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

            sage: # needs sage.rings.finite_rings
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
        return tuple([self(x) for x in generators(self._backend,
                                                  backend_parent(self._base))])

    def ngens(self, base=None):
        r"""
        Return the number of generators of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

    def gen(self):
        r"""
        Return the first generator of this extension.

        EXAMPLES::

            sage: K = GF(5^2).over()   # over GF(5)                                     # needs sage.rings.finite_rings
            sage: x = K.gen(); x                                                        # needs sage.rings.finite_rings
            z2

        Observe that the generator lives in the extension::

            sage: x.parent()                                                            # needs sage.rings.finite_rings
            Field in z2 with defining polynomial x^2 + 4*x + 2 over its base
            sage: x.parent() is K                                                       # needs sage.rings.finite_rings
            True
        """
        return self.gens()[0]

    def random_element(self):
        r"""
        Return a random element in this extension.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

    def degree_over(self, base=None):
        r"""
        Return the degree of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

            sage: K.degree_over()                                                       # needs sage.rings.finite_rings
            2
            sage: L.degree_over()                                                       # needs sage.rings.finite_rings
            3

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: K.degree_over(GF(5))                                                  # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field of size 5
        """
        base = self._check_base(base)
        return self._degree_over(base)

    cpdef _degree_over(self, Parent base):
        r"""
        Return the degree of this extension over ``base``.

        Should be implemented in subclasses.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        TESTS::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
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

            sage: # needs sage.rings.finite_rings
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

            sage: A.degree(GF(5))                                                       # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field of size 5

        .. SEEALSO::

            :meth:`relative_degree`, :meth:`absolute_degree`
        """
        return self.degree_over(base)

    def relative_degree(self):
        r"""
        Return the degree of this extension over its base.

        EXAMPLES::

            sage: A = GF(5^4).over(GF(5^2))                                             # needs sage.rings.finite_rings
            sage: A.relative_degree()                                                   # needs sage.rings.finite_rings
            2

        .. SEEALSO::

            :meth:`degree`, :meth:`absolute_degree`
        """
        return self._degree_over(self._base)

    def absolute_degree(self):
        r"""
        Return the degree of this extension over its absolute base.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

            sage: # needs sage.rings.finite_rings
            sage: K = GF(5^2).over()  # over GF(5)
            sage: L = GF(5^4).over(K)
            sage: L.is_finite_over(K)
            True
            sage: L.is_finite_over(GF(5))
            True

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: L.is_finite_over()                                                    # needs sage.rings.finite_rings
            True
        """
        cdef Parent b
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
            b = (<RingExtension_generic?>b)._base
        raise NotImplementedError

    cpdef _is_finite_over(self, Parent base):
        r"""
        Return whether or not this extension is finite over ``base``.

        Should be implemented in subclasses.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: K = GF(5^2).over()  # over GF(5)                                      # needs sage.rings.finite_rings
            sage: K.is_finite_over()  # indirect doctest                                # needs sage.rings.finite_rings
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

            sage: # needs sage.rings.finite_rings
            sage: K = GF(5^2).over()  # over GF(5)
            sage: L = GF(5^4).over(K)
            sage: L.is_free_over(K)
            True
            sage: L.is_free_over(GF(5))
            True

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: L.is_free_over()                                                      # needs sage.rings.finite_rings
            True
        """
        cdef Parent b
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
            b = (<RingExtension_generic?>b)._base
        raise NotImplementedError

    cpdef _is_free_over(self, Parent base):
        r"""
        Return whether or not this extension is finite over ``base``.

        Should be implemented in subclasses.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: K = GF(5^2).over()  # over GF(5)                                      # needs sage.rings.finite_rings
            sage: K.is_free_over()  # indirect doctest                                  # needs sage.rings.finite_rings
            True
        """
        raise NotImplementedError

    def is_field(self, proof=True):
        r"""
        Return whether or not this extension is a field.

        INPUT:

        - ``proof`` -- boolean (default: ``False``)

        EXAMPLES::

            sage: K = GF(5^5).over()  # over GF(5)                                      # needs sage.rings.finite_rings
            sage: K.is_field()                                                          # needs sage.rings.finite_rings
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

        - ``extend_base`` -- boolean (default: ``False``)

        If ``extend_base`` is ``False``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/L/K`, except
        if `L` is already a field in which base the fraction field
        of `L/K` is `L/K` itself.

        If ``extend_base`` is ``True``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/\textrm{Frac}(K)`
        (provided that the defining morphism extends to the fraction
        fields, i.e. is injective).

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = ZZ.extension(x^2 - 5)
            sage: OK = A.over()   # over ZZ
            sage: OK
            Order of conductor 2 generated by a in Number Field in a
             with defining polynomial x^2 - 5 over its base
            sage: K1 = OK.fraction_field(); K1
            Fraction Field of Order of conductor 2 generated by a
             in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K1.bases()
            [Fraction Field of Order of conductor 2 generated by a
             in Number Field in a with defining polynomial x^2 - 5 over its base,
             Order of conductor 2 generated by a in Number Field in a
             with defining polynomial x^2 - 5 over its base,
             Integer Ring]
            sage: K2 = OK.fraction_field(extend_base=True); K2
            Fraction Field of Order of conductor 2 generated by a
             in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K2.bases()
            [Fraction Field of Order of conductor 2 generated by a
             in Number Field in a with defining polynomial x^2 - 5 over its base,
             Rational Field]

        Note that there is no coercion between `K_1` and `K_2`::

            sage: K1.has_coerce_map_from(K2)                                            # needs sage.rings.number_field
            False
            sage: K2.has_coerce_map_from(K1)                                            # needs sage.rings.number_field
            False

        We check that when the extension is a field, its fraction field does not change::

            sage: K1.fraction_field() is K1                                             # needs sage.rings.number_field
            True
            sage: K2.fraction_field() is K2                                             # needs sage.rings.number_field
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
        constructor = RingExtensionFractionField, {'ring': self, 'is_backend_exposed': self._is_backend_exposed}
        return RingExtension(ring, defining_morphism, constructors=[constructor])

    cdef Map _defining_morphism_fraction_field(self, bint extend_base):
        r"""
        Return the defining morphism of the fraction field of this extension.

        This is a helper function.

        INPUT:

        - ``extend_base`` -- boolean (default: ``False``); see
          :meth:`fraction_field` for more informations

        TESTS::

            sage: K = GF(5^2).over()                                                    # needs sage.rings.finite_rings
            sage: K.fraction_field()  # indirect doctest                                # needs sage.rings.finite_rings
            Field in z2 with defining polynomial x^2 + 4*x + 2 over its base

            sage: K = QQ.over(ZZ)
            sage: K.fraction_field(extend_base=True)  # indirect doctest
            Rational Field over its base
        """
        if extend_base:
            defining_morphism = backend_morphism(self._backend_defining_morphism)
            defining_morphism = defining_morphism.extend_to_fraction_field()
            if isinstance(self._base, RingExtension_generic):
                base = self._base.fraction_field(extend_base)
                ring = defining_morphism.codomain()
                defining_morphism = RingExtensionHomomorphism(base.Hom(ring), defining_morphism)
        else:
            if self.is_field():
                defining_morphism = None
            else:
                ring = self._backend.fraction_field()
                defining_morphism = RingExtensionHomomorphism(self.Hom(ring), ring.coerce_map_from(self._backend))
        return defining_morphism

    def _Hom_(self, codomain, category):
        r"""
        Return the homset from this extension of ``codomain`` is the category ``category``.

        INPUT:

        - ``codomain`` -- a parent

        - ``category`` -- a subcategory of the category of rings

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(5^2)
            sage: K = GF(5^4).over(F)
            sage: L = GF(5^12).over(F)
            sage: K.Hom(L)  # indirect doctest
            Set of Homomorphisms
             from Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
               to Field in z12 with defining polynomial x^6 + (4*z2 + 3)*x^5 + x^4
                  + (3*z2 + 1)*x^3 + x^2 + (4*z2 + 1)*x + z2 over its base
            sage: K.Hom(L, category=Sets())
            Set of Morphisms
             from Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
               to Field in z12 with defining polynomial x^6 + (4*z2 + 3)*x^5 + x^4
                  + (3*z2 + 1)*x^3 + x^2 + (4*z2 + 1)*x + z2 over its base
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

        - ``check`` -- boolean (default: ``True``); whether to verify that the
          images of generators extend to define a map (using only canonical coercions)

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()    # over GF(5)                                # needs sage.rings.finite_rings
            sage: L.<b> = GF(5^6).over(K)                                               # needs sage.rings.finite_rings

        We define (by hand) the relative Frobenius endomorphism of the extension `L/K`::

            sage: L.hom([b^25])                                                         # needs sage.rings.finite_rings
            Ring endomorphism of
             Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: b |--> 2 + 2*a*b + (2 - a)*b^2

        Defining the absolute Frobenius of `L` is a bit more complicated
        because it is not a homomorphism of `K`-algebras.
        For this reason, the construction ``L.hom([b^5])`` fails::

            sage: L.hom([b^5])                                                          # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: images do not define a valid homomorphism

        What we need is to specify a base map::

            sage: FrobK = K.hom([a^5])                                                  # needs sage.rings.finite_rings
            sage: FrobL = L.hom([b^5], base_map=FrobK); FrobL                           # needs sage.rings.finite_rings
            Ring endomorphism of
             Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: b |--> (-1 + a) + (1 + 2*a)*b + a*b^2
                    with map on base ring:
                    a |--> 1 - a

        As a shortcut, we may use the following construction::

            sage: phi = L.hom([b^5, a^5]); phi                                          # needs sage.rings.finite_rings
            Ring endomorphism of
             Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: b |--> (-1 + a) + (1 + 2*a)*b + a*b^2
                    with map on base ring:
                    a |--> 1 - a
            sage: phi == FrobL                                                          # needs sage.rings.finite_rings
            True
        """
        if codomain is None:
            from sage.structure.sequence import Sequence
            codomain = Sequence(im_gens).universe()
        parent = self.Hom(codomain, category=category)
        return RingExtensionHomomorphism(parent, im_gens, base_map, check)

    def characteristic(self):
        r"""
        Return the characteristic of the extension as a ring.

        OUTPUT: a prime number or zero

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

        Ensure issue :issue:`34692` is fixed::

            sage: Fq = GF(11)
            sage: FqX.<X> = Fq[]
            sage: k = Frac(FqX)
            sage: K = k.over(FqX)
            sage: K.frobenius_endomorphism()
            Frobenius endomorphism x |--> x^11 of
             Fraction Field of Univariate Polynomial Ring in X over
              Finite Field of size 11 over its base
        """
        return self._backend.characteristic()


# Fraction fields
#################

cdef class RingExtensionFractionField(RingExtension_generic):
    """
    A class for ring extensions of the form `\textrm{Frac}(A)/A`.

    TESTS::

        sage: Z = ZZ.over()   # over ZZ itself
        sage: Q = Z.fraction_field()
        sage: Q
        Fraction Field of Integer Ring over its base

        sage: type(Q)
        <class 'sage.rings.ring_extension.RingExtensionFractionField'>

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

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: OK = A.over()
            sage: K = OK.fraction_field(); K
            Fraction Field of
             Maximal Order generated by a in Number Field in a
             with defining polynomial x^2 - 2 over its base
            sage: TestSuite(K).run()
        """
        RingExtension_generic.__init__(self, defining_morphism, **kwargs)
        if ring is None:
            self._ring = self._base
        else:
            self._ring = ring

    def ring(self):
        r"""
        Return the ring whose fraction field is this extension.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: OK = A.over()
            sage: K = OK.fraction_field(); K
            Fraction Field of
             Maximal Order generated by a in Number Field in a with defining polynomial x^2 - 2 over its base
            sage: K.ring()
            Maximal Order generated by a in Number Field in a with defining polynomial x^2 - 2 over its base
            sage: K.ring() is OK
            True
        """
        return self._ring

    def _repr_topring(self, **options):
        r"""
        Return a string representation of top ring of this extension.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: OK = A.over()
            sage: K = OK.fraction_field()
            sage: K._repr_topring()
            'Fraction Field of Maximal Order generated by a in Number Field in a with defining polynomial x^2 - 2'
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

cdef class RingExtensionWithBasis(RingExtension_generic):
    """
    A class for finite free ring extensions equipped
    with a basis.

    TESTS::

        sage: E = GF(5^4).over(GF(5^2)); E                                              # needs sage.rings.finite_rings
        Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base

        sage: TestSuite(E).run()                                                        # needs sage.rings.finite_rings
    """
    Element = RingExtensionWithBasisElement

    def __init__(self, defining_morphism, basis, names=None, check=True, **kwargs):
        r"""
        Initialize this ring extension.

        INPUT:

        - ``defining_morphism`` -- a ring homomorphism

        - ``basis`` -- tuple of elements in this extension

        - ``names`` -- tuple of strings or ``None`` (default: ``None``);
          the way the elements of the basis are printed

        - ``check`` -- boolean (default: ``True``); whether to check if
          ``basis`` is indeed a basis

        TESTS::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = QQ.extension(x^3 - 2)                                         # needs sage.rings.number_field
            sage: E = K.over(); E                                                       # needs sage.rings.number_field
            Field in a with defining polynomial x^3 - 2 over its base

            sage: TestSuite(E).run()                                                    # needs sage.rings.number_field
        """
        RingExtension_generic.__init__(self, defining_morphism, **kwargs)
        self._basis = [self(b) for b in basis]
        if names is None:
            names = []
            for b in self._basis:
                b = b._backend
                if b == 1:
                    names.append("")
                sb = str(b)
                if b._is_atomic() or (sb[0] == "(" and sb[-1] == ")"):
                    names.append(sb)
                else:
                    names.append("(" + sb + ")")
        else:
            if len(names) != len(self._basis):
                raise ValueError("the number of names does not match the cardinality of the basis")
        self._basis_names = names
        self._basis_latex_names = [latex_variable_name(name) for name in names]
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

            sage: # needs sage.rings.finite_rings
            sage: F = GF(5)
            sage: K = GF(5^2).over(F)
            sage: L = GF(5^4).over(K)
            sage: L._print_option_base(F) is F
            True
            sage: L._print_option_base(K) is K
            True
            sage: L._print_option_base(GF(5^2)) is K
            True
            sage: L._print_option_base(None) is K
            True
            sage: L._print_option_base(L)
            Traceback (most recent call last):
            ...
            ValueError: base must be strict
            sage: K._print_option_base(L)
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Field in z4
            with defining polynomial x^2 + (3 - z2)*x + z2 over its base
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

    cpdef _degree_over(self, Parent base):
        r"""
        Return the degree of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        TESTS::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
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

    cpdef _is_finite_over(self, Parent base):
        r"""
        Return whether or not this extension is finite over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: K = GF(5^2).over()  # over GF(5)                                      # needs sage.rings.finite_rings
            sage: K.is_finite_over()  # indirect doctest                                # needs sage.rings.finite_rings
            True
        """
        if base is self or base is self._base:
            return True
        return self._base._is_finite_over(base)

    cpdef _is_free_over(self, Parent base):
        r"""
        Return whether or not this extension is free over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: K = GF(5^2).over()  # over GF(5)                                      # needs sage.rings.finite_rings
            sage: K.is_free_over()    # indirect doctest                                # needs sage.rings.finite_rings
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

            sage: # needs sage.rings.finite_rings
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

            sage: L.basis_over()                                                        # needs sage.rings.finite_rings
            [1, c, c^2]

            sage: K.basis_over()                                                        # needs sage.rings.finite_rings
            [1, b]

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: L.degree_over(GF(5^6))                                                # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z6 of size 5^6
        """
        base = self._check_base(base)
        return self._basis_over(base)

    cpdef _basis_over(self, Parent base):
        r"""
        Return a basis of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = QQ.extension(x^3 - 2)                                         # needs sage.rings.number_field
            sage: K.<u> = A.over()                                                      # needs sage.rings.number_field
            sage: K.basis_over()                                                        # needs sage.rings.number_field
            [1, u, u^2]
        """
        if base is self:
            return [self.one()]
        if base is self._base:
            return self._basis[:]
        b = self._base._basis_over(base)
        return [x * y for x in self._basis for y in b]

    def free_module(self, base=None, map=True):
        r"""
        Return a free module V over ``base`` which is isomorphic to
        this ring

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        - ``map`` -- boolean (default: ``True``); whether to return
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
            sage: K.<a> = GF(11^2).over()                                               # needs sage.rings.finite_rings
            sage: L.<b> = GF(11^6).over(K)                                              # needs sage.rings.finite_rings

        Forgetting a part of the multiplicative structure, the field L
        can be viewed as a vector space of dimension 3 over K, equipped
        with a distinguished basis, namely `(1, b, b^2)`::

            sage: # needs sage.rings.finite_rings
            sage: V, i, j = L.free_module(K)
            sage: V
            Vector space of dimension 3 over
             Field in a with defining polynomial x^2 + 7*x + 2 over its base
            sage: i
            Generic map:
              From: Vector space of dimension 3 over
                    Field in a with defining polynomial x^2 + 7*x + 2 over its base
              To:   Field in b with defining polynomial
                    x^3 + (7 + 2*a)*x^2 + (2 - a)*x - a over its base
            sage: j
            Generic map:
              From: Field in b with defining polynomial
                    x^3 + (7 + 2*a)*x^2 + (2 - a)*x - a over its base
              To:   Vector space of dimension 3 over
                    Field in a with defining polynomial x^2 + 7*x + 2 over its base
            sage: j(b)
            (0, 1, 0)
            sage: i((1, a, a+1))
            1 + a*b + (1 + a)*b^2

        Similarly, one can view L as a F-vector space of dimension 6::

            sage: V, i, j, = L.free_module(F)                                           # needs sage.rings.finite_rings
            sage: V                                                                     # needs sage.rings.finite_rings
            Vector space of dimension 6 over Finite Field of size 11

        In this case, the isomorphisms between `V` and `L` are given by the
        basis `(1, a, b, ab, b^2, ab^2)`:

            sage: j(a*b)                                                                # needs sage.rings.finite_rings
            (0, 0, 0, 1, 0, 0)
            sage: i((1,2,3,4,5,6))                                                      # needs sage.rings.finite_rings
            (1 + 2*a) + (3 + 4*a)*b + (5 + 6*a)*b^2

        When ``base`` is omitted, the default is the base of this extension::

            sage: L.free_module(map=False)                                              # needs sage.rings.finite_rings
            Vector space of dimension 3 over
             Field in a with defining polynomial x^2 + 7*x + 2 over its base

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: L.degree(GF(11^3))                                                    # needs sage.rings.finite_rings
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

        - ``map`` -- boolean (default: ``True``); whether to return
          isomorphisms between this ring and V

        OUTPUT:

        - A finite-rank free module V over ``base``

        - The isomorphism from V to this ring corresponding to the
          basis output by the method :meth:`basis_over`
          (only included if ``map`` is ``True``)

        - The reverse isomorphism of the isomorphism above
          (only included if ``map`` is ``True``)

        TESTS::

            sage: K = GF(7^5).over()                                                    # needs sage.rings.finite_rings
            sage: L = GF(7^15).over(K)                                                  # needs sage.rings.finite_rings
            sage: for base in L.bases():                                                # needs sage.rings.finite_rings
            ....:     V, i, j = L.free_module(base)
            ....:     assert([i(v) for v in V.basis()] == L.basis_over(base))
            ....:     assert([j(x) for x in L.basis_over(base)] == V.basis())
        """
        d = self._degree_over(base)
        if map:
            return base**d, MapFreeModuleToRelativeRing(self, base), MapRelativeRingToFreeModule(self, base)
        else:
            return base**d

    @cached_method
    def fraction_field(self, extend_base=False):
        r"""
        Return the fraction field of this extension.

        INPUT:

        - ``extend_base`` -- boolean (default: ``False``)

        If ``extend_base`` is ``False``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/L/K`, except
        is `L` is already a field in which base the fraction field
        of `L/K` is `L/K` itself.

        If ``extend_base`` is ``True``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/\textrm{Frac}(K)`
        (provided that the defining morphism extends to the fraction
        fields, i.e. is injective).

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = ZZ.extension(x^2 - 5)
            sage: OK = A.over()   # over ZZ
            sage: OK
            Order of conductor 2 generated by a in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K1 = OK.fraction_field(); K1
            Fraction Field of Order of conductor 2 generated by a in Number Field in a
             with defining polynomial x^2 - 5 over its base
            sage: K1.bases()
            [Fraction Field of Order of conductor 2 generated by a in Number Field in a
              with defining polynomial x^2 - 5 over its base,
             Order of conductor 2 generated by a in Number Field in a
              with defining polynomial x^2 - 5 over its base,
             Integer Ring]
            sage: K2 = OK.fraction_field(extend_base=True); K2
            Fraction Field of Order of conductor 2 generated by a
              in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K2.bases()
            [Fraction Field of Order of conductor 2 generated by a
              in Number Field in a with defining polynomial x^2 - 5 over its base,
             Rational Field]

        Note that there is no coercion map between `K_1` and `K_2`::

            sage: K1.has_coerce_map_from(K2)                                            # needs sage.rings.number_field
            False
            sage: K2.has_coerce_map_from(K1)                                            # needs sage.rings.number_field
            False

        We check that when the extension is a field, its fraction field does not change::

            sage: K1.fraction_field() is K1                                             # needs sage.rings.number_field
            True
            sage: K2.fraction_field() is K2                                             # needs sage.rings.number_field
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
        kwargs['is_backend_exposed'] = self._is_backend_exposed
        ring = defining_morphism.codomain()
        return RingExtension(ring, defining_morphism, gen=gen, names=names, constructors=[(constructor, kwargs)])


cdef class RingExtensionWithGen(RingExtensionWithBasis):
    """
    A class for finite free ring extensions generated by
    a single element

    TESTS::

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: A.<a> = QQ.extension(x^3 - 7)
        sage: K = A.over()
        sage: type(K)
        <class 'sage.rings.ring_extension.RingExtensionWithGen'>
        sage: TestSuite(K).run()
    """
    def __init__(self, defining_morphism, gen, names, check=True, **kwargs):
        r"""
        Initialize this ring extension.

        INPUT:

        - ``defining_morphism`` -- a ring homomorphism

        - ``gen`` -- a generator of this extension

        - ``names`` -- tuple of strings or ``None`` (default: ``None``);
          the way the elements of the basis are printed

        - ``check`` -- boolean (default: ``True``); whether to check if
          ``gen`` is indeed a generator

        TESTS::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = QQ.extension(x^3 + 3*x + 1)                                   # needs sage.rings.number_field
            sage: E = K.over(); E                                                       # needs sage.rings.number_field
            Field in a with defining polynomial x^3 + 3*x + 1 over its base

            sage: TestSuite(E).run()                                                    # needs sage.rings.number_field
        """
        self._name = names[0]
        backend_base = backend_parent(defining_morphism.domain())
        _, deg_domain, deg_codomain = common_base(backend_base, defining_morphism.codomain(), True)
        degree = deg_codomain // deg_domain
        basis_names = [""]
        basis_latex_names = [""]
        if degree == 1:
            self._name = None
        else:
            basis_names += [self._name] + ["%s^%s" % (self._name, i) for i in range(2, degree)]
            latex_name = latex_variable_name(self._name)
            basis_latex_names += [latex_name] + ["%s^{%s}" % (latex_name, i) for i in range(2, degree)]
        basis = [gen ** i for i in range(degree)]
        RingExtensionWithBasis.__init__(self, defining_morphism, basis, basis_names, check, **kwargs)
        self._gen = self._backend(gen)
        self._names = (self._name,)
        self._basis_latex_names = basis_latex_names

    def _repr_topring(self, **options):
        r"""
        Return a string representation of top ring of this extension.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(5^3).over()
            sage: K._repr_topring()
            'Field in a with defining polynomial x^3 + 3*x + 3'
            sage: L.<b> = GF(5^9).over(K)
            sage: L._repr_topring()
            'Field in b with defining polynomial x^3 + (1 + 3*a^2)*x^2 + (3 + 2*a + 2*a^2)*x - a'
        """
        if self._name is None:
            return RingExtension_generic._repr_topring(self, **options)
        return "%s in %s with defining polynomial %s" % (self._type, self._name, self.modulus())

    def _latex_topring(self):
        r"""
        Return a LaTeX representation of top ring of this extension.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(5^3).over()
            sage: K._latex_topring()
            '\\Bold{F}_{5}[a]'
            sage: L.<b> = GF(5^9).over(K)
            sage: L._latex_topring()
            '\\Bold{F}_{5}[a][b]'
        """
        if self._name is None:
            return RingExtension_generic._latex_topring(self)
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

            sage: # needs sage.rings.finite_rings
            sage: K.<u> = GF(7^10).over(GF(7^2)); K
            Field in u with defining polynomial x^5 + (6*z2 + 4)*x^4
             + (3*z2 + 5)*x^3 + (2*z2 + 2)*x^2 + 4*x + 6*z2 over its base
            sage: P = K.modulus(); P
            x^5 + (6*z2 + 4)*x^4 + (3*z2 + 5)*x^3 + (2*z2 + 2)*x^2 + 4*x + 6*z2
            sage: P(u)
            0

        We can use a different variable name::

            sage: K.modulus('y')                                                        # needs sage.rings.finite_rings
            y^5 + (6*z2 + 4)*y^4 + (3*z2 + 5)*y^3 + (2*z2 + 2)*y^2 + 4*y + 6*z2
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        _, _, j = self.free_module(map=True)
        d = self.relative_degree()
        coeffs = [-c for c in j(self._gen**d)] + [1]
        S = PolynomialRing(self._base, name=var)
        return S(coeffs)

    def gens(self, base=None):
        r"""
        Return the generators of this extension over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (default: ``None``)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: K.gens()
            (a,)
            sage: L.<b> = GF(5^4).over(K)
            sage: L.gens()
            (b,)
            sage: L.gens(GF(5))
            (b, a)
        """
        if base is None:
            return (self(self._gen),)
        base = self._check_base(base)
        gens = tuple()
        b = self
        while b is not base:
            gens += b.gens()
            b = b.base()
        return gens

    @cached_method
    def fraction_field(self, extend_base=False):
        r"""
        Return the fraction field of this extension.

        INPUT:

        - ``extend_base`` -- boolean (default: ``False``)

        If ``extend_base`` is ``False``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/L/K`, except
        is `L` is already a field in which base the fraction field
        of `L/K` is `L/K` itself.

        If ``extend_base`` is ``True``, the fraction field of the
        extension `L/K` is defined as `\textrm{Frac}(L)/\textrm{Frac}(K)`
        (provided that the defining morphism extends to the fraction
        fields, i.e. is injective).

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = ZZ.extension(x^2 - 5)
            sage: OK = A.over()   # over ZZ
            sage: OK
            Order of conductor 2 generated by a in Number Field in a
             with defining polynomial x^2 - 5 over its base
            sage: K1 = OK.fraction_field(); K1
            Fraction Field of Order of conductor 2 generated by a
             in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K1.bases()
            [Fraction Field of Order of conductor 2 generated by a
              in Number Field in a with defining polynomial x^2 - 5 over its base,
             Order of conductor 2 generated by a in Number Field in a
              with defining polynomial x^2 - 5 over its base,
             Integer Ring]
            sage: K2 = OK.fraction_field(extend_base=True); K2
            Fraction Field of Order of conductor 2 generated by a
             in Number Field in a with defining polynomial x^2 - 5 over its base
            sage: K2.bases()
            [Fraction Field of Order of conductor 2 generated by a
              in Number Field in a with defining polynomial x^2 - 5 over its base,
             Rational Field]

        Note that there is no coercion map between `K_1` and `K_2`::

            sage: K1.has_coerce_map_from(K2)                                            # needs sage.rings.number_field
            False
            sage: K2.has_coerce_map_from(K1)                                            # needs sage.rings.number_field
            False

        We check that when the extension is a field, its fraction field does not change::

            sage: K1.fraction_field() is K1                                             # needs sage.rings.number_field
            True
            sage: K2.fraction_field() is K2                                             # needs sage.rings.number_field
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
        kwargs['is_backend_exposed'] = self._is_backend_exposed
        ring = defining_morphism.codomain()
        return RingExtension(ring, defining_morphism, gen=gen, names=names, constructors=[(constructor, kwargs)])
