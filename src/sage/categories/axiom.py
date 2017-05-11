# -*- coding: utf-8 -*-
"""
Axioms

At this stage, axioms in Sage are very lightweight objects, mostly
information placeholders.  On the other hand, they are used
extensively to relate categories to each other.  For details, see
:ref:`sage.categories.category_with_axiom`.  For a catalog of axioms,
see :obj:`axioms`.
"""
#*****************************************************************************
#       Copyright (C) 2017 Nicolas M. Thiéry <nthiery at users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import re
from sage.misc.fast_methods import Singleton
from sage.misc.superseded import deprecation
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.sage_object import SageObject

# utility
def uncamelcase(s, separator=" "):
    """
    Return a copy of the string with words in lower case and separated.

    INPUT:

    - ``s`` -- a string in camel case
    - ``separator`` - a string (default: a space)

    EXAMPLES::

        sage: sage.categories.category_with_axiom.uncamelcase("FiniteDimensionalAlgebras")
        'finite dimensional algebras'
        sage: sage.categories.category_with_axiom.uncamelcase("JTrivialMonoids")
        'j trivial monoids'
        sage: sage.categories.category_with_axiom.uncamelcase("FiniteDimensionalAlgebras", "_")
        'finite_dimensional_algebras'
    """
    return re.sub("(?!^)[A-Z]", lambda match: separator+match.group()[0], s).lower()

class AxiomCatalog(Singleton):
    r"""
    A catalog of the currently defined axioms.

    EXAMPLES::

        sage: list(axioms)
        [..., Finite, ...]
        sage: axioms.Distributive
        Distributive

    .. NOTE::

        As the Sage library gets lazily loaded in memory, newer axioms
        can appear in this catalog.

    TESTS::

        sage: type(axioms)
        <class 'sage.categories.axiom.AxiomCatalog'>
    """
    def __contains__(self, obj):
        r"""
        Return whether obj is a currently defined axiom or name thereof.

        EXAMPLES::

            sage: axioms.Finite in axioms
            True
            sage: "Distributive" in axioms
            True
            sage: "Foo" in axioms
            False
            sage: 3 in axioms
            False
        """
        if isinstance(obj, str):
            return obj in self.__dict__
        elif isinstance(obj, Axiom):
            return str(obj) in self.__dict__
        else:
            return False

    def __call__(self, obj):
        r"""
        Return an axiom from itself or its name.

        INPUT:

        - ``obj`` -- an existing axiom or name thereof

        This is typically used to massage user inputs, converting
        axiom names into axioms as needed. See e.g.
        :meth:`~sage.categories.category.Category._with_axiom`.

        EXAMPLES::

            sage: axioms("Associative")
            Associative

            sage: axioms(axioms.Associative)
            Associative

            sage: axioms("Foobar")
            Traceback (most recent call last):
            ...
            ValueError: No known axiom named Foobar

        TESTS

            sage: axioms("Associative") is axioms.Associative
            True
            sage: axioms(axioms.Associative) is axioms.Associative
            True
        """
        if isinstance(obj, Axiom):
            return obj
        if isinstance(obj, str):
            try:
                return getattr(self, obj)
            except AttributeError:
                raise ValueError("No known axiom named %s"%obj)
        raise ValueError("Input should be an axiom or axiom name; got %s"%obj)

    def __iter__(self):
        r"""
        Iterate through all the currently defined axioms.

        EXAMPLES::

            sage: sorted(axioms, key=str)
            [..., Commutative, ..., Distributive, ..., Finite, ..., Smooth, ...]

        .. NOTE::

            The order in which axioms are iterated through is currently unspecified.
        """
        return iter(self.__dict__.values())

    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: axioms
            A catalog of all currently defined axioms
        """
        return "A catalog of all currently defined axioms"

    def add(self, axiom_name):
        r"""
        Add a new axiom.

        EXAMPLES::

            sage: from sage.categories.category_with_axiom import Axiom
            sage: m = max(axiom.index() for axiom in axioms)
            sage: axioms.add('Awesome')
            doctest:...
            DeprecationWarning: Creating axioms with all_axioms.add("Foo") or all_axioms+=("Foo") is deprecated.
            Please use Axiom("Foo") instead.
            See http://trac.sagemath.org/22965 for details.

            sage: axioms.Awesome.index() == m + 1
            True
        """
        deprecation(22965, 'Creating axioms with all_axioms.add("Foo") or all_axioms+=("Foo") is deprecated.\nPlease use Axiom("Foo") instead.')
        if not isinstance(axiom_name, str):
            raise ValueError("axiom_name should be a string; got %s"%axiom_name)
        Axiom(axiom_name)

    def __iadd__(self, L):
        r"""
        Inline addition, which means to add a list of axioms to the container.

        EXAMPLES::

            sage: from sage.categories.category_with_axiom import Axiom
            sage: m = max(axiom.index() for axiom in axioms)
            sage: axioms += ('Fancy', 'Amazing')
            doctest:...
            DeprecationWarning: Creating axioms with all_axioms.add("Foo") or all_axioms+=("Foo") is deprecated.
            Please use Axiom("Foo") instead.
            See http://trac.sagemath.org/22965 for details.

        TESTS::

            sage: axioms.Fancy.index()   == m + 1
            True
            sage: axioms.Amazing.index() == m + 2
            True
        """
        for axiom in L:
            self.add(axiom)
        return self

all_axioms = AxiomCatalog()

class Axiom(UniqueRepresentation, SageObject):
    def __eq__(self, other):
        r"""
        Return whether ``self == other``.

        This is a slight deviation upon the equality method by ``id``
        derived from :class:`UniqueRepresentation`, in order to
        support equality with the eponym strings.  This is mostly for
        backward compatibility.

        TESTS::

            sage: axioms.Associative == "Associative"
            True
            sage: "Associative" == axioms.Associative
            True
            sage: axioms.Associative == axioms.Associative
            True
            sage: axioms.Associative == axioms.Commutative
            False
        """
        if isinstance(other, str):
            return self._name == other
        else:
            return self is other

    def __ne__(self, other):
        r"""
        Return whether ``self != other``.

        TESTS::

            sage: axioms.Associative != "Associative"
            False
            sage: "Associative" != axioms.Associative
            False
            sage: axioms.Associative != axioms.Associative
            False
            sage: axioms.Associative != axioms.Commutative
            True
        """
        return not self == other

    def __hash__(self):
        r"""
        Return a hash value for ``self``.

        This is a slight deviation upon the default hash value derived
        from :class:`UniqueRepresentation`, in order to have a
        consistent hash with the eponym strings.  This is mostly for
        backward compatibility.

        TESTS::

            sage: hash(axioms.Associative) == hash("Associative")
            True

            sage: axioms.Associative in {"Associative", "Commutative"}
            True
            sage: "Associative" in {axioms.Associative, axioms.Commutative}
            True
        """
        return hash(self._name)

    def __init__(self, name):
        """
        Initializes this axiom.

        TESTS::

            sage: from sage.categories.axiom import Axiom
            sage: Coucou = Axiom("Coucou"); Coucou
            Coucou
            sage: Coucou is Axiom("Coucou")
            True
            sage: TestSuite(Coucou).run()
        """
        assert isinstance(name, str)
        self._name = name
        if name in all_axioms.__dict__:
            raise ValueError("There already exist an axiom with name %s"%name)
        setattr(all_axioms, name, self)
        self._index = len(all_axioms.__dict__)

    def index(self):
        """
        Return the index of this axiom.

        OUTPUT:

        ``i`` such that ``self`` is the ``i``-th axiom defined in this session.

        This is used internally in the category pretty printing
        heuristic. See
        :meth:`sage.categories.category_with_axiom.CategoryWithAxiom._repr_object_names`.

        EXAMPLES::

            sage: axioms.Flying.index()
            1
        """
        return self._index

    def __repr__(self):
        """
        Return the name of this idiom as string representation.

        EXAMPLES::

            sage: axioms.FiniteDimensional.__repr__()
            'FiniteDimensional'

            sage: axioms.FiniteDimensional
            FiniteDimensional

        ..  NOTE::

            This is also used to recover the name of the idiom which
            is used intensively internaly. Hence we do NOT want the
            renaming feature of SageObject's and implement
            ``__repr__`` rather ``_repr_``.
        """
        return self._name

    def _repr_object_names(self, base_category):
        """
        Return the name of the objects in `base_category` satisfying this axiom.

        INPUT:

        - ``base_category`` -- a category

        EXAMPLES:

        This method is used to pretty print categories::

            sage: Modules(ZZ).FiniteDimensional()
            Category of finite dimensional modules over Integer Ring

        This default implementation prefixes the name of the objects
        in `base_category` with the name of this axiom in lower case,
        with words separated by spaces:

            sage: axioms.FiniteDimensional._repr_object_names(Groups())
            'finite dimensional groups'

        Subclasses of `Axiom` can customize this method to their taste::

            sage: axioms.WithBasis._repr_object_names(Modules(ZZ))
            'modules with basis over Integer Ring'

        They can safely assume that, unlike in the second example above,
        the axiom is relevant to the `base_category`.
        """
        base_object_names = base_category._repr_object_names()
        return uncamelcase(str(self)) + " " + base_object_names
