r"""
Finite complex reflection groups
----------------------------------

Let `V` be a finite-dimensional complex vector space. A reflection of
`V` is an operator `r \in \operatorname{GL}(V)` that has finite order
and fixes pointwise a hyperplane in `V`.

For more definitions and classification types of finite complex
reflection groups, see :wikipedia:`Complex_reflection_group`.

The point of entry to work with reflection groups is :func:`~sage.combinat.root_system.reflection_group_real.ReflectionGroup`
which can be used with finite Cartan-Killing types::

    sage: ReflectionGroup(['A',2])                                      # optional - gap3
    Irreducible real reflection group of rank 2 and type A2
    sage: ReflectionGroup(['F',4])                                      # optional - gap3
    Irreducible real reflection group of rank 4 and type F4
    sage: ReflectionGroup(['H',3])                                      # optional - gap3
    Irreducible real reflection group of rank 3 and type H3

or with Shephard-Todd types::

    sage: ReflectionGroup((1,1,3))                                      # optional - gap3
    Irreducible real reflection group of rank 2 and type A2
    sage: ReflectionGroup((2,1,3))                                      # optional - gap3
    Irreducible real reflection group of rank 3 and type B3
    sage: ReflectionGroup((3,1,3))                                      # optional - gap3
    Irreducible complex reflection group of rank 3 and type G(3,1,3)
    sage: ReflectionGroup((4,2,3))                                      # optional - gap3
    Irreducible complex reflection group of rank 3 and type G(4,2,3)
    sage: ReflectionGroup(4)                                            # optional - gap3
    Irreducible complex reflection group of rank 2 and type ST4
    sage: ReflectionGroup(31)                                           # optional - gap3
    Irreducible complex reflection group of rank 4 and type ST31

Also reducible types are allowed using concatenation::

    sage: ReflectionGroup(['A',3],(4,2,3))                              # optional - gap3
    Reducible complex reflection group of rank 6 and type A3 x G(4,2,3)

Some special cases also occur, among them are::

    sage: W = ReflectionGroup((2,2,2)); W                               # optional - gap3
    Reducible real reflection group of rank 2 and type A1 x A1
    sage: W = ReflectionGroup((2,2,3)); W                               # optional - gap3
    Irreducible real reflection group of rank 3 and type A3

.. WARNING:: Uses the GAP3 package *Chevie* which is available as an
             experimental package (installed by ``sage -i gap3``) or to
             download by hand from `Jean Michel's website <http://webusers.imj-prg.fr/~jean.michel/gap3/>`_.

A guided tour
-------------

We start with the example type `B_2`::

    sage: W = ReflectionGroup(['B',2]); W                               # optional - gap3
    Irreducible real reflection group of rank 2 and type B2

Most importantly, observe that the group elements are usually represented
by permutations of the roots::

    sage: for w in W: print(w)                                          # optional - gap3
    ()
    (1,3)(2,6)(5,7)
    (1,5)(2,4)(6,8)
    (1,7,5,3)(2,4,6,8)
    (1,3,5,7)(2,8,6,4)
    (2,8)(3,7)(4,6)
    (1,7)(3,5)(4,8)
    (1,5)(2,6)(3,7)(4,8)

This has the drawback that one can hardly see anything. Usually, one
would look at elements with either of the following methods::

    sage: for w in W: w.reduced_word()                                  # optional - gap3
    []
    [2]
    [1]
    [1, 2]
    [2, 1]
    [2, 1, 2]
    [1, 2, 1]
    [1, 2, 1, 2]

    sage: for w in W: w.reduced_word_in_reflections()                   # optional - gap3
    []
    [2]
    [1]
    [1, 2]
    [1, 4]
    [3]
    [4]
    [1, 3]

    sage: for w in W: w.reduced_word(); w.to_matrix(); print("")        # optional - gap3
    []
    [1 0]
    [0 1]
    <BLANKLINE>
    [2]
    [ 1  1]
    [ 0 -1]
    <BLANKLINE>
    [1]
    [-1  0]
    [ 2  1]
    <BLANKLINE>
    [1, 2]
    [-1 -1]
    [ 2  1]
    <BLANKLINE>
    [2, 1]
    [ 1  1]
    [-2 -1]
    <BLANKLINE>
    [2, 1, 2]
    [ 1  0]
    [-2 -1]
    <BLANKLINE>
    [1, 2, 1]
    [-1 -1]
    [ 0  1]
    <BLANKLINE>
    [1, 2, 1, 2]
    [-1  0]
    [ 0 -1]
    <BLANKLINE>

The standard references for actions of complex reflection groups have
the matrices acting on the right, so::

    sage: W.simple_reflection(1).to_matrix()                            # optional - gap3
    [-1  0]
    [ 2  1]

sends the simple root `\alpha_0`, or ``(1,0)`` in vector notation, to
its negative, while sending `\alpha_1` to `2\alpha_0+\alpha_1`.

.. TODO::

    - properly provide root systems for real reflection groups
    - element class should be unique to be able to work with large groups
      without creating elements multiple times
    - ``is_shephard_group``, ``is_generalized_coxeter_group``
    - exponents and coexponents
    - coinvariant ring:

      * fake degrees from Torsten Hoge
      * operation of linear characters on all characters
      * harmonic polynomials

    - linear forms for hyperplanes
    - field of definition
    - intersection lattice and characteristic polynomial::

        X = [ alpha(t) for t in W.distinguished_reflections() ]
        X = Matrix(CF,X).transpose()
        Y = Matroid(X)

    - linear characters
    - permutation pi on irreducibles
    - hyperplane orbits (76.13 in Gap Manual)
    - improve invariant_form with a code similar to the one in
      ``reflection_group_real.py``
    - add a method ``reflection_to_root`` or
      ``distinguished_reflection_to_positive_root``
    - diagrams in ASCII-art (76.15)
    - standard (BMR) presentations
    - character table directly from Chevie
    - ``GenericOrder`` (76.20), ``TorusOrder`` (76.21)
    - correct fundamental invariants for `G_34`, check the others
    - copy hardcoded data (degrees, invariants, braid relations...) into sage
    - add other hardcoded data from the tables in chevie (location is
      SAGEDIR/local/gap3/gap-jm5-2015-02-01/gap3/pkg/chevie/tbl):
      basic derivations, discriminant, ...
    - transfer code for ``reduced_word_in_reflections`` into Gap4 or Sage
    - list of reduced words for an element
    - list of reduced words in reflections for an element
    - Hurwitz action?
    - :meth:`is_crystallographic` should be hardcoded

AUTHORS:

- Christian Stump (2015): initial version
"""
#*****************************************************************************
#       Copyright (C) 2011-2016 Christian Stump <christian.stump at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.flatten import flatten
from sage.misc.misc_c import prod
from sage.categories.category import Category
from sage.categories.permutation_groups import PermutationGroups
from sage.categories.complex_reflection_groups import ComplexReflectionGroups
from sage.categories.coxeter_groups import CoxeterGroups
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.sets.family import Family
from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.perm_gps.permgroup import PermutationGroup_generic
from sage.rings.all import ZZ, QQ
from sage.matrix.all import Matrix, identity_matrix
from sage.matrix.matrix import is_Matrix
from sage.interfaces.gap3 import gap3
from sage.rings.universal_cyclotomic_field import E
from sage.arith.misc import lcm
from sage.modules.free_module_element import vector
from sage.combinat.root_system.cartan_matrix import CartanMatrix

from sage.misc.sage_eval import sage_eval

class ComplexReflectionGroup(UniqueRepresentation, PermutationGroup_generic):
    """
    A complex reflection group given as a permutation group.

    .. SEEALSO::

        :func:`ReflectionGroup`
    """
    def __init__(self, W_types, index_set=None, hyperplane_index_set=None, reflection_index_set=None):
        r"""
        TESTS::

            sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
            sage: W = ComplexReflectionGroups().example()               # optional - gap3
            sage: TestSuite(W).run()                                    # optional - gap3
        """
        W_components = []
        reflection_type = []
        for W_type in W_types:
            if W_type == (1,1,1):
                raise ValueError("the one element group is not considered a reflection group")
            elif W_type in ZZ:
                call_str = 'ComplexReflectionGroup(%s)'%W_type
            elif isinstance(W_type,CartanMatrix):
                call_str = 'PermRootGroup(IdentityMat(%s),%s)'%(W_type._rank,str(W_type._M._gap_()))
            elif is_Matrix(W_type):
                call_str = 'PermRootGroup(IdentityMat(%s),%s)'%(W_type._rank,str(W_type._gap_()))
            elif W_type in ZZ or ( isinstance(W_type, tuple) and len(W_type) == 3 ):
                call_str = 'ComplexReflectionGroup%s'%str(W_type)
            else:
                if W_type[0] == "I":
                    call_str = 'CoxeterGroup("I",2,%s)'%W_type[1]
                else:
                    call_str = 'CoxeterGroup("%s",%s)'%W_type

            W_components.append(gap3(call_str))
            X = list(W_components[-1].ReflectionType())
            if len(X) > 1:
                raise ValueError("input data %s is invalid"%W_type)
            X = X[0]
            type_dict = dict()
            type_dict["series"] = X.series.sage()
            type_dict["rank"] = X.rank.sage()
            type_dict["indices"] = X.indices.sage()
            if hasattr(X.ST,"sage"):
                type_dict["ST"] = X.ST.sage()
            elif hasattr(X.p,"sage") and hasattr(X.q,"sage"):
                type_dict["ST"] = ( X.p.sage(), X.q.sage(), X.rank.sage() )
            elif hasattr(X.bond,"sage"):
                type_dict["bond"] = X.bond.sage()
            if type_dict["series"] == "B" and (X.cartanType.sage() == 1 or X.indices.sage() == [2,1]):
                type_dict["series"] = "C"
            reflection_type.append( type_dict )

        self._type = reflection_type
        self._gap_group = prod(W_components)
        generators = [str(x) for x in self._gap_group.generators]
        self._index_set = index_set
        self._hyperplane_index_set = hyperplane_index_set
        self._reflection_index_set = reflection_index_set

        self._conjugacy_classes = {}
        self._conjugacy_classes_representatives = None
        self._reflection_representation = None

        self._rank = self._gap_group.rank.sage()
        if len(generators) == self._rank:
            category = ComplexReflectionGroups().Finite().WellGenerated()
            if all(str(W_comp).find('CoxeterGroup') >= 0 for W_comp in W_components):
                category = Category.join([category,CoxeterGroups()])
        else:
            category = ComplexReflectionGroups().Finite()
        if len(self._type) == 1:
            category = category.Irreducible()

        category = Category.join([category,PermutationGroups()]).Finite()

        PermutationGroup_generic.__init__(self, gens=generators,
                                          canonicalize=False,
                                          category=category)

        l_set = range(1, len(self.gens())+1)
        if self._index_set is None:
            self._index_set = tuple(l_set)
        else:
            if len(self._index_set) != len(l_set):
                raise ValueError("the given index set (= %s) does not have the right size"%self._index_set.values())
        self._index_set_inverse = {i: ii for ii,i in enumerate(self._index_set)}
        Nstar_set = range(1,self.number_of_reflection_hyperplanes()+1)
        if self._hyperplane_index_set is None:
            self._hyperplane_index_set = tuple(Nstar_set)
        else:
            if len(self._hyperplane_index_set) != len(Nstar_set):
                raise ValueError("the given hyperplane index set (= %s) does not have the right size"%self._index_set.values())
        self._hyperplane_index_set_inverse = {i: ii for ii,i in enumerate(self._hyperplane_index_set)}

        N_set = range(1, self.number_of_reflections()+1)
        if self._reflection_index_set is None:
            self._reflection_index_set = tuple(N_set)
        else:
            if len(self._reflection_index_set) != len(N_set):
                raise ValueError("the given reflection index set (= %s) does not have the right size"%self._index_set.values())
        self._reflection_index_set_inverse = {i: ii for ii,i in enumerate(self._reflection_index_set)}

    def _irrcomp_repr_(self,W_type):
        r"""
        Return the string representation of an irreducible component
        of ``self``.

        TESTS::

            sage: W = ReflectionGroup(25,[4,1,4],[1,1,4],[5,5,2]); W    # optional - gap3
            Reducible complex reflection group of rank 12 and type ST25 x G(4,1,4) x A3 x I2(5)
            sage: for W_type in W._type: print(W._irrcomp_repr_(W_type))    # optional - gap3
            ST25
             G(4,1,4)
             A3
             I2(5)
        """
        type_str = ''
        if "ST" in W_type:
            if W_type["ST"] in ZZ:
                type_str += "ST" + str(W_type["ST"])
            else:
                type_str += 'G' + str(W_type["ST"]).replace(' ','')
        else:
            type_str += str(W_type["series"])
            if W_type["series"] == "I":
                type_str += '2(' + str(W_type["bond"]) + ')'
            else:
                type_str += str(W_type["rank"])
        return type_str

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(25, [4,1,4],[1,1,4],[5,5,2]); W   # optional - gap3
            Reducible complex reflection group of rank 12 and type ST25 x G(4,1,4) x A3 x I2(5)
        """
        type_str = ''
        for W_type in self._type:
            type_str += self._irrcomp_repr_(W_type)
            type_str += ' x '
        type_str = type_str[:-3]
        return 'Reducible complex reflection group of rank %s and type %s'%(self._rank,type_str)

    def __iter__(self):
        r"""
        Return an iterator going through all elements in ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: for w in W: w                                         # optional - gap3
            ()
            (1,3)(2,5)(4,6)
            (1,4)(2,3)(5,6)
            (1,6,2)(3,5,4)
            (1,2,6)(3,4,5)
            (1,5)(2,4)(3,6)
        """
        from sage.combinat.root_system.reflection_group_c import iterator_tracking_words
        for w,word in iterator_tracking_words(self):
            w._reduced_word = word
            yield w

    @cached_method
    def index_set(self):
        r"""
        Return the index set of the simple reflections of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))                          # optional - gap3
            sage: W.index_set()                                         # optional - gap3
            (1, 2, 3)
            sage: W = ReflectionGroup((1,1,4), index_set=[1,3,'asdf'])  # optional - gap3
            sage: W.index_set()                                         # optional - gap3
            (1, 3, 'asdf')
            sage: W = ReflectionGroup((1,1,4), index_set=('a', 'b', 'c'))   # optional - gap3
            sage: W.index_set()                                         # optional - gap3
            ('a', 'b', 'c')
        """
        return self._index_set

    def simple_reflection(self, i):
        r"""
        Return the ``i``-th simple reflection of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.simple_reflection(1)                                # optional - gap3
            (1,4)(2,3)(5,6)
            sage: W.simple_reflections()                                # optional - gap3
            Finite family {1: (1,4)(2,3)(5,6), 2: (1,3)(2,5)(4,6)}
        """
        return self.gens()[self._index_set_inverse[i]]

    def series(self):
        r"""
        Return the series of the classification type to which ``self``
        belongs.

        For real reflection groups, these are the Cartan-Killing
        classification types "A","B","C","D","E","F","G","H","I", and
        for complx non-real reflection groups these are the
        Shephard-Todd classification type "ST".

        EXAMPLES::

            sage: ReflectionGroup((1,1,3)).series()                     # optional - gap3
            ['A']
            sage: ReflectionGroup((3,1,3)).series()                     # optional - gap3
            ['ST']
        """
        return [self._type[i]['series'] for i in range(len(self._type))]

    @cached_method
    def hyperplane_index_set(self):
        r"""
        Return the index set of the hyperplanes of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))                          # optional - gap3
            sage: W.hyperplane_index_set()                              # optional - gap3
            (1, 2, 3, 4, 5, 6)
            sage: W = ReflectionGroup((1,1,4), hyperplane_index_set=[1,3,'asdf',7,9,11])    # optional - gap3
            sage: W.hyperplane_index_set()                              # optional - gap3
            (1, 3, 'asdf', 7, 9, 11)
            sage: W = ReflectionGroup((1,1,4),hyperplane_index_set=('a','b','c','d','e','f'))   # optional - gap3
            sage: W.hyperplane_index_set()                              # optional - gap3
            ('a', 'b', 'c', 'd', 'e', 'f')
        """
        return self._hyperplane_index_set

    @cached_method
    def distinguished_reflections(self):
        r"""
        Return a finite family containing the distinguished reflections
        of ``self`` indexed by :meth:`hyperplane_index_set`.

        These are the reflections in ``self`` acting on the complement
        of the fixed hyperplane `H` as `\operatorname{exp}(2 \pi i / n)`,
        where `n` is the order of the reflection subgroup fixing `H`.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.distinguished_reflections()                         # optional - gap3
            Finite family {1: (1,4)(2,3)(5,6), 2: (1,3)(2,5)(4,6), 3: (1,5)(2,4)(3,6)}

            sage: W = ReflectionGroup((1,1,3),hyperplane_index_set=['a','b','c'])   # optional - gap3
            sage: W.distinguished_reflections()                         # optional - gap3
            Finite family {'a': (1,4)(2,3)(5,6), 'c': (1,5)(2,4)(3,6), 'b': (1,3)(2,5)(4,6)}

            sage: W = ReflectionGroup((3,1,1))                          # optional - gap3
            sage: W.distinguished_reflections()                         # optional - gap3
            Finite family {1: (1,2,3)}

            sage: W = ReflectionGroup((1,1,3),(3,1,2))                  # optional - gap3
            sage: W.distinguished_reflections()                         # optional - gap3
            Finite family {1: (1,6)(2,5)(7,8), 2: (1,5)(2,7)(6,8),
             3: (3,9,15)(4,10,16)(12,17,23)(14,18,24)(20,25,29)(21,22,26)(27,28,30),
             4: (3,11)(4,12)(9,13)(10,14)(15,19)(16,20)(17,21)(18,22)(23,27)(24,28)(25,26)(29,30),
             5: (1,7)(2,6)(5,8),
             6: (3,19)(4,25)(9,11)(10,17)(12,28)(13,15)(14,30)(16,18)(20,27)(21,29)(22,23)(24,26),
             7: (4,21,27)(10,22,28)(11,13,19)(12,14,20)(16,26,30)(17,18,25)(23,24,29),
             8: (3,13)(4,24)(9,19)(10,29)(11,15)(12,26)(14,21)(16,23)(17,30)(18,27)(20,22)(25,28)}
        """
        # makes sure that the simple reflections come first
        gens = self.gens()
        R = [t for t in gens]
        # Then import all distinguished reflections from gap,
        #   the Set is used as every such appears multiple times.
        for r in self._gap_group.Reflections():
            t = self(str(r))
            if t not in R:
                R.append(t)
        return Family(self._hyperplane_index_set,
                      lambda i: R[self._hyperplane_index_set_inverse[i]])

    def distinguished_reflection(self, i):
        r"""
        Return the ``i``-th distinguished reflection of ``self``.

        These are the reflections in ``self`` acting on the complement
        of the fixed hyperplane `H` as `\operatorname{exp}(2 \pi i / n)`,
        where `n` is the order of the reflection subgroup fixing `H`.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.distinguished_reflection(1)                         # optional - gap3
            (1,4)(2,3)(5,6)
            sage: W.distinguished_reflection(2)                         # optional - gap3
            (1,3)(2,5)(4,6)
            sage: W.distinguished_reflection(3)                         # optional - gap3
            (1,5)(2,4)(3,6)

            sage: W = ReflectionGroup((3,1,1),hyperplane_index_set=['a'])   # optional - gap3
            sage: W.distinguished_reflection('a')                       # optional - gap3
            (1,2,3)

            sage: W = ReflectionGroup((1,1,3),(3,1,2))                  # optional - gap3
            sage: for i in range(W.number_of_reflection_hyperplanes()): # optional - gap3
            ....:     W.distinguished_reflection(i+1)                   # optional - gap3
            (1,6)(2,5)(7,8)
            (1,5)(2,7)(6,8)
            (3,9,15)(4,10,16)(12,17,23)(14,18,24)(20,25,29)(21,22,26)(27,28,30)
            (3,11)(4,12)(9,13)(10,14)(15,19)(16,20)(17,21)(18,22)(23,27)(24,28)(25,26)(29,30)
            (1,7)(2,6)(5,8)
            (3,19)(4,25)(9,11)(10,17)(12,28)(13,15)(14,30)(16,18)(20,27)(21,29)(22,23)(24,26)
            (4,21,27)(10,22,28)(11,13,19)(12,14,20)(16,26,30)(17,18,25)(23,24,29)
            (3,13)(4,24)(9,19)(10,29)(11,15)(12,26)(14,21)(16,23)(17,30)(18,27)(20,22)(25,28)
        """
        return self.distinguished_reflections()[i]

    @cached_method
    def reflection_hyperplanes(self, as_linear_functionals=False):
        r"""
        Return the list of all reflection hyperplanes of ``self``,
        either as a codimension 1 space, or as its linear functional.

        INPUT:

        - ``as_linear_functionals`` -- (default:``False``) flag whether
          to return the hyperplane or its linear functional in the basis
          dual to the given root basis

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: for H in W.reflection_hyperplanes(): H                # optional - gap3
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 2]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [  1 1/2]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -1]

            sage: for H in W.reflection_hyperplanes(as_linear_functionals=True): H  # optional - gap3
            (1, -1/2)
            (1, -2)
            (1, 1)


            sage: W = ReflectionGroup((2,1,2))                          # optional - gap3
            sage: for H in W.reflection_hyperplanes(): H                # optional - gap3
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 1]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [  1 1/2]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]

            sage: for H in W.reflection_hyperplanes(as_linear_functionals=True): H  # optional - gap3
            (1, -1)
            (1, -2)
            (0, 1)
            (1, 0)
        """
        Hs = []
        for r in self.distinguished_reflections():
            mat = (r.to_matrix().transpose() - identity_matrix(self.rank()))
            if as_linear_functionals:
                Hs.append( mat.row_space().gen() )
            else:
                Hs.append( mat.right_kernel() )
        return Family(self._hyperplane_index_set,
                      lambda i: Hs[self._hyperplane_index_set_inverse[i]])

    def reflection_hyperplane(self, i, as_linear_functional=False):
        r"""
        Return the ``i``-th reflection hyperplane of ``self``.

        The ``i``-th reflection hyperplane corresponds to the ``i``
        distinguished reflection.

        INPUT:

        - ``i`` -- an index in the index set
        - ``as_linear_functionals`` -- (default:``False``) flag whether
          to return the hyperplane or its linear functional in the basis
          dual to the given root basis

        EXAMPLES::

            sage: W = ReflectionGroup((2,1,2))                          # optional - gap3
            sage: W.reflection_hyperplane(3)                            # optional - gap3
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]

        One can ask for the result as a linear form::

            sage: W.reflection_hyperplane(3, True)                      # optional - gap3
            (0, 1)
        """
        return self.reflection_hyperplanes(as_linear_functionals=as_linear_functional)[i]

    @cached_method
    def reflection_index_set(self):
        r"""
        Return the index set of the reflections of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))                          # optional - gap3
            sage: W.reflection_index_set()                              # optional - gap3
            (1, 2, 3, 4, 5, 6)
            sage: W = ReflectionGroup((1,1,4), reflection_index_set=[1,3,'asdf',7,9,11])    # optional - gap3
            sage: W.reflection_index_set()                              # optional - gap3
            (1, 3, 'asdf', 7, 9, 11)
            sage: W = ReflectionGroup((1,1,4), reflection_index_set=('a','b','c','d','e','f'))  # optional - gap3
            sage: W.reflection_index_set()                              # optional - gap3
            ('a', 'b', 'c', 'd', 'e', 'f')
        """
        return self._reflection_index_set

    @cached_method
    def reflections(self):
        r"""
        Return a finite family containing the reflections of ``self``,
        indexed by :meth:`self.reflection_index_set`.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.reflections()                                       # optional - gap3
            Finite family {1: (1,4)(2,3)(5,6), 2: (1,3)(2,5)(4,6), 3: (1,5)(2,4)(3,6)}

            sage: W = ReflectionGroup((1,1,3),reflection_index_set=['a','b','c'])   # optional - gap3
            sage: W.reflections()                                       # optional - gap3
            Finite family {'a': (1,4)(2,3)(5,6), 'c': (1,5)(2,4)(3,6), 'b': (1,3)(2,5)(4,6)}

            sage: W = ReflectionGroup((3,1,1))                          # optional - gap3
            sage: W.reflections()                                       # optional - gap3
            Finite family {1: (1,2,3), 2: (1,3,2)}

            sage: W = ReflectionGroup((1,1,3),(3,1,2))                  # optional - gap3
            sage: W.reflections()                                       # optional - gap3
            Finite family {1: (1,6)(2,5)(7,8), 2: (1,5)(2,7)(6,8),
                           3: (3,9,15)(4,10,16)(12,17,23)(14,18,24)(20,25,29)(21,22,26)(27,28,30),
                           4: (3,11)(4,12)(9,13)(10,14)(15,19)(16,20)(17,21)(18,22)(23,27)(24,28)(25,26)(29,30),
                           5: (1,7)(2,6)(5,8),
                           6: (3,19)(4,25)(9,11)(10,17)(12,28)(13,15)(14,30)(16,18)(20,27)(21,29)(22,23)(24,26),
                           7: (4,21,27)(10,22,28)(11,13,19)(12,14,20)(16,26,30)(17,18,25)(23,24,29),
                           8: (3,13)(4,24)(9,19)(10,29)(11,15)(12,26)(14,21)(16,23)(17,30)(18,27)(20,22)(25,28),
                           9: (3,15,9)(4,16,10)(12,23,17)(14,24,18)(20,29,25)(21,26,22)(27,30,28),
                           10: (4,27,21)(10,28,22)(11,19,13)(12,20,14)(16,30,26)(17,25,18)(23,29,24)}
        """
        T = self.distinguished_reflections().values()
        for i in range(self.number_of_reflection_hyperplanes()):
            for j in range(2, T[i].order()):
                T.append(T[i]**j)
        return Family(self._reflection_index_set,
                      lambda i: T[self._reflection_index_set_inverse[i]])

    def reflection(self,i):
        r"""
        Return the ``i``-th reflection of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.reflection(1)                                       # optional - gap3
            (1,4)(2,3)(5,6)
            sage: W.reflection(2)                                       # optional - gap3
            (1,3)(2,5)(4,6)
            sage: W.reflection(3)                                       # optional - gap3
            (1,5)(2,4)(3,6)

            sage: W = ReflectionGroup((3,1,1),reflection_index_set=['a','b'])   # optional - gap3
            sage: W.reflection('a')                                     # optional - gap3
            (1,2,3)
            sage: W.reflection('b')                                     # optional - gap3
            (1,3,2)
        """
        return self.reflections()[i]

    def reflection_character(self):
        r"""
        Return the reflection characters of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.reflection_character()                              # optional - gap3
            [2, 0, -1]
        """
        return self._gap_group.ReflectionCharacter().sage()

    def is_crystallographic(self):
        r"""
        Return ``True`` if self is crystallographic.

        This is, if the field of definition is the rational field.

        .. TODO::

            Make this more robust and do not use the matrix
            representation of the simple reflections.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3)); W                       # optional - gap3
            Irreducible real reflection group of rank 2 and type A2
            sage: W.is_crystallographic()                               # optional - gap3
            True

            sage: W = ReflectionGroup((2,1,3)); W                       # optional - gap3
            Irreducible real reflection group of rank 3 and type B3
            sage: W.is_crystallographic()                               # optional - gap3
            True

            sage: W = ReflectionGroup(23); W                            # optional - gap3
            Irreducible real reflection group of rank 3 and type H3
            sage: W.is_crystallographic()                               # optional - gap3
            False

            sage: W = ReflectionGroup((3,1,3)); W                       # optional - gap3
            Irreducible complex reflection group of rank 3 and type G(3,1,3)
            sage: W.is_crystallographic()                               # optional - gap3
            False

            sage: W = ReflectionGroup((4,2,2)); W                       # optional - gap3
            Irreducible complex reflection group of rank 2 and type G(4,2,2)
            sage: W.is_crystallographic()                               # optional - gap3
            False
        """
        return self.is_real() and all(t.to_matrix().base_ring() is QQ for t in self.simple_reflections())

    def _element_class(self):
        r"""
        A temporary workaround for compatibility with Sage's
        permutation groups.

        TESTS::

            sage: W = ReflectionGroup(23)                               # optional - gap3
            sage: W._element_class() is W.element_class                 # optional - gap3
            True
        """
        return self.element_class

    def number_of_irreducible_components(self):
        r"""
        Return the number of irreducible components of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.number_of_irreducible_components()                  # optional - gap3
            1

            sage: W = ReflectionGroup((1,1,3),(2,1,3))                  # optional - gap3
            sage: W.number_of_irreducible_components()                  # optional - gap3
            2
        """
        return len(self._type)

    def irreducible_components(self):
        r"""
        Return a list containing the irreducible components of ``self``
        as finite reflection groups.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.irreducible_components()                            # optional - gap3
            [Irreducible real reflection group of rank 2 and type A2]

            sage: W = ReflectionGroup((1,1,3),(2,1,3))                  # optional - gap3
            sage: W.irreducible_components()                            # optional - gap3
            [Irreducible real reflection group of rank 2 and type A2,
            Irreducible real reflection group of rank 3 and type B3]
        """
        from sage.combinat.root_system.reflection_group_real import ReflectionGroup
        irr_comps = []
        for W_type in self._type:
            if W_type["series"] in ["A","B","D","E","F","G","H","I"]:
                W_str = (W_type["series"],W_type["rank"])
            elif "ST" in W_type:
                W_str = W_type["ST"]
            irr_comps.append(ReflectionGroup(W_str))
        return irr_comps

    @cached_method
    def conjugacy_classes_representatives(self):
        r"""
        Return the shortest representatives of the conjugacy classes of
        ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: [w.reduced_word() for w in W.conjugacy_classes_representatives()] # optional - gap3
            [[], [1], [1, 2]]

            sage: W = ReflectionGroup((1,1,4))                          # optional - gap3
            sage: [w.reduced_word() for w in W.conjugacy_classes_representatives()] # optional - gap3
            [[], [1], [1, 3], [1, 2], [1, 3, 2]]

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: [w.reduced_word() for w in W.conjugacy_classes_representatives()] # optional - gap3
            [[], [1], [1, 1], [2, 1, 2, 1], [2, 1, 2, 1, 1],
             [2, 1, 1, 2, 1, 1], [2], [1, 2], [1, 1, 2]]

            sage: W = ReflectionGroup(23)                               # optional - gap3
            sage: [w.reduced_word() for w in W.conjugacy_classes_representatives()] # optional - gap3
                [[],
                 [1],
                 [1, 2],
                 [1, 3],
                 [2, 3],
                 [1, 2, 3],
                 [1, 2, 1, 2],
                 [1, 2, 1, 2, 3],
                 [1, 2, 1, 2, 3, 2, 1, 2, 3],
                 [1, 2, 1, 2, 1, 3, 2, 1, 2, 1, 3, 2, 1, 2, 3]]
        """
        S = str(gap3('List(ConjugacyClasses(%s),Representative)'%self._gap_group._name))
        exec('_conjugacy_classes_representatives=' + _gap_return(S))
        return _conjugacy_classes_representatives

    def conjugacy_classes(self):
        r"""
        Return the conjugacy classes of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: for C in W.conjugacy_classes(): sorted(C)             # optional - gap3
            [()]
            [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6), (1,5)(2,4)(3,6)]
            [(1,2,6)(3,4,5), (1,6,2)(3,5,4)]

            sage: W = ReflectionGroup((1,1,4))                          # optional - gap3
            sage: sum(len(C) for C in W.conjugacy_classes()) == W.cardinality() # optional - gap3
            True

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: sum(len(C) for C in W.conjugacy_classes()) == W.cardinality() # optional - gap3
            True

            sage: W = ReflectionGroup(23)                               # optional - gap3
            sage: sum(len(C) for C in W.conjugacy_classes()) == W.cardinality() # optional - gap3
            True
       """
        return Family(self.conjugacy_classes_representatives(),
                      lambda w: w.conjugacy_class())

    def rank(self):
        r"""
        Return the rank of ``self``.

        This is the dimension of the underlying vector space.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.rank()                                              # optional - gap3
            2
            sage: W = ReflectionGroup((2,1,3))                          # optional - gap3
            sage: W.rank()                                              # optional - gap3
            3
            sage: W = ReflectionGroup((4,1,3))                          # optional - gap3
            sage: W.rank()                                              # optional - gap3
            3
            sage: W = ReflectionGroup((4,2,3))                          # optional - gap3
            sage: W.rank()                                              # optional - gap3
            3
        """
        return self._rank

    @cached_method
    def degrees(self):
        r"""
        Return the degrees of ``self`` ordered within each irreducible
        component of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))                          # optional - gap3
            sage: W.degrees()                                           # optional - gap3
            (2, 3, 4)

            sage: W = ReflectionGroup((2,1,4))                          # optional - gap3
            sage: W.degrees()                                           # optional - gap3
            (2, 4, 6, 8)

            sage: W = ReflectionGroup((4,1,4))                          # optional - gap3
            sage: W.degrees()                                           # optional - gap3
            (4, 8, 12, 16)

            sage: W = ReflectionGroup((4,2,4))                          # optional - gap3
            sage: W.degrees()                                           # optional - gap3
            (4, 8, 8, 12)

            sage: W = ReflectionGroup((4,4,4))                          # optional - gap3
            sage: W.degrees()                                           # optional - gap3
            (4, 4, 8, 12)

        Examples of reducible types::

            sage: W = ReflectionGroup((1,1,4), (3,1,2)); W              # optional - gap3
            Reducible complex reflection group of rank 5 and type A3 x G(3,1,2)
            sage: W.degrees()                                           # optional - gap3
            (2, 3, 4, 3, 6)

            sage: W = ReflectionGroup((1,1,4), (6,1,12), 23)            # optional - gap3 # fails in GAP3
            sage: W.degrees()                                           # optional - gap3
            (2, 3, 4, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 2, 6, 10)
        """
        if self.is_irreducible():
            try:
                return tuple(sorted(self._gap_group.degrees.sage()))
            except AttributeError:
                return tuple(sorted(self._gap_group.ReflectionDegrees().sage()))
        else:
            return sum([comp.degrees() for comp in self.irreducible_components()],tuple())

    @cached_method
    def codegrees(self):
        r"""
        Return the codegrees of ``self`` ordered within each irreducible
        component of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,4))                          # optional - gap3
            sage: W.codegrees()                                         # optional - gap3
            (2, 1, 0)

            sage: W = ReflectionGroup((2,1,4))                          # optional - gap3
            sage: W.codegrees()                                         # optional - gap3
            (6, 4, 2, 0)

            sage: W = ReflectionGroup((4,1,4))                          # optional - gap3
            sage: W.codegrees()                                         # optional - gap3
            (12, 8, 4, 0)

            sage: W = ReflectionGroup((4,2,4))                          # optional - gap3
            sage: W.codegrees()                                         # optional - gap3
            (12, 8, 4, 0)

            sage: W = ReflectionGroup((4,4,4))                          # optional - gap3
            sage: W.codegrees()                                         # optional - gap3
            (8, 8, 4, 0)

            sage: W = ReflectionGroup((1,1,4), (3,1,2))                 # optional - gap3
            sage: W.codegrees()                                         # optional - gap3
            (2, 1, 0, 3, 0)

            sage: W = ReflectionGroup((1,1,4), (6,1,12), 23)            # optional - gap3 # fails in GAP3
            sage: W.codegrees()                                         # optional - gap3
            (2, 1, 0, 66, 60, 54, 48, 42, 36, 30, 24, 18, 12, 6, 0, 8, 4, 0)
        """
        if self.is_irreducible():
            if self.is_well_generated():
                h = self.coxeter_number()
                return tuple([h-d for d in self.degrees()])
            else:
                return tuple(sorted(self._gap_group.ReflectionCoDegrees().sage(), reverse=True))
        else:
            return sum([comp.codegrees() for comp in self.irreducible_components()],tuple())

    @cached_method
    def reflection_eigenvalues_family(self):
        r"""
        Return the reflection eigenvalues of ``self`` as a finite family
        indexed by the class representatives of ``self``.

        OUTPUT:

        - list with entries `k/n` representing the eigenvalue `\zeta_n^k`.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.reflection_eigenvalues_family()                     # optional - gap3
            Finite family {(): [0, 0], (1,4)(2,3)(5,6): [1/2, 0], (1,6,2)(3,5,4): [1/3, 2/3]}

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: reflection_eigenvalues = W.reflection_eigenvalues_family()    # optional - gap3
            sage: for elt in sorted(reflection_eigenvalues.keys()):     # optional - gap3
            ....:     print('%s %s'%(elt, reflection_eigenvalues[elt])) # optional - gap3
            () [0, 0]
            (1,3,9)(2,4,10)(6,11,17)(8,12,18)(14,19,23)(15,16,20)(21,22,24) [1/3, 0]
            (1,3,9)(2,16,24)(4,20,21)(5,7,13)(6,12,23)(8,19,17)(10,15,22)(11,18,14) [1/3, 1/3]
            (1,5)(2,6)(3,7)(4,8)(9,13)(10,14)(11,15)(12,16)(17,21)(18,22)(19,20)(23,24) [1/2, 0]
            (1,7,3,13,9,5)(2,8,16,19,24,17)(4,14,20,11,21,18)(6,15,12,22,23,10) [1/6, 2/3]
            (1,9,3)(2,10,4)(6,17,11)(8,18,12)(14,23,19)(15,20,16)(21,24,22) [2/3, 0]
            (1,9,3)(2,20,22)(4,15,24)(5,7,13)(6,18,19)(8,23,11)(10,16,21)(12,14,17) [1/3, 2/3]
            (1,9,3)(2,24,16)(4,21,20)(5,13,7)(6,23,12)(8,17,19)(10,22,15)(11,14,18) [2/3, 2/3]
            (1,13,9,7,3,5)(2,14,24,18,16,11)(4,6,21,23,20,12)(8,22,17,15,19,10) [1/3, 5/6]

            sage: W = ReflectionGroup(23)                               # optional - gap3
            sage: reflection_eigenvalues = W.reflection_eigenvalues_family()    # optional - gap3
            sage: for elt in sorted(reflection_eigenvalues.keys()):     # optional - gap3
            ....:     print('%s %s'%(elt, reflection_eigenvalues[elt])) # optional - gap3
            () [0, 0, 0]
            (1,8,4)(2,21,3)(5,10,11)(6,18,17)(7,9,12)(13,14,15)(16,23,19)(20,25,26)(22,24,27)(28,29,30) [1/3, 2/3, 0]
            (1,16)(2,5)(4,7)(6,9)(8,10)(11,13)(12,14)(17,20)(19,22)(21,24)(23,25)(26,28)(27,29) [1/2, 0, 0]
            (1,16)(2,9)(3,18)(4,10)(5,6)(7,8)(11,14)(12,13)(17,24)(19,25)(20,21)(22,23)(26,29)(27,28) [1/2, 1/2, 0]
            (1,16)(2,17)(3,18)(4,19)(5,20)(6,21)(7,22)(8,23)(9,24)(10,25)(11,26)(12,27)(13,28)(14,29)(15,30) [1/2, 1/2, 1/2]
            (1,19,20,2,7)(3,6,11,13,9)(4,5,17,22,16)(8,12,15,14,10)(18,21,26,28,24)(23,27,30,29,25) [1/5, 4/5, 0]
            (1,20,7,19,2)(3,11,9,6,13)(4,17,16,5,22)(8,15,10,12,14)(18,26,24,21,28)(23,30,25,27,29) [2/5, 3/5, 0]
            (1,23,26,29,22,16,8,11,14,7)(2,10,4,9,18,17,25,19,24,3)(5,21,27,30,28,20,6,12,15,13) [1/10, 1/2, 9/10]
            (1,24,17,16,9,2)(3,12,13,18,27,28)(4,21,29,19,6,14)(5,25,26,20,10,11)(7,23,30,22,8,15) [1/6, 1/2, 5/6]
            (1,29,8,7,26,16,14,23,22,11)(2,9,25,3,4,17,24,10,18,19)(5,30,6,13,27,20,15,21,28,12) [3/10, 1/2, 7/10]
            """
        class_representatives = self.conjugacy_classes_representatives()
        Ev_list = self._gap_group.ReflectionEigenvalues().sage()
        return Family(class_representatives,
                      lambda w: Ev_list[class_representatives.index(w)])

    @cached_method
    def reflection_eigenvalues(self, w, is_class_representative=False):
        r"""
        Return the reflection eigenvalue of ``w`` in ``self``.

        INPUT:

        - ``is_class_representative`` -- boolean (default ``True``) whether to
          compute instead on the conjugacy class representative.

        .. SEEALSO:: :meth:`reflection_eigenvalues_family`

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: for w in W:                                           # optional - gap3
            ....:     print('%s %s'%(w.reduced_word(), W.reflection_eigenvalues(w)))    # optional - gap3
            [] [0, 0]
            [2] [1/2, 0]
            [1] [1/2, 0]
            [1, 2] [1/3, 2/3]
            [2, 1] [1/3, 2/3]
            [1, 2, 1] [1/2, 0]
        """
        if not is_class_representative:
            w = w.conjugacy_class_representative()
        return self.reflection_eigenvalues_family()[w]

    @cached_method
    def simple_roots(self):
        r"""
        Return the simple roots of ``self``.

        These are the roots corresponding to the simple reflections.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.simple_roots()                                      # optional - gap3
            Finite family {1: (1, 0), 2: (0, 1)}

            sage: W = ReflectionGroup((1,1,4), (2,1,2))                 # optional - gap3
            sage: W.simple_roots()                                      # optional - gap3
            Finite family {1: (1, 0, 0, 0, 0), 2: (0, 1, 0, 0, 0), 3: (0, 0, 1, 0, 0), 4: (0, 0, 0, 1, 0), 5: (0, 0, 0, 0, 1)}

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: W.simple_roots()                                      # optional - gap3
            Finite family {1: (1, 0), 2: (-1, 1)}

            sage: W = ReflectionGroup((1,1,4), (3,1,2))                 # optional - gap3
            sage: W.simple_roots()                                      # optional - gap3
            Finite family {1: (1, 0, 0, 0, 0), 2: (0, 1, 0, 0, 0), 3: (0, 0, 1, 0, 0), 4: (0, 0, 0, 1, 0), 5: (0, 0, 0, -1, 1)}
        """
        from sage.sets.family import Family
        return Family({ind:self.roots()[i] for i,ind in enumerate(self._index_set)})

    def simple_root(self, i):
        r"""
        Return the simple root with index ``i``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: W.simple_root(1)                                      # optional - gap3
            (1, 0, 0)
            sage: W.simple_root(2)                                      # optional - gap3
            (0, 1, 0)
            sage: W.simple_root(3)                                      # optional - gap3
            (0, 0, 1)

        TESTS::

            sage: W.simple_root(0)                                      # optional - gap3
            Traceback (most recent call last):
            ...
            KeyError: 0
        """
        return self.simple_roots()[i]

    @cached_method
    def simple_coroots(self):
        r"""
        Return the simple coroots of ``self``.

        These are the coroots corresponding to the simple reflections.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.simple_coroots()                                    # optional - gap3
            Finite family {1: (2, -1), 2: (-1, 2)}

            sage: W = ReflectionGroup((1,1,4), (2,1,2))                 # optional - gap3
            sage: W.simple_coroots()                                    # optional - gap3
            Finite family {1: (2, -1, 0, 0, 0), 2: (-1, 2, -1, 0, 0), 3: (0, -1, 2, 0, 0), 4: (0, 0, 0, 2, -2), 5: (0, 0, 0, -1, 2)}

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: W.simple_coroots()                                    # optional - gap3
            Finite family {1: (-2*E(3) - E(3)^2, 0), 2: (-1, 1)}

            sage: W = ReflectionGroup((1,1,4), (3,1,2))                 # optional - gap3
            sage: W.simple_coroots()                                    # optional - gap3
            Finite family {1: (2, -1, 0, 0, 0), 2: (-1, 2, -1, 0, 0), 3: (0, -1, 2, 0, 0), 4: (0, 0, 0, -2*E(3) - E(3)^2, 0), 5: (0, 0, 0, -1, 1)}
        """
        from sage.sets.family import Family
        coroots = self._gap_group.simpleCoroots.sage()
        for i,coroot in enumerate(coroots):
            coroot = vector(coroot)
            coroot.set_immutable()
            coroots[i] = coroot
        return Family({ind:coroots[i] for i,ind in enumerate(self.index_set())})

    def simple_coroot(self, i):
        r"""
        Return the simple root with index ``i``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: W.simple_coroot(1)                                    # optional - gap3
            (2, -1, 0)
        """
        return self.simple_coroots()[i]

    @cached_method
    def independent_roots(self):
        r"""
        Return a collection of simple roots generating the underlying
        vector space of ``self``.

        For well-generated groups, these are all simple roots.
        Otherwise, a linearly independent subset of the simple roots is
        chosen.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.independent_roots()                                 # optional - gap3
            Finite family {1: (1, 0), 2: (0, 1)}

            sage: W = ReflectionGroup((4,2,3))                          # optional - gap3
            sage: W.simple_roots()                                      # optional - gap3
            Finite family {1: (1, 0, 0), 2: (-E(4), 1, 0), 3: (-1, 1, 0), 4: (0, -1, 1)}
            sage: W.independent_roots()                                 # optional - gap3
            Finite family {1: (1, 0, 0), 2: (-E(4), 1, 0), 4: (0, -1, 1)}
        """
        Delta = self.simple_roots()
        if self.is_well_generated():
            return Delta

        from sage.sets.family import Family
        basis = dict()
        for i,ind in enumerate(self._index_set):
            vec = Delta[ind]
            if Matrix(basis.values()+[vec]).rank() == len(basis) + 1:
                basis[ind] = vec
        return Family(basis)

    @cached_method
    def base_change_matrix(self):
        r"""
        Return the base change from the standard basis of the vector
        space of ``self`` to the basis given by the independent roots of
        ``self``.

        .. TODO::

            For non-well-generated groups there is a conflict with
            construction of the matrix for an element.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.base_change_matrix()                                # optional - gap3
            [1 0]
            [0 1]

            sage: W = ReflectionGroup(23)                               # optional - gap3
            sage: W.base_change_matrix()                                # optional - gap3
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: W.base_change_matrix()                                # optional - gap3
            [1 0]
            [1 1]

            sage: W = ReflectionGroup((4,2,2))                          # optional - gap3
            sage: W.base_change_matrix()                                # optional - gap3
            [   1    0]
            [E(4)    1]
        """
        return Matrix( list(self.independent_roots()) ).inverse()

    @cached_method
    def roots(self):
        r"""
        Return all roots corresponding to all reflections of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.roots()                                             # optional - gap3
            [(1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (-1, -1)]

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: W.roots()                                             # optional - gap3
            [(1, 0), (-1, 1), (E(3), 0), (-E(3), 1), (0, 1), (1, -1),
             (0, E(3)), (1, -E(3)), (E(3)^2, 0), (-E(3)^2, 1),
             (E(3), -1), (E(3), -E(3)), (0, E(3)^2), (1, -E(3)^2),
             (-1, E(3)), (-E(3), E(3)), (E(3)^2, -1), (E(3)^2, -E(3)),
             (E(3), -E(3)^2), (-E(3)^2, E(3)), (-1, E(3)^2),
             (-E(3), E(3)^2), (E(3)^2, -E(3)^2), (-E(3)^2, E(3)^2)]

            sage: W = ReflectionGroup((4,2,2))                          # optional - gap3
            sage: W.roots()                                             # optional - gap3
            [(1, 0), (-E(4), 1), (-1, 1), (-1, 0), (E(4), 1), (1, 1),
             (0, -E(4)), (E(4), -1), (E(4), E(4)), (0, E(4)),
             (E(4), -E(4)), (0, 1), (1, -E(4)), (1, -1), (0, -1),
             (1, E(4)), (-E(4), 0), (-1, E(4)), (E(4), 0), (-E(4), E(4)),
             (-E(4), -1), (-E(4), -E(4)), (-1, -E(4)), (-1, -1)]

            sage: W = ReflectionGroup((1,1,4), (3,1,2))                 # optional - gap3
            sage: W.roots()                                             # optional - gap3
            [(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0),
             (0, 0, 0, 1, 0), (0, 0, 0, -1, 1), (1, 1, 0, 0, 0),
             (0, 1, 1, 0, 0), (1, 1, 1, 0, 0), (-1, 0, 0, 0, 0),
             (0, -1, 0, 0, 0), (0, 0, -1, 0, 0), (-1, -1, 0, 0, 0),
             (0, -1, -1, 0, 0), (-1, -1, -1, 0, 0), (0, 0, 0, E(3), 0),
             (0, 0, 0, -E(3), 1), (0, 0, 0, 0, 1), (0, 0, 0, 1, -1),
             (0, 0, 0, 0, E(3)), (0, 0, 0, 1, -E(3)), (0, 0, 0, E(3)^2, 0),
             (0, 0, 0, -E(3)^2, 1), (0, 0, 0, E(3), -1), (0, 0, 0, E(3), -E(3)),
             (0, 0, 0, 0, E(3)^2), (0, 0, 0, 1, -E(3)^2), (0, 0, 0, -1, E(3)),
             (0, 0, 0, -E(3), E(3)), (0, 0, 0, E(3)^2, -1),
             (0, 0, 0, E(3)^2, -E(3)), (0, 0, 0, E(3), -E(3)^2),
             (0, 0, 0, -E(3)^2, E(3)), (0, 0, 0, -1, E(3)^2),
             (0, 0, 0, -E(3), E(3)^2), (0, 0, 0, E(3)^2, -E(3)^2),
             (0, 0, 0, -E(3)^2, E(3)^2)]
        """
        roots = [vector(sage_eval(str(root).replace("^","**")))
                 for root in self._gap_group.roots]
        for v in roots:
            v.set_immutable()
        return roots

    @cached_method
    def braid_relations(self):
        r"""
        Return the braid relations of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.braid_relations()                                   # optional - gap3
            [[[2, 1, 2], [1, 2, 1]]]

            sage: W = ReflectionGroup((2,1,3))                          # optional - gap3
            sage: W.braid_relations()                                   # optional - gap3
            [[[2, 1, 2, 1], [1, 2, 1, 2]], [[3, 1], [1, 3]], [[3, 2, 3], [2, 3, 2]]]

            sage: W = ReflectionGroup((2,2,3))                          # optional - gap3
            sage: W.braid_relations()                                   # optional - gap3
            [[[2, 1, 2], [1, 2, 1]], [[3, 1], [1, 3]], [[3, 2, 3], [2, 3, 2]]]
        """
        return self._gap_group.BraidRelations().sage()

    def fundamental_invariants(self):
        r"""
        Return the fundamental invariants of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: W.fundamental_invariants()                            # optional - gap3
            (-2*x0^2 + 2*x0*x1 - 2*x1^2, 6*x0^2*x1 - 6*x0*x1^2)

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: W.fundamental_invariants()                            # optional - gap3
            (x0^3 + x1^3, x0^3*x1^3)
        """
        import re
        from sage.rings.polynomial.all import PolynomialRing

        if not self.is_irreducible():
            return sum([W.fundamental_invariants() for W in self.irreducible_components() ],tuple())

        I = [ str(p) for p in gap3('List(Invariants(%s),x->ApplyFunc(x,List([0..%s],i->Mvp(SPrint("x",i)))))'%(self._gap_group._name,self.rank()-1)) ]
        P = PolynomialRing(QQ,['x%s'%i for i in range(0,self.rank())])
        x = P.gens()
        for i in range(len(I)):
            I[i] = I[i].replace('^','**')
            I[i] = re.compile('E(\d\d*)').sub(r'E(\1)', I[i])
            I[i] = re.compile('(\d)E\(').sub(r'\1*E(', I[i])
            for j in range(len(x)):
                I[i] = I[i].replace('x%s'%j,'*x[%s]'%j)
            I[i] = I[i].replace("+*","+").replace("-*","-").replace("ER(5)","*(E(5)-E(5)**2-E(5)**3+E(5)**4)").lstrip("*")
        # sage_eval is used since eval kills the rational entries!
        I = [ sage_eval(p, locals={'x':x}) for p in I ]
        return tuple(sorted(I,lambda f,g: cmp(f.degree(),g.degree())))

    def cartan_matrix(self):
        r"""
        Return the Cartan matrix associated with ``self``.

        If ``self`` is crystallographic, the returned Cartan matrix is
        an instance of :class:`CartanMatrix`, and a normal matrix
        otherwise.

        Let `s_1, \ldots, s_n` be a set of reflections which generate
        ``self`` with associated simple roots `s_1,\ldots,s_n` and
        simple coroots `s^\vee_i`. Then the Cartan matrix `C = (c_{ij})`
        is given by `s^\vee_i(s_j)`. The Cartan matrix completely
        determines the reflection representation if the `s_i` are
        linearly independent.

        EXAMPLES::

            sage: ReflectionGroup(['A',4]).cartan_matrix()              # optional - gap3
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]

            sage: ReflectionGroup(['H',4]).cartan_matrix()              # optional - gap3
            [              2 E(5)^2 + E(5)^3               0               0]
            [E(5)^2 + E(5)^3               2              -1               0]
            [              0              -1               2              -1]
            [              0               0              -1               2]

            sage: ReflectionGroup(4).cartan_matrix()                    # optional - gap3
            [-2*E(3) - E(3)^2           E(3)^2]
            [         -E(3)^2 -2*E(3) - E(3)^2]

            sage: ReflectionGroup((4,2,2)).cartan_matrix()              # optional - gap3
            [       2  -2*E(4)       -2]
            [    E(4)        2 1 - E(4)]
            [      -1 1 + E(4)        2]
        """
        # an alternative implementation is
        # Matrix(tuple(W.simple_coroots()))*Matrix(tuple(W.simple_roots())).transpose()
        # this should be implemented once we get the simple roots in an easy way
        if self.is_crystallographic():
            from sage.combinat.root_system.cartan_matrix import CartanMatrix as CartanMat
        else:
            from sage.matrix.all import Matrix as CartanMat
        return CartanMat(self._gap_group.CartanMat().sage())

    def invariant_form(self, brute_force=False):
        r"""
        Return the form that is invariant under the action of ``self``.

        This is unique only up to a global scalar on the irreducible
        components.

        INPUT:

        - ``brute_force`` -- if ``True``, the computation is done by
          applying the Reynolds operator; this is, the invariant form
          of `e_i` and `e_j` is computed as the sum
          `\langle w(e_i), w(e_j)\rangle`, where
          `\langle \cdot, \cdot\rangle` is the standard scalar product

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])                          # optional - gap3
            sage: F = W.invariant_form(); F                             # optional - gap3
            [   1 -1/2    0]
            [-1/2    1 -1/2]
            [   0 -1/2    1]

        To check that this is indeed the invariant form, see::

            sage: S = W.simple_reflections()                            # optional - gap3
            sage: all( F == S[i].matrix()*F*S[i].matrix().transpose() for i in W.index_set() )  # optional - gap3
            True

            sage: W = ReflectionGroup(['B',3])                          # optional - gap3
            sage: F = W.invariant_form(); F                             # optional - gap3
            [ 1 -1  0]
            [-1  2 -1]
            [ 0 -1  2]
            sage: w = W.an_element().to_matrix()                        # optional - gap3
            sage: w * F * w.transpose().conjugate() == F                # optional - gap3
            True

            sage: S = W.simple_reflections()                            # optional - gap3
            sage: all( F == S[i].matrix()*F*S[i].matrix().transpose() for i in W.index_set() )  # optional - gap3
            True

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: F = W.invariant_form(); F                             # optional - gap3
            [1 0]
            [0 1]

            sage: S = W.simple_reflections()                            # optional - gap3
            sage: all( F == S[i].matrix()*F*S[i].matrix().transpose().conjugate() for i in W.index_set() )  # optional - gap3
            True

        It also worked for badly generated groups::

            sage: W = ReflectionGroup(7)                                # optional - gap3
            sage: W.is_well_generated()                                 # optional - gap3
            False

            sage: F = W.invariant_form(); F                             # optional - gap3
            [1 0]
            [0 1]
            sage: S = W.simple_reflections()                            # optional - gap3
            sage: all( F == S[i].matrix()*F*S[i].matrix().transpose().conjugate() for i in W.index_set() )  # optional - gap3
            True

        And also for reducible types::

            sage: W = ReflectionGroup(['B',3],(4,2,3),4,7); W           # optional - gap3
            Reducible complex reflection group of rank 10 and type B3 x G(4,2,3) x ST4 x ST7
            sage: F = W.invariant_form(); S = W.simple_reflections()    # optional - gap3
            sage: all( F == S[i].matrix()*F*S[i].matrix().transpose().conjugate() for i in W.index_set() )  # optional - gap3
            True

        TESTS::

            sage: tests = [['A',3],['B',3],['F',4],(4,2,2),4,7]         # optional - gap3
            sage: for ty in tests:                                      # optional - gap3
            ....:     W = ReflectionGroup(ty)                           # optional - gap3
            ....:     A = W.invariant_form()                            # optional - gap3
            ....:     B = W.invariant_form(brute_force=True)            # optional - gap3
            ....:     print("{} {}".format(ty, A == B/B[0,0]))          # optional - gap3
            ['A', 3] True
            ['B', 3] True
            ['F', 4] True
            (4, 2, 2) True
            4 True
            7 True
        """
        if brute_force:
            return self._invariant_form_brute_force()

        n = self.rank()
        from sage.matrix.constructor import zero_matrix

        if self.is_crystallographic():
            ring = QQ
        else:
            from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
            ring = UniversalCyclotomicField()

        form = zero_matrix(ring, n, n)

        C = self.cartan_matrix()
        if not self.is_well_generated():
            indep_inds = sorted(self._index_set_inverse[key] for key in self.independent_roots().keys())
            C = C.matrix_from_rows_and_columns(indep_inds,indep_inds)

        for j in range(n):
            for i in range(j):
                if C[j,i] != 0:
                    form[j,j] = form[i,i] * \
                                ( C[i,j] * C[j,j].conjugate() ) / \
                                ( C[j,i].conjugate() * C[i,i] )
            if form[j,j] == 0:
                form[j,j] = ring.one()
        for j in range(n):
            for i in range(j):
                form[j, i] = C[i, j] * form[i, i] / C[i,i]
                form[i, j] = form[j, i].conjugate()

        B = self.base_change_matrix()
        form = B*form*B.conjugate().transpose()
        form /= form[0,0]
        form.set_immutable()
        return form

    def _invariant_form_brute_force(self):
        r"""
        Return the form that is invariant under the action of ``self``.

        This brute force algorithm is only kept for possible testing.

        EXAMPLES::

            sage: W = ReflectionGroup((3,1,2))                          # optional - gap3
            sage: W._invariant_form_brute_force()                       # optional - gap3
            [1 0]
            [0 1]
        """
        from sage.misc.cachefunc import cached_function

        Phi = self.roots()

        base_change = self.base_change_matrix()
        Delta = tuple(self.independent_roots())
        basis_is_Delta = base_change.is_one()
        if not basis_is_Delta:
            Delta = [beta * base_change for beta in Delta]

        S = self.simple_reflections()
        n = self.rank()

        def action_on_root(w, beta):
            if basis_is_Delta:
                return w.action_on_root(beta)
            else:
                return beta * w.to_matrix()

        @cached_function
        def invariant_value(i,j):
            if i > j:
                return invariant_value(j,i).conjugate()
            val = sum(action_on_root(w, Delta[i]) * action_on_root(w, Delta[j]).conjugate()
                      for w in self)
            if val in QQ:
                val = QQ(val)
            return val

        coeffs = []
        for i in self.index_set():
            coeff = 1-E(S[i].order())
            if coeff in QQ:
                coeff = QQ(coeff)
            coeffs.append(coeff)

        return Matrix([ [ invariant_value(i,j)/self.cardinality() for j in range(n) ] for i in range(n) ])

    def set_reflection_representation(self,refl_repr=None):
        r"""
        Set the reflection representation of ``self``.

        INPUT:

        - ``refl_repr`` -- a dictionary representing the matrices of the
          generators of ``self`` with keys given by the index set, or
          ``None`` to reset to the default reflection representation

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3
            sage: for w in W: w.to_matrix(); print("-----")             # optional - gap3
            [1 0]
            [0 1]
            -----
            [ 1  1]
            [ 0 -1]
            -----
            [-1  0]
            [ 1  1]
            -----
            [-1 -1]
            [ 1  0]
            -----
            [ 0  1]
            [-1 -1]
            -----
            [ 0 -1]
            [-1  0]
            -----

            sage: W.set_reflection_representation({1: matrix([[0,1,0],[1,0,0],[0,0,1]]), 2: matrix([[1,0,0],[0,0,1],[0,1,0]])}) # optional - gap3
            sage: for w in W: w.to_matrix(); print("-----")             # optional - gap3
            [1 0 0]
            [0 1 0]
            [0 0 1]
            -----
            [1 0 0]
            [0 0 1]
            [0 1 0]
            -----
            [0 1 0]
            [1 0 0]
            [0 0 1]
            -----
            [0 0 1]
            [1 0 0]
            [0 1 0]
            -----
            [0 1 0]
            [0 0 1]
            [1 0 0]
            -----
            [0 0 1]
            [0 1 0]
            [1 0 0]
            -----
            sage: W.set_reflection_representation()                     # optional - gap3
        """
        if refl_repr is None or set(refl_repr.keys()) == set(self.index_set()):
            self._reflection_representation = refl_repr
        else:
            raise ValueError("the reflection representation must be defined for the complete index set")

    class Element(PermutationGroupElement):

        #@cached_in_parent_method
        def conjugacy_class_representative(self):
            r"""
            Return a representative of the conjugacy class of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     print('%s %s'%(w.reduced_word(), w.conjugacy_class_representative().reduced_word()))  # optional - gap3
                [] []
                [2] [1]
                [1] [1]
                [1, 2] [1, 2]
                [2, 1] [1, 2]
                [1, 2, 1] [1]
            """
            W = self.parent()
            for w in W._conjugacy_classes.keys():
                if self in W._conjugacy_classes[w]:
                    return w
            return W.conjugacy_classes_representatives()[ gap3("PositionClass(%s,%s)"%(W._gap_group._name,self)).sage()-1 ]

        def conjugacy_class(self):
            r"""
            Return the conjugacy class of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
                sage: for w in W: sorted(w.conjugacy_class())           # optional - gap3
                [()]
                [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6), (1,5)(2,4)(3,6)]
                [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6), (1,5)(2,4)(3,6)]
                [(1,2,6)(3,4,5), (1,6,2)(3,5,4)]
                [(1,2,6)(3,4,5), (1,6,2)(3,5,4)]
                [(1,3)(2,5)(4,6), (1,4)(2,3)(5,6), (1,5)(2,4)(3,6)]
            """
            W = self.parent()
            if self not in W.conjugacy_classes_representatives():
                self = self.conjugacy_class_representative()
            if self in W._conjugacy_classes.keys():
                return W._conjugacy_classes[self]
            gens = W.simple_reflections()
            count = 0
            orbit = [self]
            orbit_set = set(orbit)
            while count < len(orbit):
                w = orbit[count]
                count += 1
                for s in gens:
                    w_new = s*w*s**-1
                    if w_new not in orbit_set:
                        orbit.append(w_new)
                        orbit_set.add(w_new)
            orbit_set = frozenset(orbit_set)
            W._conjugacy_classes[self] = orbit_set
            return orbit_set

        def reduced_word(self):
            r"""
            Return a word in the simple reflections to obtain ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((5,1,1), index_set=['a'], hyperplane_index_set=['x'], reflection_index_set=['A','B','C','D']) # optional - gap3
                sage: [w.reduced_word() for w in W]                     # optional - gap3
                [[], ['a'], ['a', 'a'], ['a', 'a', 'a'], ['a', 'a', 'a', 'a']]

            .. SEEALSO:: :meth:`reduced_word_in_reflections`
            """
            I = self.parent()._index_set
            return [I[i] for i in self._reduced_word]

        @lazy_attribute
        def _reduced_word(self):
            r"""
            Computes a reduced word and stores it into ``self._reduced_word``.

            TESTS::

                sage: W = ReflectionGroup((5,1,1))                      # optional - gap3
                sage: w = W.an_element()                                # optional - gap3
                sage: w._reduced_word                                   # optional - gap3
                [0]
            """
            W = self.parent()
            gens = [W.simple_reflection(j) for j in W._index_set]
            return _gap_factorization(self, gens)

        #@cached_in_parent_method
        def reduced_word_in_reflections(self):
            r"""
            Return a word in the reflections to obtain ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((5,1,1), index_set=['a'], reflection_index_set=['A','B','C','D']) # optional - gap3
                sage: [w.reduced_word_in_reflections() for w in W]      # optional - gap3
                [[], ['A'], ['B'], ['C'], ['D']]

            .. SEEALSO:: :meth:`reduced_word`
            """
            if self.is_one():
                return []

            W = self.parent()
            gens = [W.reflection(j) for j in W._reflection_index_set]
            word = _gap_factorization(self, gens)
            return [self.parent()._reflection_index_set[i] for i in word]

        def length(self):
            r"""
            Return the length of ``self`` in generating reflections.

            This is the minimal numbers of generating reflections needed
            to obtain ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(4)                            # optional - gap3
                sage: for w in W:         # optional - gap3
                ....:     print("{} {}".format(w.reduced_word(), w.length()))
                [] 0
                [1] 1
                [2] 1
                [1, 1] 2
                [1, 2] 2
                [2, 1] 2
                [2, 2] 2
                [1, 1, 2] 3
                [1, 2, 1] 3
                [1, 2, 2] 3
                [2, 1, 1] 3
                [2, 2, 1] 3
                [1, 1, 2, 1] 4
                [1, 1, 2, 2] 4
                [1, 2, 1, 1] 4
                [1, 2, 2, 1] 4
                [2, 1, 1, 2] 4
                [2, 2, 1, 1] 4
                [1, 1, 2, 1, 1] 5
                [1, 1, 2, 2, 1] 5
                [1, 2, 1, 1, 2] 5
                [1, 2, 2, 1, 1] 5
                [1, 1, 2, 1, 1, 2] 6
                [1, 1, 2, 2, 1, 1] 6
            """
            return len(self.reduced_word())

        #@cached_in_parent_method
        def reflection_length(self, in_unitary_group=False):
            r"""
            Return the reflection length of ``self``.

            This is the minimal numbers of reflections needed to obtain
            ``self``.

            INPUT:

            - ``in_unitary_group`` -- (default: ``False``) if ``True``,
              the reflection length is computed in the unitary group
              which is the dimension of the move space of ``self``

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
                sage: sorted([t.reflection_length() for t in W])        # optional - gap3
                [0, 1, 1, 1, 2, 2]

                sage: W = ReflectionGroup((2,1,2))                      # optional - gap3
                sage: sorted([t.reflection_length() for t in W])        # optional - gap3
                [0, 1, 1, 1, 1, 2, 2, 2]

                sage: W = ReflectionGroup((2,2,2))                      # optional - gap3
                sage: sorted([t.reflection_length() for t in W])        # optional - gap3
                [0, 1, 1, 2]

                sage: W = ReflectionGroup((3,1,2))                      # optional - gap3
                sage: sorted([t.reflection_length() for t in W])        # optional - gap3
                [0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            """
            W = self.parent()
            if self in W.conjugacy_classes_representatives():
                if in_unitary_group or W.is_real():
                    return W.rank() - self.reflection_eigenvalues(is_class_representative=True).count(0)
                else:
                    return len(self.reduced_word_in_reflections())
            else:
                w = self.conjugacy_class_representative()
                # the following assert a possible implementation bug and
                # is hopefully never needed
                assert w in self.parent().conjugacy_classes_representatives()
                return w.reflection_length(in_unitary_group=in_unitary_group)

        #@cached_in_parent_method
        def to_matrix(self, on_space="primal"):
            r"""
            Return ``self`` as a matrix acting on the underlying vector
            space.

            - ``on_space`` -- optional (default: ``"primal"``) whether
              to act as the reflection representation on the given
              basis, or to act on the dual reflection representation
              on the dual basis

            EXAMPLES::

                sage: W = ReflectionGroup((3,1,2))          # optional - gap3
                sage: for w in W:                           # optional - gap3
                ....:     w.reduced_word()                  # optional - gap3
                ....:     [w.to_matrix(), w.to_matrix(on_space="dual")] # optional - gap3
                []
                [
                [1 0]  [1 0]
                [0 1], [0 1]
                ]
                [1]
                [
                [E(3)    0]  [E(3)^2      0]
                [   0    1], [     0      1]
                ]
                [2]
                [
                [0 1]  [0 1]
                [1 0], [1 0]
                ]
                [1, 1]
                [
                [E(3)^2      0]  [E(3)    0]
                [     0      1], [   0    1]
                ]
                [1, 2]
                [
                [   0 E(3)]  [     0 E(3)^2]
                [   1    0], [     1      0]
                ]
                [2, 1]
                [
                [   0    1]  [     0      1]
                [E(3)    0], [E(3)^2      0]
                ]
                [1, 1, 2]
                [
                [     0 E(3)^2]  [   0 E(3)]
                [     1      0], [   1    0]
                ]
                [1, 2, 1]
                [
                [   0 E(3)]  [     0 E(3)^2]
                [E(3)    0], [E(3)^2      0]
                ]
                [2, 1, 1]
                [
                [     0      1]  [   0    1]
                [E(3)^2      0], [E(3)    0]
                ]
                [2, 1, 2]
                [
                [   1    0]  [     1      0]
                [   0 E(3)], [     0 E(3)^2]
                ]
                [1, 1, 2, 1]
                [
                [     0 E(3)^2]  [     0   E(3)]
                [  E(3)      0], [E(3)^2      0]
                ]
                [1, 2, 1, 1]
                [
                [     0   E(3)]  [     0 E(3)^2]
                [E(3)^2      0], [  E(3)      0]
                ]
                [1, 2, 1, 2]
                [
                [E(3)    0]  [E(3)^2      0]
                [   0 E(3)], [     0 E(3)^2]
                ]
                [2, 1, 1, 2]
                [
                [     1      0]  [   1    0]
                [     0 E(3)^2], [   0 E(3)]
                ]
                [1, 1, 2, 1, 1]
                [
                [     0 E(3)^2]  [   0 E(3)]
                [E(3)^2      0], [E(3)    0]
                ]
                [1, 1, 2, 1, 2]
                [
                [E(3)^2      0]  [  E(3)      0]
                [     0   E(3)], [     0 E(3)^2]
                ]
                [1, 2, 1, 1, 2]
                [
                [  E(3)      0]  [E(3)^2      0]
                [     0 E(3)^2], [     0   E(3)]
                ]
                [1, 1, 2, 1, 1, 2]
                [
                [E(3)^2      0]  [E(3)    0]
                [     0 E(3)^2], [   0 E(3)]
                ]
            """
            W = self.parent()
            if W._reflection_representation is None:
                Delta = W.independent_roots()
                Phi = W.roots()
                inds = [ W._index_set_inverse[i] for i in W.independent_roots().keys() ]
                M = Matrix([Phi[self(i+1)-1] for i in inds])
                mat = W.base_change_matrix() * M
            else:
                refl_repr = W._reflection_representation
                id_mat = identity_matrix(QQ, refl_repr[W.index_set()[0]].nrows())
                mat = prod([refl_repr[i] for i in self.reduced_word()], id_mat)

            if on_space == "primal":
                pass
            elif on_space == "dual":
                mat = mat.inverse().transpose()
            else:
                raise ValueError('on_space must be "primal" or "dual"')

            mat.set_immutable()
            return mat

        matrix = to_matrix

        def action(self, vec, on_space="primal"):
            r"""
            Return the image of ``vec`` under the action of ``self``.

            INPUT:

            - ``vec`` -- vector in the basis given by the simple root

            - ``on_space`` -- optional (default: ``"primal"``) whether
              to act as the reflection representation on the given
              basis, or to act on the dual reflection representation
              on the dual basis

            EXAMPLES::

                sage: W = ReflectionGroup((3,1,2))                      # optional - gap3
                sage: w = W.from_reduced_word([1, 2, 1, 1, 2])          # optional - gap3
                sage: for alpha in W.independent_roots():               # optional - gap3
                ....:     print("%s -> %s"%(alpha,w.action(alpha)))     # optional - gap3
                (1, 0) -> (E(3), 0)
                (-1, 1) -> (-E(3), E(3)^2)
            """
            mat = self.matrix(on_space=on_space)
            return vec * mat

        def _act_on_(self, vec, self_on_left):
            r"""
            Defines the action of ``self`` as a linear transformation
            on the vector space, in the basis given by the simple
            roots.

            - ``vec`` -- the vector (an iterable) to act on

            - ``self_on_left`` -- whether the action of ``self`` is on
              the left or on the right

            .. WARNING::

                This coercion is currently only defined for
                ``self_on_left`` being ``False``.

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])          # optional - gap3
                sage: w = W.from_reduced_word([1,2])        # optional - gap3
                sage: for root in W.positive_roots():       # optional - gap3
                ....:     print("%s -> %s"%(root, w*root))  # optional - gap3
                (1, 0) -> (0, 1)
                (0, 1) -> (-1, -1)
                (1, 1) -> (-1, 0)

                sage: for root in W.positive_roots():       # optional - gap3
                ....:     print("%s -> %s"%(root, root*w))  # optional - gap3
                (1, 0) -> (-1, -1)
                (0, 1) -> (1, 0)
                (1, 1) -> (0, -1)
            """
            if not self_on_left:
                return self.action(vec, side="right")

        def action_on_root_indices(self, i):
            """
            Return the right action on the set of roots.

            INPUT:

            - ``i`` -- index of the root to act on

            EXAMPLES::

                sage: W = ReflectionGroup(['A',3])           # optional - gap3
                sage: w = W.w0                               # optional - gap3
                sage: N = len(W.roots())                     # optional - gap3
                sage: [w.action_on_root_indices(i) for i in range(N)]    # optional - gap3
                [8, 7, 6, 10, 9, 11, 2, 1, 0, 4, 3, 5]

                sage: W = ReflectionGroup(['A',2], reflection_index_set=['A','B','C'])   # optional - gap3
                sage: w = W.w0                               # optional - gap3
                sage: N = len(W.roots())                     # optional - gap3
                sage: [w.action_on_root_indices(i) for i in range(N)]    # optional - gap3
                [4, 3, 5, 1, 0, 2]

            TESTS::

                sage: W = ReflectionGroup(4)                 # optional - gap3
                sage: N = len(W.roots())                     # optional - gap3
                sage: all(sorted([w.action_on_root_indices(i) for i in range(N)]) == range(N) for w in W)   # optional - gap3
                True
            """
            return self(i + 1) - 1

        def action_on_root(self, root):
            r"""
            Return the root obtained by applying ``self`` to ``root``
            on the right.

            INPUT:

            - ``root`` -- the root to act on

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     print("%s %s"%(w.reduced_word(),              # optional - gap3
                ....:           [w.action_on_root(beta,side="left") for beta in W.positive_roots()]))  # optional - gap3
                [] [(1, 0), (0, 1), (1, 1)]
                [2] [(1, 1), (0, -1), (1, 0)]
                [1] [(-1, 0), (1, 1), (0, 1)]
                [1, 2] [(0, 1), (-1, -1), (-1, 0)]
                [2, 1] [(-1, -1), (1, 0), (0, -1)]
                [1, 2, 1] [(0, -1), (-1, 0), (-1, -1)]

            TESTS::

                sage: W = ReflectionGroup(4); Phi = sorted(W.roots())   # optional - gap3
                sage: all(sorted([w.action_on_root(beta) for beta in Phi]) == Phi for w in W)   # optional - gap3
                True
            """
            Phi = self.parent().roots()
            return Phi[self.action_on_root_indices(Phi.index(root))]

        def to_permutation_of_roots(self):
            r"""
            Return ``self`` as a permutation of the roots with indexing
            starting at `1`.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
                sage: for w in W:       # optional - gap3
                ....:     perm = w.to_permutation_of_roots()
                ....:     print("{} {}".format(perm, perm == w))
                () True
                (1,3)(2,5)(4,6) True
                (1,4)(2,3)(5,6) True
                (1,6,2)(3,5,4) True
                (1,2,6)(3,4,5) True
                (1,5)(2,4)(3,6) True
            """
            return PermutationGroupElement(self)

        #@cached_in_parent_method
        def fix_space(self):
            r"""
            Return the fix space of ``self``.

            This is the sub vector space of the underlying vector space
            on which ``self`` acts trivially.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     w.reduced_word()                              # optional - gap3
                ....:     w.fix_space()                                 # optional - gap3
                []
                Vector space of degree 2 and dimension 2 over Rational Field
                Basis matrix:
                [1 0]
                [0 1]
                [2]
                Vector space of degree 2 and dimension 1 over Rational Field
                Basis matrix:
                [1 0]
                [1]
                Vector space of degree 2 and dimension 1 over Rational Field
                Basis matrix:
                [0 1]
                [1, 2]
                Vector space of degree 2 and dimension 0 over Rational Field
                Basis matrix:
                []
                [2, 1]
                Vector space of degree 2 and dimension 0 over Rational Field
                Basis matrix:
                []
                [1, 2, 1]
                Vector space of degree 2 and dimension 1 over Rational Field
                Basis matrix:
                [ 1 -1]

                sage: W = ReflectionGroup(23)                           # optional - gap3
                sage: W.an_element().fix_space()                        # optional - gap3
                Vector space of degree 3 and dimension 2 over Universal Cyclotomic Field
                Basis matrix:
                [0 1 0]
                [0 0 1]
            """
            I = identity_matrix(QQ, self.parent().rank())
            return (self.to_matrix() - I).right_kernel()

        #@cached_in_parent_method
        def reflection_eigenvalues(self, is_class_representative=False):
            r"""
            Return the reflection eigenvalues of ``self``.

            INPUT:

            - ``is_class_representative`` -- (default: ``False``) whether
              to first replace ``self`` by the representative of its
              conjugacy class

            EXAMPLES::

                sage: W = ReflectionGroup(4)                            # optional - gap3
                sage: for w in W: w.reflection_eigenvalues()            # optional - gap3
                [0, 0]
                [1/3, 0]
                [1/3, 0]
                [2/3, 0]
                [1/6, 1/2]
                [1/6, 1/2]
                [2/3, 0]
                [1/4, 3/4]
                [1/4, 3/4]
                [1/4, 3/4]
                [1/4, 3/4]
                [1/4, 3/4]
                [1/3, 0]
                [1/2, 5/6]
                [1/3, 0]
                [1/2, 5/6]
                [1/2, 5/6]
                [1/2, 5/6]
                [1/6, 1/2]
                [2/3, 0]
                [1/6, 1/2]
                [2/3, 0]
                [1/2, 1/2]
                [1/4, 3/4]
            """
            return self.parent().reflection_eigenvalues(self, is_class_representative=is_class_representative)

        #@cached_in_parent_method
        def galois_conjugates(self):
            r"""
            Return all Galois conjugates of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(4)                            # optional - gap3
                sage: for w in W: print(w.galois_conjugates())          # optional - gap3
                [[1 0]
                 [0 1]]
                [[   1    0]
                 [   0 E(3)], [     1      0]
                 [     0 E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2],
                 [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]]
                [[     1      0]
                 [     0 E(3)^2], [   1    0]
                 [   0 E(3)]]
                [[ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2],
                 [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2],
                 [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2],
                 [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2],
                 [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2],
                 [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2],
                 [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2],
                 [-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2],
                 [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2],
                 [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2],
                 [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2],
                 [-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2],
                 [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]]
                [[   -1     0]
                 [    0 -E(3)], [     -1       0]
                 [      0 -E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2],
                 [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [ 4/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]]
                [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2],
                 [-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]]
                [[-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2],
                 [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]]
                [[     -1       0]
                 [      0 -E(3)^2], [   -1     0]
                 [    0 -E(3)]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [-4/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2],
                 [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]]
                [[-1  0]
                 [ 0 -1]]
                [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
                 [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2],
                 [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
                 [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]]
            """
            rk = self.parent().rank()
            M = self.to_matrix().list()
            m = lcm([x.conductor() if hasattr(x,"conductor") else 1 for x in M])
            M_gals = [x.galois_conjugates(m) if hasattr(x,"galois_conjugates") else [x] for x in M]
            conjugates = []
            for i in range(len(M_gals[0])):
                conjugates.append(Matrix(rk, [X[i] for X in M_gals]))
            return conjugates

        def __cmp__(self, other):
            r"""
            Compare ``self`` with ``other``.

            Without this comparison method, the initialization of this
            permutation group fails :-(.

            TESTS::

                sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
                sage: a,b = W.gens()                                    # optional - gap3
                sage: a.__cmp__(b)                                      # optional - gap3
                1
                sage: b.__cmp__(a)                                      # optional - gap3
                -1
            """
            return super(ComplexReflectionGroup.Element, self).__cmp__(other)

class IrreducibleComplexReflectionGroup(ComplexReflectionGroup):

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3)); W                       # optional - gap3
            Irreducible real reflection group of rank 2 and type A2
            sage: W = ReflectionGroup((3,1,4)); W                       # optional - gap3
            Irreducible complex reflection group of rank 4 and type G(3,1,4)
        """
        type_str = self._irrcomp_repr_(self._type[0])
        return 'Irreducible complex reflection group of rank %s and type %s'%(self._rank,type_str)

    def elements_below_coxeter_element(self, c=None):
        r"""
        Return all elements in ``self`` in the interval `[1,c]` in the
        absolute order of ``self``.

        This order is defined by

        .. MATH::

            \omega \leq_R \tau \Leftrightarrow \ell_R(\omega) +
            \ell_R(\omega^{-1} \tau) = \ell_R(\tau),

        where `\ell_R` denotes the reflection length.

        .. NOTE::

            ``self`` is assumed to be well-generated.

        INPUT:

        - ``c`` -- (default: ``None``) if an element ``c`` is given, it
          is used as the maximal element in the interval; if a list is
          given, the union of the various maximal elements is computed

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3

            sage: sorted( w.reduced_word() for w in W.elements_below_coxeter_element() )    # optional - gap3
            [[], [1], [1, 2], [1, 2, 1], [2]]

            sage: sorted( w.reduced_word() for w in W.elements_below_coxeter_element(W.from_reduced_word([2,1])) )  # optional - gap3
            [[], [1], [1, 2, 1], [2], [2, 1]]

            sage: sorted( w.reduced_word() for w in W.elements_below_coxeter_element(W.from_reduced_word([2])) )    # optional - gap3
            [[], [2]]
        """
        if c in self:
            cs = [c]
        elif c is None:
            cs = [self.coxeter_element()]
        else:
            cs = list(c)
        l = cs[0].reflection_length(in_unitary_group=True)
        f = lambda pi: any(pi.reflection_length(in_unitary_group=True)
                           + (c*pi**-1).reflection_length(in_unitary_group=True) == l
                           for c in cs)
        # first computing the conjugacy classes only needed if the interaction with gap3 is slow due to a bug
        #self.conjugacy_classes()
        return filter(f, self)

    # TODO: lift to ComplexReflectionGroups.Finite
    #       this method can be defined for well-generated, finite,
    #       irreducible complex reflection group. The current
    #       implementation uses this particular connection to chevie.
    # TODO: have a cached and an uncached version
    @cached_method
    def noncrossing_partition_lattice(self, c=None, L=None, in_unitary_group=False):
        r"""
        Return the interval `[1,c]` in the absolute order of
        ``self`` as a finite lattice.

        .. SEEALSO:: :meth:`elements_below_coxeter_element`

        INPUT:

        - ``c`` -- (default: ``None``) if an element ``c`` in ``self`` is
          given, it is used as the maximal element in the interval

        - ``L`` -- (default: ``None``) if a subset ``L`` (must be hashable!)
          of ``self`` is given, it is used as the underlying set (only
          cover relations are checked)

        - ``in_unitary_group`` -- (default: ``False``) if ``False``, the
          relation is given by `\sigma \leq \tau` if
          `l_R(\sigma) + l_R(\sigma^{-1}\tau) = l_R(\tau)`;
          if ``True``, the relation is given by `\sigma \leq \tau` if
          `\dim(\mathrm{Fix}(\sigma)) + \dim(\mathrm{Fix}(\sigma^{-1}\tau))
          = \dim(\mathrm{Fix}(\tau))`

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3

            sage: sorted( w.reduced_word() for w in W.noncrossing_partition_lattice() ) # optional - gap3
            [[], [1], [1, 2], [1, 2, 1], [2]]

            sage: sorted( w.reduced_word() for w in W.noncrossing_partition_lattice(W.from_reduced_word([2,1])) )   # optional - gap3
            [[], [1], [1, 2, 1], [2], [2, 1]]

            sage: sorted( w.reduced_word() for w in W.noncrossing_partition_lattice(W.from_reduced_word([2])) ) # optional - gap3
            [[], [2]]
        """
        from sage.combinat.posets.all import Poset, LatticePoset
        if c is None:
            c = self.coxeter_element()

        smart_covers = not in_unitary_group

        if self.is_real():
            smart_covers = in_unitary_group = True

        R = self.reflections()
        if L is None:
            L = self.elements_below_coxeter_element(c=c)
            if c.is_coxeter_element():
                smart_covers = in_unitary_group = True
        rels = []
        ref_lens = {w: w.reflection_length(in_unitary_group=in_unitary_group)
                    for w in L}
        if smart_covers:
            for pi in L:
                for t in R:
                    tau = pi*t
                    if tau in L and ref_lens[pi] + 1 == ref_lens[tau]:
                        rels.append((pi,tau))
        else:
            rels = [(pi,tau) for pi in L for tau in L
                    if ref_lens[pi] + ref_lens[pi.inverse()*tau] == ref_lens[tau]]
        P = Poset((L,rels), cover_relations=smart_covers, facade=True)
        if P.is_lattice():
            return LatticePoset(P)
        else:
            return P

    # TODO: lift to ComplexReflectionGroups.Finite
    #       this method can be defined for well-generated, finite,
    #       irreducible complex reflection group. The current
    #       implementation uses this particular connection to chevie.
    def generalized_noncrossing_partitions(self, m, c=None, positive=False):
        r"""
        Return the set of all chains of length ``m`` in the noncrossing
        partition lattice of ``self``, see
        :meth:`noncrossing_partition_lattice`.

        .. NOTE::

            ``self`` is assumed to be well-generated

        INPUT:

        - ``c`` -- (default: ``None``) if an element ``c`` in ``self``
          is given, it is used as the maximal element in the interval

        - ``positive`` -- (default: ``False``) if ``True``, only those
          generalized noncrossing partitions of full support are returned

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                          # optional - gap3

            sage: sorted([w.reduced_word() for w in chain]              # optional - gap3
            ....:        for chain in W.generalized_noncrossing_partitions(2))  # optional - gap3
            [[[], [], [1, 2]],
             [[], [1], [2]],
             [[], [1, 2], []],
             [[], [1, 2, 1], [1]],
             [[], [2], [1, 2, 1]],
             [[1], [], [2]],
             [[1], [2], []],
             [[1, 2], [], []],
             [[1, 2, 1], [], [1]],
             [[1, 2, 1], [1], []],
             [[2], [], [1, 2, 1]],
             [[2], [1, 2, 1], []]]

            sage: sorted([w.reduced_word() for w in chain]              # optional - gap3
            ....:        for chain in W.generalized_noncrossing_partitions(2, positive=True))   # optional - gap3
            [[[], [1, 2], []],
             [[], [1, 2, 1], [1]],
             [[1], [2], []],
             [[1, 2], [], []],
             [[1, 2, 1], [], [1]],
             [[1, 2, 1], [1], []],
             [[2], [1, 2, 1], []]]
        """
        from sage.combinat.combination import Combinations
        NC = self.noncrossing_partition_lattice(c=c)
        one = self.one()
        if c is None:
            c = self.coxeter_element()
        chains = NC.chains()
        NCm = set()
        iter = chains.breadth_first_search_iterator()
        chain = next(iter)
        chain = next(iter)
        while len(chain) <= m:
            chain.append( c )
            for i in range(len(chain)-1, 0, -1):
                chain[i] = chain[i-1]**-1 * chain[i]
            k = m + 1 - len(chain)
            for positions in Combinations(range(m+1),k):
                ncm = []
                for l in range(m+1):
                    if l in positions:
                        ncm.append(one)
                    else:
                        l_prime = l - len([i for i in positions if i <= l])
                        ncm.append(chain[l_prime])
                if not positive or prod(ncm[:-1]).has_full_support():
                    NCm.add(tuple(ncm))
            try:
                chain = next(iter)
            except StopIteration:
                chain = range(m + 1)
        return NCm

    # TODO: lift to ComplexReflectionGroups.Finite
    #       this method can be defined for well-generated, finite,
    #       irreducible complex reflection group. The current
    #       implementation uses this particular connection to chevie.
    # TODO: have a cached and an uncached version
    def absolute_poset(self, in_unitary_group=False):
        r"""
        Return the poset induced by the absolute order of ``self`` as a
        finite lattice.

        INPUT:

        - ``in_unitary_group`` -- (default: ``False``) if ``False``, the
          relation is given by ``\sigma \leq \tau`` if
          `l_R(\sigma) + l_R(\sigma^{-1}\tau) = l_R(\tau)`
          If ``True``, the relation is given by `\sigma \leq \tau` if
          `\dim(\mathrm{Fix}(\sigma)) + \dim(\mathrm{Fix}(\sigma^{-1}\tau))
          = \dim(\mathrm{Fix}(\tau))`.

        .. SEEALSO:: :meth:`noncrossing_partition_lattice`

        EXAMPLES::

            sage: P = ReflectionGroup((1,1,3)).absolute_poset(); P      # optional - gap3
            Finite poset containing 6 elements

            sage: sorted(w.reduced_word() for w in P)                   # optional - gap3
            [[], [1], [1, 2], [1, 2, 1], [2], [2, 1]]

            sage: W = ReflectionGroup(4); W                             # optional - gap3
            Irreducible complex reflection group of rank 2 and type ST4
            sage: W.absolute_poset()                                    # optional - gap3
            Finite poset containing 24 elements
        """
        return self.noncrossing_partition_lattice(L=self)

    class Element(ComplexReflectionGroup.Element):

        # TODO: lift to ComplexReflectionGroups.Finite
        #       this method can be defined for well-generated, finite,
        #       irreducible complex reflection group. The current
        #       implementation uses this particular connection to chevie.
        #@cached_in_parent_method
        def is_coxeter_element(self, which_primitive=1, is_class_representative=False):
            r"""
            Return ``True`` if ``self`` is a Coxeter element.

            INPUT:

            - ``which_primitive`` -- (default:``1``) for which power of
              the first primitive ``h``-th root of unity to look as a
              reflection eigenvalue for a regular element

            - ``is_class_representative`` -- boolean (default ``True``) whether
              to compute instead on the conjugacy class representative

            .. SEEALSO::

                :meth:`~IrreducibleComplexReflectionGroup.coxeter_element`

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     print('%s %s'%(w.reduced_word(), w.is_coxeter_element())) # optional - gap3
                [] False
                [2] False
                [1] False
                [1, 2] True
                [2, 1] True
                [1, 2, 1] False
            """
            if not self.parent().is_irreducible() or not self.parent().is_well_generated():
                raise ValueError("this method is available for elements in irreducible, well-generated complex reflection groups")
            h = self.parent().coxeter_number()
            # to check regularity for a Coxeter number h, we get that an eigenvector is regular for free
            return any(QQ(ev).denom() == h and QQ(ev).numer() == which_primitive
                       for ev in self.reflection_eigenvalues(is_class_representative=is_class_representative))

        #@cached_in_parent_method
        def is_h_regular(self, is_class_representative=False):
            r"""
            Return whether ``self`` is regular.

            This is if ``self`` has an eigenvector with eigenvalue `h`
            and which does not lie in any reflection hyperplane.
            Here, `h` denotes the Coxeter number.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     print('%s %s'%(w.reduced_word(), w.is_h_regular()))   # optional - gap3
                [] False
                [2] False
                [1] False
                [1, 2] True
                [2, 1] True
                [1, 2, 1] False
            """
            if not self.parent().is_irreducible() or not self.parent().is_well_generated():
                raise ValueError("This method is available for elements in irreducible, well-generated complex reflection groups")
            h = self.parent().coxeter_number()
            # to check regularity for a Coxeter number h, we get that an eigenvector is regular for free
            return any(QQ(ev).denom() == h
                       for ev in self.reflection_eigenvalues(is_class_representative=is_class_representative))

        #@cached_in_parent_method
        def is_regular(self, h, is_class_representative=False):
            r"""
            Return whether ``self`` is regular.

            This is, if ``self`` has an eigenvector with eigenvalue
            ``h`` and which does not lie in any reflection hyperplane.

            - ``is_class_representative`` -- boolean (default ``True``) whether
              to compute instead on the conjugacy class representative

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3)); h = W.coxeter_number()  # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     print("{} {}".format(w.reduced_word(), w.is_regular(h)))
                [] False
                [2] False
                [1] False
                [1, 2] True
                [2, 1] True
                [1, 2, 1] False

                sage: W = ReflectionGroup(23); h = W.coxeter_number()   # optional - gap3
                sage: for w in W:                                       # optional - gap3
                ....:     if w.is_regular(h):                           # optional - gap3
                ....:         w.reduced_word()                          # optional - gap3
                [1, 2, 3]
                [2, 1, 3]
                [1, 3, 2]
                [3, 2, 1]
                [2, 1, 2, 3, 2]
                [2, 3, 2, 1, 2]
                [1, 2, 1, 2, 3, 2, 1]
                [1, 2, 3, 2, 1, 2, 1]
                [1, 2, 1, 2, 3, 2, 1, 2, 3]
                [2, 1, 2, 1, 3, 2, 1, 2, 3]
                [2, 1, 2, 3, 2, 1, 2, 1, 3]
                [1, 2, 3, 2, 1, 2, 1, 3, 2]
                [3, 2, 1, 2, 1, 3, 2, 1, 2]
                [1, 2, 1, 2, 1, 3, 2, 1, 2]
                [2, 3, 2, 1, 2, 1, 3, 2, 1]
                [2, 1, 2, 1, 3, 2, 1, 2, 1]
                [2, 3, 2, 1, 2, 1, 3, 2, 1, 2, 3]
                [1, 3, 2, 1, 2, 1, 3, 2, 1, 2, 3]
                [1, 2, 1, 2, 1, 3, 2, 1, 2, 1, 3]
                [1, 2, 1, 2, 3, 2, 1, 2, 1, 3, 2]
                [1, 2, 3, 2, 1, 2, 1, 3, 2, 1, 2]
                [2, 1, 2, 3, 2, 1, 2, 1, 3, 2, 1]
                [2, 1, 2, 3, 2, 1, 2, 1, 3, 2, 1, 2, 3]
                [1, 2, 1, 3, 2, 1, 2, 1, 3, 2, 1, 2, 3]
            """
            evs = self.reflection_eigenvalues(is_class_representative=is_class_representative)
            I = identity_matrix(self.parent().rank())
            for ev in evs:
                ev = QQ(ev)
                if h == ev.denom():
                    M = self.to_matrix() - E(ev.denom(),ev.numer()) * I
                    V = M.right_kernel()
                    if all(not V.is_subspace(H) for H in self.parent().reflection_hyperplanes()):
                        return True
            return False

def _gap_factorization(w, gens):
    r"""
    Return a factorization of ``w`` using the generators ``gens``.

    .. WARNING::

        This is only available through GAP3 and Chevie.

    EXAMPLES::

        sage: from sage.combinat.root_system.reflection_group_complex import _gap_factorization
        sage: W = ReflectionGroup((1,1,3))                              # optional - gap3
        sage: gens = [ W.simple_reflection(i) for i in W.index_set() ]  # optional - gap3
        sage: for w in W: _gap_factorization(w,gens)                    # optional - gap3
        []
        [1]
        [0]
        [0, 1]
        [1, 0]
        [0, 1, 0]
    """
    gap3.execute('W := GroupWithGenerators(%s)'%str(gens))
    gap3.execute(_gap_factorization_code)
    fac = gap3('MinimalWord(W,%s)'%str(w)).sage()
    return [i-1 for i in fac]

_gap_factorization_code = """
# MinimalWord(G,w)
# given a permutation group G find some expression of minimal length in the
# generators of G and their inverses of the element w (an inverse is
# representated by a negative index).
# To speed up  later calls to  the same function  the fields G.base, G.words,
# G.nbwordslength are kept.
MinimalWord:=function(G,w)
  local decode,i,p,g,h,n,bag,nbe,nbf,new,gens,inds;
  # to save space elements of G are represented as image of the base, and
  # words are represented as: index of previous elt, last generator applied;
  if not IsBound(G.base) then
    StabChain(G);g:=G; G.base:=[];
    while IsBound(g.orbit) do Add(G.base,g.orbit[1]); g:=g.stabilizer; od;
  fi;
  w:=OnTuples(G.base,w);
  if not IsBound(G.words) then
    G.words:=[G.base]; G.lastmult:=[[0,0]];
    G.nbwordslength:=[1];
  fi;
  gens:=ShallowCopy(G.generators);inds:=[1..Length(gens)];
  #  for g in G.generators do
  #    if g<>g^-1 then Add(gens,g^-1);Add(inds,-Position(gens,g));fi;
  #  od;
  bag:=Set(G.words);
  nbe:=0;nbf:=0;
  decode:=function(i)local w;w:=[];
    while i<>1 do Add(w,G.lastmult[i][2]); i:=G.lastmult[i][1];od;
    return Reversed(w);
  end;
  while true do
    if w in bag then return decode(Position(G.words,w));fi;
    new:=Length(G.words);
    for g in [1..Length(gens)] do
      for h in [1+Sum(G.nbwordslength{[1..Length(G.nbwordslength)-1]})..new] do
         n:=OnTuples(G.words[h],gens[g]);
         if n in bag then
           nbe:=nbe+1;# if nbe mod 500=1 then Print(".\c");fi;
         else
           nbf:=nbf+1;# if nbf mod 500=1 then Print("*\c");fi;
       Add(G.words,n);Add(G.lastmult,[h,inds[g]]);AddSet(bag,n);
         fi;
       od;
    od;
    Add(G.nbwordslength,Length(G.words)-new);
    Print("\n",G.nbwordslength[Length(G.nbwordslength)]," elements of length ",
      Length(G.nbwordslength)-1);
  od;
end;"""

def _gap_return(S, coerce_obj='self'):
    r"""
    Return the string ``S`` after a few modifications are done.

    This is a stupid internal function to take GAP output as a string,
    replace a few things, to then turn it into a Sage object.

    TESTS::

        sage: from sage.combinat.root_system.reflection_group_complex import _gap_return
        sage: _gap_return("[ (), (1,4)(2,3)(5,6), (1,6,2)(3,5,4) ]")    # optional - gap3
        "[self('()',check=False),self('(1,4)(2,3)(5,6)',check=False),self('(1,6,2)(3,5,4)',check=False)]"
    """
    S = S.replace(' ','').replace('\n','')
    S = S.replace(',(','\',check=False),%s(\'('%coerce_obj).replace('[','[%s(\''%coerce_obj).replace(']','\',check=False)]')
    return S

