r"""
Implementation of Flag, elements of :class:`CombinatorialTheory`

Cython class for flags and types. Types will be called ftype 
(short for flag type, to distinguish from the type keyword in python)
They are elements of :class:`CombinatorialTheory`. They also behave as
basis elements for :class:`FlagAlgebra`, hence basic operations
like addition and multiplication are defined.


The class :class:`CombinatorialTheory` acts as a parent for flag elements.
Its members are theories, they come with signature, ways to generate 
elements, ways to check if they are equal. This example uses the 
GraphTheory object. To find out more about combinatorial theories,
other pre-implemented members, and how to solve problems inside them
using flag algebras, see the documentation of 
:mod:`sage.algebras.flag_algebras`. This file is about flags mainly.

The examples will use GraphTheory

    sage: from sage.algebras.flag_algebras import GraphTheory

To create flags from a theory we can call for example ::
    
    sage: g = GraphTheory(3, edges=[[0, 1]])

This creates a graph on `3` vertices. They are called `[0, 1, 2]` for 
simplicity, and it defines the single edge `[0, 1]`. The result is a Flag
`g`. The ftype is a collection of marked vertices (a constant in the theory). 
For example, we can define the same graph but mark one vertex ::
    
    sage: g0 = GraphTheory(3, edges=[[0, 1]], ftype=[0])
    sage: g1 = GraphTheory(3, edges=[[0, 1]], ftype=[1])
    sage: g2 = GraphTheory(3, edges=[[0, 1]], ftype=[2])

Note that marking `2` defines a different flag than the other two since the
only nontrivial automorphism of the graph is the [0, 1] pair transposition ::
    
    sage: g0==g2
    False
    sage: g0==g1
    True

On flags we can do various computations. For example we can define an edge ::

    sage: e = GraphTheory(2, edges=[[0, 1]])

And calculate the sum of the two elements resulting in a 
FlagAlgebraElement ::

    sage: e+g
    Flag Algebra Element over Rational Field
    0   - Flag on 3 points, ftype from [] with edges=[]
    4/3 - Flag on 3 points, ftype from [] with edges=[[0, 2]]
    2/3 - Flag on 3 points, ftype from [] with edges=[[0, 2], [1, 2]]
    1   - Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]]

We can also calculate the product of the elements ::
    
    sage: e*g
    Flag Algebra Element over Rational Field
    1/5  - Flag on 5 points, ftype from [] with edges=[[0, 3], [1, 4]]
    ...
    3/10 - Flag on 5 points, ftype from [] with edges=[[0, 2], [0, 3], [0, 4], [1, 4], [2, 3], [2, 4], [3, 4]]
    1/10 - Flag on 5 points, ftype from [] with edges=[[0, 2], [0, 3], [0, 4], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]

In order to add or multiply elements, it is required to have the same ftype
and the result will have the same ftype. In the examples above 
`e, g, e+g, e*g` all have the same empty ftype. It is the default ftype when 
nothing is provided. `g0, g1, g2` have the same point
ftype ::
    
    sage: g.ftype()
    Ftype on 0 points with edges=[]
    sage: g1.ftype()
    Ftype on 1 points with edges=[]

This implementation defines ftypes as flags. The only distinction is that
an ftype has all it's points part of the ftype, while flags have points
outside the ftype.

To change between ftypes, one can project the flag to a smaller ftype. 
A projection is defined by an injective map from a smaller ftype into a larger
one. For example we can define the following flag with a largeish ftype ::

    sage: gg = GraphTheory(5, edges=[[0, 1], [2, 3]], ftype=[0, 2, 3, 4])
    sage: gg
    Flag on 5 points, ftype from [0, 2, 3, 4] with edges=[[0, 1], [2, 3]]
    sage: gg.ftype()
    Ftype on 4 points with edges=[[1, 2]]
    
The injection is simply a list of the smaller ftype's points inside
the larger ftype's points. For example picking `[0, 1, 2]` as the
injection, would result in an ftype on `[0, 2, 3]` points with the same 
edge `[2, 3]` which after renaming to `[0, 1, 2]` would be `[1, 2]`::
    
    sage: gg.project([0, 1, 2])
    Flag Algebra Element over Rational Field
    1/2 - Flag on 5 points, ftype from [0, 1, 4] with edges=[[0, 3], [1, 4]]
    sage: gg.project([0, 1, 2]).ftype()
    Ftype on 3 points with edges=[[1, 2]]
    
The result is normalized, but note that it is isomorphic to the one claimed
before. Three points `[0, 1, 4]` and one edge between, `[1, 4]`. If nothing
is provided then the flag is projected to the empty ftype ::
    
    sage: gg.project()
    Flag Algebra Element over Rational Field
    1/15 - Flag on 5 points, ftype from [] with edges=[[0, 3], [1, 4]]

For a little extra speed, multiplication and then projection is implemented::
    
    sage: gp = GraphTheory(2, ftype=[0])
    sage: (gp*gp).project()
    Flag Algebra Element over Rational Field
    1   - Flag on 3 points, ftype from [] with edges=[]
    1/3 - Flag on 3 points, ftype from [] with edges=[[0, 2]]
    0   - Flag on 3 points, ftype from [] with edges=[[0, 2], [1, 2]]
    0   - Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]]
    sage: gp.mul_project(gp) == (gp*gp).project()
    True

.. SEEALSO::
    :func:`Flag.__init__`
    :func:`Flag.ftype`
    :func:`Flag._add_`
    :func:`Flag._mul_`
    :func:`Flag.project`
    :func:`Flag.mul_project`
    :class:`CombinatorialTheory`
    :class:`FlagAlgebra`
    :class:`FlagAlgebraElement`

AUTHORS:

- Levente Bodnar (Dec 2023): Initial version

"""

# ****************************************************************************
#       Copyright (C) 2023 LEVENTE BODNAR <bodnalev at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools
from sage.rings.rational_field import QQ
from cysignals.signals cimport sig_check
from sage.structure.element cimport Element

cpdef _subblock_helper(list points, list block):
    r"""
    Helper to find induced substructures
    """
    cdef bint gd = 0
    ret = []
    if len(block)==0:
        return ret
    for xx in block:
        gd = 1
        for yy in xx:
            if yy not in points:
                gd = 0
                break
        if gd:
            ret.append([points.index(ii) for ii in xx])
    return ret

cdef class Flag(Element):
    
    cdef int _n
    cdef int _ftype_size
    
    cdef list _ftype_points
    cdef list _not_ftype_points
    cdef dict _blocks
    cdef tuple _unique
    
    cdef Flag _ftype
    
    def __init__(self, theory, n, **params):
        r"""
        Initialize a :class:`Flag` element

        INPUT:

        - ``theory`` -- :class:`CombinatorialTheory`; the underlying theory, 
            which is also the parent
        - ``n`` -- integer; the size (number of vertices) of the flag
        - ``**params`` -- optional parameters; can contain the points of the
            ftype as a list of points under the name "ftype" of "ftype_points", 
            if not provided then empty ftype is assumed. Can also contain the 
            blocks for each element in the signature, if not provided then an 
            empty block set is assumed.

        OUTPUT: The resulting flag

        EXAMPLES::

        Create a simple GraphTheory triangle ::
            
            sage: from sage.algebras.flag_algebras import *
            sage: from sage.algebras.flag import Flag
            sage: Flag(GraphTheory, 3, edges=[[0, 1], [0, 2], [1, 2]])
            Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]]
        
        Create a DiGraphTheory edge out from a pointed ftype ::

            sage: Flag(DiGraphTheory, 2, edges=[[0, 1]], ftype=[0])
            Flag on 2 points, ftype from [0] with edges=[[0, 1]]
        
        .. NOTE::

            It is recommended to create flags using the parent's element constructor

        .. SEEALSO::

            :func:`CombinatorialTheory._element_constructor_`
        """
        self._n = int(n)
        
        if 'ftype_points' in params:
            ftype_points = params['ftype_points']
        elif 'ftype' in params:
            ftype_points = params['ftype']
        else:
            ftype_points = []
        
        self._ftype_size = len(ftype_points)
        self._ftype_points = list(ftype_points)
        self._not_ftype_points = None
        self._blocks = {}
        for xx in theory._signature.keys():
            if xx in params:
                self._blocks[xx] = [list(yy) for yy in params[xx]]
            else:
                self._blocks[xx] = []
        self._unique = None
        
        self._ftype = None
        Element.__init__(self, theory)
    
    def _repr_(self):
        r"""
        Return a nice representation.
        
        If it is an ftype (all the points are part of the ftype), then
        the text changes to indicate this fact.

        EXAMPLES::

        For flags ::

            sage: from sage.algebras.flag_algebras import *
            sage: GraphTheory(3) 
            Flag on 3 points, ftype from [] with edges=[]
            
        For ftypes ::
        
            sage: ThreeGraphTheory(4, ftype_points=[0, 1, 2, 3], edges=[[0, 1, 2]])
            Ftype on 4 points with edges=[[0, 1, 2]]
        """
        blocks = self.blocks()
        strblocks = ', '.join([xx+'='+str(blocks[xx]) for xx in blocks.keys()])
        if self.is_ftype():
            return 'Ftype on {} points with {}'.format(self.size(), strblocks)
        return 'Flag on {} points, ftype from {} with {}'.format(self.size(), self.ftype_points(), strblocks)
    
    cpdef subflag(self, points=None, ftype_points=None):
        r"""
        Returns the induced subflag.
        
        The resulting sublaf contains the union of points and ftype_points
        and has ftype constructed from ftype_points. 

        INPUT:

        - ``points`` -- list (default: `None`); the points inducing the subflag.
            If not provided (or `None`) then this is the entire vertex set, so
            only the ftype changes
        - ``ftype_points`` list (default: `None`); the points inducing the ftype
            of the subflag. If not provided (or `None`) then the original ftype
            point set is used, so the result of the ftype will be the same

        OUTPUT: The induced sub Flag

        EXAMPLES::

        Same ftype ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3, edges=[[0, 1]], ftype=[0])
            sage: g.subflag([0, 2])
            Flag on 2 points, ftype from [0] with edges=[]
            
        Only change ftype ::
            
            sage: g.subflag(ftype_points=[0, 1])
            Flag on 3 points, ftype from [0, 1] with edges=[[0, 1]]

        .. NOTE::

            As the ftype points can be chosen, the result can have different
            ftype as self.

        TESTS::

            sage: g.subflag()==g
            True
        """
        if ftype_points==None:
            ftype_points = self._ftype_points
        if points==None:
            points = list(range(self._n))
        else:
            points = [ii for ii in range(self._n) if (ii in points or ii in ftype_points)]
        if len(points)==self._n and ftype_points==self._ftype_points:
            return self
        blocks = {xx: _subblock_helper(points, self._blocks[xx]) for xx in self._blocks.keys()}
        new_ftype_points = [points.index(ii) for ii in ftype_points]
        return self.__class__(self.parent(), len(points), ftype=new_ftype_points, **blocks)
    
    def combinatorial_theory(self):
        r"""
        Returns the combinatorial theory this flag is a member of
        
        This is the same as the parent.

        .. SEEALSO::

            :func:`theory`
            :func:`parent`
        """
        return self.parent()
    
    theory = combinatorial_theory
    
    def as_flag_algebra_element(self, basis=QQ):
        r"""
        Transforms this `Flag` to a `FlagAlgebraElement` over a given basis

        INPUT:

        - ``basis`` -- Ring (default: `QQ`); the base of
            the FlagAlgebra where the target will live.

        OUTPUT: A `FlagAlgebraElement` representing this `Flag`

        .. SEEALSO::

            :class:`FlagAlgebra`
            :func:`FlagAlgebra._element_constructor_`
            :class:`FlagAlgebraElement`
        """
        from sage.algebras.flag_algebras import FlagAlgebra
        targ_alg = FlagAlgebra(basis, self.theory(), self.ftype())
        return targ_alg(self)
    
    afae = as_flag_algebra_element
    
    def as_operand(self):
        r"""
        Turns this `Flag` into a `FlagAlgebraElement` so operations can be performed on it

        .. SEEALSO::

            :func:`as_flag_algebra_element`
        """
        return self.afae(QQ)
    
    cpdef size(self):
        r"""
        Returns the size of the vertex set of this Flag.

        OUTPUT: integer, the number of vertices.

        EXAMPLES::

        This is the size parameter in the `Flag` initialization ::

            sage: from sage.algebras.flag_algebras import *
            sage: GraphTheory(4).size()
            4
        """
        return self._n
    
    vertex_number = size
    
    cpdef blocks(self, as_tuple=False):
        r"""
        Returns the blocks

        INPUT:

        - ``as_tuple`` -- boolean (default: `False`); if the result should
            contain the blocks as a tuple

        OUTPUT: A dictionary, one entry for each element in the signature
            and list (or tuple) of the blocks for that signature.
        """
        if as_tuple:
            ret = {}
            for xx in self._blocks:
                ret[xx] = tuple([tuple(yy) for yy in self._blocks[xx]])
            return ret
        return self._blocks
    
    cpdef ftype(self):
        r"""
        Returns the ftype of this `Flag`

        EXAMPLES::

        Ftype of a pointed triangle is just a point ::

            sage: from sage.algebras.flag_algebras import *
            sage: pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0])
            sage: pointed_triangle.ftype()
            Ftype on 1 points with edges=[]
        
        And with two points it is ::
        
            sage: two_pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0, 1])
            sage: two_pointed_triangle.ftype()
            Ftype on 2 points with edges=[[0, 1]]

        .. NOTE::

            This is essentially the subflag, but the order of points matter. The result is saved
            for speed.

        .. SEEALSO::

            :func:`subflag`
        """
        if self._ftype==None:
            if self.is_ftype():
                self._ftype = self
            self._ftype = self.subflag([])
        return self._ftype
    
    cpdef ftype_points(self):
        r"""
        The points of the ftype inside self.
        
        This gives an injection of ftype into self

        OUTPUT: list of integers

        EXAMPLES::

            sage: from sage.algebras.flag_algebras import *
            sage: two_pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0, 1])
            sage: two_pointed_triangle.ftype_points()
            [0, 1]

        .. SEEALSO::

            :func:`__init__`
        """
        return self._ftype_points
    
    cpdef not_ftype_points(self):
        r"""
        This is a helper function, caches the points that are not
        part of the ftype.
        """
        if self._not_ftype_points != None:
            return self._not_ftype_points
        self._not_ftype_points = [ii for ii in range(self.size()) if ii not in self._ftype_points]
        return self._not_ftype_points
    
    def unique(self):
        r"""
        This returns a unique identifier that can equate isomorphic
        objects

        EXAMPLES::

        Isomorphic graphs have the same :func:`unique` value ::

            sage: from sage.algebras.flag_algebras import *
            sage: b1 = [[0, 1], [0, 2], [0, 4], [1, 3], [2, 4]]
            sage: b2 = [[0, 4], [1, 2], [1, 3], [2, 3], [3, 4]]
            sage: g1 = GraphTheory(5, edges=b1)
            sage: g2 = GraphTheory(5, edges=b2)
            sage: g1.unique() == g2.unique()
            True
            
        .. NOTE::

            The value returned here depends on the values of
            the parent `CombinatorialTheory`

        .. SEEALSO::

            :func:`theory`
            :func:`CombinatorialTheory.identify`
            :func:`__eq__`
        """
        if self._unique==None:
            self._unique = self.theory().identify(
                self._n, self._ftype_points, **self._blocks)
        return self._unique
    
    cpdef is_ftype(self):
        r"""
        Returns `True` if this flag is an ftype.

        .. SEEALSO::

            :func:`_repr_`
        """
        return self._n == self._ftype_size
    
    def _add_(self, other):
        r"""
        Add two Flags together
        
        The flags must have the same ftype. Different sizes are 
            all shifted to the larger one.

        OUTPUT: The :class:`FlagAlgebraElement` object, 
            which is the sum of the two parameters

        EXAMPLES::

        Adding to self is 2*self ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: g+g==2*g
            True
        
        Adding two distinct elements with the same size gives a vector 
        with exactly two `1` entries ::

            sage: h = GraphTheory(3, edges=[[0, 1]])
            sage: (g+h).values()
            (1, 1, 0, 0)
        
        Adding with different size the smaller flag
        is shifted to have the same size ::
        
            sage: e = GraphTheory(2)
            sage: (e+h).values()
            (1, 5/3, 1/3, 0)

        .. SEEALSO::

            :func:`FlagAlgebraElement._add_`
            :func:`__lshift__`

        """
        if self.ftype()!=other.ftype():
            raise TypeError("The terms must have the same ftype")
        return self.afae()._add_(other.afae())
    
    def _sub_(self, other):
        r"""
        Subtract a Flag from `self`
        
        The flags must have the same ftype. Different sizes are 
            all shifted to the larger one.

        EXAMPLES::
            
            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(2)
            sage: h = GraphTheory(3, edges=[[0, 1]])
            sage: (g-h).values()
            (1, -1/3, 1/3, 0)

        .. SEEALSO::

            :func:`_add_`
            :func:`__lshift__`
            :func:`FlagAlgebraElement._sub_`

        """
        if self.ftype()!=other.ftype():
            raise TypeError("The terms must have the same ftype")
        return self.afae()._sub_(other.afae())
    
    def _mul_(self, other):
        r"""
        Multiply two flags together.
        
        The flags must have the same ftype. The result
        will have the same ftype and size 
        `self.size() + other.size() - self.ftype().size()`

        OUTPUT: The :class:`FlagAlgebraElement` object, 
            which is the product of the two parameters

        EXAMPLES::

        Pointed edge multiplied by itself ::

            sage: from sage.algebras.flag_algebras import *
            sage: pe = GraphTheory(2, edges=[[0, 1]], ftype=[0])
            sage: (pe*pe).values()
            (0, 0, 0, 0, 1, 1)

        .. SEEALSO::

            :func:`FlagAlgebraElement._mul_`
            :func:`mul_project`
            :func:`CombinatorialTheory.mul_project_table`

        TESTS::

            sage: sum((pe*pe*pe*pe).values())
            11
            sage: e = GraphTheory(2)
            sage: (e*e).values()
            (1, 2/3, 1/3, 0, 2/3, 1/3, 0, 0, 1/3, 0, 0)
        """
        if self.ftype()!=other.ftype():
            raise TypeError("The terms must have the same ftype")
        return self.afae()._mul_(other.afae())
    
    def __lshift__(self, amount):
        r"""
        `FlagAlgebraElement`, equal to this, with size is shifted by the amount

        EXAMPLES::

        Edge shifted to size `3` ::

            sage: from sage.algebras.flag_algebras import *
            sage: edge = GraphTheory(2, edges=[[0, 1]])
            sage: (edge<<1).values()
            (0, 1/3, 2/3, 1)

        .. SEEALSO::

            :func:`FlagAlgebraElement.__lshift__`
        """
        return self.afae().__lshift__(amount)
    
    def __truediv__(self, other):
        r"""
        Divide by a scalar

        INPUT:

        - ``other`` -- number; any number such that `1` can be divided with that

        OUTPUT: The `FlagAlgebraElement` resulting from the division

        EXAMPLES::

        Divide by `2` ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: (g/2).values()
            (1/2, 0, 0, 0)
            
        Even for `x` symbolic `1/x` is defined, so the division is understood ::
            sage: var('x')
            x
            sage: g = GraphTheory(2)
            sage: g/x
            Flag Algebra Element over Symbolic Ring
            1/x - Flag on 2 points, ftype from [] with edges=[]
            0   - Flag on 2 points, ftype from [] with edges=[[0, 1]]
        
        .. NOTE::

            Dividing by `Flag` or `FlagAlgebraElement` is not allowed, only
            numbers such that the division is defined in some extension
            of the rationals.

        .. SEEALSO::

            :func:`FlagAlgebraElement.__truediv__`
        """
        return self.afae().__truediv__(other)
    
    def __eq__(self, other):
        r"""
        Compare two flags for == (equality)

        .. SEEALSO::

            :func:`unique`
            :func:`theory`
            :func:`CombinatorialTheory.identify`
        """
        if type(other)!=type(self):
            return False
        if self.parent()!=other.parent():
            return False
        return self.unique() == other.unique()
    
    def __lt__(self, other):
        r"""
        Compare two flags for < (proper induced inclusion)
        
        Returns true if self appears as a proper induced structure 
        inside other.

        .. SEEALSO::

            :func:`__le__`
        """
        if type(other)!=type(self):
            return False
        if self.parent()!=other.parent():
            return False
        if self.size()>=other.size():
            return False
        if self.ftype() != other.ftype():
            return False
        for subp in itertools.combinations(other.not_ftype_points(), self.size()-self.ftype().size()):
            if other.subflag(subp).__eq__(self):
                return True
        return False
    
    def __le__(self, other):
        r"""
        Compare two flags for <= (induced inclusion)
        
        Returns true if self appears as an induced structure inside
        other.

        EXAMPLES::

        Edge appears in a 4 star ::

            sage: from sage.algebras.flag_algebras import *
            sage: star = GraphTheory(4, edges=[[0, 1], [0, 2], [0, 3]])
            sage: edge = GraphTheory(2, edges=[[0, 1]])
            sage: edge <= star
            True
            
        The ftypes must agree ::
        
            sage: p_edge = GraphTheory(2, edges=[[0, 1]], ftype_points=[0])
            sage: p_edge <= star
            False
        
        But when ftypes agree, the inclusion must respect it ::
            
            sage: pstar = star.subflag(ftype_points=[0])
            sage: sub1 = GraphTheory(3, ftype=[0], edges=[[0, 1], [0, 2]])
            sage: sub1 <= pstar
            True
            sage: sub2 = GraphTheory(3, ftype=[1], edges=[[0, 1], [0, 2]])
            sage: sub2 <= pstar
            False

        .. SEEALSO::

            :func:`__lt__`
            :func:`__eq__`
            :func:`unique`
        """
        return self==other or self<other
    
    def __hash__(self):
        r"""
        A hash based on the unique identifier
        so this is compatible with `__eq__`.
        """
        return hash(self.unique())
    
    def __getstate__(self):
        r"""
        Saves this flag to a dictionary
        """
        dd = {'theory': self.theory(),
              'n': self._n, 
              'ftype_points': self._ftype_points, 
              'blocks':self._blocks, 
              'unique':self._unique}
        return dd
    
    def __setstate__(self, dd):
        r"""
        Loads this flag from a dictionary
        """
        self._set_parent(dd['theory'])
        self._n = dd['n']
        self._ftype_points = dd['ftype_points']
        self._ftype_size = len(self._ftype_points)
        self._not_ftype_points = None
        self._blocks = dd['blocks']
        self._unique = dd['unique']
    
    def project(self, ftype_inj=tuple()):
        r"""
        Project this `Flag` to a smaller ftype
        

        INPUT:

        - ``ftype_inj`` -- tuple (default: (, )); the injection of the
            projected ftype inside the larger ftype

        OUTPUT: the `FlagAlgebraElement` resulting from the projection

        EXAMPLES::

        If the center of a cherry is flagged, then the projection has
        coefficient 1/3 ::

            sage: from sage.algebras.flag_algebras import *
            sage: p_cherry = GraphTheory(3, edges=[[0, 1], [0, 2]], ftype_points=[0])
            sage: p_cherry.project().values()
            (0, 0, 1/3, 0)

        .. NOTE::

            If `ftype_inj==tuple(range(self.ftype().size()))` then this
            does nothing.

        .. SEEALSO::

            :func:`FlagAlgebraElement.project`
        """
        return self.afae().project(ftype_inj)
    
    def mul_project(self, other, ftype_inj=tuple()):
        r"""
        Multiply self with other, and the project the result.

        INPUT:

        - ``ftype_inj`` -- tuple (default: (, )); the injection of the
            projected ftype inside the larger ftype

        OUTPUT: the `FlagAlgebraElement` resulting from the multiplication
            and projection

        EXAMPLES::

        Pointed edge multiplied with itself and projected ::

            sage: from sage.algebras.flag_algebras import *
            sage: p_edge = GraphTheory(2, edges=[[0, 1]], ftype_points=[0])
            sage: p_edge.mul_project(p_edge).values()
            (0, 0, 1/3, 1)

        .. NOTE::

            If `ftype_inj==tuple(range(self.ftype().size()))` then this
            is the same as usual multiplication.

        .. SEEALSO::

            :func:`_mul_`
            :func:`project`
            :func:`FlagAlgebraElement.mul_project`
        """
        return self.afae().mul_project(other, ftype_inj)
    
    def density(self, other):
        r"""
        The density of self in other.
        
        Randomly choosing self.size() points in other, the
        probability of getting self.

        EXAMPLES::

        Density of an edge in the cherry graph is 2/3 ::

            sage: from sage.algebras.flag_algebras import *
            sage: cherry = GraphTheory(3, edges=[[0, 1], [0, 2]])
            sage: edge = GraphTheory(2, edges=[[0, 1]])
            sage: cherry.density(edge)
            2/3
        
        .. SEEALSO::
        
            :func:`FlagAlgebraElement.density`
        """
        safae = self.afae()
        oafae = safae.parent(other)
        return self.afae().density(other)
    
    def _ftypes_inside(self, target):
        r"""
        Returns the possible ways self ftype appears in target

        INPUT:

        - ``target`` -- Flag; the flag where we are looking for copies of self

        OUTPUT: list of Flags with ftype matching as self, not necessarily unique
        """
        ret = []
        lrp = list(range(target.size()))
        for ftype_points in itertools.permutations(range(target.size()), self._n):
            if target.subflag(ftype_points, ftype_points)==self:
                ret.append(target.subflag(lrp, ftype_points))
        return ret
    
    cpdef densities(self, n1, n1flgs, n2, n2flgs, ftype_remap, large_ftype, small_ftype):
        r"""
        Returns the density matrix, indexed by the entries of `n1flgs` and `n2flgs`
        
        The matrix returned has entry `(i, j)` corresponding to the possibilities of
        `n1flgs[i]` and `n2flgs[j]` inside self, projected to the small ftype. 
        
        This is the same as counting the ways we can choose `n1` and `n2` points
        inside `self`, such that the two sets cover the entire `self.size()` point
        set, and calculating the probability that the overlap induces an ftype
        isomorphic to `large_ftype` and the points sets are isomorphic to 
        `n1flgs[i]` and `n2flgs[j]`.

        INPUT:

        - ``n1`` -- integer; the size of the first flag list
        - ``n1flgs`` -- list of flags; the first flag list (each of size `n1`)
        - ``n2`` -- integer; the size of the second flag list
        - ``n2flgs`` -- list of flags; the second flag list (each of size `n2`)
        - ``ftype_remap`` -- list; shows how to remap `small_ftype` into `large_ftype`
        - ``large_ftype`` -- ftype; the ftype of the overlap
        - ``small_ftype`` -- ftype; the ftype of self

        OUTPUT: a sparse matrix corresponding with the counts

        .. SEEALSO::

            :func:`CombinatorialTheory.mul_project_table`
            :func:`FlagAlgebra.mul_project_table`
            :func:`FlagAlgebraElement.mul_project`
        """
        cdef int N = self._n
        cdef int small_size = small_ftype.size()
        cdef int large_size = large_ftype.size()
        cdef int ctr = 0
        cdef bint chk = 0
        
        ret = {}
        small_points = self._ftype_points
        for difference in itertools.permutations(self.not_ftype_points(), large_size - small_size):
            sig_check()
            large_points = [0]*len(ftype_remap)
            for ii in range(len(ftype_remap)):
                vii = ftype_remap[ii]
                if vii<small_size:
                    large_points[ii] = small_points[vii]
                else:
                    large_points[ii] = difference[vii-small_size]
            ind_large_ftype = self.subflag([], ftype_points=large_points)
            if ind_large_ftype==large_ftype:
                not_large_points = [ii for ii in range(N) if ii not in large_points]
                for n1_extra_points in itertools.combinations(not_large_points, n1 - large_size):
                    n2_extra_points = [ii for ii in not_large_points if ii not in n1_extra_points]
                    try:
                        n1_ind = n1flgs.index(self.subflag(n1_extra_points, ftype_points=large_points))
                    except ValueError:
                        subf = self.subflag(n1_extra_points, ftype_points=large_points)
                        raise ValueError("Could not find \n", subf, "\nin the list of ", \
                                         n1, " sized flags with ", large_ftype, \
                                         ".\nThis can happen if the generator and identifier ",\
                                         "(from the current CombinatorialTheory) is incompatible, ",\
                                         "or if the theory is not heredetary")
                    try:
                        n2_ind = n2flgs.index(self.subflag(n2_extra_points, ftype_points=large_points))
                    except:
                        subf = self.subflag(n2_extra_points, ftype_points=large_points)
                        raise ValueError("Could not find \n", subf, "\nin the list of ", \
                                         n2, " sized flags with ", large_ftype, \
                                         ".\nThis can happen if the generator and identifier ",\
                                         "(from the current CombinatorialTheory) is incompatible, ",\
                                         "or if the theory is not heredetary")
                    try:
                        ret[(n1_ind, n2_ind)] += 1
                    except:
                        ret[(n1_ind, n2_ind)] = 1
        return (len(n1flgs), len(n2flgs), ret)