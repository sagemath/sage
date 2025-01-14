r"""

.. SEEALSO::
    :func:`CombinatorialTheory.__init__`
    :func:`CombinatorialTheory.exclude`
    :func:`CombinatorialTheory.optimize_problem`
    :func:`CombinatorialTheory.generate_flags`

AUTHORS:

- Levente Bodnar (2023-2025): Main development

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

from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element, get_coercion_model
from sage.all import QQ, Integer
from sage.algebras.flag import Flag, Pattern

from sage.all import vector


class FlagAlgebraElement(Element):
    def __init__(self, parent, n, values):
        r"""
        Initialize a Flag Algebra Element
        
        INPUT:

        - ``parent`` -- FlagAlgebra; The parent :class:`FlagAlgebra`
        - ``n`` -- integer; the size of the flags
        - ``values`` -- vector; the values representing the
            linear combination of :class:`Flag` elements

        OUTPUT: A FlagAlgebraElement object with the given initial parameters

        .. NOTE::

            It is recommended to use the FlagAlgebra element constructor

        .. SEEALSO::

            :func:`FlagAlgebra._element_constructor_`
        """
        if len(values)!=parent.get_size(n):
            raise ValueError("The coefficients must have the same length " + 
                             "as the number of flags")
        self._n = n
        base = parent.base()
        self._values = vector(base, values, sparse=True)
        Element.__init__(self, parent)
    
    def ftype(self):
        r"""
        Return the ftype of the parent FlagAlgebra
        
        The algebra is defined over elements of a CombinatorialTheory, all
        having the same ftype.

        OUTPUT: The ftype of the parent FlagAlgebra. A :class:`Flag` element

        EXAMPLES::

        The ftype of a :class:`Flag` is the same as the ftype of the 
        :class:`FlagAlgebraElement` we can construct from it ::
        
            
            sage: g = GraphTheory(3)
            sage: g.ftype()
            Ftype on 0 points with edges=()
            sage: g.ftype()==g.afae().ftype()
            True
        
        .. SEEALSO::

            :func:`ftype` in :class:`FlagAlgebra`
            :func:`ftype` in :class:`Flag`
        """
        return self.parent().ftype()
    
    def size(self):
        r"""
        Return the size of the vertex set of flags in this element
        
        OUTPUT: The size of each flag is :func:`flags`. 

        TESTS::

            
            sage: FG = FlagAlgebra(GraphTheory, QQ)
            sage: FGElem = FG._an_element_()
            sage: FGElem.size() == FGElem.flags()[0].size()
            True
        """
        return self._n
    
    vertex_number = size
    
    def flags(self):
        r"""
        Returns the list of flags, corresponding to each base element
        
        The flags returned are the list of flags with size equal to 
        `self.size()` and ftype to `self.ftype()`. Their number is 
        the same as the length of `self.values()`. 

        OUTPUT: The list of flags

        EXAMPLES::

        3 vertex graphs with empty ftype ::

            
            sage: g = GraphTheory(3)
            sage: g.afae()
            Flag Algebra Element over Rational Field
            1 - Flag on 3 points, ftype from () with edges=()
            0 - Flag on 3 points, ftype from () with edges=(01)
            0 - Flag on 3 points, ftype from () with edges=(01 02)
            0 - Flag on 3 points, ftype from () with edges=(01 02 12)
            sage: g.afae().flags()
            (Flag on 3 points, ftype from () with edges=(),
             Flag on 3 points, ftype from () with edges=(01),
             Flag on 3 points, ftype from () with edges=(01 02),
             Flag on 3 points, ftype from () with edges=(01 02 12))

        .. NOTE::

            This is the same as 
            `self.theory().generate_flags(self.size(), self.ftype())`

        .. SEEALSO::

            :func:`CombinatorialTheory.generate_flags`
            :func:`size`
            :func:`ftype`
            :func:`values`
            :func:`Flag.afae`

        TESTS::

            
            sage: g.afae().flags() == g.theory().generate_flags(g.size(), g.ftype())
            True
        """
        return self.parent().generate_flags(self._n)
    
    def values(self):
        r"""
        Returns the vector of values, corresponding to each element 
        in :func:`flags`
        
        OUTPUT: A vector

        EXAMPLES::

        A flag transformed to a flag algebra element has 
        all zeroes except one entry, itself ::

            
            sage: g = GraphTheory(3)
            sage: g.afae().values()
            (1, 0, 0, 0)

        .. SEEALSO::

            :func:`flags`
            :func:`_vector_`
            :func:`Flag.afae`
        """
        return self._values
    
    def _vector_(self, R):
        r"""
        Returns the vector of values, corresponding to each element 
        in :func:`flags` in a given base.
        
        OUTPUT: A vector

        EXAMPLES::

        A flag transformed to a flag algebra element has 
        all zeroes except one entry, itself ::

            
            sage: g = GraphTheory(3)
            sage: g.afae()._vector_(QQ['x'])
            (1, 0, 0, 0)

        .. SEEALSO::

            :func:`values()`
            :func:`Flag.afae`
        """
        return vector(R, self._values, sparse=True)
    
    def __len__(self):
        r"""
        Returns the length, the number of elements
        in :func:`flags` or :func:`values` (which is the same)

        EXAMPLES::

            sage: g = GraphTheory(3)
            sage: len(g.afae())
            4

        .. SEEALSO::

            :func:`flags`
            :func:`values`
            :func:`Flag.afae`
            :func:`__iter__`
        """
        return len(self._values)
    
    def __iter__(self):
        r"""
        Iterates over the elements of self
        
        It yields (number, flag) tuples, indicating the
        coefficient of the flag.

        EXAMPLES::

            sage: g = GraphTheory(3)
            sage: for x in g.afae():
            ....:   print(x)
            (1, Flag on 3 points, ftype from () with edges=())
            (0, Flag on 3 points, ftype from () with edges=(01))
            (0, Flag on 3 points, ftype from () with edges=(01 02))
            (0, Flag on 3 points, ftype from () with edges=(01 02 12))


        .. SEEALSO::

            :func:`flags`
            :func:`values`
            :func:`Flag.afae`
            :func:`__len__`
        """
        for ii, fl in enumerate(self.parent().generate_flags(self.size())):
            yield (self._values[ii], fl)
    
    def _repr_(self):
        r"""
        Give a string representation
        
        Lists the flags and the corresponding coefficients,
        each on a separate line. If the list is too long 
        then only shows nonzero entries.


        EXAMPLES::

        Short list, so display all ::

            
            sage: gf = GraphTheory(3).afae()
            sage: gf
            Flag Algebra Element over Rational Field
            1 - Flag on 3 points, ftype from () with edges=()
            0 - Flag on 3 points, ftype from () with edges=(01)
            0 - Flag on 3 points, ftype from () with edges=(01 02)
            0 - Flag on 3 points, ftype from () with edges=(01 02 12)
            
        Long list, only the nonzero entries are displayed::
        
            sage: g1 = GraphTheory(5)
            sage: g2 = GraphTheory(5, edges=[[0, 1], [3, 4]])
            sage: g1+g2
            Flag Algebra Element over Rational Field
            1 - Flag on 5 points, ftype from () with edges=()
            1 - Flag on 5 points, ftype from () with edges=(02 14)
            
        .. SEEALSO::

            :func:`Flag._repr_`
        """
        sttrl = ['Flag Algebra Element over {}'.format(
            self.parent().base()
            )]
        strs = [str(xx) for xx in self.values()]
        maxstrlen = max([len(xx) for xx in strs])
        for ii, fl in enumerate(self.parent().generate(self.size())):
            if len(self)<10:
                sttrl.append(('{:<'+str(maxstrlen)+'} - {}').format(
                    strs[ii], str(fl)
                    ))
            else:
                include = True
                try: 
                    include = abs(float(self.values()[ii]))>=1e-8
                except: 
                    include = self.values()[ii]!=0
                if include:
                    sttrl.append(('{:<'+str(maxstrlen)+'} - {}').format(
                        strs[ii], str(fl)
                        ))
        return "\n".join(sttrl)
    
    def base(self):
        return self.parent().base()

    def theory(self):
        return self.parent().theory()

    def custom_coerce(self, other):
        if isinstance(other, Flag) or isinstance(other, Pattern):
            if self.ftype()!=other.ftype():
                raise ValueError("The ftypes must agree.")
            alg = self.parent()
            return (self, alg(other))
        elif isinstance(other, FlagAlgebraElement):
            if self.ftype()!=other.ftype():
                raise ValueError("The ftypes must agree.")
            sbase = self.base()
            obase = other.base()
            base = get_coercion_model().common_parent(sbase, obase)
            alg = FlagAlgebra(self.theory(), base, self.ftype())
            return (alg(self), alg(other))
        else:
            base = get_coercion_model().common_parent(
                self.base(), other.parent()
                )
            alg = FlagAlgebra(self.theory(), base, self.ftype())
            return (alg(self), alg(base(other)))

    def as_flag_algebra_element(self):
        r"""
        Returns self.
        
        Only here to allow calling this function on
        both flags and flag algebra elements

        .. SEEALSO::

            :func:`Flag.afae`
        """
        return self
    
    afae = as_flag_algebra_element
    
    def _add_(self, other):
        r"""
        Adds to FlagAlgebraElements together

        OUTPUT: The sum

        EXAMPLES::

        The smaller size is shifted to match the larger ::

            
            sage: g = GraphTheory(3).afae()
            sage: e = GraphTheory(2).afae()
            sage: e+g
            Flag Algebra Element over Rational Field
            2   - Flag on 3 points, ftype from () with edges=()
            2/3 - Flag on 3 points, ftype from () with edges=(01)
            1/3 - Flag on 3 points, ftype from () with edges=(01 02)
            0   - Flag on 3 points, ftype from () with edges=(01 02 12)

        .. NOTE::

            The result's size will match the size of the larger component

        .. SEEALSO::

            :func:`Flag._add_`
            :func:`__lshift__`
            :func:`_sub_`
        """
        nm = max(self.size(), other.size())
        vals = (self<<(nm-self.size())).values() + \
            (other<<(nm-other.size())).values()
        return self.__class__(self.parent(), nm, vals)
    
    def _sub_(self, other):
        r"""
        Subtract a FlagAlgebraElement from this

        EXAMPLES::

        This also shifts the smaller flag to match the larger ::

            
            sage: g = GraphTheory(3).afae()
            sage: e = GraphTheory(2).afae()
            sage: e-g
            Flag Algebra Element over Rational Field
            0   - Flag on 3 points, ftype from () with edges=()
            2/3 - Flag on 3 points, ftype from () with edges=(01)
            1/3 - Flag on 3 points, ftype from () with edges=(01 02)
            0   - Flag on 3 points, ftype from () with edges=(01 02 12)

        .. SEEALSO::

            :func:`Flag._sub_`
            :func:`__lshift__`
            :func:`_add_`
        """
        nm = max(self.size(), other.size())
        vals = (self<<(nm-self.size())).values() - \
            (other<<(nm-other.size())).values()
        return self.__class__(self.parent(), nm, vals)
    
    def _neg_(self):
        return self.__class__(
            self.parent(), 
            self.size(), 
            self.values()*(-1))

    def _mul_(self, other):
        r"""
        Multiplies two elements together
        
        The result will have size 
        `self.size() + other.size() - self.ftype().size()`

        EXAMPLES::

        Two empty edges multiplied together has size 4 ::

            
            sage: e = GraphTheory(2).afae()
            sage: (e*e).size()
            4
            
        But if pointed (size of ftype is 1), then the size is 3 ::
            
            sage: pe = GraphTheory(2, ftype=[0])
            sage: pe*pe
            Flag Algebra Element over Rational Field
            1 - Flag on 3 points, ftype from (0,) with edges=()
            0 - Flag on 3 points, ftype from (0,) with edges=(01)
            1 - Flag on 3 points, ftype from (2,) with edges=(01)
            0 - Flag on 3 points, ftype from (0,) with edges=(01 02)
            0 - Flag on 3 points, ftype from (1,) with edges=(01 02)
            0 - Flag on 3 points, ftype from (0,) with edges=(01 02 12)
            
        Can also multiply with constants:
            
            sage: pe*3
            Flag Algebra Element over Rational Field
            3 - Flag on 2 points, ftype from (0,) with edges=()
            0 - Flag on 2 points, ftype from (0,) with edges=(01)
        
        .. NOTE::
            
            Multiplying and then projecting to a smaller ftype
            can be performed more efficiently with :func:`mul_project`
        
        .. SEEALSO::

            :func:`Flag._mul_`
            :func:`mul_project`
        """
        if self.size()==self.ftype().size():
            vals = other.values() * self.values()[0]
            return self.__class__(self.parent(), other.size(), vals)
        if other.size()==self.ftype().size():
            vals = self.values() * other.values()[0]
            return self.__class__(self.parent(), self.size(), vals)
        table = self.parent().mpt(self.size(), other.size())
        N = -self.ftype().size() + self.size() + other.size()
        vals = [self.values() * mat * other.values() for mat in table]
        return self.__class__(self.parent(), N, vals)
    
    def __truediv__(self, other):
        r"""
        Divide by a scalar

        INPUT:

        - ``other`` -- number; any number such that `1` can be divided with that

        OUTPUT: The `FlagAlgebraElement` resulting from the division

        EXAMPLES::

        If 1 can be divided by that, then the division is allowed ::

            
            sage: var('x')
            x
            sage: g = GraphTheory(2)
            sage: g.afae()/x
            Flag Algebra Element over Symbolic Ring
            1/x - Flag on 2 points, ftype from () with edges=()
            0   - Flag on 2 points, ftype from () with edges=(01)
        
        .. NOTE::
            
            This is the linear extension of :func:`Flag.__truediv__`
        
        .. SEEALSO::

            :func:`Flag.afae`
            :func:`Flag.__truediv__`
        """
        return self * (1/other)
    
    def __lshift__(self, amount):
        r"""
        `FlagAlgebraElement`, equal to this, with size is 
        shifted by the amount
        
        The result will have size equal to 
        `self.size() + amount`, but the elements will be equal
        
        EXAMPLES::

        Edge shifted to size `3` ::

            
            sage: edge = GraphTheory(2, edges=[[0, 1]])
            sage: (edge.afae()<<1).values()
            (0, 1/3, 2/3, 1)

        .. NOTE::
            
            This is the linear extension of :func:`Flag.__lshift__`

        .. SEEALSO::

            :func:`Flag.__lshift__`
            :func:`Flag.afae()`
        """
        if amount<0:
            raise ValueError('Can not have negative shifts')
        if amount==0:
            return self
        ressize = amount + self.size()
        table = self.parent().mpt(self.size(), self.ftype().size(), 
                                  target_size=ressize)
        vals = [sum(self.values() * mat) for mat in table]
        return self.__class__(self.parent(), ressize, vals)
    
    def __getitem__(self, flag):
        if isinstance(flag, Flag):
            ind = self.parent().get_index(flag)
        elif isinstance(flag, Integer) and \
            0 <= flag and \
                flag < self.parent().get_size(self.size()):
            ind = flag
        if ind == -1:
            raise TypeError("Indecies must be Flags with matching " + 
                            "ftype and size, or integers. " + 
                            "Not {}".format(str(type(flag))))
        return self._values[ind]
    
    def __setitem__(self, flag, value):
        if isinstance(flag, Flag):
            ind = self.parent().get_index(flag)
        elif isinstance(flag, Integer) and \
            0 <= flag and \
                flag < self.parent().get_size(self.size()):
            ind = flag
        if ind == -1:
            raise TypeError("Indecies must be Flags with matching " + 
                            "ftype and size, or integers. " + 
                            "Not {}".format(str(type(flag))))
        self.values()[self.flags().index(flag)] = value
    
    def project(self, ftype_inj=tuple()):
        r"""
        Project this to a smaller ftype
        

        INPUT:

        - ``ftype_inj`` -- tuple (default: (, )); the injection of the
            projected ftype inside the larger ftype

        OUTPUT: the `FlagAlgebraElement` resulting from the projection

        EXAMPLES::

        If the center of a cherry is flagged, then the projection has
        coefficient 1/3 ::

            
            sage: p_cherry = GraphTheory(3, edges=[[0, 1], [0, 2]], ftype_points=[0])
            sage: p_cherry.afae().project().values()
            (0, 0, 1/3, 0)

        .. NOTE::
            
            This is the linear extension of :func:`Flag.project`

            If `ftype_inj==tuple(range(self.ftype().size()))` then this
            does nothing.

        .. SEEALSO::

            :func:`Flag.project`
        """
        return self.mul_project(1, ftype_inj)
    
    def mul_project(self, other, ftype_inj=tuple(), target_size=None):
        r"""
        Multiply self with other, and the project the result.

        INPUT:

        - ``ftype_inj`` -- tuple (default: (, )); the injection of the
            projected ftype inside the larger ftype

        OUTPUT: the `FlagAlgebraElement` resulting from the multiplication
            and projection

        EXAMPLES::

        Pointed edge multiplied with itself and projected ::

            
            sage: felem = GraphTheory(2, edges=[[0, 1]], ftype_points=[0]).afae()
            sage: felem.mul_project(felem).values()
            (0, 0, 1/3, 1)

        .. NOTE::
            
            This is the bilinear extension of :func:`Flag.mul_project`

            If `ftype_inj==tuple(range(self.ftype().size()))` then this
            is the same as usual multiplication.

        .. SEEALSO::

            :func:`_mul_`
            :func:`project`
            :func:`Flag.mul_project`
        """
        other = self.parent(other)
        ftype_inj = tuple(ftype_inj)
        new_ftype = self.ftype().subflag([], ftype_points=ftype_inj)
        if new_ftype==None or new_ftype.unique()==None:
            raise ValueError("The ftype injection maps to an invalid ftype.")
        N = -self.ftype().size() + self.size() + other.size()
        if target_size!=None:
            if target_size<N:
                raise ValueError(
                    "Target size is smaller than minimum allowed size " + 
                    "for this operation."
                    )
            N = target_size
        table = self.parent().mpt(
            self.size(), other.size(), 
            ftype_inj=ftype_inj, target_size=N
            )
        vals = [self.values() * mat * other.values() for mat in table]
        
        TargetAlgebra = FlagAlgebra(
            self.parent().combinatorial_theory(), self.parent().base(), new_ftype
            )
        return TargetAlgebra(N, vals)
    
    def density(self, other):
        r"""
        The density of self in other.
        
        Randomly choosing self.size() points in other, the
        probability of getting self.

        EXAMPLES::

        Density of an edge in the cherry graph is 2/3 ::

            
            sage: cherry = GraphTheory(3, edges=[[0, 1], [0, 2]]).afae()
            sage: edge = GraphTheory(2, edges=[[0, 1]]).afae()
            sage: cherry.density(edge)
            2/3

        .. NOTE::
            
            This is the bilinear extension of :func:`Flag.mul_project`
        
        .. SEEALSO::
        
            :func:`Flag.density`
        """
        diff = self.size() - other.size()
        if diff < 0:
            raise ValueError('The target can not be larger')
        return self.values() * (other<<diff).values()
    
    def set_sum(self, target=1):
        r"""
        Make the symbolic variables sum to the given target.
        """
        valvec = self.values()
        if len(valvec)==0:
            return self
        par = valvec[0].parent()
        try:
            gs = par.gens()
        except:
            return self
        repl = {gs[-1]: target - sum(gs[:-1])}
        nvec = vector([xx.subs(repl) for xx in valvec])
        return self.parent()(self.size(), nvec)

    def subs(self, args, ring=QQ):
        r"""
        Symbolic substitution.
        """
        valvec = self.values()
        if len(valvec)==0:
            return self
        par = valvec[0].parent()
        try:
            gs = par.gens()
        except:
            return self
        repl = {gs[ii]:args[ii] for ii in range(min(len(args), len(gs)))}
        nvec = vector([xx.subs(repl) for xx in valvec])
        retalg = FlagAlgebra(self.parent().theory(), ring)
        return retalg(self.size(), nvec)

    def derivative(self, times):
        r"""
        Symbolic differentiation.
        """
        valvec = self.values()
        if len(valvec)==0:
            return self
        par = valvec[0].parent()
        try:
            gs = par.gens()
        except:
            return self
        rvec = []
        for ff, vv in enumerate(valvec):
            if vv==0:
                rvec.append(0)
                continue
            aval = vv
            for ii in range(min(len(times), len(gs))):
                aval = aval.derivative(gs[ii], times[ii])
            rvec.append(aval)
        return self.parent()(self.size(), vector(rvec))

    def derivatives(self, point, ring=QQ):
        r"""
        Returns all symbolic derivatives evaluated at a given point.
        """
        times = []
        valvec = self.values()
        if len(valvec)==0:
            return self
        par = valvec[0].parent()
        try:
            gs = par.gens()
        except:
            return self
        for xx in self.values():
            if xx==0:
                continue
            dd = xx.dict()
            for kk in dd.keys():
                kk = tuple(kk)
                if kk in times:
                    continue
                for ll in itertools.product(*[range(ii+1) for ii in kk]):
                    if ll not in times:
                        times.append(ll)
        res = []
        for xx in times:
            der = self.derivative(xx).subs(point, ring=ring)
            minnz = None
            for yy in der.values():
                if ring(yy)!=0:
                    if minnz==None:
                        minnz = abs(yy)
                    else:
                        minnz = min(abs(yy), minnz)
            if minnz != None:
                res.append(der/minnz)
        return res
    
    def _richcmp_(self, other, op):
        r"""
        Compares the elements.
        
        Since the parent agrees, the ftype too. They are shifted to the
        same size and the values compared elementwise.

        EXAMPLES::

        Trivial example `g <= 2*g` ::

            
            sage: g = GraphTheory(3).afae()
            sage: g <= 2*g
            True
            
        Since `g` has zero coefficients ::
        
            sage: g < 2*g
            False
            
        Shifting the size gives equal elements ::
        
            sage: g == (g<<1)
            True
            
        More sophisticated example, Mantel's theorem.
        Using the fact that for large enough structures
        `x.mul_project(x)` is always positive, we get that 
        `1/2 >= GraphTheory(2)`, so K_3 free graphs can not
        contain more than 1/2 density of K_2 ::
            
            sage: GraphTheory.exclude(GraphTheory(3))
            sage: f = GraphTheory(2, ftype=[0]) - 1/2
            sage: 1/2 >= f.mul_project(f)*2 + GraphTheory(2)
            True

        .. NOTE::
            
            When comparing Flags with FlagAlgebraElements, it
            will return False, unless the Flag is transformed
            to a FlagAlgebraElement.

        .. SEEALSO::

            :func:`mul_project`
        """
        nm = max(self.size(), other.size())
        v1 = (self<<(nm-self.size())).values()
        v2 = (other<<(nm-other.size())).values()
        return all([richcmp(v1[ii], v2[ii], op) for ii in range(len(v1))])

class FlagAlgebra(Parent, UniqueRepresentation):
    
    def __init__(self, theory, base=QQ, ftype=None):
        r"""
        Initialize a FlagAlgebra

        INPUT:

        - ``base`` -- Ring; The base ring, this FlagAlgebra is constructed over
            This must contain the rationals `QQ`
        - ``theory`` -- CombinatorialTheory; The combinatorial theory
            this flag algebra is based on
        - ``ftype`` -- Flag (default=`None`); The ftype of the elements
            in this FlagAlgebra. The default `None` gives the empty type

        OUTPUT: The resulting FlagAlgebra

        EXAMPLES::

        Create the FlagAlgebra for GraphTheory (without any ftype) ::

            sage: GraphFlagAlgebra = FlagAlgebra(GraphTheory, QQ)
            sage: GraphFlagAlgebra
            Flag Algebra with Ftype on 0 points with edges=() over Rational Field
        """
        if ftype==None:
            ftype = theory.empty_element()
        else:
            if not ftype.is_ftype():
                raise ValueError('{} is not an Ftype'.format(ftype))
            if ftype.theory() != theory:
                raise ValueError('{} is not a part of {}'.format(ftype, theory))
        if not base.has_coerce_map_from(QQ):
            raise ValueError('The base must contain the rationals')
        self._theory = theory
        self._ftype = ftype
        self._index_set = {}
        self._size_set = {}
        self._curr_excluded = theory._excluded
        Parent.__init__(self, base)
    
    Element = FlagAlgebraElement
    
    def get_index_set(self, n):
        if self._theory._excluded != self._curr_excluded:
            self._curr_excluded = self._theory._excluded
            self._size_set = {}
            self._index_set = {}
        if n not in self._index_set:
            fls = self.generate_flags(n)
            fldict = dict(zip(fls, range(len(fls))))
            self._index_set[n] = fldict
        return self._index_set[n]

    def get_index(self, flag):
        if not isinstance(flag, Flag):
            return -1
        if flag.ftype()!=self.ftype():
            return -1
        indn = self.get_index_set(flag.size())
        if flag not in indn:
            return -1
        return indn[flag]

    def get_size(self, n):
        if self._theory._excluded != self._curr_excluded:
            self._curr_excluded = self._theory._excluded
            self._size_set = {}
            self._index_set = {}
        if n not in self._size_set:
            self._size_set[n] = len(self.generate_flags(n))
        return self._size_set[n]

    def _element_constructor_(self, *args, **kwds):
        r"""
        Constructs a FlagAlgebraElement with the given parameters
        
        If a single value is provided it constructs the constant in the Algebra
        If a Flag is provided it constructs the element whose only `1` value is
        that Flag.
        If a different FlagAlgebraElement is provided, then checks if the
        `theory` and the `ftype` agrees, then tries to coerce the values to the
        base of self.
        Otherwise uses the constructor of FlagAlgebraElement, which accepts a
        size value and a list of coefficients, whose length must be precisely the
        number of flags.

        EXAMPLES::

        Construct from a constant ::

            
            sage: FA = FlagAlgebra(GraphTheory, QQ)
            sage: FA(3)
            Flag Algebra Element over Rational Field
            3 - Ftype on 0 points with edges=()
        
        Construct from a flag ::
        
            sage: g = GraphTheory(2)
            sage: el = FA(g)
            sage: el
            Flag Algebra Element over Rational Field
            1 - Flag on 2 points, ftype from () with edges=()
            0 - Flag on 2 points, ftype from () with edges=(01)
            
        Construct from a FlagAlgebraElement with smaller base ::
        
            sage: FAX = FlagAlgebra(GraphTheory, QQ['x'])
            sage: FAX(el)
            Flag Algebra Element over Univariate Polynomial Ring in x over Rational Field
            1 - Flag on 2 points, ftype from () with edges=()
            0 - Flag on 2 points, ftype from () with edges=(01)
            
        Constructing the element directly from coefficients ::
            
            sage: FA(2, [3, 4])
            Flag Algebra Element over Rational Field
            3 - Flag on 2 points, ftype from () with edges=()
            4 - Flag on 2 points, ftype from () with edges=(01)

        .. SEEALSO::

            :func:`FlagAlgebraElement.__init__`
        """
        if len(args)==1:
            v = args[0]
            base = self.base()
            if isinstance(v, Flag):
                ind = self.get_index(v)
                size = self.get_size(v.size())
                if ind!=-1:
                    return self.element_class(
                        self, v.size(), vector(self.base(), size, {ind:1})
                        )
            elif isinstance(v, Pattern):
                if self.ftype() == v.ftype():
                    dvec = {self.get_index(xx):1 for xx in v.compatible_flags()}
                    size = self.get_size(v.size())
                    return self.element_class(
                        self, v.size(), vector(self.base(), size, dvec)
                        )
            elif isinstance(v, FlagAlgebraElement):
                if v.ftype()==self.ftype():
                    if self.base()==v.parent().base():
                        return v
                    elif self.base().has_coerce_map_from(v.parent().base()):
                        vals = vector(self.base(), v.values(), sparse=True)
                        return self.element_class(self, v.size(), vals)
            elif v in base:
                return self.element_class(self, self.ftype().size(), vector(base, [v], sparse=True))
            raise ValueError('Can\'t construct an element from {} for the theory\n{}'.format(v, self))
        return self.element_class(self, *args, **kwds)
    
    def _coerce_map_from_(self, S):
        r"""
        Checks if it can be coerced from S
        """
        if self.base().has_coerce_map_from(S):
            return True
        if S==self.theory():
            return True
        if isinstance(S, FlagAlgebra):
            if S.ftype()==self.ftype() and self.base().has_coerce_map_from(S.base()):
                return True
        return False
    
    def _pushout_(self, S):
        r"""
        Constructs the pushout FlagAlgebra
        """
        if S.has_coerce_map_from(self.base()):
            return FlagAlgebra(self.theory(), S, self.ftype())
        return None
    
    def _repr_(self):
        r"""
        Returns a short text representation

        EXAMPLES::

            
            sage: FlagAlgebra(GraphTheory, QQ)
            Flag Algebra with Ftype on 0 points with edges=() over Rational Field

        .. SEEALSO::

            :func:`Flag._repr_`
        """
        return 'Flag Algebra with {} over {}'.format(self.ftype(), self.base())
    
    def ftype(self):
        r"""
        Returns the ftype of this FlagAlgebra.

        EXAMPLES::

        Without specifying anything in the constructor, the ftype
        is empty ::

            
            sage: FA = FlagAlgebra(GraphTheory, QQ)
            sage: FA.ftype()
            Ftype on 0 points with edges=()

        .. NOTE::

            This is the same ftype as the ftype of the elements, 
            and the ftype of the flags in those elements. 

        .. SEEALSO::

            :func:`Flag.ftype`
            :func:`FlagAlgebraElement.ftype`
        """
        return self._ftype
    
    def combinatorial_theory(self):
        r"""
        Returns the :class:`CombinatorialTheory` object, whose
        flags form the basis of this FlagAlgebra

        EXAMPLES::

        This is the same as provided in the constructor ::

            
            sage: FA = FlagAlgebra(GraphTheory, QQ)
            sage: FA.theory()
            Theory for Graph

        .. SEEALSO::

            :func:`__init__`
            :class:`CombinatorialTheory`
        """
        return self._theory
    
    theory = combinatorial_theory
    
    def base_ring(self):
        r"""
        Returns the base_ring
        
        Same as `base_ring` of the `base` provided in the constructor

        EXAMPLES::

            
            sage: FA = FlagAlgebra(GraphTheory, QQ)
            sage: FA.base_ring()
            Rational Field

        .. SEEALSO::

            :func:`base`
        """
        return self.base().base_ring()
    
    def characteristic(self):
        r"""
        Returns the characteristic
        
        Same as `characteristic` of the `base` provided in the constructor

        EXAMPLES::

            
            sage: FA = FlagAlgebra(GraphTheory, QQ)
            sage: FA.characteristic()
            0

        .. SEEALSO::

            :func:`base`
        """
        return self.base().characteristic()
    
    def generate_flags(self, n):
        r"""
        Generates flags of a given size, and `ftype` of `self`

        .. NOTE::

            Same as `CombinatorialTheory.generate_flags` with `ftype`
            from `self`

        .. SEEALSO::

            :func:`CombinatorialTheory.generate_flags`
        """
        return self.theory().generate_flags(n, self.ftype())
    
    generate = generate_flags

    def _an_element_(self):
        r"""
        Returns an element
        """
        a = self.base().an_element()
        f = self.combinatorial_theory()._an_element_(n=self.ftype().size(), ftype=self.ftype())
        return self(f)*a
    
    def some_elements(self):
        r"""
        Returns a small list of elements
        """
        return [self.an_element(),self(self.base().an_element())]
    
    def mul_project_table(self, n1, n2, ftype_inj=None, target_size=None):
        r"""
        Returns the multiplication projection table
        
        This is the same as :func:`CombinatorialTheory.mul_project_table`
        with self.ftype() substituted in.

        .. SEEALSO::

            :func:`CombinatorialTheory.mul_project_table`
        """
        return self.theory().mul_project_table(n1, n2, self.ftype(), ftype_inj, target_size)
    
    mpt = mul_project_table
