r"""
Quotients of finite rank free modules over a field.
"""

####################################################################################
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
####################################################################################

from .free_module import FreeModule_ambient_field


class FreeModule_ambient_field_quotient(FreeModule_ambient_field):
    """
    A quotient `V/W` of two vector spaces as a vector space.

    To obtain `V` or `W` use ``self.V()`` and ``self.W()``.

    EXAMPLES::

        sage: k.<i> = QuadraticField(-1)
        sage: A = k^3; V = A.span([[1,0,i], [2,i,0]])
        sage: W = A.span([[3,i,i]])
        sage: U = V/W; U
        Vector space quotient V/W of dimension 1 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I where
        V: Vector space of degree 3 and dimension 2 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
        Basis matrix:
        [ 1  0  i]
        [ 0  1 -2]
        W: Vector space of degree 3 and dimension 1 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
        Basis matrix:
        [    1 1/3*i 1/3*i]
        sage: U.V()
        Vector space of degree 3 and dimension 2 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
        Basis matrix:
        [ 1  0  i]
        [ 0  1 -2]
        sage: U.W()
        Vector space of degree 3 and dimension 1 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
        Basis matrix:
        [    1 1/3*i 1/3*i]
        sage: U.quotient_map()
        Vector space morphism represented by the matrix:
        [  1]
        [3*i]
        Domain: Vector space of degree 3 and dimension 2 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
        Basis matrix:
        [ 1  0  i]
        [ 0  1 -2]
        Codomain: Vector space quotient V/W of dimension 1 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I where
        V: Vector space of degree 3 and dimension 2 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
        Basis matrix:
        [ 1  0  i]
        [ 0  1 -2]
        W: Vector space of degree 3 and dimension 1 over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
        Basis matrix:
        [    1 1/3*i 1/3*i]
        sage: Z = V.quotient(W)
        sage: Z == U
        True

    We create three quotient spaces and compare them::

        sage: A = QQ^2
        sage: V = A.span_of_basis([[1,0], [1,1]])
        sage: W0 = V.span([V.1, V.0])
        sage: W1 = V.span([V.1])
        sage: W2 = V.span([V.1])
        sage: Q0 = V/W0
        sage: Q1 = V/W1
        sage: Q2 = V/W2

        sage: Q0 == Q1
        False
        sage: Q1 == Q2
        True

    TESTS::

        sage: A = QQ^0; V = A.span([]) # corner case
        sage: W = A.span([])
        sage: U = V/W

        sage: loads(dumps(U)) == U
        True
        sage: type(loads(dumps(U)) )
        <class 'sage.modules.quotient_module.FreeModule_ambient_field_quotient_with_category'>
    """
    def __init__(self, domain, sub, quotient_matrix, lift_matrix, inner_product_matrix=None):
        """
        Create this quotient space, from the given domain, submodule,
        and quotient_matrix.

        EXAMPLES::

            sage: A = QQ^5; V = A.span_of_basis([[1,0,-1,1,1], [1,-1,0,2/3,3/4]]); V
            Vector space of degree 5 and dimension 2 over Rational Field
            User basis matrix:
            [  1   0  -1   1   1]
            [  1  -1   0 2/3 3/4]
            sage: W = V.span_of_basis([V.0 - 2/3*V.1]); W
            Vector space of degree 5 and dimension 1 over Rational Field
            User basis matrix:
            [1/3 2/3  -1 5/9 1/2]

        This creates a quotient vector space::

            sage: Q = V / W

        Behold the type of Q::

            sage: type(Q)
            <class 'sage.modules.quotient_module.FreeModule_ambient_field_quotient_with_category'>

        We do some consistency checks on the extra quotient and
        lifting structure of Q::

            sage: Q(V.0)
            (1)
            sage: Q( V.0 - 2/3*V.1 )
            (0)
            sage: v = Q.lift(Q.0); v
            (1, 0, -1, 1, 1)
            sage: Q( v )
            (1)
        """
        base_field = domain.base_field()
        dimension = quotient_matrix.ncols()
        sparse = domain.is_sparse()
        self.__sub = sub
        self.__domain = domain
        self.__hash = hash((domain, sub))
        FreeModule_ambient_field.__init__(self, base_field, dimension, sparse)
        self.__quo_map = domain.Hom(self)(quotient_matrix)
        self.__quo_map.register_as_coercion()
        self.__lift_map = self.Hom(domain)(lift_matrix)

    def _repr_(self):
        r"""
        Return the rather verbose string representation of this quotient space V/W.

        EXAMPLES:

        We create a quotient vector space over a finite field::

            sage: k.<a> = GF(9); A = k^3; V = A.span_of_basis([[1,0,a], [a,a,1]]); W = V.span([V.1])
            sage: Q = V/W

        Note the type::

            sage: type(Q)
            <class 'sage.modules.quotient_module.FreeModule_ambient_field_quotient_with_category'>

        The string representation mentions that this is a quotient
        `V/W`, that the quotient has dimension 1 and is over a finite
        field, and also describes `V` and `W`::

            sage: Q._repr_()
            'Vector space quotient V/W of dimension 1 over Finite Field in a of size 3^2 where\nV: Vector space of degree 3 and dimension 2 over Finite Field in a of size 3^2\nUser basis matrix:\n[1 0 a]\n[a a 1]\nW: Vector space of degree 3 and dimension 1 over Finite Field in a of size 3^2\nBasis matrix:\n[    1     1 a + 2]'
        """
        return "%s space quotient V/W of dimension %s over %s where\nV: %s\nW: %s"%(
            "Sparse vector" if self.is_sparse() else "Vector",
            self.dimension(), self.base_ring(),
            self.V(), self.W())

    def __hash__(self):
        """
        Return hash of this quotient space `V/W`, which is, by definition,
        the hash of the tuple `(V, W)`.

        EXAMPLES:

        We compute the hash of a certain 0-dimension quotient vector
        space::

            sage: A = QQ^2; V = A.span_of_basis([[1,0], [1,1]]); W = V.span([V.1, V.0])
            sage: Q = V/W; Q.dimension()
            0
            sage: hash(Q) == hash((V,W))
            True
        """
        return self.__hash

    def _element_constructor_(self, x):
        """
        Convert an element into this quotient space `V/W` if there is
        a way to make sense of it.

        An element converts into self if it can be converted into `V`,
        or if not at least if it can be made sense of as a list of
        length the dimension of self.

        EXAMPLES:

        We create a 2-dimensional quotient of a 3-dimension ambient
        vector space::

            sage: M = QQ^3 / [[1,2,3]]

        A list of length 3 converts into ``QQ^3``, so it converts into
        `M`::

            sage: M([1,2,4])  #indirect doctest
            (-1/3, -2/3)
            sage: M([1,2,3])
            (0, 0)

        A list of length 2 converts into M, where here it just gives
        the corresponding linear combination of the basis for `M`::

            sage: M([1,2])
            (1, 2)
            sage: M.0 + 2*M.1
            (1, 2)

        Of course, elements of ``QQ^3`` convert into the quotient
        module as well. Here is a different example::

            sage: V = QQ^3; W = V.span([[1,0,0]]); Q = V/W
            sage: Q(V.0)
            (0, 0)
            sage: Q(V.1)
            (1, 0)
            sage: Q.0
            (1, 0)
            sage: Q.0 + V.1
            (2, 0)

        Here we start with something that is over ZZ, so it
        canonically coerces into ``QQ^3``, hence into ``self``::

            sage: Q((ZZ^3)([1,2,3]))
            (2, 3)

        """
        if isinstance(x, self.element_class) and x.parent() is self:
            return x
        if isinstance(x, (list, tuple)) and len(x) == self.__domain.rank():
            return self.__quo_map(self.__domain(x))
        return FreeModule_ambient_field._element_constructor_(self, x)

    def _coerce_map_from_(self, M):
        """
        Return a coercion map from `M` to ``self``, or None.

        EXAMPLES::

            sage: V = QQ^2 / [[1, 2]]
            sage: V.coerce_map_from(ZZ^2)
            Composite map:
              From: Ambient free module of rank 2 over the principal ideal domain Integer Ring
              To:   Vector space quotient V/W of dimension 1 over Rational Field where
            V: Vector space of dimension 2 over Rational Field
            W: Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 2]
              Defn:   Coercion map:
                      From: Ambient free module of rank 2 over the principal ideal domain Integer Ring
                      To:   Vector space of dimension 2 over Rational Field
                    then
                      Vector space morphism represented by the matrix:
                    [   1]
                    [-1/2]
                    Domain: Vector space of dimension 2 over Rational Field
                    Codomain: Vector space quotient V/W of dimension 1 over Rational Field where
                    V: Vector space of dimension 2 over Rational Field
                    W: Vector space of degree 2 and dimension 1 over Rational Field
                    Basis matrix:
                    [1 2]

        Make sure :trac:`10513` is fixed (no coercion from an abstract
        vector space to an isomorphic quotient vector space)::

            sage: V = QQ^3 / [[1,2,3]]
            sage: V.coerce_map_from(QQ^2)

        """
        from sage.modules.free_module import FreeModule_ambient
        if (isinstance(M, FreeModule_ambient)
            and not (isinstance(M, FreeModule_ambient_field_quotient)
                     and self.W() == M.W())):
            # No map between different quotients.
            # No map from quotient to abstract module.
            return None
        f = super(FreeModule_ambient_field, self)._coerce_map_from_(M)
        if f is not None:
            return f
        f = self.__domain.coerce_map_from(M)
        if f is not None:
            return self.__quo_map * f
        return None

    def quotient_map(self):
        """
        Given this quotient space `Q = V / W`, return the natural quotient
        map from `V` to `Q`.

        EXAMPLES::

            sage: M = QQ^3 / [[1,2,3]]
            sage: M.quotient_map()
            Vector space morphism represented by the matrix:
            [   1    0]
            [   0    1]
            [-1/3 -2/3]
            Domain: Vector space of dimension 3 over Rational Field
            Codomain: Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of dimension 3 over Rational Field
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]

            sage: M.quotient_map()( (QQ^3)([1,2,3]) )
            (0, 0)
        """
        return self.__quo_map

    def lift_map(self):
        r"""
        Given this quotient space `Q = V / W`, return a fixed choice of
        linear homomorphism (a section) from `Q` to `V`.

        EXAMPLES::

            sage: M = QQ^3 / [[1,2,3]]
            sage: M.lift_map()
            Vector space morphism represented by the matrix:
            [1 0 0]
            [0 1 0]
            Domain: Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of dimension 3 over Rational Field
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
            Codomain: Vector space of dimension 3 over Rational Field
        """
        return self.__lift_map

    def lift(self, x):
        r"""
        Lift element of this quotient `V / W` to `V` by applying
        the fixed lift homomorphism.

        The lift is a fixed homomorphism.

        EXAMPLES::

            sage: M = QQ^3 / [[1,2,3]]
            sage: M.lift(M.0)
            (1, 0, 0)
            sage: M.lift(M.1)
            (0, 1, 0)
            sage: M.lift(M.0 - 2*M.1)
            (1, -2, 0)
        """
        return self.__lift_map(x)

    def W(self):
        """
        Given this quotient space `Q = V/W`, return `W`.

        EXAMPLES::

            sage: M = QQ^10 / [list(range(10)), list(range(2,12))]
            sage: M.W()
            Vector space of degree 10 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1 -2 -3 -4 -5 -6 -7 -8]
            [ 0  1  2  3  4  5  6  7  8  9]
        """
        return self.__sub

    def V(self):
        """
        Given this quotient space `Q = V/W`, return `V`.

        EXAMPLES::

            sage: M = QQ^10 / [list(range(10)), list(range(2,12))]
            sage: M.V()
            Vector space of dimension 10 over Rational Field
        """
        return self.__domain

    def cover(self):
        """
        Given this quotient space `Q = V/W`, return `V`.

        This is the same as :meth:`V`.

        EXAMPLES::

            sage: M = QQ^10 / [list(range(10)), list(range(2,12))]
            sage: M.cover()
            Vector space of dimension 10 over Rational Field
        """
        return self.V()

    def relations(self):
        """
        Given this quotient space `Q = V/W`, return `W`.

        This is the same as :meth:`W`.

        EXAMPLES::

            sage: M = QQ^10 / [list(range(10)), list(range(2,12))]
            sage: M.relations()
            Vector space of degree 10 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1 -2 -3 -4 -5 -6 -7 -8]
            [ 0  1  2  3  4  5  6  7  8  9]
        """
        return self.W()

