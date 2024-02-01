r"""
Generalized monomial orders on `\ZZ^n` defined in [PU1999]_. 

EXAMPLES::

    sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder

    sage: order = GeneralizedMonomialOrder(2); order
    Generalized monomial order in 2 variables using (lex, min)

    sage: order.n_cones()
    3
    sage: order.cone(1)
    [-1  0]
    [-1  1]
    sage: order.greatest_tuple((1,2),(-2,3))
    (-2, 3)
    sage: order.greatest_tuple([(1,2),(-3,2),(2,-1)])
    (-3, 2)
    sage: order.greatest_tuple_for_cone(2,[(-1,2), (-1,2), (-2,-3)])
    (-2, -3)

AUTHORS:

- Legrand Lucas (2024-26-01): initial version

"""
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.rings.polynomial.term_order import TermOrder
from sage.matrix.constructor import matrix
from sage.structure.sage_object import SageObject

class GeneralizedMonomialOrder(SageObject):
    def __init__(self,n,group_order="lex",score_function="min"):
        r"""
        Create a generalized monomial order.

        INPUTS:

            - ``n`` -- The number of variables.
            - ``group_order`` (default: ``lex``) -- The name of a group order on `\ZZ^n`. Choices are: "lex".
            - ``score_function`` (default: ``min``) -- The name of a score function. Choices are: "min", "degmin".

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: GeneralizedMonomialOrder(2)
            Generalized monomial order in 2 variables using (lex, min)

            sage: GeneralizedMonomialOrder(8, group_order="lex")
            Generalized monomial order in 8 variables using (lex, min)

            sage: GeneralizedMonomialOrder(3, score_function="min")
            Generalized monomial order in 3 variables using (lex, min)

            sage: GeneralizedMonomialOrder(3, group_order="lex", score_function="degmin")
            Generalized monomial order in 3 variables using (lex, degmin)
        """
        self._n = n
        self._n_cones = self._n + 1
        self._cones = build_cones(self._n) 
        self._group_order = get_group_order(group_order)
        self._group_order_name = group_order
        self._score_function = get_score_function(score_function)
        self._score_function_name = score_function

    def __hash__(self):
        r"""
        Return the hash of self. It depends on the number of variables, the group_order and the score function
        """
        return hash(self._group_order_name + self._score_function_name + str(self._n))


    def _repr_(self):
        r"""

        TESTS::
            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: GeneralizedMonomialOrder(2)._repr_()
            'Generalized monomial order in 2 variables using (lex, min)'
        """
        group = self._group_order_name
        function = self._score_function_name
        n = str(self._n)
        return "Generalized monomial order in %s variables using (%s, %s)" % (n, group, function)

    def __hash__(self):
        r"""
        Return the hash of self. It depends on the number of variables, the group_order and the score function.
        """
        return hash(self._group_order_name + self._score_function_name + str(self._n))

    def n_cones(self):
        r"""
        Return the number of cones (which is the number of variables plus one).

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2)
            sage: order.n_cones()
            3
        """
        return self._n_cones

    def cones(self):
        r"""
        Return the list of matrices containing the generators of the cones.

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2)
            sage: order.cones()
            [
            [1 0]  [-1  0]  [ 1 -1]
            [0 1], [-1  1], [ 0 -1]
            ]
        """
        return self._cones

    def cone(self,i):
        r"""
        Return the matrix whose columns are the generators of the `i`-th cone.
    
        INPUTS:
            - `i` -- a cone index
        
        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2)
            sage: order.cone(1)
            [-1 0]
            [-1 1]
        """
        if i < 0 or i > self._n:
            raise ValueError("Cone index out of range")
        return self._cones[i]

    def group_order(self):
        r"""
        Return the underlying group order.
        
        EXAMPLES::
           
            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)
            sage: G.group_order()
            Lexicographic term order
        """
        return self._group_order

    def group_order_name(self):
        r"""
        Return the name of the underlying group order.

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2,group_order="lex")
            sage: order.group_order_name()
            'lex'
        """
        return self._group_order_name
    
    def score_function_name(self):
        r"""
        Return the name of the underlying score function.

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2,score_function="degmin")
            sage: order.score_function_name()
            'degmin'
        """
        return self._score_function_name

    def score(self,t):
        r"""
        Compute the score of a tuple ``t``.

        INPUTS:

            - `t` -- a tuple

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2)
            sage: order.score((-2,3))
            2
        """
        return self._score_function(t)
                                
    def compare(self,a,b):
        r"""
        Return 1, 0 or -1 whether tuple ``a`` is greater than, equal or less than to tuple ``b`` respectively.

        INPUTS:

            - `a` and `b` -- two tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(3)
            sage: order.compare((1,2,-2), (3,-4,5))
            -1
            sage: order.compare((3,-4,5), (1,2,-2))
            1
            sage: order.compare((3,-4,5), (3,-4,5))
            0
        """

        if not isinstance(a, tuple) or not isinstance(b, tuple):
            raise ValueError("Args must be tuples")
        if a == b: return 0
        diff = self._score_function(a) - self._score_function(b)
        if diff != 0: 
            return 1 if diff > 0 else -1
        else: 
            return 1 if self._group_order.greater_tuple(a,b) == a else -1

    def greatest_tuple(self,a,b=None):
        r"""
        Return the greatest of the two tuples ``a`` and ``b`` if ``b`` is not None
        Return the greatest tuple in the list ``L`` if ``b`` is not a tuple

        INPUTS:

            - ``a``,``b`` -- two tuples or ``a`` a list of tuples and ``b`` anything except a tuple

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(3)
            sage: order.greatest_tuple((0,0,1),(2,3,-2))
            (2, 3, -2)
            sage: L = [(1,2,-1),(3,-3,0),(4,-5,-6)]
            sage: order.greatest_tuple(L)
            (4, -5, -6)
        """
        
        if (not isinstance(a,tuple) or not isinstance(b, tuple)) and not isinstance(a,list):
            raise TypeError("Args must be either two tuples, or first argument must be a list of tuples")
        if isinstance(a,tuple):
            return a if self.compare(a,b) == 1 else b
        else:
            if len(a) == 1: return a[0]
            else: return self.greatest_tuple(a[0],self.greatest_tuple(a[1:]))

    def translate_to_cone(self,i,L):
        r"""
        Return a tuple ``t`` such that `t + L` is contained in the ``i``-th cone.

        INPUTS:

            - ``i`` -- a cone index
            - ``L`` -- a list of tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2)
            sage: order.translate_to_cone(0, [(1,2),(-2,-3),(1,-4)])
            (2, 4)
        """
        cone_matrix = self._cones[i]
        T = matrix(ZZ, [cone_matrix*vector(v) for v in L])
        return tuple(cone_matrix*vector([-min(0,*c) for c in T.columns()]))
        
    def greatest_tuple_for_cone(self, i, L):
        r"""
        Return the greatest tuple of the list of tuple `L` with respect to the `i`-th cone.

        This is the unique tuple `t` in `L`  such that each time the greatest tuple of `s + L`
        for a tuple `s` is contained in the `i`-th cone, it is equal to `s + t`.

        INPUTS:

            - ``i`` -- a cone index
            - ``L`` -- a list of tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2)
            sage: L = [(1,2), (-2,-2), (-4,5), (5,-6)]
            sage: t = order.greatest_tuple_for_cone(1, L);t
            (-4, 5)

        We can check the result::

            sage: s = order.translate_to_cone(1, L)
            sage: sL = [tuple(vector(s) + vector(l)) for l in L]
            sage: tuple(vector(s) + vector(t)) == order.greatest_tuple(sL)
            True

        """
        t = vector(self.translate_to_cone(i,L))
        L = [vector(l) for l in L]
        return tuple(vector(self.greatest_tuple([tuple(t + l) for l in L])) - t)
              
    def is_in_cone(self,i,v): 
        r"""
        Test whether the tuple `t` is contained in the `i`-th cone or not.

        INPUTS:
            - ``i`` -- a cone index
            - ``t`` -- a tuple

        EXAMPLES::
            
            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2)
            sage: order.is_in_cone(0, (1,2)) 
            True
            sage: order.is_in_cone(0,(-2,3))
            False
        """
        return all(c >=0 for c in self._cones[i]*vector(v))

    def generator(self,i,L):
        r"""
        Return the generator of the module over the ``i``-th cone for ``L``.

        This is the monoïd of elements `t \in \ZZ^n` such that the greatest 
        tuple of ``t + L`` is contained in the ``i``-th cone.

        INPUTS:
            
            - ``i`` -- a cone index
            - ``L`` -- a list of tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(3)
            sage: L = [(1,2,-2), (3,4,-5), (-6,2,-7), (3,-6,1)]
            sage: order.generator(0,L)
            (6, 6, 7)

            sage: order.generator(2,L)
            (5, 5, 6)
        """
    
        cone_matrix = self._cones[i]
        t = vector(self.translate_to_cone(i,L))
        L = [vector(l) for l in L]
        for c in cone_matrix.columns():
            while self.is_in_cone(i,self.greatest_tuple([tuple(t + l) for l in L])):
                t = t - c
            t = t + c
        return tuple(t)

    def generator_for_pair(self,i,L1,L2):
        r"""
        Return the generator of the module over the ``i``-th cone for ``L1`` and ``L2``. 

        INPUTS:

            - ``i`` -- A cone index.
            - ``L1`` -- A list of tuples.
            - ``L2`` - -A list of tuples.

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: order = GeneralizedMonomialOrder(2)
            sage: L1 = [(2,3),(-4,2),(-1,-2)]
            sage: L2 = [(1,-6),(5,-2),(3,4)]
            sage: order.generator_for_pair(1,L1,L2)
            (1, 5)

        """
        cone_matrix = self._cones[i]
        lm1 = vector(self.greatest_tuple_for_cone(i,L1))
        lm2 = vector(self.greatest_tuple_for_cone(i,L2))
        g1 = vector(self.generator(i,L1))
        g2 = vector(self.generator(i,L2))
        m = vector([max(a,b) for a,b in zip(lm1 + g1, lm2 + g2)])
        return tuple(cone_matrix*m)

def build_cones(n):
    r"""
    Return a list of integer matrices such that the colums of the `j`-th matrix in the list
    are the generators of the `j`-th cone.
    
    INPUTS:

        - ``n`` -- The number of variables (which is the number of cones minus one).

    TESTS::

        sage: from sage.rings.polynomial.generalized_monomial_order import build_cones
        sage: build_cones(3)
        [
        [1 0 0]  [-1  0  0]  [ 1 -1  0]  [ 1  0 -1]
        [0 1 0]  [-1  1  0]  [ 0 -1  0]  [ 0  1 -1]
        [0 0 1], [-1  0  1], [ 0 -1  1], [ 0  0 -1]
        ]
    """
    L = [matrix.identity(ZZ,n)]
    for i in range(0,n):
        mat = matrix.identity(ZZ,n)
        mat.set_column(i,vector(ZZ,[-1]*n))
        L.append(mat)
    return L

def get_score_function(name):
    r"""
    Return the score function specified by ``name``.
    
    INPUTS:
        
        - ``name`` -- Name of a score function within "min", "degmin"
    """
    
    def min_score_function(t):
        return -min(0,*t)

    def degmin_score_function(t):
        return sum(t) -(len(t)+1)*min(0,*t)

    if name == "min":
        return min_score_function
    elif name == "degmin":
        return degmin_score_function
    else:
        raise ValueError("Available score functions are: 'min', 'degmin'")
        
def get_group_order(name):
    r"""
    Return the group order specified by ``name``.

    INPUTS:

        - ``name`` -- Name of a group order within "lex"

    TESTS::

        sage: from sage.rings.polynomial.generalized_monomial_order import get_group_order
        sage: get_group_order("lex")
        Lexicographic term order
    """
    
    if name in ["lex"]:
        return TermOrder(name=name)
    else:
        raise ValueError("Available group order are: 'lex', 'invlex'")


