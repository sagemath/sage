r"""
Generalized monomial orders on the monoïd `\ZZ^n`. 

Generalized monomial orders on `\ZZ^n` have been introduced by Pauer and Unterkircher.
These are total orders on `\ZZ^Unterkircher 


EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Legrand Lucas (2024-26-01): initial version

"""

from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.rings.polynomial.term_order import TermOrder
from sage.matrix.constructor import matrix

class GeneralizedMonomialOrder:
    def __init__(self,n,group_order="lex",score_function="min"):
        r"""
        Create a generalized monomial order.

        INPUTS:

        - ``n`` -- the number of variables
        - ``group_order`` (default: ``lex``) -- a group order on `\ZZ^n`. A class with a method 'greater_tuple'
          taking as inputs two tuples in `\ZZ^n` and returning the greatest with respect to a group order
        - ``score_function`` (default: ``min``) -- a function `\ZZ^n \to \QQ_{\le 0}` satisfying propeties (1) and (2) 

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)

        """
        self._n = n
        self._cones = self.__build_cones() 
        self.__set_group_order(group_order)
        self.__set_score_function(score_function)

    def __build_cones(self):
        r"""
        Return a list of integer matrices such that the colums of the `j`-th matrix in the list
        are the generators of the `j`th cone in self._cones
        
        TESTS::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(3)
            sage: G.get_cones()
            [
            [1 0 0]  [-1  0  0]  [ 1 -1  0]  [ 1  0 -1]
            [0 1 0]  [-1  1  0]  [ 0 -1  0]  [ 0  1 -1]
            [0 0 1], [-1  0  1], [ 0 -1  1], [ 0  0 -1]
            ]
        """
        L = [matrix.identity(ZZ,self._n)]
        for i in range(0,self._n):
            mat = matrix.identity(ZZ,self._n)
            mat.set_column(i,vector(ZZ,[-1]*self._n))
            L.append(mat)
        return L

    def get_cones(self):
        r"""
        Return the list of matrices containing the generators of the cones

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)
            sage: G.get_cones()
            [
            [1 0]  [-1  0]  [ 1 -1]
            [0 1], [-1  1], [ 0 -1]
            ]
        """
        return self._cones

    def get_cone(self,i):
        r"""
        Return the matrix whose columns are the generators of the `i`-th cone.
    
        INPUTS:
            - `i` -- a cone index
        
        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)
            sage: G.get_cone(1)
            [-1 0]
            [-1 1]
        """
        if i < 0 or i > self._n:
            raise ValueError("Cone index out of range")
        return self._cones[i]

    def __set_group_order(self,group_order):
        r"""
        Set the underlying group order

        INPUTS:

            - ``group_order`` -- either one of the strings ('lex', ...), or any class with a public method
              greater_tuple comparing two tuple with respect to a group order on `\ZZ^n`

        TESTS::
            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)
            sage: G.get_group_order()
            Lexicographic term order
        """
        if isinstance(group_order,str):
            if group_order == "lex":
                self._group_order = TermOrder(name="lex")
            else:
                raise ValueError("Currently, only 'lex' is available as built-in group order. You can also provide your own defined group order.")
        else:
            self._group_order = group_order

    def get_group_order(self):
        r"""
        Return the underlying group order.
        
        EXAMPLES::
           
            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)
            sage: G.get_group_order()
            Lexicographic term order
        """
        return self._group_order

    def __set_score_function(self,score_function):
        r"""
        Set the underlying score function

        INPUTS:

            - score_functionction -- A function from tuples to the positive rationals
              satisfying properties (1) and (2)

        TESTS::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)
            sage: G.score((1,-6))
            6
        """
        if isinstance(score_function,str):
            if score_function == "min":
                self._score_function = min_score_function
            else:
                raise ValueError("Currently, only 'min' is available as built-in score function. You can also provide your own defined score function.")
        else:
            self._score_function = score_function

    def score(self,t):
        r"""
        Compute the score of a tuple ``t``

        INPUTS:

            - `t` -- a tuple

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)
            sage: G.score((-2,3))
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
            sage: G = GeneralizedMonomialOrder(3)
            sage: G.compare((1,2,-2), (3,-4,5))
            -1
            sage: G.compare((3,-4,5), (1,2,-2))
            1
            sage: G.compare((3,-4,5), (3,-4,5))
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

    def greatest_tuple(self,a,b=0):
        r"""
        Return the greatest of the two tuples ``a`` and ``b`` if ``b`` is not zero
        Return the greatest tuple in the list ``L`` if ``b`` is not a tuple

        INPUTS:

            - ``a``,``b`` -- two tuples or ``a`` a list of tuples and ``b`` anything except a tuple

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(3)
            sage: G.greatest_tuple((0,0,1),(2,3,-2))
            (2, 3, -2)
            sage: L = [(1,2,-1),(3,-3,0),(4,-5,-6)]
            sage: G.greatest_tuple(L)
            sage: (4, -5, -6)
        """
        
        if (not isinstance(a,tuple) or not isinstance(b, tuple)) and not isinstance(a,list):
            raise TypeError("Args must be either two tuples, or first argument must be a list of tuples")
        if b != 0:
            return a if self.compare(a,b) == 1 else b
        else:
            if len(a) == 1: return a[0]
            else: return self.greatest_tuple(a[0],self.greatest_tuple(a[1:]))

    def translate_to_cone(self,i,L):
        r"""
        Return a tuple `t` such that `t + L` is included in the ``i``-th cone.

        INPUTS:

            - ``i`` -- a cone index
            - ``L`` -- a list of tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)
            sage: G.translate_to_cone(0, [(1,2),(-2,-3),(1,-4)])
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
            sage: G = GeneralizedMonomialOrder(2)
            sage: L = [(1,2), (-2,-2), (-4,5), (5,-6)]
            sage: t = G.greatest_tuple_for_cone(1, L);t
            (-4, 5)

        We can check the result::

            sage: s = G.translate_to_cone(1, L)
            sage: sL = [tuple(vector(s) + vector(l)) for l in L]
            sage: tuple(vector(s) + vector(t)) == G.greatest_tuple(sL)
            True

        """
        t = vector(self.translate_to_cone(i,L))
        L = [vector(l) for l in L]
        return tuple(vector(self.greatest_tuple([tuple(t + l) for l in L])) - t)
              
    def is_in_cone(self,i,v): 
        r"""
        Test whether the tuple `t` is contained in the `i`-th cone or not

        INPUTS:
            - ``i`` -- a cone index
            - ``t`` -- a tuple

        EXAMPLES::
            
            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(2)
            sage: G.is_in_cone(0, (1,2)) 
            True
            sage: G.is_in_cone(0,(-2,3))
            False
        """
        return all(c >=0 for c in self._cones[i]*vector(v))

    def generator(self,i,L):
        r"""
        Return the generator of the module over the `i`-th cone of elements `t \in \ZZ^n`
        such that the greatest_tuple of `t + L` is contained in the `i`-th cone.

        INPUTS:
            
            - ``i`` -- a cone index
            - ``L`` -- a list of tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_monomial_order import GeneralizedMonomialOrder
            sage: G = GeneralizedMonomialOrder(3)
            sage: L = [(1,2,-2), (3,4,-5), (-6,2,-7), (3,-6,1)]
            sage: G.generator(0,L)
            (6, 6, 7)

            sage: G.generator(2,L)
            (5, 5, 6)
        """

        cone_matrix = self.get_cone(i)
        t = vector(self.translate_to_cone(i,L))
        L = [vector(l) for l in L]
        for c in cone_matrix.columns():
            while self.is_in_cone(i,self.greatest_tuple([tuple(t + l) for l in L])):
                t = t - c
            t = t + c
        return tuple(t)

    def generator_for_pair(self,i,L1,L2):
        cone_matrix = self.get_cone(i)
        lm1 = vector(self.greatest_tuple_for_cone(i,L1))
        lm2 = vector(self.greatest_tuple_for_cone(i,L2))
        g1 = vector(self.generator(i,L1))
        g2 = vector(self.generator(i,L2))
        v1 = lm1 + g1
        v2 = lm2 + g2
        # print(lm1 + g1)
        # print(lm2 + g2)
        # print(cone_matrix)
        m = vector([max(a,b) for a,b in zip(v1,v2)])
        # m2 = max(lm1 + g1, lm2 + g2)
        # print(m)
        # print(m2)
        return tuple(cone_matrix*m)

def min_score_function(v):
    return -min(0,*v)



