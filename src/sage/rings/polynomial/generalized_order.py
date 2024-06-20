r"""
Generalized monomial orders on `\ZZ^n` are total orders on `\ZZ` introduced in
[PU1999]_ and defined by a triple:

- Conic decomposition: a finite covering of the additive monoïd `\ZZ^n` by submonoids.
  In the current implementation, the conic decomposition of Example 2.2 in 
  [PU1999]_ is imposed.

- Group order: this can be any class with a 'greatest_tuple' method returning
  the greatest of two tuples based on a group order on `\ZZ^n`. 
  Available choices include 'lex'.

- Score function: a function from the set of integer tuples to the 
  positive rationals, satisfying additional conditions (see Lemma 2.1 in [PU1999]_).
  Choices for this function include 'min' and 'degmin'.

When creating an instance of this class, you must specify at least the number of variables.
Optionally, you can also provide the names of the group order and the score function. 
Defaults are 'lex' and 'min' respectively.

These orders are useful for defining Gröbner bases for ideals in Laurent polynomial rings.

EXAMPLES::

    sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder

    sage: order = GeneralizedOrder(2); order
    Generalized monomial order in 2 variables using (lex, min)
    sage: order.n_cones()
    3
    sage: order.cone(1)
    [-1  0]
    [-1  1]
    sage: order.greatest_tuple((1,2),(-2,3))
    (-2, 3)
    sage: order.greatest_tuple((1,2),(-3,2),(2,-1))
    (-3, 2)
    sage: order.greatest_tuple_for_cone(2, (-1,2), (-1,2), (-2,-3))
    (-2, -3)

Using the ``degmin`` score function::

    sage: order2 = GeneralizedOrder(3, score_function="degmin"); order2
    Generalized monomial order in 3 variables using (lex, degmin)
    sage: order2.greatest_tuple((1, -1, 1), (2, -3, 4), (5, -6, 2), (-1, -2, -3))
    (5, -6, 2)

AUTHORS:

- Legrand Lucas (2024-26-01): initial version

"""
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.rings.polynomial.term_order import TermOrder
from sage.matrix.constructor import matrix
from sage.structure.sage_object import SageObject

class GeneralizedOrder(SageObject):
    def __init__(self,n,group_order="lex",score_function="min"):
        r"""
        Create a generalized monomial order in ``n`` varaibles with optional group order and score function.

        INPUT:

        - ``n`` -- the number of variables
        - ``group_order`` (default: ``lex``) -- the name of a group order on `\ZZ^n`, choices are: "lex"
        - ``score_function`` (default: ``min``) -- the name of a score function, choices are: "min", "degmin"

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: GeneralizedOrder(2)
            Generalized monomial order in 2 variables using (lex, min)

            sage: GeneralizedOrder(8, group_order="lex")
            Generalized monomial order in 8 variables using (lex, min)

            sage: GeneralizedOrder(3, score_function="min")
            Generalized monomial order in 3 variables using (lex, min)

            sage: GeneralizedOrder(3, group_order="lex", score_function="degmin")
            Generalized monomial order in 3 variables using (lex, degmin)
        """
        self._n = n
        self._n_cones = self._n + 1
        # Build cones        
        self._cones = [matrix.identity(ZZ,n)]
        for i in range(0,n):
            mat = matrix.identity(ZZ,n)
            mat.set_column(i,vector(ZZ,[-1]*n))
            self._cones.append(mat)
        # Set group order. Add more here.
        if group_order in ["lex"]:
            self._group_order = TermOrder(name=group_order)
        else:
            raise ValueError("Available group order are: 'lex'")
        # Set score function. Add more here.
        if score_function == "min":
            self._score_function = min_score_function
        elif score_function == "degmin":
            self._score_function = degmin_score_function
        else:
            raise ValueError("Available score function are: 'min', 'degmin'")
        # Store names
        self._group_order_name = group_order
        self._score_function_name = score_function

    def __hash__(self):
        r"""
        Return the hash of self. It depends on the number of variables, the group_order and the score function
        """
        return hash(self._group_order_name + self._score_function_name + str(self._n))

    def _repr_(self):
        r"""
        TESTS::
            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: GeneralizedOrder(2)._repr_()
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

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2)
            sage: order.n_cones()
            3
        """
        return self._n_cones

    def cones(self):
        r"""
        Return the list of matrices containing the generators of the cones.

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2)
            sage: order.cones()
            [
            [1 0]  [-1  0]  [ 1 -1]
            [0 1], [-1  1], [ 0 -1]
            ]
        """
        return self._cones

    def cone(self,i):
        r"""
        Return the matrix whose columns are the generators of the ``i``-th cone.
    
        INPUT:
        
        - `i` -- an integer, a cone index
        
        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2)
            sage: order.cone(1)
            [-1 0]
            [-1 1]
        """
        if i < 0 or i > self._n:
            raise IndexError("cone index out of range")
        return self._cones[i]

    def group_order(self):
        r"""
        Return the underlying group order.
        
        EXAMPLES::
           
            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: G = GeneralizedOrder(2)
            sage: G.group_order()
            Lexicographic term order
        """
        return self._group_order

    def group_order_name(self):
        r"""
        Return the name of the underlying group order.

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2,group_order="lex")
            sage: order.group_order_name()
            'lex'
        """
        return self._group_order_name
    
    def score_function_name(self):
        r"""
        Return the name of the underlying score function.

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2,score_function="degmin")
            sage: order.score_function_name()
            'degmin'
        """
        return self._score_function_name

    def score(self,t):
        r"""
        Compute the score of a tuple ``t``.

        INPUT:

        - `t` -- a tuple of integers

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2)
            sage: order.score((-2,3))
            2
        """
        return self._score_function(t)
                                
    def compare(self,a,b):
        r"""
        Return 1, 0 or -1 whether tuple ``a`` is greater than, equal or less than to tuple ``b`` respectively.

        INPUT:

        - `a` and `b` -- two tuples of integers

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(3)
            sage: order.compare((1,2,-2), (3,-4,5))
            -1
            sage: order.compare((3,-4,5), (1,2,-2))
            1
            sage: order.compare((3,-4,5), (3,-4,5))
            0
        """
        if a == b: return 0
        diff = self._score_function(a) - self._score_function(b)
        if diff != 0: 
            return 1 if diff > 0 else -1
        else: 
            return 1 if self._group_order.greater_tuple(a,b) == a else -1

    def greatest_tuple(self,*L):
        r"""
        Return the greatest tuple in ``L``.

        INPUT:

        - *L -- integer tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(3)
            sage: order.greatest_tuple((0,0,1),(2,3,-2))
            (2, 3, -2)
            sage: L = [(1,2,-1),(3,-3,0),(4,-5,-6)]
            sage: order.greatest_tuple(*L)
            (4, -5, -6)
        """
        n = len(L) 
        if n == 0:
           raise ValueError("empty list of tuples")

        if n == 1:
            return L[0]
        else:
            a = L[0]
            b = self.greatest_tuple(*L[1:])
            return a if self.compare(a, b) == 1 else b

    def translate_to_cone(self,i,L):
        r"""
        Return a tuple ``t`` such that `t + L` is contained in the ``i``-th cone.

        INPUT:

        - ``i`` -- an integer, a cone index
        - ``L`` -- a list of integer tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2)
            sage: order.translate_to_cone(0, [(1,2),(-2,-3),(1,-4)])
            (2, 4)
        """
        cone_matrix = self._cones[i]
        T = matrix(ZZ, [cone_matrix*vector(v) for v in L])
        return tuple(cone_matrix*vector([-min(0,*c) for c in T.columns()]))
        
    def greatest_tuple_for_cone(self, i, *L):
        r"""
        Return the greatest tuple of the list of tuple `L` with respect to the `i`-th cone.

        This is the unique tuple `t` in `L`  such that each time the greatest tuple of `s + L`
        for a tuple `s` is contained in the `i`-th cone, it is equal to `s + t`.

        INPUT:

        - ``i`` -- an integer, a cone index
        - ``L`` -- a list of integer tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2)
            sage: L = [(1,2), (-2,-2), (-4,5), (5,-6)]
            sage: t = order.greatest_tuple_for_cone(1, *L);t
            (-4, 5)

        We can check the result::

            sage: s = order.translate_to_cone(1, L)
            sage: sL = [tuple(vector(s) + vector(l)) for l in L]
            sage: tuple(vector(s) + vector(t)) == order.greatest_tuple(*sL)
            True

        """
        t = vector(self.translate_to_cone(i,L))
        L = [vector(l) for l in L]
        return tuple(vector(self.greatest_tuple(*[tuple(t + l) for l in L])) - t)
              
    def is_in_cone(self,i,t): 
        r"""
        Test whether the tuple ``t`` is contained in the ``i``-th cone or not.

        INPUT:

        - ``i`` -- an integer, a cone index
        - ``t`` -- a tuple of integers

        EXAMPLES::
            
            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2)
            sage: order.is_in_cone(0, (1,2)) 
            True
            sage: order.is_in_cone(0,(-2,3))
            False
        """
        return all(c >=0 for c in self.cone(i)*vector(t))

    def generator(self,i,L):
        r"""
        Return the generator of the module over the ``i``-th cone for ``L``.

        This is the monoïd of elements `t \in \ZZ^n` such that the greatest 
        tuple of `t + L` is contained in the ``i``-th cone.

        INPUT:
            
        - ``i`` -- an integer, a cone index
        - ``L`` -- a list of iinteger tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(3)
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
            while self.is_in_cone(i,self.greatest_tuple(*[tuple(t + l) for l in L])):
                t = t - c
            t = t + c
        return tuple(t)

    def generator_for_pair(self,i,L1,L2):
        r"""
        Return the generator of the module over the ``i``-th cone for ``L1`` and ``L2``. 

        INPUT:

        - ``i`` -- an integer, a cone index
        - ``L1`` -- a list of integer tuples
        - ``L2`` -- a list of integer tuples

        EXAMPLES::

            sage: from sage.rings.polynomial.generalized_order import GeneralizedOrder
            sage: order = GeneralizedOrder(2)
            sage: L1 = [(2,3),(-4,2),(-1,-2)]
            sage: L2 = [(1,-6),(5,-2),(3,4)]
            sage: order.generator_for_pair(1,L1,L2)
            (1, 5)

        """
        cone_matrix = self._cones[i]
        lm1 = vector(self.greatest_tuple_for_cone(i,*L1))
        lm2 = vector(self.greatest_tuple_for_cone(i,*L2))
        g1 = vector(self.generator(i,L1))
        g2 = vector(self.generator(i,L2))
        m = vector([max(a,b) for a,b in zip(lm1 + g1, lm2 + g2)])
        return tuple(cone_matrix*m)
   
# Score functions. Add more here and update the __init__ method.
def min_score_function(t):
        return -min(0,*t)

def degmin_score_function(t):
        return sum(t) -(len(t)+1)*min(0,*t)



