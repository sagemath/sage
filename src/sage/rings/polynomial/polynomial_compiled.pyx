# sage_setup: distribution = sagemath-categories
"""
Polynomial Compilers

AUTHORS:

- Tom Boothby, initial design & implementation
- Robert Bradshaw, bug fixes / suggested & assisted with significant design improvements
"""

################################################################################
#       Copyright (C) 2007 Tom Boothby <boothby@u.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
################################################################################

from sage.misc.binary_tree cimport BinaryTree


cdef class CompiledPolynomialFunction:
    """
    Build a reasonably optimized directed acyclic graph representation
    for a given polynomial.  A ``CompiledPolynomialFunction`` is callable from
    python, though it is a little faster to call the eval function from
    pyrex.

    This class is not intended to be called by a user, rather, it is
    intended to improve the performance of immutable polynomial objects.

    .. TODO::

        - Recursive calling
        - Faster casting of coefficients / argument
        - Multivariate polynomials
        - Cython implementation of Pippenger's Algorithm that doesn't
          depend heavily upon dicts.
        - Computation of parameter sequence suggested by Pippenger
        - Univariate exponentiation can use Brauer's method to improve
          extremely sparse polynomials of very high degree
    """

    def __init__(self, coeffs, algorithm='binary'):
        """
        Compiles a polynomial into an evaluation DAG representation which
        is at least as optimal as using Horner's Rule.  For polynomials
        which have a relatively large number of zero coefficients, the
        improvement over Horner's Rule grows significantly.

        Here is a rough description of the algorithm.  Actual
        implementation differs slightly, for sake of speed.  Specifically,
        steps 1 and 3 are done at the same time.

        Step 1: Collect Coefficient Gaps.
                Scan through coefficient list, and record the lengths of
                sequences of zero coefficients.  This corresponds to
                collapsing Horner's Form into a reduced representation.
                For example, ::

                  x^8 + x^4 + x^2 + 1
                   = ((((((((1)*x + 0)*x+0)*x+0)*x+1)*x+0)*x+0)*x+1)*x+0)*x+1
                   = ((((1)*x^4 + 1)*x^2 + 1)*x^2 + 1

                gives a list of "gaps": [2,4]

        Step 2: Fill in Gap Structure.
                Given the list of gaps, find a reasonable sequence of
                of multiplications / squarings of x that will result in
                the computation of all desired exponents.  Record this
                sequence of steps as an evaluation DAG, and retain
                references to the nodes representing the desired
                exponents.  For the example above, we would have::

                  x^2 = x*x
                  x^4 = (x^2) * (x^2)

        Step 3: Construct Evaluation Dag.
                Rescan the coefficient list, and build an evaluation DAG
                representation for the reduced Horner Form as above.
                Whenever a sequence of zeros is encountered, multiply by
                the appropriate "gap" node.  Retain a reference to the
                result node.

        Implementation considerations:

        * By combining steps 1 and 3, we greatly improve the speed of
          this construction, but some complexity is introduced.  The
          solution to this is that "dummy" nodes are created to represent
          the desired gaps.  As the structure of the gaps is filled in,
          these dummies get references to usable DAG nodes.  After all gaps
          are filled in, we strip out dummy nodes, and are left with a
          complete representation of our polynomial.

        * The "binary" algorithm (currently the only algorithm; others are
          forthcoming) requires the gaps to considered in order, and adds
          additional dummies as it goes.  Hence, the gaps are put into a
          binary tree.
        """
        cdef generic_pd max_gap, dag
        cdef BinaryTree gaps

        self._coeffs = coeffs
        gaps, dag = self._parse_structure()
        max_gap = <generic_pd>(gaps.get_max())
        if max_gap.label > 1:
            if algorithm == 'binary':
                self._fill_gaps_binary(gaps)
            elif algorithm == 'pippenger':
                raise NotImplementedError("Implementation of Pippenger's Algorithm is not ready for prime time.")
            else:
                raise RuntimeError("Method '%s' not supported.")

        self._dag = dag.nodummies()

    def __repr__(self):
        return "CompiledPolynomialFunction(%s)" % self._dag

    def __call__(self, x):
        return self.eval(x)

    cdef object eval(CompiledPolynomialFunction self, object x):
        cdef object temp
        try:
            pd_eval(self._dag, x, self._coeffs)  #see further down
            temp = self._dag.value               #for an explanation
            pd_clean(self._dag)                  #of these 3 lines
            return temp
        except TypeError as msg:
            self._dag.reset()
            raise TypeError(msg)

    cdef object _parse_structure(CompiledPolynomialFunction self):
        """
        Loop through the coefficients of the polynomial, and collect
        coefficient gap widths.  Meanwhile, construct the evaluation
        DAG; inserting dummy nodes wherever there are gaps of width
        greater than 1.

        Return the resultant head DAG node, and a binary tree
        containing the dummy nodes.
        """
        cdef BinaryTree gaps
        cdef int d
        cdef generic_pd s

        s = univar_pd()
        gaps = BinaryTree()
        gaps.insert(1, s)

        d = len(self._coeffs)-1

        s = coeff_pd(d)
        gap_width = 0

        d -= 1
        while d > 0:
            gap_width += 1
            if self._coeffs[d]:
                s = abc_pd(s, self._get_gap(gaps, gap_width), d)
                gap_width = 0
            d-=1

        gap_width += 1
        if self._coeffs[0]:
            s = abc_pd(s, self._get_gap(gaps, gap_width), 0)
        else:
            s = mul_pd(s, self._get_gap(gaps, gap_width))

        return gaps, s

    cdef generic_pd _get_gap(CompiledPolynomialFunction self, BinaryTree gaps, int gap):
        """
        Find an entry in the BinaryTree gaps, identified by the int gap.
        If such an entry does not exist, create it and put it in the tree.
        """

        cdef generic_pd g
        g = gaps.get(gap)
        if g is not None:
            return g
        else:
            g = dummy_pd(gap)
            gaps.insert(gap,g)
            return g

    cdef void _fill_gaps_binary(CompiledPolynomialFunction self, BinaryTree gaps) noexcept:
        """
        The desired gaps come in a tree, filled with dummy nodes (with the
        exception of the var node, which is not a dummy).  The nodes are
        labeled with their width.  When we refer to multiplication or
        squaring nodes, we give the current dummy node a usable dag of the
        appropriate type, which treats the nodes being operated on as
        argument.  That is, if we want to multiply nodes A and B, and
        give the result to M, we would call ::

            M.fill(mul_pd(A,B))

        and refer to this operation as ::

            M = A*B

        Sometimes we want a node that isn't in the tree.  In that case, we
        create a new dummy node, and put it into the tree.  The phrase
        "get the node for..." means that we obtain one with the _get_gap
        method, which may create a new node.


        Fill in the gaps with the following algorithm:

        | Step 1:
        |     Remove max node from the tree, denote its width m.
        |     If m == 1, halt.
        |
        | Step 2:
        |     If m is even, check if m/2 is already in our tree.
        |     If yes, square the node corresponding to m/2, and go to step 1.
        |
        | Step 3:
        |     Peek at the next-largest, and denote its width n.
        |     Write m = qn+r.
        |     If r != 0, get the node R, whose gap-width is r.
        |         Also, get the node QN, whose width is qn.
        |         Then, set M = R*QN and go to step 1.
        |     If r == 0, we have two cases:
        |         q is even:
        |             get/create the node for q/n, denote it T.
        |             Set M = T**2 and go to step 1.
        |         q is odd:
        |             get/create the node for (Q-1)*N, denote it T.
        |             Set M = T*N and go to step 1.

        The r == 0 case in step 3 is equivalent to binary exponentiation.
        """
        cdef int m,n,k,r,half
        cdef generic_pd T,N,H
        cdef dummy_pd M
        T = gaps.pop_max()
        while T is not None:
            N = gaps.get_max()
            if N is None:
                return
            else:
                M = <dummy_pd>T
            m = M.label
            n = N.label
            k = m/n
            r = m%n
            half = m/2

            found = False
            if m % 2 == 0 and n >= half:
                H = gaps.get(half)
                if H is not None:
                    M.fill(sqr_pd(H))
                    found = True

            if not found:
                if r > 0:
                    M.fill(mul_pd(self._get_gap(gaps, n*k), self._get_gap(gaps, r)))
                elif k % 2 == 0:
                    M.fill(sqr_pd(self._get_gap(gaps, n*k/2)))
                else:
                    M.fill(mul_pd(self._get_gap(gaps, n*(k-1)), N))

            T = gaps.pop_max()


########################################################
#
#  Polynomial DAG Code
#
# An "evaluation DAG" is a rather generic concept.  What
# follows is a set of DAG classes which are suitable for
# evaluating polynomials in Horner's form, and very
# little else.  Let's call them polydags from here down.
# DAG in all-caps is a whole graph. However,
# polydag or pd := polynomial DAG node
#
# 5 operations, and 3 types of constants are implemented:
#
#  univar_pd: a variable of a univariate polynomial.
#  var_pd: a variable of a multivariate polynomial.
#          parametrized by an index number.
#  const_pd: a coefficient of a polynomial.  Also
#            parametrized by an index number.
#
# unary operations:
#
#  sqr_pd: takes another polydag as an argument, squares
#          by multiplying argument's value by itself.
#  pow_pd: takes another polydag, and an int, as arguments.
#          exponentiates with the ** operator.
#
# binary operations:
#
#  add_pd: adds values of two other polydags.
#  mul_pd: like add_pd, but we with multiplication
#  abc_pd: not technically binary, but it only depends
#          upon two other polydags, so it fits into the
#          class structure this way.  Multiplies the
#          values of two other polydags, and adds the
#          coefficient corresponding to the given index.
#
# Polydag structure:
#  properties:
#    value: some python object.  Should be None before
#           and after the DAG is evaluated.
#    refs: number of polydag nodes that reference self
#    hits: number of times, during current computation,
#          that self.value has been accessed
#    label: used to identify nodes from the outside.
#           specifically, when they're in a binary tree.
#
#  methods:
#    eval: runs the computation for which the node is
#          intended, and stores the result in self.value.
#          Takes 2 python objects as arguments.  The first
#          is the value of the variable, for univariate
#          polynomials, or a tuple of values for multi-
#          variate.  The second is the list of coefficients.
#          The reason for taking the list of coefficients,
#          rather than holding references to them, is that
#          we can implement recursive calling at some later
#          date.
#    nodummies: recursively evict dummies, replacing them
#               with the non-dummy nodes that they
#               reference.


# These inline functions are called wherever a node gets
# evaluated.  First, pd_eval is called to ensure that the
# target DAG node will have its .value property set. It
# serves two purposes:
#   1) Count the number of times this node has been asked
#      for its value.
#   2) Call the node's eval() function if and only if it
#      has not been called yet.  These functions are
#      inlined to prevent excessive function calls.
# Then, pd_clean is called.  It checks to see if all of
# the polydag's dependents have phoned in, and resets
# self.value to None if they have.  Thus, we don't hold
# on to intermediate values any longer than we have to.

cdef inline int pd_eval(generic_pd pd, object vars, object coeffs) except -2:
    if pd.value is None:
        pd.eval(vars, coeffs)
    pd.hits += 1

cdef inline void pd_clean(generic_pd pd) noexcept:
    if pd.hits >= pd.refs:
        pd.value = None
        pd.hits = 0

cdef class generic_pd:
    def __init__(generic_pd self):
        self.value = None
        self.hits = 0
        self.refs = 0
        self.label = -1

    cdef int eval(generic_pd self, object vars, object coeffs) except -2:
        raise NotImplementedError

    cdef generic_pd nodummies(generic_pd self):
        return self

    cdef void reset(generic_pd self) noexcept:
        self.hits = 0
        self.value = None

cdef class dummy_pd(generic_pd):
    def __init__(dummy_pd self, int label):
        self.label = label

    cdef void fill(dummy_pd self, generic_pd link) noexcept:
        self.link = link

    cdef generic_pd nodummies(dummy_pd self):
        #sorry guys, this is my stop
        self.link.refs = self.refs
        return self.link.nodummies()

cdef class var_pd(generic_pd):
    def __init__(var_pd self, int index):
        generic_pd.__init__(self)
        self.index = index
    cdef int eval(var_pd self, object vars, object coeffs) except -2:
        self.value = vars[self.index]

    def __repr__(self):
        return "x[%s]" % (self.index)


cdef class univar_pd(generic_pd):
    def __init__(univar_pd self):
        generic_pd.__init__(self)
        self.label = 1
    cdef int eval(univar_pd self, object var, object coeffs) except -2:
        self.value = var

    def __repr__(self):
        return "x"

cdef class coeff_pd(generic_pd):
    def __init__(coeff_pd self, int index):
        generic_pd.__init__(self)
        self.index = index
    cdef int eval(coeff_pd self, object vars, object coeffs) except -2:
        self.value = coeffs[self.index]

    def __repr__(self):
        return "a%s" % (self.index)

    cdef void reset(self) noexcept:
        pass

cdef class unary_pd(generic_pd):
    def __init__(unary_pd self, generic_pd operand):
        generic_pd.__init__(self)
        self.operand = operand
        self.operand.refs += 1

    cdef generic_pd nodummies(self):
        self.operand = self.operand.nodummies()
        return self

    cdef void reset(self) noexcept:
        generic_pd.reset(self)
        self.operand.reset()


cdef class sqr_pd(unary_pd):
    cdef int eval(sqr_pd self, object vars, object coeffs) except -2:
        pd_eval(self.operand, vars, coeffs)
        self.value = self.operand.value * self.operand.value
        pd_clean(self.operand)

    def __repr__(self):
        return "(%s)^2" % (self.operand)

cdef class pow_pd(unary_pd):
    def __init__(unary_pd self, generic_pd base, object exponent):
        unary_pd.__init__(self, base)
        self.exponent = exponent

    cdef int eval(pow_pd self, object vars, object coeffs) except -2:
        pd_eval(self.operand, vars, coeffs)
        self.value = self.operand.value ** self.exponent
        pd_clean(self.operand)

    def __repr__(self):
        return "(%s^%s)" % (self.left, self.exponent)


cdef class binary_pd(generic_pd):
    def __init__(binary_pd self, generic_pd left, generic_pd right):
        generic_pd.__init__(self)
        self.left = left
        self.right = right
        self.left.refs+= 1
        self.right.refs+= 1

    cdef generic_pd nodummies(self):
        self.left = self.left.nodummies()
        self.right = self.right.nodummies()
        return self

    cdef void reset(self) noexcept:
        generic_pd.reset(self)
        self.left.reset()
        self.right.reset()

cdef class add_pd(binary_pd):
    cdef int eval(add_pd self, object vars, object coeffs) except -2:
        pd_eval(self.left, vars, coeffs)
        pd_eval(self.right, vars, coeffs)
        self.value = self.left.value + self.right.value
        pd_clean(self.left)
        pd_clean(self.right)

    def __repr__(self):
        return "(%s+%s)" % (self.left, self.right)

cdef class mul_pd(binary_pd):
    cdef int eval(mul_pd self, object vars, object coeffs) except -2:
        pd_eval(self.left, vars, coeffs)
        pd_eval(self.right, vars, coeffs)
        self.value = self.left.value * self.right.value
        pd_clean(self.left)
        pd_clean(self.right)

    def __repr__(self):
        return "(%s*%s)" % (self.left, self.right)

cdef class abc_pd(binary_pd):
    def __init__(abc_pd self, generic_pd left, generic_pd right, int index):
        binary_pd.__init__(self, left, right)
        self.index = index

    def __repr__(self):
        return "(%s*%s+a%s)" % (self.left, self.right, self.index)

    cdef int eval(abc_pd self, object vars, object coeffs) except -2:
        pd_eval(self.left, vars, coeffs)
        pd_eval(self.right, vars, coeffs)
        self.value = self.left.value * self.right.value + coeffs[self.index]
        pd_clean(self.left)
        pd_clean(self.right)
