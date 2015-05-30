r"""
ClusterSeed

A *cluster seed* is a pair `(B,\mathbf{x})` with `B` being a *skew-symmetrizable* `(n+m \times n)` *-matrix*
and with `\mathbf{x}` being an `n`-tuple of *independent elements* in the field of rational functions in `n` variables.

For the compendium on the cluster algebra and quiver package see
:arxiv:`1102.4844`.

AUTHORS:

- Gregg Musiker
- Christian Stump

.. seealso:: For mutation types of cluster seeds, see :meth:`sage.combinat.cluster_algebra_quiver.quiver_mutation_type.QuiverMutationType`. Cluster seeds are closely related to :meth:`sage.combinat.cluster_algebra_quiver.quiver.ClusterQuiver`.
"""

#*****************************************************************************
#       Copyright (C) 2011 Gregg Musiker <musiker@math.mit.edu>
#                          Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import time
from sage.structure.sage_object import SageObject
from copy import copy
from sage.rings.all import QQ, infinity
from sage.rings.integer_ring import ZZ
from sage.rings.all import FractionField, PolynomialRing
from sage.rings.fraction_field_element import FractionFieldElement
from sage.sets.all import Set
from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import  QuiverMutationType_Irreducible, QuiverMutationType_Reducible
from sage.combinat.cluster_algebra_quiver.mutation_type import is_mutation_finite

class ClusterSeed(SageObject):
    r"""
    The *cluster seed* associated to an *exchange matrix*.

    INPUT:

    - ``data`` -- can be any of the following::

        * QuiverMutationType
        * str - a string representing a QuiverMutationType or a common quiver type (see Examples)
        * ClusterQuiver
        * Matrix - a skew-symmetrizable matrix
        * DiGraph - must be the input data for a quiver
        * List of edges - must be the edge list of a digraph for a quiver

    EXAMPLES::

        sage: S = ClusterSeed(['A',5]); S
        A seed for a cluster algebra of rank 5 of type ['A', 5]

        sage: S = ClusterSeed(['A',[2,5],1]); S
        A seed for a cluster algebra of rank 7 of type ['A', [2, 5], 1]

        sage: T = ClusterSeed( S ); T
        A seed for a cluster algebra of rank 7 of type ['A', [2, 5], 1]

        sage: T = ClusterSeed( S._M ); T
        A seed for a cluster algebra of rank 7

        sage: T = ClusterSeed( S.quiver()._digraph ); T
        A seed for a cluster algebra of rank 7

        sage: T = ClusterSeed( S.quiver()._digraph.edges() ); T
        A seed for a cluster algebra of rank 7

        sage: S = ClusterSeed(['B',2]); S
        A seed for a cluster algebra of rank 2 of type ['B', 2]

        sage: S = ClusterSeed(['C',2]); S
        A seed for a cluster algebra of rank 2 of type ['B', 2]

        sage: S = ClusterSeed(['A', [5,0],1]); S
        A seed for a cluster algebra of rank 5 of type ['D', 5]

        sage: S = ClusterSeed(['GR',[3,7]]); S
        A seed for a cluster algebra of rank 6 of type ['E', 6]

        sage: S = ClusterSeed(['F', 4, [2,1]]); S
        A seed for a cluster algebra of rank 6 of type ['F', 4, [1, 2]]
    """
    def __init__(self, data, frozen=None, is_principal=None):
        r"""
        TESTS::

            sage: S = ClusterSeed(['A',4])
            sage: TestSuite(S).run()
        """
        from quiver import ClusterQuiver

        # constructs a cluster seed from a cluster seed
        if isinstance(data, ClusterSeed):
            if frozen:
                print "The input \'frozen\' is ignored"
            self._M = copy( data._M )
            self._cluster = copy(data._cluster)
            self._n = data._n
            self._m = data._m
            self._R = data._R
            self._quiver = ClusterQuiver( data._quiver ) if data._quiver else None
            self._mutation_type = copy( data._mutation_type )
            self._description = copy( data._description )
            self._is_principal = data._is_principal

        # constructs a cluster seed from a quiver
        elif isinstance(data, ClusterQuiver):
            if frozen:
                print "The input \'frozen\' is ignored"

            quiver = ClusterQuiver( data )
            self._M = copy(quiver._M)
            self._n = quiver._n
            self._m = quiver._m
            self._quiver = quiver
            self._mutation_type = quiver._mutation_type
            self._description = 'A seed for a cluster algebra of rank %d' %(self._n)
            self._R = FractionField(PolynomialRing(QQ,['x%s'%i for i in range(0,self._n)]+['y%s'%i for i in range(0,self._m)]))
            self._cluster = list(self._R.gens())
            self._is_principal = None

        # in all other cases, we construct the corresponding ClusterQuiver first
        else:
            quiver = ClusterQuiver( data, frozen=frozen )
            self.__init__( quiver )

        if is_principal is not None:
            self._is_principal = is_principal

    def __eq__(self, other):
        r"""
        Returns True iff ``self`` represent the same cluster seed as ``other``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',5])
            sage: T = S.mutate( 2, inplace=False )
            sage: S.__eq__( T )
            False

            sage: T.mutate( 2 )
            sage: S.__eq__( T )
            True
        """
        return isinstance(other, ClusterSeed) and self._M == other._M and self._cluster == other._cluster

    def _repr_(self):
        r"""
        Returns the description of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',5])
            sage: S._repr_()
            "A seed for a cluster algebra of rank 5 of type ['A', 5]"

            sage: S=ClusterSeed(['B',2])
            sage: T=S.principal_extension()
            sage: T._repr_()
            "A seed for a cluster algebra of rank 2 of type ['B', 2] with principal coefficients"
        """
        name = self._description
        if self._mutation_type:
            if type( self._mutation_type ) in [QuiverMutationType_Irreducible,QuiverMutationType_Reducible]:
                name += ' of type ' + str(self._mutation_type)
            # the following case allows description of 'undetermined finite mutation type'
            else:
                name += ' of ' + self._mutation_type
        if self._is_principal:
            name += ' with principal coefficients'
        elif self._m == 1:
            name += ' with %s frozen variable'%self._m
        elif self._m > 1:
            name += ' with %s frozen variables'%self._m
        return name

    def plot(self, circular=False, mark=None, save_pos=False):
        r"""
        Returns the plot of the quiver of ``self``.

        INPUT:

        - ``circular`` -- (default:False) if True, the circular plot is chosen, otherwise >>spring<< is used.
        - ``mark`` -- (default: None) if set to i, the vertex i is highlighted.
        - ``save_pos`` -- (default:False) if True, the positions of the vertices are saved.

        EXAMPLES::

            sage: S = ClusterSeed(['A',5])
            sage: pl = S.plot()
            sage: pl = S.plot(circular=True)
        """
        return self.quiver().plot(circular=circular,mark=mark,save_pos=save_pos)

    def show(self, fig_size=1, circular=False, mark=None, save_pos=False):
        r"""
        Shows the plot of the quiver of ``self``.

        INPUT:

        - ``fig_size`` -- (default: 1) factor by which the size of the plot is multiplied.
        - ``circular`` -- (default: False) if True, the circular plot is chosen, otherwise >>spring<< is used.
        - ``mark`` -- (default: None) if set to i, the vertex i is highlighted.
        - ``save_pos`` -- (default:False) if True, the positions of the vertices are saved.

        TESTS::

            sage: S = ClusterSeed(['A',5])
            sage: S.show() # long time
        """
        self.quiver().show(fig_size=fig_size, circular=circular,mark=mark,save_pos=save_pos)

    def interact(self, fig_size=1, circular=True):
        r"""
        Only in *notebook mode*. Starts an interactive window for cluster seed mutations.

        INPUT:

        - ``fig_size`` -- (default: 1) factor by which the size of the plot is multiplied.
        - ``circular`` -- (default: True) if True, the circular plot is chosen, otherwise >>spring<< is used.

        TESTS::

            sage: S = ClusterSeed(['A',4])
            sage: S.interact() # long time
            'The interactive mode only runs in the Sage notebook.'
        """
        from sage.plot.plot import EMBEDDED_MODE
        from sagenb.notebook.interact import interact, selector
        from sage.misc.all import html,latex

        if not EMBEDDED_MODE:
            return "The interactive mode only runs in the Sage notebook."
        else:
            seq = []
            sft = [True]
            sss = [True]
            ssv = [True]
            ssm = [True]
            ssl = [True]
            @interact
            def player(k=selector(values=range(self._n),nrows = 1,label='Mutate at: '), show_seq=("Mutation sequence:", True), show_vars=("Cluster variables:", True), show_matrix=("B-Matrix:", True), show_lastmutation=("Show last mutation:", True) ):
                ft,ss,sv,sm,sl = sft.pop(), sss.pop(), ssv.pop(), ssm.pop(), ssl.pop()
                if ft:
                    self.show(fig_size=fig_size, circular=circular)
                elif show_seq is not ss or show_vars is not sv or show_matrix is not sm or show_lastmutation is not sl:
                    if seq and show_lastmutation:
                        self.show(fig_size=fig_size, circular=circular, mark=seq[len(seq)-1])
                    else:
                        self.show(fig_size=fig_size, circular=circular )
                else:
                    self.mutate(k)
                    seq.append(k)
                    if not show_lastmutation:
                        self.show(fig_size=fig_size, circular=circular)
                    else:
                        self.show(fig_size=fig_size, circular=circular,mark=k)
                sft.append(False)
                sss.append(show_seq)
                ssv.append(show_vars)
                ssm.append(show_matrix)
                ssl.append(show_lastmutation)
                if show_seq: html( "Mutation sequence: $" + str( [ seq[i] for i in xrange(len(seq)) ] ).strip('[]') + "$" )
                if show_vars:
                    html( "Cluster variables:" )
                    table = "$\\begin{align*}\n"
                    for i in xrange(self._n):
                        table += "\tv_{%s} &= "%i + latex( self._cluster[i] ) + "\\\\ \\\\\n"
                    table += "\\end{align*}$"
                    html( "$ $" )
                    html( table )
                    html( "$ $" )
                if show_matrix:
                    html( "B-Matrix:" )
                    m = self._M
                    #m = matrix(range(1,self._n+1),sparse=True).stack(m)
                    m = latex(m)
                    m = m.split('(')[1].split('\\right')[0]
                    html( "$ $" )
                    html( "$\\begin{align*} " + m + "\\end{align*}$" )
                    #html( "$" + m + "$" )
                    html( "$ $" )

    def save_image(self, filename, circular=False, mark=None, save_pos=False):
        r"""
        Saves the plot of the underlying digraph of the quiver of ``self``.

        INPUT:

        - ``filename`` -- the filename the image is saved to.
        - ``circular`` -- (default: False) if True, the circular plot is chosen, otherwise >>spring<< is used.
        - ``mark`` -- (default: None) if set to i, the vertex i is highlighted.
        - ``save_pos`` -- (default:False) if True, the positions of the vertices are saved.

        EXAMPLES::

            sage: S = ClusterSeed(['F',4,[1,2]])
            sage: S.save_image(os.path.join(SAGE_TMP, 'sage.png'))
        """
        graph_plot = self.plot( circular=circular, mark=mark, save_pos=save_pos)
        graph_plot.save( filename=filename )

    def b_matrix(self):
        r"""
        Returns the `B` *-matrix* of ``self``.

        EXAMPLES::

            sage: ClusterSeed(['A',4]).b_matrix()
            [ 0  1  0  0]
            [-1  0 -1  0]
            [ 0  1  0  1]
            [ 0  0 -1  0]

            sage: ClusterSeed(['B',4]).b_matrix()
            [ 0  1  0  0]
            [-1  0 -1  0]
            [ 0  1  0  1]
            [ 0  0 -2  0]

            sage: ClusterSeed(['D',4]).b_matrix()
            [ 0  1  0  0]
            [-1  0 -1 -1]
            [ 0  1  0  0]
            [ 0  1  0  0]

            sage: ClusterSeed(QuiverMutationType([['A',2],['B',2]])).b_matrix()
            [ 0  1  0  0]
            [-1  0  0  0]
            [ 0  0  0  1]
            [ 0  0 -2  0]
        """
        return copy( self._M )

    def ground_field(self):
        r"""
        Returns the *ground field* of the cluster of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.ground_field()
            Fraction Field of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        """
        return self._R

    def x(self,k):
        r"""
        Returns the `k` *-th initial cluster variable* for the associated cluster seed.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.mutate([2,1])
            sage: S.x(0)
            x0

            sage: S.x(1)
            x1

            sage: S.x(2)
            x2
        """
        if k in range(self._n):
            x = self._R.gens()[k]
            return ClusterVariable( x.parent(), x.numerator(), x.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable' )
        else:
            raise ValueError("The input is not in an index of a cluster variable.")

    def y(self,k):
        r"""
        Returns the `k` *-th initial coefficient (frozen variable)* for the associated cluster seed.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1])
            sage: S.y(0)
            y0

            sage: S.y(1)
            y1

            sage: S.y(2)
            y2
        """
        if k in range(self._m):
            x = self._R.gens()[self._n+k]
            return ClusterVariable( x.parent(), x.numerator(), x.denominator(), mutation_type=self._mutation_type, variable_type='frozen variable' )
        else:
            raise ValueError("The input is not in an index of a frozen variable.")

    def n(self):
        r"""
        Returns the number of *exchangeable variables* of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.n()
            3
        """
        return self._n

    def m(self):
        r"""
        Returns the number of *frozen variables* of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.n()
            3

            sage: S.m()
            0

            sage: S = S.principal_extension()
            sage: S.m()
            3
        """
        return self._m

    def cluster_variable(self,k):
        r"""
        Returns the `k`-th *cluster variable* of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.mutate([1,2])

            sage: [S.cluster_variable(k) for k in range(3)]
            [x0, (x0*x2 + 1)/x1, (x0*x2 + x1 + 1)/(x1*x2)]
        """
        if k not in range(self._n):
            raise ValueError("The cluster seed does not have a cluster variable of index %s."%k)
        f = self._cluster[k]
        return ClusterVariable( f.parent(), f.numerator(), f.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable' )

    def cluster(self):
        r"""
        Returns the *cluster* of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.cluster()
            [x0, x1, x2]

            sage: S.mutate(1)
            sage: S.cluster()
            [x0, (x0*x2 + 1)/x1, x2]

            sage: S.mutate(2)
            sage: S.cluster()
            [x0, (x0*x2 + 1)/x1, (x0*x2 + x1 + 1)/(x1*x2)]

            sage: S.mutate([2,1])
            sage: S.cluster()
            [x0, x1, x2]
        """
        return [ self.cluster_variable(k) for k in range(self._n) ]

    def f_polynomial(self,k,ignore_coefficients=False):
        r"""
        Returns the ``k``-th *F-polynomial* of ``self``. It is obtained from the
        ``k``-th cluster variable by setting all `x_i` to `1`.

        Requires principal coefficients, initialized by using principal_extension(),
        or the user can set 'ignore_coefficients=True' to bypass this restriction.

        Warning: this method assumes the sign-coherence conjecture and that the
        input seed is sign-coherent (has an exchange matrix with columns of like signs).
        Otherwise, computational errors might arise.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: [S.f_polynomial(k) for k in range(3)]
            [1, y1*y2 + y2 + 1, y1 + 1]

            sage: S = ClusterSeed(Matrix([[0,1],[-1,0],[1,0],[-1,1]])); S
            A seed for a cluster algebra of rank 2 with 2 frozen variables
            sage: T = ClusterSeed(Matrix([[0,1],[-1,0]])).principal_extension(); T
            A seed for a cluster algebra of rank 2 with principal coefficients
            sage: S.mutate(0)
            sage: T.mutate(0)
            sage: S.f_polynomials()
            Traceback (most recent call last):
            ...
            ValueError: No principal coefficients initialized. Use principal_extension, or ignore_coefficients to ignore this.
            sage: S.f_polynomials(ignore_coefficients=True)
            [y0 + y1, 1]
            sage: T.f_polynomials()
            [y0 + 1, 1]
        """
        if not (ignore_coefficients or self._is_principal):
            raise ValueError("No principal coefficients initialized. Use principal_extension, or ignore_coefficients to ignore this.")

        if k not in range(self._n):
            raise ValueError("The cluster seed does not have a cluster variable of index %s."%k)
        eval_dict = dict( [ ( self.x(i), 1 ) for i in range(self._n) ] )
        return self.cluster_variable(k).subs(eval_dict)

    def f_polynomials(self,ignore_coefficients=False):
        r"""
        Returns all *F-polynomials* of ``self``. These are obtained from the
        cluster variables by setting all `x_i`'s to `1`.

        Requires principal coefficients, initialized by using principal_extension(),
        or the user can set 'ignore_coefficients=True' to bypass this restriction.

        Warning: this method assumes the sign-coherence conjecture and that the
        input seed is sign-coherent (has an exchange matrix with columns of like signs).
        Otherwise, computational errors might arise.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: S.f_polynomials()
            [1, y1*y2 + y2 + 1, y1 + 1]
        """
        if not (ignore_coefficients or self._is_principal):
            raise ValueError("No principal coefficients initialized. Use principal_extension, or ignore_coefficients to ignore this.")

        return [ self.f_polynomial(k,ignore_coefficients=ignore_coefficients) for k in range(self._n) ]

    def g_vector(self,k,ignore_coefficients=False):
        r"""
        Returns the ``k``-th *g-vector* of ``self``. This is the degree vector
        of the ``k``-th cluster variable after setting all `y_i`'s to `0`.

        Requires principal coefficients, initialized by using principal_extension(),
        or the user can set 'ignore_coefficients=True' to bypass this restriction.

        Warning: this method assumes the sign-coherence conjecture and that the
        input seed is sign-coherent (has an exchange matrix with columns of like signs).
        Otherwise, computational errors might arise.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: [ S.g_vector(k) for k in range(3) ]
            [(1, 0, 0), (0, 0, -1), (0, -1, 0)]
        """
        if not (ignore_coefficients or self._is_principal):
            raise ValueError("No principal coefficients initialized. Use principal_extension, or ignore_coefficients to ignore this.")
        if k not in range(self._n):
            raise ValueError("The cluster seed does not have a cluster variable of index %s."%k)
        f = self.cluster_variable(k)
        eval_dict = dict( [ ( self.y(i), 0 ) for i in range(self._m) ] )
        f0 = f.subs(eval_dict)
        d1 = f0.numerator().degrees()
        d2 = f0.denominator().degrees()
        return tuple( d1[i] - d2[i] for i in range(self._n) )

    def g_matrix(self,ignore_coefficients=False):
        r"""
        Returns the matrix of all *g-vectors* of ``self``. This are the degree vectors
        of the cluster variables after setting all `y_i`'s to `0`.

        Requires principal coefficients, initialized by using principal_extension(),
        or the user can set 'ignore_coefficients=True' to bypass this restriction.

        Warning: this method assumes the sign-coherence conjecture and that the
        input seed is sign-coherent (has an exchange matrix with columns of like signs).
        Otherwise, computational errors might arise.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: S.g_matrix()
            [ 1  0  0]
            [ 0  0 -1]
            [ 0 -1  0]

            sage: S = ClusterSeed(['A',3])
            sage: S2 = S.principal_extension()
            sage: S.mutate([0,1])
            sage: S2.mutate([0,1])
            sage: S.g_matrix()
            Traceback (most recent call last):
            ...
            ValueError: No principal coefficients initialized. Use
            principal_extension, or ignore_coefficients to ignore this.
            sage: S.g_matrix(ignore_coefficients=True)
            [-1  0  0]
            [ 1  0  0]
            [ 0  1  1]
            sage: S2.g_matrix()
            [-1 -1  0]
            [ 1  0  0]
            [ 0  0  1]
        """
        from sage.matrix.all import matrix
        if not (ignore_coefficients or self._is_principal):
            raise ValueError("No principal coefficients initialized. Use principal_extension, or ignore_coefficients to ignore this.")
        return matrix( [ self.g_vector(k,ignore_coefficients=ignore_coefficients) for k in range(self._n) ] ).transpose()

    def c_vector(self,k,ignore_coefficients=False):
        r"""
        Returns the ``k``-th *c-vector* of ``self``. It is obtained as the
        ``k``-th column vector of the bottom part of the ``B``-matrix of ``self``.

        Requires principal coefficients, initialized by using principal_extension(),
        or the user can set 'ignore_coefficients=True' to bypass this restriction.

        Warning: this method assumes the sign-coherence conjecture and that the
        input seed is sign-coherent (has an exchange matrix with columns of like signs).
        Otherwise, computational errors might arise.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: [ S.c_vector(k) for k in range(3) ]
            [(1, 0, 0), (0, 0, -1), (0, -1, 0)]

            sage: S = ClusterSeed(Matrix([[0,1],[-1,0],[1,0],[-1,1]])); S
            A seed for a cluster algebra of rank 2 with 2 frozen variables
            sage: S.c_vector(0)
            Traceback (most recent call last):
            ...
            ValueError: No principal coefficients initialized. Use principal_extension, or ignore_coefficients to ignore this.
            sage: S.c_vector(0,ignore_coefficients=True)
            (1, -1)
        """
        if k not in range(self._n):
            raise ValueError("The cluster seed does not have a c-vector of index %s."%k)
        if not (ignore_coefficients or self._is_principal):
            raise ValueError("No principal coefficients initialized. Use principal_extension, or ignore_coefficients to ignore this.")
        return tuple( self._M[i,k] for i in range(self._n,self._n+self._m) )

    def c_matrix(self,ignore_coefficients=False):
        r"""
        Returns all *c-vectors* of ``self``.

        Requires principal coefficients, initialized by using principal_extension(),
        or the user can set 'ignore_coefficients=True' to bypass this restriction.

        Warning: this method assumes the sign-coherence conjecture and that the
        input seed is sign-coherent (has an exchange matrix with columns of like signs).
        Otherwise, computational errors might arise.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: S.c_matrix()
            [ 1  0  0]
            [ 0  0 -1]
            [ 0 -1  0]
        """
        if not (ignore_coefficients or self._is_principal):
            raise ValueError("No principal coefficients initialized. Use principal_extension, or ignore_coefficients to ignore this.")

        return self._M.submatrix(self._n,0)

    def coefficient(self,k):
        r"""
        Returns the *coefficient* of ``self`` at index ``k``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: [ S.coefficient(k) for k in range(3) ]
            [y0, 1/y2, 1/y1]
        """
        from sage.misc.all import prod
        if k not in range(self._n):
            raise ValueError("The cluster seed does not have a coefficient of index %s."%k)
        if self._m == 0:
            return self.x(0)**0
        #### Note: this special case m = 0 no longer needed except if we want type(answer) to be a cluster variable rather than an integer.
        else:
            exp = self.c_vector(k,ignore_coefficients=True)
            return prod( self.y(i)**exp[i] for i in xrange(self._m) )

    def coefficients(self):
        r"""
        Returns all *coefficients* of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: S.coefficients()
            [y0, 1/y2, 1/y1]
        """
        return [ self.coefficient(k) for k in range(self._n) ]

    def quiver(self):
        r"""
        Returns the *quiver* associated to ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.quiver()
            Quiver on 3 vertices of type ['A', 3]
        """
        from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        if self._quiver is None:
            self._quiver = ClusterQuiver( self._M )
        return self._quiver

    def is_acyclic(self):
        r"""
        Returns True iff self is acyclic (i.e., if the underlying quiver is acyclic).

        EXAMPLES::

            sage: ClusterSeed(['A',4]).is_acyclic()
            True

            sage: ClusterSeed(['A',[2,1],1]).is_acyclic()
            True

            sage: ClusterSeed([[0,1],[1,2],[2,0]]).is_acyclic()
            False
        """
        return self.quiver()._digraph.is_directed_acyclic()

    def is_bipartite(self,return_bipartition=False):
        r"""
        Returns True iff self is bipartite (i.e., if the underlying quiver is bipartite).

        INPUT:

        - return_bipartition -- (default:False) if True, the bipartition is returned in the case of ``self`` being bipartite.

        EXAMPLES::

            sage: ClusterSeed(['A',[3,3],1]).is_bipartite()
            True

            sage: ClusterSeed(['A',[4,3],1]).is_bipartite()
            False
        """
        return self.quiver().is_bipartite(return_bipartition=return_bipartition)

    def mutate(self, sequence, inplace=True):
        r"""
        Mutates ``self`` at a vertex or a sequence of vertices.

        INPUT:

        - ``sequence`` -- a vertex of self or an iterator of vertices of self.
        - ``inplace`` -- (default: True) if False, the result is returned, otherwise ``self`` is modified.

        EXAMPLES::

            sage: S = ClusterSeed(['A',4]); S.b_matrix()
            [ 0  1  0  0]
            [-1  0 -1  0]
            [ 0  1  0  1]
            [ 0  0 -1  0]

            sage: S.mutate(0); S.b_matrix()
            [ 0 -1  0  0]
            [ 1  0 -1  0]
            [ 0  1  0  1]
            [ 0  0 -1  0]

            sage: T = S.mutate(0, inplace=False); T
            A seed for a cluster algebra of rank 4 of type ['A', 4]

            sage: S.mutate(0)
            sage: S == T
            True

            sage: S.mutate([0,1,0])
            sage: S.b_matrix()
            [ 0 -1  1  0]
            [ 1  0  0  0]
            [-1  0  0  1]
            [ 0  0 -1  0]

            sage: S = ClusterSeed(QuiverMutationType([['A',1],['A',3]]))
            sage: S.b_matrix()
            [ 0  0  0  0]
            [ 0  0  1  0]
            [ 0 -1  0 -1]
            [ 0  0  1  0]

            sage: T = S.mutate(0,inplace=False)
            sage: S == T
            False
        """
        if inplace:
            seed = self
        else:
            seed = ClusterSeed( self )

        n, m = seed._n, seed._m
        V = range(n)

        if sequence in V:
            seq = [sequence]
        else:
            seq = sequence
        if isinstance(seq, tuple):
            seq = list( seq )
        if not isinstance(seq, list):
            raise ValueError('The quiver can only be mutated at a vertex or at a sequence of vertices')
        if not isinstance(inplace, bool):
            raise ValueError('The second parameter must be boolean.  To mutate at a sequence of length 2, input it as a list.')
        if any( v not in V for v in seq ):
            v = filter( lambda v: v not in V, seq )[0]
            raise ValueError('The quiver cannot be mutated at the vertex ' + str( v ))

        for k in seq:
            M = seed._M
            cluster = seed._cluster
            mon_p = seed._R(1)
            mon_n = seed._R(1)

            for j in range(n+m):
                if M[j,k] > 0:
                    mon_p = mon_p*cluster[j]**M[j,k]
                elif M[j,k] < 0:
                    mon_n = mon_n*cluster[j]**(-M[j,k])

            cluster[k] = (mon_p+mon_n)*cluster[k]**(-1)
            seed._M.mutate(k)
            #seed._M = _matrix_mutate( seed._M, k )

        seed._quiver = None
        if not inplace:
            return seed

    def mutation_sequence(self, sequence, show_sequence=False, fig_size=1.2,return_output='seed'):
        r"""
        Returns the seeds obtained by mutating ``self`` at all vertices in ``sequence``.

        INPUT:

        - ``sequence`` -- an iterable of vertices of self.
        - ``show_sequence`` -- (default: False) if True, a png containing the associated quivers is shown.
        - ``fig_size`` -- (default: 1.2) factor by which the size of the plot is multiplied.
        - ``return_output`` -- (default: 'seed') determines what output is to be returned::

            * if 'seed', outputs all the cluster seeds obtained by the ``sequence`` of mutations.
            * if 'matrix', outputs a list of exchange matrices.
            * if 'var', outputs a list of new cluster variables obtained at each step.

        EXAMPLES::

            sage: S = ClusterSeed(['A',2])
            sage: for T in S.mutation_sequence([0,1,0]):
            ...     print T.b_matrix()
            [ 0 -1]
            [ 1  0]
            [ 0  1]
            [-1  0]
            [ 0 -1]
            [ 1  0]

            sage: S=ClusterSeed(['A',2])
            sage: S.mutation_sequence([0,1,0,1],return_output='var')
            [(x1 + 1)/x0, (x0 + x1 + 1)/(x0*x1), (x0 + 1)/x1, x0]
        """
        seed = ClusterSeed( self )

        new_clust_var = []
        seed_sequence = []

        for v in sequence:
            seed = seed.mutate(v,inplace=False)
            new_clust_var.append( seed._cluster[v])
            seed_sequence.append( seed )

        if show_sequence:
            self.quiver().mutation_sequence(sequence=sequence, show_sequence=True, fig_size=fig_size )

        if return_output=='seed':
            return seed_sequence
        elif return_output=='matrix':
            return [ seed._M for seed in seed_sequence ]
        elif return_output=='var':
            return new_clust_var
        else:
            raise ValueError('The parameter `return_output` can only be `seed`, `matrix`, or `var`.')

    def exchangeable_part(self):
        r"""
        Returns the restriction to the principal part (i.e. the exchangeable variables) of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',4])
            sage: T = ClusterSeed( S.quiver().digraph().edges(), frozen=1 )
            sage: T.quiver().digraph().edges()
            [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]

            sage: T.exchangeable_part().quiver().digraph().edges()
            [(0, 1, (1, -1)), (2, 1, (1, -1))]

            sage: S2 = S.principal_extension()
            sage: S3 = S2.principal_extension(ignore_coefficients=True)
            sage: S2.exchangeable_part() == S3.exchangeable_part()
            True
        """
        from sage.combinat.cluster_algebra_quiver.mutation_class import _principal_part
        eval_dict = dict( [ ( self.y(i), 1 ) for i in xrange(self._m) ] )
        seed = ClusterSeed( _principal_part( self._M ) )
        seed._cluster = [ self._cluster[k].subs(eval_dict) for k in xrange(self._n) ]
        seed._mutation_type = self._mutation_type
        return seed

    def universal_extension(self):
        r"""
        Returns the universal extension of ``self``.

        This is the initial seed of the associated cluster algebra
        with universal coefficients, as defined in section 12 of
        :arxiv:`math/0602259`.

        This method works only if ``self`` is a bipartite, finite-type seed.

        Due to some limitations in the current implementation of
        ``CartanType``, we need to construct the set of almost positive
        coroots by hand. As a consequence their ordering is not the
        standard one (the rows of the bottom part of the exchange
        matrix might be a shuffling of those you would expect).

        EXAMPLES::

            sage: S = ClusterSeed(['A',2])
            sage: T = S.universal_extension()
            sage: T.b_matrix()
            [ 0  1]
            [-1  0]
            [-1  0]
            [ 1  0]
            [ 1 -1]
            [ 0  1]
            [ 0 -1]

            sage: S = ClusterSeed(['A',3])
            sage: T = S.universal_extension()
            sage: T.b_matrix()
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]
            [-1  0  0]
            [ 1  0  0]
            [ 1 -1  0]
            [ 1 -1  1]
            [ 0  1  0]
            [ 0 -1  0]
            [ 0 -1  1]
            [ 0  0 -1]
            [ 0  0  1]

            sage: S = ClusterSeed(['B',2])
            sage: T = S.universal_extension()
            sage: T.b_matrix()
            [ 0  1]
            [-2  0]
            [-1  0]
            [ 1  0]
            [ 1 -1]
            [ 2 -1]
            [ 0  1]
            [ 0 -1]

        """
        if self._m != 0:
            raise ValueError("To have universal coefficients we need "
                             "to start from a coefficient-free seed")
        if not self.is_bipartite() or not self.is_finite():
            raise ValueError("Universal coefficients are defined only "
                             "for finite type cluster algebras at a "
                             "bipartite initial cluster")

        from sage.matrix.all import matrix
        from sage.combinat.root_system.cartan_matrix import CartanMatrix

        A = 2 - self.b_matrix().apply_map(abs).transpose()

        rs = CartanMatrix(A).root_space()
        almost_positive_coroots = rs.almost_positive_roots()

        sign = [-1 if all(x <= 0 for x in self.b_matrix()[i]) else 1
                for i in range(self._n)]
        C = matrix([[sign[j] * alpha[j + 1] for j in range(self._n)]
                    for alpha in almost_positive_coroots])

        M = self._M.stack(C)
        seed = ClusterSeed(M, is_principal=False)
        seed._mutation_type = self._mutation_type
        return seed

    def principal_extension(self,ignore_coefficients=False):
        r"""
        Returns the principal extension of self, yielding a 2n-by-n matrix.  Raises an error if the input seed has a non-square exchange matrix,
        unless 'ignore_coefficients=True' is set.  In this case, the method instead adds n frozen variables to any previously frozen variables.
        I.e., the seed obtained by adding a frozen variable to every exchangeable variable of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed([[0,1],[1,2],[2,3],[2,4]]); S
            A seed for a cluster algebra of rank 5

            sage: T = S.principal_extension(); T
            A seed for a cluster algebra of rank 5 with principal coefficients

            sage: T.b_matrix()
            [ 0  1  0  0  0]
            [-1  0  1  0  0]
            [ 0 -1  0  1  1]
            [ 0  0 -1  0  0]
            [ 0  0 -1  0  0]
            [ 1  0  0  0  0]
            [ 0  1  0  0  0]
            [ 0  0  1  0  0]
            [ 0  0  0  1  0]
            [ 0  0  0  0  1]

            sage: T2 = T.principal_extension()
            Traceback (most recent call last):
            ...
            ValueError: The b-matrix is not square. Use ignore_coefficients to ignore this.

            sage: T2 = T.principal_extension(ignore_coefficients=True); T2.b_matrix()
            [ 0  1  0  0  0]
            [-1  0  1  0  0]
            [ 0 -1  0  1  1]
            [ 0  0 -1  0  0]
            [ 0  0 -1  0  0]
            [ 1  0  0  0  0]
            [ 0  1  0  0  0]
            [ 0  0  1  0  0]
            [ 0  0  0  1  0]
            [ 0  0  0  0  1]
            [ 1  0  0  0  0]
            [ 0  1  0  0  0]
            [ 0  0  1  0  0]
            [ 0  0  0  1  0]
            [ 0  0  0  0  1]
            """
        from sage.matrix.all import identity_matrix
        if not ignore_coefficients and self._m != 0:
            raise ValueError("The b-matrix is not square. Use ignore_coefficients to ignore this.")
        M = self._M.stack(identity_matrix(self._n))
        is_principal = (self._m == 0)
        seed = ClusterSeed( M, is_principal=is_principal )
        seed._mutation_type = self._mutation_type
        return seed

    def reorient( self, data ):
        r"""
        Reorients ``self`` with respect to the given total order,
        or with respect to an iterator of ordered pairs.

        WARNING:

        - This operation might change the mutation type of ``self``.
        - Ignores ordered pairs `(i,j)` for which neither `(i,j)` nor `(j,i)` is an edge of ``self``.

        INPUT:

        - ``data`` -- an iterator defining a total order on ``self.vertices()``, or an iterator of ordered pairs in ``self`` defining the new orientation of these edges.

        EXAMPLES::

            sage: S = ClusterSeed(['A',[2,3],1])
            sage: S.mutation_type()
            ['A', [2, 3], 1]

            sage: S.reorient([(0,1),(2,3)])
            sage: S.mutation_type()
            ['D', 5]

            sage: S.reorient([(1,0),(2,3)])
            sage: S.mutation_type()
            ['A', [1, 4], 1]

            sage: S.reorient([0,1,2,3,4])
            sage: S.mutation_type()
            ['A', [1, 4], 1]
        """
        if not self._quiver:
            self.quiver()
        self._quiver.reorient( data )
        self._M = self._quiver._M
        self.reset_cluster()
        self._mutation_type = None

    def set_cluster( self, cluster ):
        r"""
        Sets the cluster for ``self`` to ``cluster``.

        INPUT:

        - ``cluster`` -- an iterable defining a cluster for ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: cluster = S.cluster()
            sage: S.mutate([1,2,1])
            sage: S.cluster()
            [x0, (x1 + 1)/x2, (x0*x2 + x1 + 1)/(x1*x2)]

            sage: S.set_cluster(cluster)
            sage: S.cluster()
            [x0, x1, x2]
        """
        if not len(cluster) == self._n+self._m:
            raise ValueError('The number of given cluster variables is wrong')
        if any(c not in self._R for c in cluster):
            raise ValueError('The cluster variables are not all contained in %s'%self._R)
        self._cluster = [ self._R(x) for x in cluster ]
        self._is_principal = None

    def reset_cluster( self ):
        r"""
        Resets the cluster of ``self`` to the initial cluster.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.mutate([1,2,1])
            sage: S.cluster()
            [x0, (x1 + 1)/x2, (x0*x2 + x1 + 1)/(x1*x2)]

            sage: S.reset_cluster()
            sage: S.cluster()
            [x0, x1, x2]

            sage: T = S.principal_extension()
            sage: T.cluster()
            [x0, x1, x2]
            sage: T.mutate([1,2,1])
            sage: T.cluster()
            [x0, (x1*y2 + x0)/x2, (x1*y1*y2 + x0*y1 + x2)/(x1*x2)]

            sage: T.reset_cluster()
            sage: T.cluster()
            [x0, x1, x2]
        """
        self.set_cluster(self._R.gens())

    def reset_coefficients( self ):
        r"""
        Resets the coefficients of ``self`` to the frozen variables but keeps the current cluster.
        Raises an error if the number of frozen variables is different than the number of exchangeable variables.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.b_matrix()
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]
            [ 1  0  0]
            [ 0  1  0]
            [ 0  0  1]
            sage: S.mutate([1,2,1])
            sage: S.b_matrix()
            [ 0  1 -1]
            [-1  0  1]
            [ 1 -1  0]
            [ 1  0  0]
            [ 0  1 -1]
            [ 0  0 -1]
            sage: S.reset_coefficients()
            sage: S.b_matrix()
            [ 0  1 -1]
            [-1  0  1]
            [ 1 -1  0]
            [ 1  0  0]
            [ 0  1  0]
            [ 0  0  1]
        """
        n,m = self._n, self._m
        if not n == m:
            raise ValueError("The numbers of cluster variables and of frozen variables do not coincide.")
        for i in xrange(m):
            for j in xrange(n):
                if i == j:
                    self._M[i+n,j] = 1
                else:
                    self._M[i+n,j] = 0
        self._quiver = None
        self._is_principal = None

    def mutation_class_iter( self, depth=infinity, show_depth=False, return_paths=False, up_to_equivalence=True, only_sink_source=False ):
        r"""
        Returns an iterator for the mutation class of ``self`` with respect to certain constrains.

        INPUT:

        - ``depth`` -- (default: infinity) integer or infinity, only seeds with distance at most ``depth`` from ``self`` are returned.
        - ``show_depth`` -- (default: False) if True, the current depth of the mutation is shown while computing.
        - ``return_paths`` -- (default: False) if True, a shortest path of mutations from ``self`` to the given quiver is returned as well.
        - ``up_to_equivalence`` -- (default: True) if True, only one seed up to simultaneous permutation of rows and columns of the exchange matrix is recorded.
        - ``sink_source`` -- (default: False) if True, only mutations at sinks and sources are applied.

        EXAMPLES:

        A standard finite type example::

            sage: S = ClusterSeed(['A',3])
            sage: it = S.mutation_class_iter()
            sage: for T in it: print T
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]

        A finite type example with given depth::

            sage: it = S.mutation_class_iter(depth=1)
            sage: for T in it: print T
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]
            A seed for a cluster algebra of rank 3 of type ['A', 3]

        A finite type example where the depth is shown while computing::

            sage: it = S.mutation_class_iter(show_depth=True)
            sage: for T in it: pass
            Depth: 0     found: 1          Time: ... s
            Depth: 1     found: 4          Time: ... s
            Depth: 2     found: 9          Time: ... s
            Depth: 3     found: 13         Time: ... s
            Depth: 4     found: 14         Time: ... s

        A finite type example with shortest paths returned::

            sage: it = S.mutation_class_iter(return_paths=True)
            sage: for T in it: print T
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [2])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [1])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [0])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [2, 1])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [0, 2])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [0, 1])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [1, 2])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [1, 0])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [0, 2, 1])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [0, 1, 2])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [2, 1, 0])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [1, 0, 2])
            (A seed for a cluster algebra of rank 3 of type ['A', 3], [0, 1, 2, 0])

        Finite type examples not considered up to equivalence::

            sage: it = S.mutation_class_iter(up_to_equivalence=False)
            sage: len( [ T for T in it ] )
            84

            sage: it = ClusterSeed(['A',2]).mutation_class_iter(return_paths=True,up_to_equivalence=False)
            sage: for T in it: print T
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [])
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [1])
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [0])
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [0, 1])
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [1, 0])
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [1, 0, 1])
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [0, 1, 0])
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [1, 0, 1, 0])
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [0, 1, 0, 1])
            (A seed for a cluster algebra of rank 2 of type ['A', 2], [1, 0, 1, 0, 1])

        Check that :trac:`14638` is fixed::

            sage: S = ClusterSeed(['E',6])
            sage: MC = S.mutation_class(depth=7); len(MC)
            534

        Infinite type examples::

            sage: S = ClusterSeed(['A',[1,1],1])
            sage: it = S.mutation_class_iter()
            sage: next(it)
            A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1]
            sage: next(it)
            A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1]
            sage: next(it)
            A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1]
            sage: next(it)
            A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1]

            sage: it = S.mutation_class_iter(depth=3, return_paths=True)
            sage: for T in it: print T
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [1])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [0])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [1, 0])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [0, 1])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [1, 0, 1])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [0, 1, 0])
        """
        depth_counter = 0
        n = self._n
        timer = time.time()
        if return_paths:
            yield (self,[])
        else:
            yield self
        if up_to_equivalence:
            cl = Set( self._cluster )
        else:
            cl = tuple( self._cluster )
        clusters = {}
        clusters[ cl ] = [ self, range(n), [] ]
        gets_bigger = True
        if show_depth:
            timer2 = time.time()
            dc = str(depth_counter)
            dc += ' ' * (5-len(dc))
            nr = str(len(clusters))
            nr += ' ' * (10-len(nr))
            print "Depth: %s found: %s Time: %.2f s"%(dc,nr,timer2-timer)
        while gets_bigger and depth_counter < depth:
            gets_bigger = False
            keys = clusters.keys()
            for key in keys:
                sd = clusters[key]
                while sd[1]:
                    i = sd[1].pop()
                    if not only_sink_source or all( entry >= 0 for entry in sd[0]._M.row( i ) ) or all( entry <= 0 for entry in sd[0]._M.row( i ) ):
                        sd2  = sd[0].mutate( i, inplace=False )
                        if up_to_equivalence:
                            cl2 = Set(sd2._cluster)
                        else:
                            cl2 = tuple(sd2._cluster)
                        if cl2 in clusters:
                            if not up_to_equivalence and i in clusters[cl2][1]:
                                clusters[cl2][1].remove(i)
                        else:
                            gets_bigger = True
                            if only_sink_source:
                                orbits = range(n)
                            else:
                                orbits = [ index for index in xrange(n) if index > i or sd2._M[index,i] != 0 ]

                            clusters[ cl2 ] = [ sd2, orbits, clusters[key][2]+[i] ]
                            if return_paths:
                                yield (sd2,clusters[cl2][2])
                            else:
                                yield sd2
            depth_counter += 1
            if show_depth and gets_bigger:
                timer2 = time.time()
                dc = str(depth_counter)
                dc += ' ' * (5-len(dc))
                nr = str(len(clusters))
                nr += ' ' * (10-len(nr))
                print "Depth: %s found: %s Time: %.2f s"%(dc,nr,timer2-timer)

    def mutation_class( self, depth=infinity, show_depth=False, return_paths=False, up_to_equivalence=True, only_sink_source=False ):
        r"""
        Returns the mutation class of ``self`` with respect to certain constraints.

        INPUT:

        - ``depth`` -- (default: infinity) integer, only seeds with distance at most depth from self are returned.
        - ``show_depth`` -- (default: False) if True, the actual depth of the mutation is shown.
        - ``return_paths`` -- (default: False) if True, a shortest path of mutation sequences from self to the given quiver is returned as well.
        - ``up_to_equivalence`` -- (default: True) if True, only seeds up to equivalence are considered.
        - ``sink_source`` -- (default: False) if True, only mutations at sinks and sources are applied.

        EXAMPLES:

        - for examples see :meth:`mutation_class_iter`

        TESTS::

            sage: A = ClusterSeed(['A',3]).mutation_class()
        """
        if depth is infinity and not self.is_finite():
            raise ValueError('The mutation class can - for infinite types - only be computed up to a given depth')
        return list( S for S in self.mutation_class_iter( depth=depth, show_depth=show_depth, return_paths=return_paths, up_to_equivalence=up_to_equivalence, only_sink_source=only_sink_source ) )

    def cluster_class_iter(self, depth=infinity, show_depth=False, up_to_equivalence=True):
        r"""
        Returns an iterator through all clusters in the mutation class of ``self``.

        INPUT:

        - ``depth`` -- (default: infinity) integer or infinity, only seeds with distance at most depth from self are returned
        - ``show_depth`` -- (default False) - if True, ignored if depth is set; returns the depth of the mutation class, i.e., the maximal distance from self of an element in the mutation class
        - ``up_to_equivalence`` -- (default: True) if True, only clusters up to equivalence are considered.

        EXAMPLES:

        A standard finite type example::

            sage: S = ClusterSeed(['A',3])
            sage: it = S.cluster_class_iter()
            sage: for T in it: print T
            [x0, x1, x2]
            [x0, x1, (x1 + 1)/x2]
            [x0, (x0*x2 + 1)/x1, x2]
            [(x1 + 1)/x0, x1, x2]
            [x0, (x0*x2 + x1 + 1)/(x1*x2), (x1 + 1)/x2]
            [(x1 + 1)/x0, x1, (x1 + 1)/x2]
            [(x1 + 1)/x0, (x0*x2 + x1 + 1)/(x0*x1), x2]
            [x0, (x0*x2 + 1)/x1, (x0*x2 + x1 + 1)/(x1*x2)]
            [(x0*x2 + x1 + 1)/(x0*x1), (x0*x2 + 1)/x1, x2]
            [(x1 + 1)/x0, (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2), (x1 + 1)/x2]
            [(x1 + 1)/x0, (x0*x2 + x1 + 1)/(x0*x1), (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2)]
            [(x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2), (x0*x2 + x1 + 1)/(x1*x2), (x1 + 1)/x2]
            [(x0*x2 + x1 + 1)/(x0*x1), (x0*x2 + 1)/x1, (x0*x2 + x1 + 1)/(x1*x2)]
            [(x0*x2 + x1 + 1)/(x1*x2), (x0*x2 + x1 + 1)/(x0*x1), (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2)]

        A finite type example with given depth::

            sage: it = S.cluster_class_iter(depth=1)
            sage: for T in it: print T
            [x0, x1, x2]
            [x0, x1, (x1 + 1)/x2]
            [x0, (x0*x2 + 1)/x1, x2]
            [(x1 + 1)/x0, x1, x2]

        A finite type example where the depth is returned while computing::

            sage: it = S.cluster_class_iter(show_depth=True)
            sage: for T in it: print T
            [x0, x1, x2]
            Depth: 0     found: 1          Time: ... s
            [x0, x1, (x1 + 1)/x2]
            [x0, (x0*x2 + 1)/x1, x2]
            [(x1 + 1)/x0, x1, x2]
            Depth: 1     found: 4          Time: ... s
            [x0, (x0*x2 + x1 + 1)/(x1*x2), (x1 + 1)/x2]
            [(x1 + 1)/x0, x1, (x1 + 1)/x2]
            [(x1 + 1)/x0, (x0*x2 + x1 + 1)/(x0*x1), x2]
            [x0, (x0*x2 + 1)/x1, (x0*x2 + x1 + 1)/(x1*x2)]
            [(x0*x2 + x1 + 1)/(x0*x1), (x0*x2 + 1)/x1, x2]
            Depth: 2     found: 9          Time: ... s
            [(x1 + 1)/x0, (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2), (x1 + 1)/x2]
            [(x1 + 1)/x0, (x0*x2 + x1 + 1)/(x0*x1), (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2)]
            [(x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2), (x0*x2 + x1 + 1)/(x1*x2), (x1 + 1)/x2]
            [(x0*x2 + x1 + 1)/(x0*x1), (x0*x2 + 1)/x1, (x0*x2 + x1 + 1)/(x1*x2)]
            Depth: 3     found: 13         Time: ... s
            [(x0*x2 + x1 + 1)/(x1*x2), (x0*x2 + x1 + 1)/(x0*x1), (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2)]
            Depth: 4     found: 14         Time: ... s

        Finite type examples not considered up to equivalence::

            sage: it = S.cluster_class_iter(up_to_equivalence=False)
            sage: len( [ T for T in it ] )
            84

            sage: it = ClusterSeed(['A',2]).cluster_class_iter(up_to_equivalence=False)
            sage: for T in it: print T
            [x0, x1]
            [x0, (x0 + 1)/x1]
            [(x1 + 1)/x0, x1]
            [(x1 + 1)/x0, (x0 + x1 + 1)/(x0*x1)]
            [(x0 + x1 + 1)/(x0*x1), (x0 + 1)/x1]
            [(x0 + x1 + 1)/(x0*x1), (x1 + 1)/x0]
            [(x0 + 1)/x1, (x0 + x1 + 1)/(x0*x1)]
            [x1, (x1 + 1)/x0]
            [(x0 + 1)/x1, x0]
            [x1, x0]

        Infinite type examples::

            sage: S = ClusterSeed(['A',[1,1],1])
            sage: it = S.cluster_class_iter()
            sage: next(it)
            [x0, x1]
            sage: next(it)
            [x0, (x0^2 + 1)/x1]
            sage: next(it)
            [(x1^2 + 1)/x0, x1]
            sage: next(it)
            [(x0^4 + 2*x0^2 + x1^2 + 1)/(x0*x1^2), (x0^2 + 1)/x1]
            sage: next(it)
            [(x1^2 + 1)/x0, (x1^4 + x0^2 + 2*x1^2 + 1)/(x0^2*x1)]

            sage: it = S.cluster_class_iter(depth=3)
            sage: for T in it: print T
            [x0, x1]
            [x0, (x0^2 + 1)/x1]
            [(x1^2 + 1)/x0, x1]
            [(x0^4 + 2*x0^2 + x1^2 + 1)/(x0*x1^2), (x0^2 + 1)/x1]
            [(x1^2 + 1)/x0, (x1^4 + x0^2 + 2*x1^2 + 1)/(x0^2*x1)]
            [(x0^4 + 2*x0^2 + x1^2 + 1)/(x0*x1^2), (x0^6 + 3*x0^4 + 2*x0^2*x1^2 + x1^4 + 3*x0^2 + 2*x1^2 + 1)/(x0^2*x1^3)]
            [(x1^6 + x0^4 + 2*x0^2*x1^2 + 3*x1^4 + 2*x0^2 + 3*x1^2 + 1)/(x0^3*x1^2), (x1^4 + x0^2 + 2*x1^2 + 1)/(x0^2*x1)]
        """
        mc_iter = self.mutation_class_iter( depth=depth, show_depth=show_depth, up_to_equivalence=up_to_equivalence )
        for c in mc_iter:
            yield c.cluster()

    def cluster_class(self, depth=infinity, show_depth=False, up_to_equivalence=True):
        r"""
        Returns the cluster class of ``self`` with respect to certain constraints.

        INPUT:

        - ``depth`` -- (default: infinity) integer, only seeds with distance at most depth from self are returned
        - ``return_depth`` -- (default False) - if True, ignored if depth is set; returns the depth of the mutation class, i.e., the maximal distance from self of an element in the mutation class
        - ``up_to_equivalence`` -- (default: True) if True, only clusters up to equivalence are considered.

        EXAMPLES:

        - for examples see :meth:`cluster_class_iter`

        TESTS::

            sage: A = ClusterSeed(['A',3]).cluster_class()
        """
        if depth is infinity and not self.is_finite():
            raise ValueError('The variable class can - for infinite types - only be computed up to a given depth')

        return [ c for c in self.cluster_class_iter(depth=depth, show_depth=show_depth, up_to_equivalence=up_to_equivalence) ]

    def b_matrix_class_iter(self, depth=infinity, up_to_equivalence=True):
        r"""
        Returns an iterator through all `B`-matrices in the mutation class of ``self``.

        INPUT:

        - ``depth`` -- (default:infinity) integer or infinity, only seeds with distance at most depth from self are returned
        - ``up_to_equivalence`` -- (default: True) if True, only 'B'-matrices up to equivalence are considered.

        EXAMPLES:

        A standard finite type example::

            sage: S = ClusterSeed(['A',4])
            sage: it = S.b_matrix_class_iter()
            sage: for T in it: print T
            [ 0  0  0  1]
            [ 0  0  1  1]
            [ 0 -1  0  0]
            [-1 -1  0  0]
            [ 0  0  0  1]
            [ 0  0  1  0]
            [ 0 -1  0  1]
            [-1  0 -1  0]
            [ 0  0  1  1]
            [ 0  0  0 -1]
            [-1  0  0  0]
            [-1  1  0  0]
            [ 0  0  0  1]
            [ 0  0 -1  1]
            [ 0  1  0 -1]
            [-1 -1  1  0]
            [ 0  0  0  1]
            [ 0  0 -1  0]
            [ 0  1  0 -1]
            [-1  0  1  0]
            [ 0  0  0 -1]
            [ 0  0 -1  1]
            [ 0  1  0 -1]
            [ 1 -1  1  0]

        A finite type example with given depth::

            sage: it = S.b_matrix_class_iter(depth=1)
            sage: for T in it: print T
            [ 0  0  0  1]
            [ 0  0  1  1]
            [ 0 -1  0  0]
            [-1 -1  0  0]
            [ 0  0  0  1]
            [ 0  0  1  0]
            [ 0 -1  0  1]
            [-1  0 -1  0]
            [ 0  0  1  1]
            [ 0  0  0 -1]
            [-1  0  0  0]
            [-1  1  0  0]

        Finite type example not considered up to equivalence::

            sage: S = ClusterSeed(['A',3])
            sage: it = S.b_matrix_class_iter(up_to_equivalence=False)
            sage: for T in it: print T
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]
            [ 0  1  0]
            [-1  0  1]
            [ 0 -1  0]
            [ 0 -1  0]
            [ 1  0  1]
            [ 0 -1  0]
            [ 0 -1  0]
            [ 1  0 -1]
            [ 0  1  0]
            [ 0 -1  1]
            [ 1  0 -1]
            [-1  1  0]
            [ 0  1 -1]
            [-1  0  1]
            [ 1 -1  0]
            [ 0  0  1]
            [ 0  0 -1]
            [-1  1  0]
            [ 0 -1  1]
            [ 1  0  0]
            [-1  0  0]
            [ 0  0 -1]
            [ 0  0  1]
            [ 1 -1  0]
            [ 0  1 -1]
            [-1  0  0]
            [ 1  0  0]
            [ 0  1  1]
            [-1  0  0]
            [-1  0  0]
            [ 0 -1 -1]
            [ 1  0  0]
            [ 1  0  0]
            [ 0  0 -1]
            [ 0  0 -1]
            [ 1  1  0]
            [ 0  0  1]
            [ 0  0  1]
            [-1 -1  0]

        Infinite (but finite mutation) type example::

            sage: S = ClusterSeed(['A',[1,2],1])
            sage: it = S.b_matrix_class_iter()
            sage: for T in it: print T
            [ 0  1  1]
            [-1  0  1]
            [-1 -1  0]
            [ 0 -2  1]
            [ 2  0 -1]
            [-1  1  0]

        Infinite mutation type example::

            sage: S = ClusterSeed(['E',10])
            sage: it = S.b_matrix_class_iter(depth=3)
            sage: len ( [T for T in it] )
            266
        """
        Q = self.quiver()
        for M in Q.mutation_class_iter( depth=depth, up_to_equivalence=up_to_equivalence, data_type='matrix' ):
            yield M

    def b_matrix_class(self, depth=infinity, up_to_equivalence=True):
        r"""
        Returns all `B`-matrices in the mutation class of ``self``.

        INPUT:

        - ``depth`` -- (default:infinity) integer or infinity, only seeds with distance at most depth from self are returned
        - ``up_to_equivalence`` -- (default: True) if True, only 'B'-matrices up to equivalence are considered.

        EXAMPLES:

        - for examples see :meth:`b_matrix_class_iter`

        TESTS::

            sage: A = ClusterSeed(['A',3]).b_matrix_class()
            sage: A = ClusterSeed(['A',[2,1],1]).b_matrix_class()
        """
        if depth is infinity and not self.is_mutation_finite():
            raise ValueError('The B-matrix class can - for infinite mutation types - only be computed up to a given depth')

        return [ M for M in self.b_matrix_class_iter( depth=depth, up_to_equivalence=up_to_equivalence ) ]

    def variable_class_iter(self, depth=infinity, ignore_bipartite_belt=False):
        r"""
        Returns an iterator for all cluster variables in the mutation class of ``self``.

        INPUT:

            - ``depth`` -- (default:infinity) integer, only seeds with distance at most depth from self are returned
            - ``ignore_bipartite_belt`` -- (default:False) if True, the algorithms does not use the bipartite belt

        EXAMPLES:

        A standard finite type example::

            sage: S = ClusterSeed(['A',3])
            sage: it = S.variable_class_iter()
            sage: for T in it: print T
            x0
            x1
            x2
            (x1 + 1)/x0
            (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2)
            (x1 + 1)/x2
            (x0*x2 + x1 + 1)/(x0*x1)
            (x0*x2 + 1)/x1
            (x0*x2 + x1 + 1)/(x1*x2)

        Finite type examples with given depth::

            sage: it = S.variable_class_iter(depth=1)
            sage: for T in it: print T
            Found a bipartite seed - restarting the depth counter at zero and constructing the variable class using its bipartite belt.
            x0
            x1
            x2
            (x1 + 1)/x0
            (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2)
            (x1 + 1)/x2
            (x0*x2 + x1 + 1)/(x0*x1)
            (x0*x2 + 1)/x1
            (x0*x2 + x1 + 1)/(x1*x2)

        Note that the notion of *depth* depends on whether a bipartite seed is found or not, or if it is manually ignored::

            sage: it = S.variable_class_iter(depth=1,ignore_bipartite_belt=True)
            sage: for T in it: print T
            x0
            x1
            x2
            (x1 + 1)/x2
            (x0*x2 + 1)/x1
            (x1 + 1)/x0

            sage: S.mutate([0,1])
            sage: it2 = S.variable_class_iter(depth=1)
            sage: for T in it2: print T
            (x1 + 1)/x0
            (x0*x2 + x1 + 1)/(x0*x1)
            x2
            (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2)
            x1
            (x0*x2 + 1)/x1

        Infinite type examples::

            sage: S = ClusterSeed(['A',[1,1],1])
            sage: it = S.variable_class_iter(depth=2)
            sage: for T in it: print T
            Found a bipartite seed - restarting the depth counter at zero and constructing the variable class using its bipartite belt.
            x0
            x1
            (x1^2 + 1)/x0
            (x1^4 + x0^2 + 2*x1^2 + 1)/(x0^2*x1)
            (x0^4 + 2*x0^2 + x1^2 + 1)/(x0*x1^2)
            (x0^2 + 1)/x1
            (x1^6 + x0^4 + 2*x0^2*x1^2 + 3*x1^4 + 2*x0^2 + 3*x1^2 + 1)/(x0^3*x1^2)
            (x1^8 + x0^6 + 2*x0^4*x1^2 + 3*x0^2*x1^4 + 4*x1^6 + 3*x0^4 + 6*x0^2*x1^2 + 6*x1^4 + 3*x0^2 + 4*x1^2 + 1)/(x0^4*x1^3)
            (x0^8 + 4*x0^6 + 3*x0^4*x1^2 + 2*x0^2*x1^4 + x1^6 + 6*x0^4 + 6*x0^2*x1^2 + 3*x1^4 + 4*x0^2 + 3*x1^2 + 1)/(x0^3*x1^4)
            (x0^6 + 3*x0^4 + 2*x0^2*x1^2 + x1^4 + 3*x0^2 + 2*x1^2 + 1)/(x0^2*x1^3)
        """
        mut_iter = self.mutation_class_iter( depth=depth,show_depth=False )
        var_class = set()

        for seed in mut_iter:
            if seed is self:
                seed = ClusterSeed(seed)
            if not ignore_bipartite_belt and seed.is_bipartite():
                bipartition = seed.is_bipartite(return_bipartition=True)
                bipartition = (list(bipartition[0]),list(bipartition[1]))
                if depth is not infinity:
                    print "Found a bipartite seed - restarting the depth counter at zero and constructing the variable class using its bipartite belt."
                depth_counter = 0
                end = False
                seed2 = ClusterSeed(seed)
                for c in seed._cluster:
                    if c not in var_class:
                        yield ClusterVariable( c.parent(), c.numerator(), c.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable' )
                var_class = var_class.union( seed._cluster )

                init_cluster = set(seed._cluster)
                while not end and depth_counter < depth:
                    depth_counter += 1
                    seed.mutate(bipartition[0])
                    seed.mutate(bipartition[1])
                    if set(seed._cluster) in [set(seed2._cluster),init_cluster]:
                        end = True
                    if not end:
                        for c in seed._cluster:
                            if c not in var_class:
                                yield ClusterVariable( c.parent(), c.numerator(), c.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable' )
                        var_class = var_class.union( seed._cluster )
                        seed2.mutate(bipartition[1])
                        seed2.mutate(bipartition[0])
                        if set(seed2._cluster) in [set(seed._cluster),init_cluster]:
                            end = True
                        if not end:
                            for c in seed2._cluster:
                                if c not in var_class:
                                    yield ClusterVariable( c.parent(), c.numerator(), c.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable' )
                            var_class = var_class.union(seed2._cluster)
                return
            else:
                for c in seed._cluster:
                    if c not in var_class:
                        yield ClusterVariable( c.parent(), c.numerator(), c.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable' )
                var_class = var_class.union(seed._cluster)

    def variable_class(self, depth=infinity, ignore_bipartite_belt=False):
        r"""
        Returns all cluster variables in the mutation class of ``self``.

        INPUT:

            - ``depth`` -- (default:infinity) integer, only seeds with distance at most depth from self are returned
            - ``ignore_bipartite_belt`` -- (default:False) if True, the algorithms does not use the bipartite belt


        EXAMPLES:

        - for examples see :meth:`variable_class_iter`

        TESTS::

            sage: A = ClusterSeed(['A',3]).variable_class()
        """
        if depth is infinity and not self.is_finite():
            raise ValueError('The variable class can - for infinite types - only be computed up to a given depth')

        var_iter = self.variable_class_iter( depth=depth, ignore_bipartite_belt=ignore_bipartite_belt )
        return sorted(var_iter)

    def is_finite( self ):
        r"""
        Returns True if ``self`` is of finite type.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.is_finite()
            True

            sage: S = ClusterSeed(['A',[2,2],1])
            sage: S.is_finite()
            False
        """
        mt = self.mutation_type()
        if isinstance(mt, str):
            return False
        else:
            return mt.is_finite()

    def is_mutation_finite( self, nr_of_checks=None, return_path=False ):
        r"""
        Returns True if ``self`` is of finite mutation type.

        INPUT:

        - ``nr_of_checks`` -- (default: None) number of mutations applied. Standard is 500*(number of vertices of self).
        - ``return_path`` -- (default: False) if True, in case of self not being mutation finite, a path from self to a quiver with an edge label (a,-b) and a*b > 4 is returned.

        ALGORITHM:

        - A cluster seed is mutation infinite if and only if every `b_{ij}*b_{ji} > -4`. Thus, we apply random mutations in random directions

        WARNING:

        - Uses a non-deterministic method by random mutations in various directions.
        - In theory, it can return a wrong True.

        EXAMPLES::

            sage: S = ClusterSeed(['A',10])
            sage: S._mutation_type = None
            sage: S.is_mutation_finite()
            True

            sage: S = ClusterSeed([(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(2,9)])
            sage: S.is_mutation_finite()
            False
        """
        is_finite, path = is_mutation_finite(copy(self._M),nr_of_checks=nr_of_checks)
        if return_path:
            return is_finite, path
        else:
            return is_finite

    def mutation_type(self):
        r"""
        Returns the mutation_type of each connected component of ``self``, if it can be determined.
        Otherwise, the mutation type of this component is set to be unknown.

        The mutation types of the components are ordered by vertex labels.

        WARNING:

        - All finite types can be detected,
        - All affine types can be detected, EXCEPT affine type D (the algorithm is not yet implemented)
        - All exceptional types can be detected.

        - Might fail to work if it is used within different Sage processes simultaneously (that happend in the doctesting).

        EXAMPLES:

        - finite types::

            sage: S = ClusterSeed(['A',5])
            sage: S._mutation_type = S._quiver._mutation_type = None
            sage: S.mutation_type()
            ['A', 5]

            sage: S = ClusterSeed([(0,1),(1,2),(2,3),(3,4)])
            sage: S.mutation_type()
            ['A', 5]

        - affine types::

            sage: S = ClusterSeed(['E',8,[1,1]]); S
            A seed for a cluster algebra of rank 10 of type ['E', 8, [1, 1]]
            sage: S._mutation_type = S._quiver._mutation_type = None; S
            A seed for a cluster algebra of rank 10
            sage: S.mutation_type() # long time
            ['E', 8, [1, 1]]

        - the not yet working affine type D::

            sage: S = ClusterSeed(['D',4,1])
            sage: S._mutation_type = S._quiver._mutation_type = None
            sage: S.mutation_type() # todo: not implemented
            ['D', 4, 1]

        - the exceptional types::

            sage: S = ClusterSeed(['X',6])
            sage: S._mutation_type = S._quiver._mutation_type = None
            sage: S.mutation_type() # long time
            ['X', 6]

        -  infinite types::

            sage: S = ClusterSeed(['GR',[4,9]])
            sage: S._mutation_type = S._quiver._mutation_type = None
            sage: S.mutation_type()
            'undetermined infinite mutation type'
        """
        if self._mutation_type is None:
            if self._quiver is None:
                self.quiver()
            self._mutation_type = self._quiver.mutation_type()
        return self._mutation_type

    def greedy(self, a1, a2, method='by_recursion'):
        r"""
        Returns the greedy element `x[a_1,a_2]` assuming that self is rank two.

        The third input can be 'by_recursion', 'by_combinatorics', or
        'just_numbers' to specify if the user wants the element
        computed by the recurrence, combinatorial formula, or wants to
        set `x_1` and `x_2` to be one.

        See [LeeLiZe]_ for more details.

        EXAMPLES::

            sage: S = ClusterSeed(['R2', [3, 3]])
            sage: S.greedy(4, 4)
            (x0^12 + x1^12 + 4*x0^9 + 4*x1^9 + 6*x0^6 + 4*x0^3*x1^3 + 6*x1^6 + 4*x0^3 + 4*x1^3 + 1)/(x0^4*x1^4)
            sage: S.greedy(4, 4, 'by_combinatorics')
            (x0^12 + x1^12 + 4*x0^9 + 4*x1^9 + 6*x0^6 + 4*x0^3*x1^3 + 6*x1^6 + 4*x0^3 + 4*x1^3 + 1)/(x0^4*x1^4)
            sage: S.greedy(4, 4, 'just_numbers')
            35
            sage: S = ClusterSeed(['R2', [2, 2]])
            sage: S.greedy(1, 2)
            (x0^4 + 2*x0^2 + x1^2 + 1)/(x0*x1^2)
            sage: S.greedy(1, 2, 'by_combinatorics')
            (x0^4 + 2*x0^2 + x1^2 + 1)/(x0*x1^2)

        REFERENCES:

        .. [LeeLiZe] Lee-Li-Zelevinsky, Greedy elements in rank 2
           cluster algebras, :arxiv:`1208.2391`
        """
        if self.b_matrix().dimensions() == (2, 2):
            b = abs(self.b_matrix()[0, 1])
            c = abs(self.b_matrix()[1, 0])
            if method == 'by_recursion':
                ans = self.x(0)**(-a1)*self.x(1)**(-a2)
                for p in range(max(a2, 0)+1):
                    for q in range(max(a1, 0)+1):
                        if p != 0 or q != 0:
                            ans += self._R(coeff_recurs(p, q, a1, a2, b, c))*self.x(0)**(b*p-a1)*self.x(1)**(c*q-a2)
                return(ans)
            elif method == 'by_combinatorics':
                if b == 0:
                    S = ClusterSeed([['A', 1], ['A', 1]])
                else:
                    S = ClusterSeed(['R2', [b, b]])
                ans = 0
                if a1 >= a2:
                    PS = PathSubset(a1, a2)
                elif a1 < a2:
                    PS = PathSubset(a2, a1)
                from sage.combinat.subset import Subsets
                for T in Subsets(PS):
                    if a1 >= a2:
                        if is_LeeLiZel_allowable(T, a1, a2, b, c):
                            oddT = set(T).intersection(PathSubset(a1, 0))
                            evenT = set(T).symmetric_difference(oddT)
                            ans = ans + S.x(0)**(b*len(evenT)) * S.x(1)**(c*len(oddT))
                    elif a1 < a2:
                        if is_LeeLiZel_allowable(T, a2, a1, b, c):
                            oddT = set(T).intersection(PathSubset(a2, 0))
                            evenT = set(T).symmetric_difference(oddT)
                            ans = ans + S.x(0)**(b*len(oddT)) * S.x(1)**(c*len(evenT))
                ans = ans*S.x(0)**(-a1)*S.x(1)**(-a2)
                return ans
            elif method == 'just_numbers':
                ans = 1
                for p in range(max(a2, 0)+1):
                    for q in range(max(a1, 0)+1):
                        if p != 0 or q != 0:
                            ans += coeff_recurs(p, q, a1, a2, b, c)
                return(ans)
            else:
                raise ValueError("The third input should be 'by_recursion', "
                                 "'by_combinatorics', or 'just_numbers'.")
        else:
            raise ValueError("Greedy elements are only currently "
                             "defined for cluster seeds of rank two.")

def _bino(n, k):
    """
    Binomial coefficient which we define as zero for negative n.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.cluster_seed import _bino
        sage: _bino(3, 2)
        3
        sage: _bino(-3, 2)
        0
    """
    if n >= 0:
        from sage.rings.arith import binomial
        return binomial(n, k)
    else:
        return 0

def coeff_recurs(p, q, a1, a2, b, c):
    """
    Coefficients in Laurent expansion of greedy element, as defined by recursion.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.cluster_seed import coeff_recurs
        sage: coeff_recurs(1, 1, 5, 5, 3, 3)
        10
    """
    if p == 0 and q == 0:
        return 1
    elif p < 0 or q < 0:
        return 0
    else:
        if c*a1*q <= b*a2*p:
            return sum((-1)**(k-1)*coeff_recurs(p-k, q, a1, a2, b, c)*_bino(a2-c*q+k-1, k)
                       for k in range(1, p+1))
        else:
            return sum((-1)**(k-1)*coeff_recurs(p, q-k, a1, a2, b, c)*_bino(a1-b*p+k-1, k)
                       for k in range(1, q+1))

def PathSubset(n,m):
    r"""
    Encodes a *maximal* Dyck path from (0,0) to (n,m) (for n >= m >= 0) as a subset of {0,1,2,..., 2n-1}.
    The encoding is given by indexing horizontal edges by odd numbers and vertical edges by evens.

    The horizontal between (i,j) and (i+1,j) is indexed by the odd number 2*i+1.
    The vertical between (i,j) and (i,j+1) is indexed by the even number 2*j.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.cluster_seed import PathSubset
        sage: PathSubset(4,0)
        {1, 3, 5, 7}
        sage: PathSubset(4,1)
        {1, 3, 5, 6, 7}
        sage: PathSubset(4,2)
        {1, 2, 3, 5, 6, 7}
        sage: PathSubset(4,3)
        {1, 2, 3, 4, 5, 6, 7}
        sage: PathSubset(4,4)
        {0, 1, 2, 3, 4, 5, 6, 7}
    """
    from sage.misc.misc import union
    from sage.functions.other import floor
    S = [ ]
    for i in range(n):
        S = union(S, [2*i+1])
    if m > 0:
        for j in range(n):
            if floor((j+1)*m/n) - floor(j*m/n) == 1:
                S = union(S, [2*j])
    return set(S)

def SetToPath(T):
    r"""
    Rearranges the encoding for a *maximal* Dyck path (as a set) so that it is a list in the proper order of the edges.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.cluster_seed import PathSubset
        sage: from sage.combinat.cluster_algebra_quiver.cluster_seed import SetToPath
        sage: SetToPath(PathSubset(4,0))
        [1, 3, 5, 7]
        sage: SetToPath(PathSubset(4,1))
        [1, 3, 5, 7, 6]
        sage: SetToPath(PathSubset(4,2))
        [1, 3, 2, 5, 7, 6]
        sage: SetToPath(PathSubset(4,3))
        [1, 3, 2, 5, 4, 7, 6]
        sage: SetToPath(PathSubset(4,4))
        [1, 0, 3, 2, 5, 4, 7, 6]
    """
    n = (max(T)+1) // 2
    ans = [1]
    for i in range(n-1):
        if 2*i in T:
            ans.append(2*i)
        ans.append(2*i+3)
    if 2*n-2 in T:
        ans.append(2*n-2)
    return ans

def is_LeeLiZel_allowable(T,n,m,b,c):
    """
    Check if the subset T contributes to the computation of the greedy
    element x[m,n] in the rank two (b,c)-cluster algebra.

    This uses the conditions of Lee-Li-Zelevinsky's paper [LeeLiZe]_.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.cluster_seed import is_LeeLiZel_allowable
        sage: is_LeeLiZel_allowable({1,3,2,5,7,6},4,2,6,6)
        False
        sage: is_LeeLiZel_allowable({1,2,5},3,3,1,1)
        True
    """
    horiz = set(T).intersection( PathSubset(n, 0))
    vert = set(T).symmetric_difference(horiz)
    if len(horiz) == 0 or len(vert) == 0:
        return True
    else:
        Latt = SetToPath(PathSubset(n, m))
        for u in horiz:
            from sage.combinat.words.word import Word
            from sage.modules.free_module_element import vector
            WW = Word(Latt)
            LattCycled = vector(WW.conjugate(Latt.index(u))).list()
            for v in vert:
                uv_okay = False
                for A in range(LattCycled.index(v)):
                    EA = []
                    AF = copy(LattCycled)
                    for i in range(LattCycled.index(v), len(LattCycled)-1):
                        AF.pop()
                    AF.reverse()
                    for i in range(A+1):
                        EA.append(LattCycled[i])
                        AF.pop()
                    AF.reverse()
                    nAF1 = 0
                    for i in range(len(AF)):
                        if AF[i] % 2 == 1:
                            nAF1 += 1
                    nAF2 = 0
                    for i in range(len(AF)):
                        if AF[i] % 2 == 0 and AF[i] in vert:
                            nAF2 += 1
                    nEA2 = 0
                    for i in range(len(EA)):
                        if EA[i] % 2 == 0:
                            nEA2 += 1
                    nEA1 = 0
                    for i in range(len(EA)):
                        if EA[i] % 2 == 1 and EA[i] in horiz:
                            nEA1 += 1
                    if nAF1 == b*nAF2 or nEA2 == c*nEA1:
                        uv_okay = True
                if uv_okay == False:
                        return False
        return True


class ClusterVariable(FractionFieldElement):
    r"""
    This class is a thin wrapper for cluster variables in cluster seeds.

    It provides the extra feature to store if a variable is frozen or not.

    - the associated positive root::

        sage: S = ClusterSeed(['A',3])
        sage: for T in S.variable_class_iter(): print T, T.almost_positive_root()
        x0 -alpha[1]
        x1 -alpha[2]
        x2 -alpha[3]
        (x1 + 1)/x0 alpha[1]
        (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2) alpha[1] + alpha[2] + alpha[3]
        (x1 + 1)/x2 alpha[3]
        (x0*x2 + x1 + 1)/(x0*x1) alpha[1] + alpha[2]
        (x0*x2 + 1)/x1 alpha[2]
        (x0*x2 + x1 + 1)/(x1*x2) alpha[2] + alpha[3]
    """
    def __init__( self, parent, numerator, denominator, coerce=True, reduce=True, mutation_type=None, variable_type=None ):
        r"""
        Initializes a cluster variable in the same way that elements in the field of rational functions are initialized.

        .. see also:: :class:`Fraction Field of Multivariate Polynomial Ring`

        TESTS::

            sage: S = ClusterSeed(['A',2])
            sage: for f in S.cluster():
            ...     print type(f)
            <class 'sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterVariable'>
            <class 'sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterVariable'>

            sage: S.variable_class()
            [(x0 + x1 + 1)/(x0*x1), (x1 + 1)/x0, (x0 + 1)/x1, x1, x0]
        """
        FractionFieldElement.__init__( self, parent, numerator, denominator, coerce=coerce, reduce=reduce )
        self._mutation_type = mutation_type
        self._variable_type = variable_type

    def almost_positive_root( self ):
        r"""
        Returns the *almost positive root* associated to ``self`` if ``self`` is of finite type.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: for T in S.variable_class_iter(): print T, T.almost_positive_root()
            x0 -alpha[1]
            x1 -alpha[2]
            x2 -alpha[3]
            (x1 + 1)/x0 alpha[1]
            (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2) alpha[1] + alpha[2] + alpha[3]
            (x1 + 1)/x2 alpha[3]
            (x0*x2 + x1 + 1)/(x0*x1) alpha[1] + alpha[2]
            (x0*x2 + 1)/x1 alpha[2]
            (x0*x2 + x1 + 1)/(x1*x2) alpha[2] + alpha[3]
        """
        if self._variable_type == 'frozen variable':
            raise ValueError('The variable is frozen.')
        if isinstance(self._mutation_type, str):
            raise ValueError('The cluster algebra for %s is not of finite type.'%self._repr_())
        else:
            if self._mutation_type is None:
                self._mutation_type = self.parent().mutation_type()
            if self._mutation_type.is_finite():
                from sage.combinat.root_system.root_system import RootSystem
                # the import above is used in the line below
                mt = self._mutation_type._repr_()
                # mt is a string of the shape "['A', 15]"
                # where A is a single letter and 15 is an integer
                Phi = RootSystem([mt[2: 3], ZZ(mt[6: -1])])
                Phiplus = Phi.root_lattice().simple_roots()
                if self.denominator() == 1:
                    return -Phiplus[ self.numerator().degrees().index(1) + 1 ]
                else:
                    root = self.denominator().degrees()
                    return sum( [ root[i]*Phiplus[ i+1 ] for i in range(len(root)) ] )
            else:
                raise ValueError('The cluster algebra for %s is not of finite type.'%self._repr_())
