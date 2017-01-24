# -*- coding: utf-8 -*-
r"""
ClusterSeed

A *cluster seed* is a pair `(B,\mathbf{x})` with `B` being a *skew-symmetrizable* `(n+m \times n)` *-matrix*
and with `\mathbf{x}` being an `n`-tuple of *independent elements* in the field of rational functions in `n` variables.

For the compendium on the cluster algebra and quiver package see
:arxiv:`1102.4844`.

AUTHORS:

- Gregg Musiker: Initial Version
- Christian Stump: Initial Version
- Aram Dermenjian (2015-07-01): Updating ability to not rely solely on clusters
- Jesse Levitt (2015-07-01): Updating ability to not rely solely on clusters

REFERENCES:

.. [BDP2013] Thomas Brüstle, Grégoire Dupont, Matthieu Pérotin
   *On Maximal Green Sequences*
   :arxiv:`1205.2050`

.. [FZ2007] Sergey Fomin, Andrei Zelevinsky
   "Cluster Algebras IV: coefficients"
   :arxiv:`0602259`

.. SEEALSO:: For mutation types of cluster seeds, see :meth:`sage.combinat.cluster_algebra_quiver.quiver_mutation_type.QuiverMutationType`. Cluster seeds are closely related to :meth:`sage.combinat.cluster_algebra_quiver.quiver.ClusterQuiver`.
"""

#*****************************************************************************
#       Copyright (C) 2011 Gregg Musiker <musiker@math.mit.edu>
#                          Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from six.moves import range

import time
from operator import pos
from sage.structure.sage_object import SageObject
from copy import copy
from sage.rings.all import QQ, infinity
from sage.rings.integer_ring import ZZ
from sage.rings.all import FractionField, PolynomialRing
from sage.rings.fraction_field_element import FractionFieldElement
from sage.sets.all import Set
from sage.graphs.digraph import DiGraph
from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import QuiverMutationType_Irreducible, QuiverMutationType_Reducible
from sage.combinat.cluster_algebra_quiver.mutation_type import is_mutation_finite
from sage.misc.misc import exists
from random import randint
from sage.misc.all import prod
from sage.matrix.all import identity_matrix
from sage.matrix.constructor import matrix
from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
from sage.rings.integer import Integer

from sage.misc.decorators import rename_keyword

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

        sage: S = ClusterSeed(['A',4]); S._use_fpolys
        True

        sage: S._use_d_vec
        True

        sage: S._use_g_vec
        True

        sage: S._use_c_vec
        True

        sage: S = ClusterSeed(['A',4]); S.use_fpolys(False); S._use_fpolys
        False

    """
    def __init__(self, data, frozen=None, is_principal=False, user_labels=None, user_labels_prefix='x'):
        r"""

        Initializes the ClusterSeed ``self`` with the following range of possible attributes:

            * self._n - the number of mutable elements of the cluster seed.
            * self._m - the number of immutable elements of the cluster seed.
            * self._M - the 'n + m' x 'n' exchange matrix associated to the cluster seed.
            * self._B - the mutable part of self._M.
            * self._b_initial - the initial exchange matrix
            * self._description - the description of the ClusterSeed
            * self._frozen - the number of frozen vertices, assumed to be the later vertices, when the input is a ClusterQuiver.
            * self._use_fpolys - a boolean tracking whether F-polynomials and cluster variables will be tracked as part of every mutation.
            * self._cluster - a list tracking the current names of cluster elements.
            * self._user_labels_prefix - the prefix for every named cluster element. Defaults to 'x'.
            * self._user_labels - an optional dictionary or list of user defined names for all cluster elements. Defaults to 'x_i' for mutable elements and 'y_i' for immutable elements.
            * self._init_vars - an internal list for defining ambient the algebraic setting and naming quiver vertices.
            * self._init_exch - the dictionary storing the initial mutable cluster variable names.
            * self._U - the coefficient tuple of the initial cluster seed.
            * self._F - the dictionary of F-polynomials.
            * self._R - the ambient polynomial ring.
            * self._y - the coefficient tuple for the current cluster seed.
            * self._yhat - the mixed coefficient tuple appearing in Proposition 3.9 of [FZ2007]
            * self._use_g_vec - a boolean stating if g-vectors for the cluster seed are being tracked. User input overridden as needed.
            * self._G - the matrix containing all g-vectors.
            * self._use_d_vec - a boolean stating if d-vectors for the cluster seed are being tracked.
            * self._D - the matrix containing all d-vectors.
            * self._bot_is_c - a boolean stating if the c-vectors are stored on the bottom of the exchange matrix M.
            * self._use_c_vec - a boolean stating if c-vectors for the cluster seed are being tracked. User input overridden as needed.
            * self._C - the matrix containing all c-vectors.
            * self._BC - an extended matrix involving the B and C matrices used for simplifying mutation calculations.
            * self._is_principal - a boolean tracking whether the ClusterSeed contains immutable elements coming from a principal extension of the mutable vertices. To be deprecated in future versions.

            * self._quiver - the ClusterQuiver corresponding to the exchange matrix self._M .
            * self._mutation_type - the mutation type of self._quiver .

            * self._track_mut - a boolean tracking whether the a ClusterSeed's mutation path is being recorded.
            * self._mut_path - the list of integers recording the mutation path of a seed - with consecutive repeats deleted since mutations is an involution.

        TESTS::

            sage: S = ClusterSeed(['A',4])
            sage: TestSuite(S).run()
        """
        from sage.matrix.matrix import Matrix

        #initialize a null state ClusterSeed object so all tests run and fail as appropriate.
        # numerous doctests if this null state is not first initialized.
        self._n = 0
        self._m = 0
        self._M = None
        self._B = None
        self._b_initial = None
        self._description = None
        self._frozen = 0
        self._use_fpolys = None
        self._cluster = None
        self._user_labels_prefix = None
        self._user_labels = None
        self._init_vars = None
        self._init_exch = None
        self._U = None
        self._F = None
        self._R = None
        self._y = None
        self._yhat = None

        self._use_g_vec = None
        self._G = None

        self._use_d_vec = None
        self._D = None

        self._bot_is_c = None
        self._use_c_vec = None
        self._C = None
        self._BC = None
        self._is_principal = None

        self._quiver = None
        self._mutation_type = None

        self._track_mut = None
        self._mut_path = None

        # constructs a cluster seed from a cluster seed
        if isinstance(data, ClusterSeed):
            if frozen:
                print("The input \'frozen\' is ignored")

            # Copy the following attributes from data
            self._M = copy( data._M )
            self._M.set_immutable()
            self._B = copy( data._B )
            self._n = data._n
            self._m = data._m

            # initialize matrix of g-vectors if desired and possible
            if data._use_g_vec and (data._G or data._cluster or (data._B.is_skew_symmetric() and data._C) or data._track_mut):
                self._G = data.g_matrix()

            # initialize matrix of c-vectors if desired and possible
            if data._use_c_vec and (data._C or (data._B.is_skew_symmetric() and (data._cluster or (data._use_g_vec and data._G)) or data._track_mut)):
                self._C = data.c_matrix()
                self._BC = copy(self._M).stack(copy(self._C))
            else:
                self._BC = copy(self._M)

            # initialize matrix of d-vectors if desired and possible
            if data._use_d_vec and (data._D or data._cluster or data._track_mut):
                self._D = data.d_matrix()

            self._cluster = copy( data._cluster)

            self._b_initial = copy( data._b_initial)

            self._mutation_type = copy( data._mutation_type)
            self._description = copy( data._description)
            self._quiver = ClusterQuiver( data._quiver ) if data._quiver else None

            # copy all previous booleans
            self._use_fpolys = data._use_fpolys
            self._use_g_vec = data._use_g_vec
            self._use_d_vec = data._use_d_vec
            self._bot_is_c = data._bot_is_c
            self._use_c_vec = data._use_c_vec
            self._track_mut = data._track_mut
            self._is_principal = data._is_principal

            # copy all previous dictionaries, names and data
            self._user_labels = copy( data._user_labels)
            self._user_labels_prefix = copy( data._user_labels_prefix)
            self._init_vars = copy( data._init_vars)
            self._init_exch = copy( data._init_exch)
            self._U = copy( data._U)
            self._F = copy( data._F)
            self._R = copy( data._R)
            self._y = copy( data._y)
            self._yhat = copy( data._yhat)
            self._mut_path = copy( data._mut_path)

        # constructs a cluster seed from a quiver
        elif isinstance(data, ClusterQuiver):
            if frozen:
                print("The input \'frozen\' is ignored")

            quiver = ClusterQuiver( data )

            self._M = copy(quiver._M)    # B-tilde exchange matrix
            self._M.set_immutable()
            self._n = quiver._n
            self._m = quiver._m
            self._B = copy(self._M[:self._n,:self._n])  # Square Part of the B_matrix

            # If initializing from a ClusterQuiver rather than a ClusterSeed, the initial B-matrix is reset to be the input B-matrix.
            self._b_initial = copy(self._M)
            self._mutation_type = copy(quiver._mutation_type)
            self._description = 'A seed for a cluster algebra of rank %d' %(self._n)
            self._quiver = quiver

            # We are now updating labels from user's most recent choice.
            self._is_principal = is_principal
            self._user_labels = user_labels
            self._user_labels_prefix = user_labels_prefix

            # initialize the rest
 
            self._C = matrix.identity(self._n)
            self._use_c_vec = True

            self._G = matrix.identity(self._n)
            self._use_g_vec = True

            self._BC = copy(self._M).stack(self.c_matrix())
            self._bot_is_c=False

            self._D = -matrix.identity(self._n)
            self._use_d_vec = True

            self._mut_path = [ ]
            self._track_mut = True

            if user_labels:
                self._sanitize_init_vars(user_labels, user_labels_prefix)
            else:
                xs = {i:'x%s'%i for i in range(self._n)}
                ys = {(i+self._n):'y%s'%i for i in range(self._n+self._m)}
                self._init_vars = copy(xs)
                self._init_vars.update(ys)

            self._init_exch = dict(self._init_vars.items()[:self._n])
            self._U = PolynomialRing(QQ,['y%s'%i for i in range(self._n)])
            self._F = dict([(i,self._U(1)) for i in self._init_exch.values()])
            self._R = PolynomialRing(QQ,[val for val in self._init_vars.values()])
            self._y = dict([ (self._U.gen(j),prod([self._R.gen(i)**self._M[i,j] for i in range(self._n,self._n+self._m)])) for j in range(self._n)])
            self._yhat = dict([ (self._U.gen(j),prod([self._R.gen(i)**self._M[i,j] for i in range(self._n+self._m)])) for j in range(self._n)])
            #self._cluster = None
            self._use_fpolys = True

        # in all other cases, we construct the corresponding ClusterQuiver first
        else:
            quiver = ClusterQuiver( data, frozen=frozen )
            self.__init__( quiver, is_principal=is_principal, user_labels=user_labels, user_labels_prefix=user_labels_prefix)

    def use_c_vectors(self, use=True, bot_is_c=False, force=False):
        r"""
        Reconstruct c vectors from other data or initialize if no usable data exists.

        Warning: Initialization may lead to inconsistent data.

        INPUT:

        - ``use`` -- (default:True) If True, will use c vectors
        - ``bot_is_c`` -- (default:False) If True and ClusterSeed self has self._m == self._n, then will assume bottom half of the extended exchange matrix is the c-matrix. If true, lets the ClusterSeed know c-vectors can be calculated.

        EXAMPLES::

            sage: S = ClusterSeed(['A',4]);
            sage: S.use_c_vectors(False); S.use_g_vectors(False); S.use_fpolys(False); S.track_mutations(False)
            sage: S.use_c_vectors(True)
            Warning: Initializing c-vectors at this point could lead to inconsistent seed data.

            sage: S.use_c_vectors(True, force=True)
            sage: S.c_matrix()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: S = ClusterSeed(['A',4]);
            sage: S.use_c_vectors(False); S.use_g_vectors(False); S.use_fpolys(False); S.track_mutations(False)
            sage: S.mutate(1);
            sage: S.use_c_vectors(True, force=True)
            sage: S.c_matrix()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

        """
        if self._use_c_vec != use:
            self._use_c_vec = use
            if self._use_c_vec:
                #self._C = matrix.identity(self._n)
                try:
                    self._use_c_vec = False   # temporarily turns off c-vectors to see if they can be recovered.
                    self._C = self.c_matrix()    # if not just sets it to be identity matrix, i.e. reinitialized.
                    self._BC = copy(self._M).stack(self.c_matrix())
                    self._use_c_vec = True
                except ValueError:
                    if not force:
                        print("Warning: Initializing c-vectors at this point could lead to inconsistent seed data.")
                    else:
                        self._use_c_vec = True
                        self._C = matrix.identity(self._n)
                        self._BC = copy(self._M).stack(self.c_matrix())
                except AttributeError:
                    if not force:
                        print("Warning: Initializing c-vectors at this point could lead to inconsistent seed data.")
                    else:
                        self._use_c_vec = True
                        self._C = matrix.identity(self._n)
                        self._BC = copy(self._M).stack(self.c_matrix())
            else:
                self._C = None
                self._BC = copy(self._M)
        if self._bot_is_c != bot_is_c: # If we need to do this. It overrides the previous designations.
            self._bot_is_c = bot_is_c
            if self._bot_is_c:
                self._use_c_vec = True
                if self._m == self._n: # in this case, the second half of a 2n x n matrix is a c-matrix.
                    self._C = copy(self._M[self._n:(self._n+self._m),:self._n])
                    self._BC = copy(self._M)
                else: # self._n != self._m
                    raise ValueError('There are immutable elements not in the c-matrix. Storing the c-matrix separately.')
                    self._C = copy(self._M[self._m:(self._n+self._m),:self._n])
                    self._BC = copy(self._M)
                    self._M = self._M[:self._m:self._n]
                    self._M.set_immutable()
                    self._bot_is_c = False

    def use_g_vectors(self, use=True, force=False):
        r"""
        Reconstruct g vectors from other data or initialize if no usable data exists.

        Warning: Initialization may lead to inconsistent data.

        INPUT:

        - ``use`` -- (default:True) If True, will use g vectors

        EXAMPLES::

            sage: S = ClusterSeed(['A',4]);
            sage: S.use_g_vectors(False); S.use_fpolys(False)
            sage: S.use_g_vectors(True)
            sage: S.g_matrix()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: S = ClusterSeed(['A',4]);
            sage: S.use_g_vectors(False); S.use_fpolys(False)
            sage: S.mutate(1);
            sage: S.use_g_vectors(True)
            sage: S.g_matrix()
            [ 1  0  0  0]
            [ 0 -1  0  0]
            [ 0  0  1  0]
            [ 0  0  0  1]

            sage: S = ClusterSeed(['A',4]);
            sage: S.use_g_vectors(False); S.use_fpolys(False); S.track_mutations(False)
            sage: S.mutate(1);
            sage: S.use_c_vectors(False)
            sage: S.g_matrix()
            Traceback (most recent call last):
            ...
            ValueError: Unable to calculate g-vectors. Need to use g vectors.

            sage: S = ClusterSeed(['A',4]);
            sage: S.use_g_vectors(False); S.use_fpolys(False); S.track_mutations(False)
            sage: S.mutate(1);
            sage: S.use_c_vectors(False)
            sage: S.use_g_vectors(True)
            Warning: Initializing g-vectors at this point could lead to inconsistent seed data.

            sage: S.use_g_vectors(True, force=True)
            sage: S.g_matrix()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        if self._use_g_vec != use:
            self._use_g_vec = use
            if self._use_g_vec:
                #self._G = matrix.identity(self._n) if self._use_g_vec else None
                try:
                    self._use_g_vec = False   # temporarily turns off g-vectors to see if they can be recovered.
                    self._G = self.g_matrix()    # if not just sets it to be identity matrix, i.e. reinitialized.
                    self._use_g_vec = True
                except ValueError:
                    if not force:
                        print("Warning: Initializing g-vectors at this point could lead to inconsistent seed data.")
                    else:
                        self._use_g_vec = True
                        self._G = matrix.identity(self._n)
                except AttributeError:
                    if not force:
                        print("Warning: Initializing g-vectors at this point could lead to inconsistent seed data.")
                    else:
                        self._use_g_vec = True
                        self._G = matrix.identity(self._n)
            else:
                self._G = None

            # Initially coded so c_vectors would be turned back on but now each of these boolean flags are independent
            #if self._use_g_vec and not self._use_c_vec:
            #    self.use_c_vectors(True)

    def use_d_vectors(self, use=True, force=False):
        r"""
        Reconstruct d vectors from other data or initialize if no usable data exists.

        Warning: Initialization may lead to inconsistent data.

        INPUT:

        - ``use`` -- (default:True) If True, will use d vectors

        EXAMPLES::

            sage: S = ClusterSeed(['A',4]);
            sage: S.use_d_vectors(True)
            sage: S.d_matrix()
            [-1  0  0  0]
            [ 0 -1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]

            sage: S = ClusterSeed(['A',4]); S.use_d_vectors(False); S.track_mutations(False); S.mutate(1); S.d_matrix()
            [-1  0  0  0]
            [ 0  1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]
            sage: S.use_fpolys(False)
            sage: S.d_matrix()
            Traceback (most recent call last):
            ...
            ValueError: Unable to calculate d-vectors. Need to use d vectors.

            sage: S = ClusterSeed(['A',4]); S.use_d_vectors(False); S.track_mutations(False); S.mutate(1); S.d_matrix()
            [-1  0  0  0]
            [ 0  1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]
            sage: S.use_fpolys(False)
            sage: S.use_d_vectors(True)
            Warning: Initializing d-vectors at this point could lead to inconsistent seed data.

            sage: S.use_d_vectors(True, force=True)
            sage: S.d_matrix()
            [-1  0  0  0]
            [ 0 -1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]

            sage: S = ClusterSeed(['A',4]); S.mutate(1); S.d_matrix()
            [-1  0  0  0]
            [ 0  1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]
            sage: S = ClusterSeed(['A',4]); S.use_d_vectors(True); S.mutate(1); S.d_matrix()
            [-1  0  0  0]
            [ 0  1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]
        """
        if self._use_d_vec != use:
            self._use_d_vec = use
            if self._use_d_vec:
                try:
                    self._use_d_vec = False # temporarily turns off d-vectors to see if they can be recovered.
                    self._D = self.d_matrix()
                    self._use_d_vec = True
                except ValueError:
                    if not force:
                        print("Warning: Initializing d-vectors at this point could lead to inconsistent seed data.")
                    else:
                        self._use_d_vec = True  # if not just sets it to be negative identity matrix, i.e. reinitialized.
                        self._D = -matrix.identity(self._n)
                except AttributeError:
                    if not force:
                        print("Warning: Initializing d-vectors at this point could lead to inconsistent seed data.")
                    else:
                        self._use_d_vec = True  # if not just sets it to be negative identity matrix, i.e. reinitialized.
                        self._D = -matrix.identity(self._n)
            else:
                self._D = None

    def use_fpolys(self, use=True, user_labels=None, user_labels_prefix=None):
        r"""
        Use F-polynomials in our Cluster Seed

        Note: This will automatically try to recompute the cluster variables if possible

        INPUT:

        - ``use`` -- (default:True) If True, will use F-polynomials
        - ``user_labels`` -- (default:None) If set will overwrite the default cluster variable labels
        - ``user_labels_prefix`` -- (default:None) If set will overwrite the default

        EXAMPLES::

            sage: S = ClusterSeed(['A',4]); S.use_fpolys(False); S._cluster
            sage: S.use_fpolys(True)
            sage: S.cluster()
            [x0, x1, x2, x3]

            sage: S = ClusterSeed(['A',4]); S.use_fpolys(False); S.track_mutations(False); S.mutate(1)
            sage: S.use_fpolys(True)
            Traceback (most recent call last):
            ...
            ValueError: F-polynomials and Cluster Variables cannot be reconstructed from given data.
            sage: S.cluster()
            Traceback (most recent call last):
            ...
            ValueError: Clusters not being tracked

        """
        if user_labels:
            self._user_labels = user_labels
        if user_labels_prefix:
            self._user_labels_prefix = user_labels_prefix
        if self._use_fpolys != use:
            self._use_fpolys = use

            if self._use_fpolys:

                if user_labels:
                    self._sanitize_init_vars(user_labels, user_labels_prefix)
                else:
                    xs = {i:'x%s'%i for i in range(self._n)}
                    ys = {(i+self._n):'y%s'%i for i in range(self._n+self._m)}
                    self._init_vars = copy(xs)
                    self._init_vars.update(ys)

                if self._G == matrix.identity(self._n): # If we are at the root
                    if not self._use_g_vec:
                        self.use_g_vectors(True)
                    self._init_exch = dict(self._init_vars.items()[:self._n])
                    self._U = PolynomialRing(QQ,['y%s'%i for i in range(self._n)])
                    self._F = dict([(i,self._U(1)) for i in self._init_exch.values()])
                    self._R = PolynomialRing(QQ,[val for val in self._init_vars.values()])
                    self._y = dict([ (self._U.gen(j),prod([self._R.gen(i)**self._M[i,j] for i in range(self._n,self._n+self._m)])) for j in range(self._n)])
                    self._yhat = dict([ (self._U.gen(j),prod([self._R.gen(i)**self._M[i,j] for i in range(self._n+self._m)])) for j in range(self._n)])
                elif self._cluster:
                    raise ValueError("should not be possible to have cluster variables without f-polynomials")    # added this as a sanity check.  This error should never appear however.
                elif self._track_mut: # If we can navigate from the root to where we are
                    if not self._use_g_vec:
                        self.use_g_vectors(True)
                    catchup = ClusterSeed(self._b_initial, user_labels=user_labels, user_labels_prefix=user_labels_prefix)
                    catchup.use_c_vectors(use=self._use_c_vec,bot_is_c=self._bot_is_c)
                    catchup.mutate(self.mutations())

                    self._init_exch = catchup._init_exch
                    self._U = catchup._U
                    self._F = catchup._F
                    self._R = catchup._R
                    self._y = catchup._y
                    self._yhat = catchup._yhat
                else:
                    self._use_fpolys = False
                    self._cluster = None
                    raise ValueError("F-polynomials and Cluster Variables cannot be reconstructed from given data.")

                # since we have F polynomials, set up clusters properly
                self._cluster = None
                self.cluster()
            else:
                if user_labels:
                    print("Warning: since 'use_fpolys' is False, the parameter 'user_labels' is ignored.")
                self._init_vars = None
                self._init_exch = None
                self._U = None
                self._F = None
                self._R = None
                self._y = None
                self._yhat = None
                self._cluster = None

    def track_mutations(self, use=True):
        r"""
        Begins tracking the mutation path.

        Warning: May initialize all other data to ensure that all c, d, and g vectors agree on the start of mutations.

        INPUT:

        - ``use`` -- (default:True) If True, will begin filling the mutation path

        EXAMPLES::

            sage: S = ClusterSeed(['A',4]); S.track_mutations(False)
            sage: S.mutate(0)
            sage: S.mutations()
            Traceback (most recent call last):
            ...
            ValueError: Not recording mutation sequence.  Need to track mutations.
            sage: S.track_mutations(True)
            sage: S.g_matrix()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: S.mutate([0,1])
            sage: S.mutations()
            [0, 1]
        """

        if self._track_mut != use:
            self._track_mut = use
            if self._track_mut:
                self._b_initial = self.b_matrix()
                self._mut_path = []
                use_c = self._use_c_vec
                use_d = self._use_d_vec
                use_g = self._use_g_vec
                # use_f = self._use_fpolys    #### also reinitialize F polynomials? -J
                # if use_f:
                #    self.use_fpolys(False)
                if use_g:
                    self.use_g_vectors(False)
                if use_c:
                    self.use_c_vectors(False)
                if use_d:
                    self.use_d_vectors(False)
                if use_g:
                    self.use_g_vectors()
                if use_c:
                    self.use_c_vectors()
                if use_d:
                    self.use_d_vectors()
                # if use_f:
                #    self.use_fpolys()
            else:
                self._mut_path = None

    def _sanitize_init_vars(self, user_labels, user_labels_prefix = 'x'):
        r"""
        Warning: This is an internal method that rewrites a user-given set of cluster variable names into a format that Sage can utilize.

        INPUT:

        - ``user_labels`` -- The labels that need sanitizing
        - ``user_labels_prefix`` -- (default:'x') The prefix to use for labels if integers given for labels
 
        EXAMPLES::

            sage: S = ClusterSeed(['A',4]); S._init_vars
            {0: 'x0', 1: 'x1', 2: 'x2', 3: 'x3', 4: 'y0', 5: 'y1', 6: 'y2', 7: 'y3'}
            sage: S._sanitize_init_vars([1,2,3,4],'z')
            sage: S._init_vars
            {0: 'z1', 1: 'z2', 2: 'z3', 3: 'z4'}

            sage: S = ClusterSeed(['A',4]); S._init_vars
            {0: 'x0', 1: 'x1', 2: 'x2', 3: 'x3', 4: 'y0', 5: 'y1', 6: 'y2', 7: 'y3'}
            sage: S._sanitize_init_vars(['a', 'b', 'c', 'd'])
            sage: S._init_vars
            {0: 'a', 1: 'b', 2: 'c', 3: 'd'}
        """
        if isinstance(user_labels,list):
            self._init_vars = {}
            for i in range(len(user_labels)):
                if isinstance(user_labels[i], Integer):
                    self._init_vars[i] = user_labels_prefix+user_labels[i].str()
                elif isinstance(user_labels[i], list):
                    self._user_labels_prefix = user_labels_prefix
                    strng = self._user_labels_prefix
                    for j in user_labels[i]:
                        if isinstance(j, Integer):
                            strng = strng+"_"+j.str()
                        else:
                            strng = strng+"_"+j
                    self._init_vars[i] = strng
                else:
                    self._init_vars[i] = user_labels[i]
        elif isinstance(user_labels,dict):
            self._init_vars = user_labels
        else:
            raise ValueError("The input 'user_labels' must be a dictionary or a list.")
        if len(self._init_vars.keys()) != self._n+self._m:
            raise ValueError("The number of user-defined labels is not the number of exchangeable and frozen variables.")

    def set_c_matrix(self, data):
        r"""
        Will force set the c matrix according to a matrix, a quiver, or a seed.

        INPUT:

        - ``data`` -- The matrix to set the c matrix to.  Also allowed
          to be a quiver or cluster seed, in which case the b_matrix
          is used.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]);
            sage: X = matrix([[0,0,1],[0,1,0],[1,0,0]])
            sage: S.set_c_matrix(X)
            sage: S.c_matrix()
            [0 0 1]
            [0 1 0]
            [1 0 0]

            sage: Y = matrix([[-1,0,1],[0,1,0],[1,0,0]])
            sage: S.set_c_matrix(Y)
            C matrix does not look to be valid - there exists a column containing positive and negative entries.
            Continuing...

            sage: Z = matrix([[1,0,1],[0,1,0],[2,0,2]])
            sage: S.set_c_matrix(Z)
            C matrix does not look to be valid - not a linearly independent set.
            Continuing...
        """
        if isinstance(data, ClusterQuiver):
            data = data.b_matrix()
        if isinstance(matrix, ClusterSeed):
            data=data.b_matrix()

        if data.determinant() == 0:
            print("C matrix does not look to be valid - not a linearly independent set.")
            print("Continuing...")

        # Do a quick check to make sure that each column is either all positive or all negative.
        # Can do this through green/red vertices
        greens = Set(get_green_vertices(data))
        reds = Set(get_red_vertices(data))
        if greens.intersection(reds).cardinality() > 0 or greens.union(reds).cardinality() < data.ncols():
            print("C matrix does not look to be valid - there exists a column containing positive and negative entries.")
            print("Continuing...")

        self._C = data
        self._BC = copy(self._M.stack(self._C))

    def __eq__(self, other):
        r"""
        Returns True iff ``self`` represent the same cluster seed as ``other`` and all tracked data agrees.

        EXAMPLES::

            sage: S = ClusterSeed(['A',5])
            sage: T = S.mutate( 2, inplace=False )
            sage: S.__eq__( T )
            False

            sage: T.mutate( 2 )
            sage: S.__eq__( T )
            True

            sage: S = ClusterSeed(['A',2])
            sage: T = ClusterSeed(S)
            sage: S.__eq__( T )
            True

            sage: S.mutate([0,1,0,1,0])
            sage: S.__eq__( T )
            False
            sage: S.cluster()
            [x1, x0]
            sage: T.cluster()
            [x0, x1]

            sage: S.mutate([0,1,0,1,0])
            sage: S.__eq__( T )
            True
            sage: S.cluster()
            [x0, x1]
        """
        if not isinstance(other, ClusterSeed):
            return False
        g_vec = True
        c_vec = True
        d_vec = True
        clusters = True
        ExMat = self._M == other._M
        if self._use_fpolys and other._use_fpolys:
            clusters = self.cluster() == other.cluster() and self.ground_field() == other.ground_field()
        elif self._use_g_vec and other._use_g_vec:
            g_vec = self.g_matrix() == other.g_matrix()
        if self._use_c_vec and other._use_c_vec:
            c_vec = self.c_matrix() == other.c_matrix()
        if self._use_d_vec and other._use_d_vec:
            d_vec = self.d_matrix() == other.d_matrix()
        return g_vec and c_vec and d_vec and clusters and ExMat

    def __hash__(self):
        """
        Return a hash of ``self``.

        EXAMPLES::

            sage: Q = ClusterSeed(['A',5])
            sage: hash(Q)  # indirect doctest
            -5649412990944896369  # 64-bit
            222337679  # 32-bit
        """
        # mat_hash = self._M.__hash__()
        if self._use_fpolys:
            return hash(tuple(self.cluster()))
        elif self._use_g_vec:
            return hash(self.g_matrix())
        elif self._use_c_vec:
            return hash(self.c_matrix())
        elif self._use_d_vec:
            return hash(self.d_matrix())

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

    def plot(self, circular=False, mark=None, save_pos=False, force_c =False, with_greens=False, add_labels = False):
        r"""
        Returns the plot of the quiver of ``self``.

        INPUT:

        - ``circular`` -- (default:False) if True, the circular plot is chosen, otherwise >>spring<< is used.
        - ``mark`` -- (default: None) if set to i, the vertex i is highlighted.
        - ``save_pos`` -- (default:False) if True, the positions of the vertices are saved.
        - ``force_c`` -- (default:False) if True, will show the frozen vertices even if they were never initialized
        - ``with_greens`` -- (default:False) if True, will display the green vertices in green
        - ``add_labels`` -- (default:False) if True, will use the initial variables as labels

        EXAMPLES::

            sage: S = ClusterSeed(['A',5])
            sage: pl = S.plot()
            sage: pl = S.plot(circular=True)
        """

        greens = []
        if with_greens:
            greens = self.green_vertices()

        if force_c:
            quiver = ClusterQuiver(self._BC)
        elif add_labels:
            # relabelling multiple times causes errors, so we need to always do it in place
            quiver = self.quiver().relabel(self._init_vars, inplace=True)
        else:
            quiver = self.quiver()

        return quiver.plot(circular=circular,mark=mark,save_pos=save_pos, greens=greens)

    def show(self, fig_size=1, circular=False, mark=None, save_pos=False, force_c = False, with_greens= False, add_labels = False):
        r"""
        Shows the plot of the quiver of ``self``.

        INPUT:

        - ``fig_size`` -- (default: 1) factor by which the size of the plot is multiplied.
        - ``circular`` -- (default: False) if True, the circular plot is chosen, otherwise >>spring<< is used.
        - ``mark`` -- (default: None) if set to i, the vertex i is highlighted.
        - ``save_pos`` -- (default:False) if True, the positions of the vertices are saved.
        - ``force_c`` -- (default:False) if True, will show the frozen vertices even if they were never initialized
        - ``with_greens`` -- (default:False) if True, will display the green vertices in green
        - ``add_labels`` -- (default:False) if True, will use the initial variables as labels

        TESTS::

            sage: S = ClusterSeed(['A',5])
            sage: S.show() # long time
        """

        greens = []
        if with_greens:
            greens = self.green_vertices()

        if force_c:
            quiver = ClusterQuiver(self._BC)
        elif add_labels:
            # relabelling multiple times causes errors, so we need to always do it in place
            quiver = self.quiver().relabel(self._init_vars, inplace=True)
        else:
            quiver = self.quiver()

        quiver.show(fig_size=fig_size, circular=circular,mark=mark,save_pos=save_pos, greens=greens)

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
        # Also update so works in cloud and not just notebook
        from sage.plot.plot import EMBEDDED_MODE
        from sagenb.notebook.interact import interact, selector
        from sage.misc.all import html,latex
        from sage.repl.rich_output.pretty_print import pretty_print

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
            def player(k=selector(values=list(range(self._n)), nrows = 1,
                                  label='Mutate at: '),
                       show_seq=("Mutation sequence:", True),
                       show_vars=("Cluster variables:", True),
                       show_matrix=("B-Matrix:", True),
                       show_lastmutation=("Show last mutation:", True)):
                ft, ss, sv, sm, sl = sft.pop(), sss.pop(), ssv.pop(), ssm.pop(), ssl.pop()
                if ft:
                    self.show(fig_size=fig_size, circular=circular)
                elif show_seq is not ss or show_vars is not sv or show_matrix is not sm or show_lastmutation is not sl:
                    if seq and show_lastmutation:
                        self.show(fig_size=fig_size, circular=circular,
                                  mark=seq[len(seq) - 1])
                    else:
                        self.show(fig_size=fig_size, circular=circular )
                else:
                    self.mutate(k)
                    seq.append(k)
                    if not show_lastmutation:
                        self.show(fig_size=fig_size, circular=circular)
                    else:
                        self.show(fig_size=fig_size, circular=circular, mark=k)
                sft.append(False)
                sss.append(show_seq)
                ssv.append(show_vars)
                ssm.append(show_matrix)
                ssl.append(show_lastmutation)
                if show_seq:
                    pretty_print(html("Mutation sequence: $" + str( [ seq[i] for i in range(len(seq)) ] ).strip('[]') + "$"))
                if show_vars:
                    pretty_print(html("Cluster variables:"))
                    table = "$\\begin{align*}\n"
                    for i in range(self._n):
                        table += "\tv_{%s} &= " % i + latex(self.cluster_variable(i)) + "\\\\ \\\\\n"
                    table += "\\end{align*}$"
                    pretty_print(html("$ $"))
                    pretty_print(html(table))
                    pretty_print(html("$ $"))
                if show_matrix:
                    pretty_print(html("B-Matrix:"))
                    pretty_print(html(self._M))
                    pretty_print(html("$ $"))

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
        return copy(self._M)

    def ground_field(self):
        r"""
        Returns the *ground field* of the cluster of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.ground_field()
            Multivariate Polynomial Ring in x0, x1, x2, y0, y1, y2 over Rational Field
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

        if self._use_fpolys and k in range(self._n):
            x = self._R.gens()[k]
            return ClusterVariable(FractionField(self._R), x.numerator(), x.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable' ,xdim=self._n)
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

        if self._use_fpolys and k in range(self._m):
            x = self._R.gens()[self._n+k]
            return ClusterVariable( FractionField(self._R), x.numerator(), x.denominator(), mutation_type=self._mutation_type, variable_type='frozen variable',xdim=self._n )
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

    def mutations(self):
        r"""
        Returns the list of mutations ``self`` has undergone if they are being tracked.

        Examples::

            sage: S = ClusterSeed(['A',3])
            sage: S.mutations()
            []

            sage: S.mutate([0,1,0,2])
            sage: S.mutations()
            [0, 1, 0, 2]

            sage: S.track_mutations(False)
            sage: S.mutations()
            Traceback (most recent call last):
            ...
            ValueError: Not recording mutation sequence.  Need to track mutations.
        """
        if self._track_mut:
            return copy(self._mut_path)
        else:
            raise ValueError("Not recording mutation sequence.  Need to track mutations.")

    def cluster_variable(self, k):
        r"""
        Generates a cluster variable using F-polynomials

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.mutate([0,1])
            sage: S.cluster_variable(0)
            (x1 + 1)/x0
            sage: S.cluster_variable(1)
            (x0*x2 + x1 + 1)/(x0*x1)
        """
        if self._use_fpolys:
            IE = self._init_exch.values()
            if (k in range(self._n)) or (k in IE):
                if k in range(self._n):
                    pass
                elif k in IE:
                    k = IE.index(k)

                g_mon = prod([self._R.gen(i)**self._G[i,k] for i in range(self._n)])
                F_num = self._F[IE[k]].subs(self._yhat)
                F_den = self._R(self._F[IE[k]].subs(self._y).denominator())
                cluster_variable = g_mon*F_num*F_den

                return ClusterVariable(FractionField(self._R), cluster_variable.numerator(), cluster_variable.denominator(), mutation_type=self._mutation_type,  variable_type='cluster variable',xdim=self._n)
            else:
                raise ValueError('No cluster variable with index or label ' + str(k) + '.')
        elif self._track_mut: # if we can recreate the clusters
            catchup = ClusterSeed(self._b_initial, user_labels=self._user_labels, user_labels_prefix=self._user_labels_prefix)
            catchup.use_c_vectors(use=self._use_c_vec, bot_is_c=self._bot_is_c)
            catchup.mutate(self.mutations())
            return catchup.cluster_variable(k)
        else:
            raise ValueError('Clusters not being tracked')
            return None

    def cluster(self):
        r"""
        Returns a copy of the *cluster* of ``self``.

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

        if not self._use_fpolys:
            if self._track_mut: # if we can recreate the clusters
                catchup = ClusterSeed(self._b_initial, user_labels=self._user_labels, user_labels_prefix=self._user_labels_prefix)
                catchup.use_c_vectors(use=self._use_c_vec, bot_is_c=self._bot_is_c)
                catchup.mutate(self.mutations())
                return catchup.cluster()
            else:
                raise ValueError('Clusters not being tracked')
        elif self._cluster is None:
            self._cluster = [self.cluster_variable(k) for k in range(self._n)]
        return copy(self._cluster)

    def _f_mutate( self, k):
        r"""
        An internal procedure that returns ``self`` with F-polynomials mutated at k.

        WARNING: This function assumes you are sending it good data

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S._f_mutate(0)
            sage: S.f_polynomial(0)
            y0 + 1
        """
        if self._use_fpolys:
            IE = self._init_exch.values()
        else:
            IE = []

        F = self._F
        B = self.b_matrix()
        C = self.c_matrix()

        # F-polynomials
        pos = self._U(1)
        neg = self._U(1)

        for j in range(self._n):
            if C[j,k] > 0:
                pos *= self._U.gen(j)**C[j,k]
            else:
                neg *= self._U.gen(j)**(-C[j,k])
            if B[j,k] > 0:
                pos *= F[IE[j]]**B[j,k]
            else:
                neg *= F[IE[j]]**(-B[j,k])

        # can the following be improved?
        self._F[IE[k]] = (pos+neg)//F[IE[k]]

    def f_polynomial(self,k):
        r"""
        Returns the ``k``-th *F-polynomial* of ``self``. It is obtained from the
        ``k``-th cluster variable by setting all `x_i` to `1`.

        Warning: this method assumes the sign-coherence conjecture and that the
        input seed is sign-coherent (has an exchange matrix with columns of like signs).
        Otherwise, computational errors might arise.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: [S.f_polynomial(k) for k in range(3)]
            [1, y1*y2 + y2 + 1, y1 + 1]

            sage: S = ClusterSeed(Matrix([[0,1],[-1,0],[1,0],[-1,1]])); S.use_c_vectors(bot_is_c=True); S
            A seed for a cluster algebra of rank 2 with 2 frozen variables
            sage: T = ClusterSeed(Matrix([[0,1],[-1,0]])).principal_extension(); T
            A seed for a cluster algebra of rank 2 with principal coefficients
            sage: S.mutate(0)
            sage: T.mutate(0)
            sage: S.f_polynomials()
            [y0 + y1, 1]
            sage: T.f_polynomials()
            [y0 + 1, 1]
        """
        if self._use_fpolys:
            IE = self._init_exch.values()
            if k in range(self._n):
                pass
            elif k in IE:
                k = IE.index(k)
            else:
                raise ValueError("The cluster seed does not have a cluster variable of index %s."%k)

            return self._F[IE[k]]
        elif self._track_mut:
            catchup = ClusterSeed(self._b_initial, user_labels=self._user_labels, user_labels_prefix=self._user_labels_prefix)
            catchup.use_c_vectors(use=self._use_c_vec, bot_is_c=self._bot_is_c)
            catchup.mutate(self.mutations())

            return catchup.f_polynomial(k)
        else:
            raise ValueError("Turn on use_fpolys to get F polynomial %s."%k)

    def f_polynomials(self):
        r"""
        Returns all *F-polynomials* of ``self``. These are obtained from the
        cluster variables by setting all `x_i`'s to `1`.

        Warning: this method assumes the sign-coherence conjecture and that the
        input seed is sign-coherent (has an exchange matrix with columns of like signs).
        Otherwise, computational errors might arise.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: S.f_polynomials()
            [1, y1*y2 + y2 + 1, y1 + 1]
        """

        return [self.f_polynomial(i) for i in range(self._n)]

    def g_vector(self,k):
        r"""
        Returns the ``k``-th *g-vector* of ``self``. This is the degree vector
        of the ``k``-th cluster variable after setting all `y_i`'s to `0`.

        Warning: this method assumes the sign-coherence conjecture and that the
        input seed is sign-coherent (has an exchange matrix with columns of like signs).
        Otherwise, computational errors might arise.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: [ S.g_vector(k) for k in range(3) ]
            [(1, 0, 0), (0, 0, -1), (0, -1, 0)]
        """

        if not (self._is_principal or self._use_g_vec or (self._use_fpolys and self._cluster)):
            raise ValueError("Unable to calculate g-vectors. Need to use g vectors.")
        if k not in range(self._n):
            raise ValueError("The cluster seed does not have a cluster variable of index %s."%k)

        if self._use_g_vec: # This implies the g-matrix is maintained by the mutate function and will always be up to date
            return copy(self._G.column(k))
        elif self._use_fpolys and self._cluster:
            f = copy(self.cluster_variable(k))
            eval_dict = dict( [ ( self.y(i), 0 ) for i in range(self._m) ] )
            f0 = f.subs(eval_dict)
            d1 = f0.numerator().degrees()
            d2 = f0.denominator().degrees()
            return tuple( d1[i] - d2[i] for i in range(self._n) )
        else: # in the is_principal=True case
            try:
                # ensure that we cannot create a loop by calling g_matrix() here by filtering out loop causing conditions in the previous if-elif sections
                return self.g_matrix().column(k)
            except ValueError:
                raise ValueError("Unable to calculate g-vectors. Need to use g vectors.")

    def g_matrix(self, show_warnings=True):
        r"""
        Returns the matrix of all *g-vectors* of ``self``. These are the degree vectors
        of the cluster variables after setting all `y_i`'s to `0`.

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
            sage: S.mutate([0,1])
            sage: S.g_matrix()
            [-1 -1  0]
            [ 1  0  0]
            [ 0  0  1]

            sage: S = ClusterSeed(['A',4]); S.use_g_vectors(False); S.use_fpolys(False); S.g_matrix()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: S = ClusterSeed(['A',4])
            sage: S.use_g_vectors(False); S.use_c_vectors(False); S.use_fpolys(False); S.track_mutations(False); S.g_matrix()
            Traceback (most recent call last):
            ...
            ValueError: Unable to calculate g-vectors. Need to use g vectors.
        """

        from sage.matrix.all import matrix
        if self._use_g_vec:
            return copy(self._G)
        elif self._use_fpolys and self._cluster: # This only calls g_vector when it will not create a loop.
            return matrix( [ self.g_vector(k) for k in range(self._n) ] ).transpose()
        elif self._use_c_vec:
            if self.b_matrix().is_skew_symmetric():
                return copy(self._C).inverse().transpose()
            elif self._track_mut:
                BC1 = copy(self._b_initial[0:self._n])
                BC1 = -BC1.transpose()
                BC1 = BC1.stack(matrix.identity(self._n))
                seq = iter(self.mutations())
                for k in seq:
                    BC1.mutate(k)
                return copy(BC1[self._n:2*self._n]).inverse().transpose()
            else:
                raise ValueError("Unable to calculate g-vectors. Need to use g vectors.")
        elif self._track_mut:
            catchup = ClusterSeed(self._b_initial)
            catchup.use_fpolys(False)
            catchup.mutate(self.mutations())
            return catchup.g_matrix()
        elif show_warnings:
            raise ValueError("Unable to calculate g-vectors. Need to use g vectors.")
        else:
            return None

    def _g_mutate(self, k):
        r"""
        An internal procedure that returns ``self`` with g-vectors mutated at k.

        .. WARNING::

            This function assumes you are sending it good data.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S._g_mutate(0)
            sage: S.g_vector(0)
            (-1, 1, 0)

        REFERENCES:

        .. [NaZe2011] Tomoki Nakanishi and Andrei Zelevinsky
           *On Tropical Dualities In Cluster Algebras*
           :arxiv:`1101.3736v3`
        """
        from sage.matrix.all import identity_matrix

        if self._use_fpolys:
            IE = self._init_exch.values()
        else:
            IE = []

        B = self.b_matrix()
        C = self.c_matrix()

        # G-matrix
        J = identity_matrix(self._n)
        if any(x > 0 for x in C.column(k)):
            eps = +1
        else:
            eps = -1
        for j in range(self._n):
            J[j,k] += max(0, -eps*B[j,k])
        J[k,k] = -1
        self._G = self._G*J

    def c_vector(self,k):
        r"""
        Returns the ``k``-th *c-vector* of ``self``. It is obtained as the
        ``k``-th column vector of the bottom part of the ``B``-matrix of ``self``.

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
            (1, 0)

            sage: S = ClusterSeed(Matrix([[0,1],[-1,0],[1,0],[-1,1]])); S.use_c_vectors(bot_is_c=True); S
            A seed for a cluster algebra of rank 2 with 2 frozen variables
            sage: S.c_vector(0)
            (1, -1)

        """
        if k not in range(self._n):
            raise ValueError("The cluster seed does not have a c-vector of index %s."%k)
        if not (self._is_principal or self._use_c_vec):
            raise ValueError("Requires C vectors to use.")
        if self._use_c_vec:
            return self.c_matrix().column(k)
        else:
            return tuple( self._M[i,k] for i in range(self._n,self._n+self._m) )

    def c_matrix(self,show_warnings=True):
        r"""
        Returns all *c-vectors* of ``self``.

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

            sage: S = ClusterSeed(['A',4]);
            sage: S.use_g_vectors(False); S.use_fpolys(False);  S.use_c_vectors(False);  S.use_d_vectors(False); S.track_mutations(False); 
            sage: S.c_matrix()
            Traceback (most recent call last):
            ...
            ValueError: Unable to calculate c-vectors. Need to use c vectors.
        """

        if self._bot_is_c:
            return copy(self._M[self._m:(self._n+self._m),:self._n])
        elif self._use_c_vec:
            return copy(self._C)
        elif self._use_g_vec or self._use_fpolys: #both of these will populate g_matrix() successfully
            if self.b_matrix().is_skew_symmetric():
                return self.g_matrix().inverse().transpose()
            elif self._track_mut:
                BC1 = copy(self._b_initial[0:self._n])
                BC1 = BC1.stack(matrix.identity(self._n))
                seq = iter(self.mutations())
                for k in seq:
                    BC1.mutate(k)
                return copy(BC1[self._n:2*self._n])
            else:
                raise ValueError("Unable to calculate c-vectors. Need to use c vectors.")
        elif self._track_mut:
            BC1 = copy(self._b_initial[0:self._n])
            BC1 = BC1.stack(matrix.identity(self._n))
            seq = iter(self.mutations())
            for k in seq:
                BC1.mutate(k)
            return copy(BC1[self._n:2*self._n])
        elif show_warnings:
            raise ValueError("Unable to calculate c-vectors. Need to use c vectors.")
        else:
            return None

    def d_vector(self, k):
        r"""
        Returns the ``k``-th *d-vector* of ``self``. This is the exponent vector
        of the denominator of the ``k``-th cluster variable.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S.mutate([2,1,2])
            sage: [ S.d_vector(k) for k in range(3) ]
            [(-1, 0, 0), (0, 1, 1), (0, 1, 0)]
        """
        from sage.modules.free_module_element import vector

        if self._use_d_vec:
            return copy(self._D).column(k)
        elif self._use_fpolys:
            f = self.cluster_variable(k)
            if f in self._R.gens():
                return -vector(f.numerator().monomials()[0].exponents()[0][:self._n])
            return vector(f.denominator().monomials()[0].exponents()[0][:self._n])
        elif self._track_mut:
            catchup = ClusterSeed(self._b_initial)
            catchup.use_fpolys(False)
            catchup.use_g_vectors(False)
            catchup.use_c_vectors(False)

            catchup.mutate(self.mutations())
            return copy(catchup._D).column(k)
        else:
            raise ValueError("Unable to calculate d-vector %s. Need to use d vectors."%k)

    def d_matrix(self, show_warnings=True):
        r"""
        Returns the matrix of *d-vectors* of ``self``.

        EXAMPLES::
            sage: S = ClusterSeed(['A',4]); S.d_matrix()
            [-1  0  0  0]
            [ 0 -1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]
            sage: S.mutate([1,2,1,0,1,3]); S.d_matrix()
            [1 1 0 1]
            [1 1 1 1]
            [1 0 1 1]
            [0 0 0 1]


        """
        if not (self._use_d_vec or self._use_fpolys or self._track_mut):
            #raise ValueError("No d-vectors initialized.")
            raise ValueError("Unable to calculate d-vectors. Need to use d vectors.")
        if self._use_d_vec:
            return copy(self._D)
        elif self._use_fpolys:
            return matrix( [ self.d_vector(k) for k in range(self._n) ] ).transpose()
        elif self._track_mut:
            catchup = ClusterSeed(self._b_initial)
            catchup.use_fpolys(False)
            catchup.use_g_vectors(False)
            catchup.use_c_vectors(False)
            catchup.track_mutations(False)

            catchup.mutate(self.mutations())
            return catchup.d_matrix()
        elif show_warnings:
            raise ValueError("No valid way to calculate d-vectors")

    def _d_mutate(self, k):
        r"""
        An internal procedure that returns ``self`` with d-vectors mutated at k.

        WARNING: This function assumes you are sending it good data (does not check for sanitized inputs)

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: S._d_mutate(0)
            sage: S.d_matrix()
            [ 1  0  0]
            [ 0 -1  0]
            [ 0  0 -1]
            sage: S.d_vector(0)
            (1, 0, 0)

        """
        if self._use_fpolys:
            IE = self._init_exch.values()
        else:
            IE = []

        B = self.b_matrix()
        D = copy(self._D)
        dnew = copy(-D.column(k))
        dp = copy( dnew.parent().zero() )
        dn = copy( dnew.parent().zero() )
        dmax = copy( dnew.parent().zero() )

        for j in range(self._n):
            if B[j,k] >0:
                dp += B[j,k]*D.column(j)
            elif B[j,k] <0:
                dn -= B[j,k]*D.column(j)
        for i in range(self._n):
            dmax[i] = max(dp[i],dn[i])
        self._D.set_column(k,dnew+dmax)

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
        else:
            try:    # are c vectors being tracked?
                exp = self.c_vector(k)
            except Exception:    # if not try and reconstruct them
                exp = self.c_matrix().column(k)

            return prod( self.y(i)**exp[i] for i in range(self._m) )

    def coefficients(self):
        r"""
        Returns all *coefficients* of ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutate([2,1,2])
            sage: S.coefficients()
            [y0, 1/y2, 1/y1]
        """
        # exceptions are caught in the subroutine.
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

    def green_vertices(self):
        r"""
        Return the list of green vertices of ``self``.

        A vertex is defined to be green if its c-vector has all non-positive
        entries. More information on green vertices can be found at [BDP2013]_

        OUTPUT:

        The green vertices as a list of integers.

        EXAMPLES::

            sage: ClusterSeed(['A',3]).principal_extension().green_vertices()
            [0, 1, 2]

            sage: ClusterSeed(['A',[3,3],1]).principal_extension().green_vertices()
            [0, 1, 2, 3, 4, 5]
        """

        # Make sure we have c vectors
        if not self._use_c_vec:
            raise ValueError("Must use c vectors to grab the vertices.")

        return get_green_vertices(self._C)

    def first_green_vertex(self):
        r"""
        Return the first green vertex of ``self``.

        A vertex is defined to be green if its c-vector has all non-positive entries.
        More information on green vertices can be found at [BDP2013]_

        EXAMPLES::

            sage: ClusterSeed(['A',3]).principal_extension().first_green_vertex()
            0

            sage: ClusterSeed(['A',[3,3],1]).principal_extension().first_green_vertex()
            0
        """
        # Make sure we have c vectors
        if not self._use_c_vec:
            raise ValueError("Must use c vectors to grab the vertices.")

        greens = self.green_vertices()
        if len(greens) > 0:
            return greens[0]

        return None

    def red_vertices(self):
        r"""
        Return the list of red vertices of ``self``.

        A vertex is defined to be red if its c-vector has all non-negative entries.
        More information on red vertices can be found at [BDP2013]_.

        OUTPUT:

        The red vertices as a list of integers.

        EXAMPLES::

            sage: ClusterSeed(['A',3]).principal_extension().red_vertices()
            []

            sage: ClusterSeed(['A',[3,3],1]).principal_extension().red_vertices()
            []

            sage: Q = ClusterSeed(['A',[3,3],1]).principal_extension();
            sage: Q.mutate(1);
            sage: Q.red_vertices()
            [1]

        """
        # Make sure we have c vectors on
        if not self._use_c_vec:
            raise ValueError("Must use c vectors to grab the vertices.")

        return get_red_vertices(self._C)

    def first_red_vertex(self):
        r"""
        Return the first red vertex of ``self``.

        A vertex is defined to be red if its c-vector has all non-negative entries.
        More information on red vertices can be found at [BDP2013]_.

        EXAMPLES::

            sage: ClusterSeed(['A',3]).principal_extension().first_red_vertex()

            sage: ClusterSeed(['A',[3,3],1]).principal_extension().first_red_vertex()

            sage: Q = ClusterSeed(['A',[3,3],1]).principal_extension();
            sage: Q.mutate(1);
            sage: Q.first_red_vertex()
            1

        """
        # Make sure we have c vectors
        if not self._use_c_vec:
            raise ValueError("Must use c vectors to grab the vertices.")

        reds = self.red_vertices()
        if len(reds) > 0:
            return reds[0]

        return None

    def urban_renewals(self, return_first=False):
        r"""
        Return the list of the urban renewal vertices of ``self``.

        An urban renewal vertex is one in which there are two arrows pointing
        toward the vertex and two arrows pointing away.

        INPUT:

        - ``return_first`` -- (default:False) if True, will return the first urban renewal

        OUTPUT:

        A list of vertices (as integers)

        EXAMPLES::

            sage: G = ClusterSeed(['GR',[4,9]]); G.urban_renewals()
            [5, 6]
        """
        vertices = []
        for i in range(self._n):
            if self.quiver().digraph().in_degree(i) == 2 and self.quiver().digraph().out_degree(i) == 2:
                if return_first:
                    return i
                vertices.append(i)

        if return_first:
            return None
        return vertices

    def first_urban_renewal(self):
        r"""
        Return the first urban renewal vertex.

        An urban renewal vertex is one in which there are two arrows pointing
        toward the vertex and two arrows pointing away.

        EXAMPLES::

            sage: G = ClusterSeed(['GR',[4,9]]); G.first_urban_renewal()
            5
        """
        return self.urban_renewals(return_first=True)

    def highest_degree_denominator(self, filter=None):
        r"""
        Return the vertex of the cluster polynomial with highest degree in the denominator.

        INPUT:

        - ``filter`` - Filter should be a list or iterable

        OUTPUT:

        Returns an integer.

        EXAMPLES::

            sage: B = matrix([[0,-1,0,-1,1,1],[1,0,1,0,-1,-1],[0,-1,0,-1,1,1],[1,0,1,0,-1,-1],[-1,1,-1,1,0,0],[-1,1,-1,1,0,0]])
            sage: C = ClusterSeed(B).principal_extension(); C.mutate([0,1,2,4,3,2,5,4,3])
            sage: C.highest_degree_denominator()
            5
        """
        if filter is None:
            filter = list(range(len(self.cluster())))
        degree = 0
        vertex_to_mutate = []

        # if we have d vectors use those, else see if we have clusters
        if self._use_d_vec:
            for i in list(enumerate(self.d_matrix().columns())):
                if i[0] not in filter:
                    continue
                col = i[1]
                vertex = i[0]
                cur_vertex_degree = sum(col)
                if degree == cur_vertex_degree:
                    vertex_to_mutate.append(vertex)
                if degree < cur_vertex_degree:
                    degree = cur_vertex_degree
                    vertex_to_mutate = [vertex]
        elif self._use_fpolys:
            for i in list(enumerate(self.cluster())):
                if i[0] not in filter:
                    continue
                vari = i[1]
                vertex = i[0]
                denom = vari.denominator()
                cur_vertex_degree = denom.degree()
                if degree == cur_vertex_degree:
                    vertex_to_mutate.append(vertex)
                if degree < cur_vertex_degree:
                    degree = cur_vertex_degree
                    vertex_to_mutate = [vertex]


        return_key = randint(0,len(vertex_to_mutate) - 1)
        return vertex_to_mutate[return_key]

    def smallest_c_vector(self):
        r"""
        Return the vertex with the smallest c vector

        OUTPUT:

        Returns an integer.

        EXAMPLES::
            sage: B = matrix([[0,2],[-2,0]])
            sage: C = ClusterSeed(B).principal_extension();
            sage: C.mutate(0)
            sage: C.smallest_c_vector()
            0

        """
        min_sum = infinity
        vertex_to_mutate = []

        for i in list(enumerate(self.c_matrix().columns())):
            col = i[1]
            vertex=i[0]
            cur_vertex_sum = abs(sum(col))
            if min_sum == cur_vertex_sum:
                vertex_to_mutate.append(vertex)
            if min_sum > cur_vertex_sum:
                min_sum = cur_vertex_sum
                vertex_to_mutate = [vertex]

        return_key = randint(0,len(vertex_to_mutate) - 1)
        return vertex_to_mutate[return_key]

    def most_decreased_edge_after_mutation(self):
        r"""

        Return the vertex that will produce the least degrees after mutation

        EXAMPLES::

            sage: S = ClusterSeed(['A',5])
            sage: S.mutate([0,2,3,1,2,3,1,2,0,2,3])
            sage: S.most_decreased_edge_after_mutation()
            2

        """
        analysis = self.mutation_analysis(['edge_diff'])
        least_edge = infinity
        least_vertex = []
        for edge,edge_analysis in analysis.items():
            if least_edge == edge_analysis['edge_diff']:
                least_vertex.append(edge)
            if least_edge > edge_analysis['edge_diff']:
                least_edge = edge_analysis['edge_diff']
                least_vertex = [edge]

        # if we have one vertex, return it
        if len(least_vertex) == 1:
            return least_vertex[0]

        # if not then do a test based on which one currently has the highest degree
        return self.highest_degree_denominator(least_vertex)

    def most_decreased_denominator_after_mutation(self):
        r"""

        Return the vertex that will produce the most decrease in denominator degrees after mutation

        EXAMPLES::

            sage: S = ClusterSeed(['A',5])
            sage: S.mutate([0,2,3,1,2,3,1,2,0,2,3])
            sage: S.most_decreased_denominator_after_mutation()
            2

        """
        analysis = self.mutation_analysis(['d_matrix'])
        least_change = infinity
        least_vertex = []
        current_columns = [sum(i) for i in self.d_matrix().columns()]
        for vertex,edge_analysis in analysis.items():
            mutated_column = sum(edge_analysis['d_matrix'].column(vertex))

            diff = mutated_column - current_columns[vertex]
            if least_change == diff:
                least_vertex.append(vertex)
            if diff < least_change:
                least_change = diff
                least_vertex = [vertex]

        return_key = randint(0,len(least_vertex) - 1)
        return least_vertex[return_key]

    def mutate(self, sequence, inplace=True):
        r"""
        Mutates ``self`` at a vertex or a sequence of vertices.

        INPUT:

        - ``sequence`` -- a vertex of ``self``, an iterator of vertices of ``self``,
          a function which takes in the ClusterSeed and returns a vertex or an iterator of vertices,
          or a string representing a type of vertices to mutate.
        - ``inplace`` -- (default: True) if False, the result is returned, otherwise ``self`` is modified.

        Possible values for vertex types in ``sequence`` are:

        - ``"first_source"``: mutates at first found source vertex,
        - ``"sources"``: mutates at all sources,
        - ``"first_sink"``: mutates at first sink,
        - ``"sinks"``: mutates at all sink vertices,
        - ``"green"``: mutates at the first green vertex,
        - ``"red"``: mutates at the first red vertex,
        - ``"urban_renewal"`` or ``"urban"``: mutates at first urban renewal vertex,
        - ``"all_urban_renewals"`` or ``"all_urban"``: mutates at all urban renewal vertices.

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

            sage: Q = ClusterSeed(['A',3]);Q.b_matrix()
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]

            sage: Q.mutate('first_sink');Q.b_matrix()
            [ 0 -1  0]
            [ 1  0  1]
            [ 0 -1  0]

            sage: def last_vertex(self): return self._n - 1
            sage: Q.mutate(last_vertex); Q.b_matrix()
            [ 0 -1  0]
            [ 1  0 -1]
            [ 0  1  0]

            sage: S = ClusterSeed(['A',4], user_labels=['a','b','c','d']);
            sage: S.mutate('a'); S.mutate('(b+1)/a')
            sage: S.cluster()
            [a, b, c, d]

            sage: S = ClusterSeed(['A',4], user_labels=['a','b','c']);
            Traceback (most recent call last):
            ...
            ValueError: The number of user-defined labels is not the number of exchangeable and frozen variables.

            sage: S = ClusterSeed(['A',4],user_labels=['x','y','w','z'])
            sage: S.mutate('x')
            sage: S.cluster()
            [(y + 1)/x, y, w, z]
            sage: S.mutate('(y+1)/x')
            sage: S.cluster()
            [x, y, w, z]
            sage: S.mutate('y')
            sage: S.cluster()
            [x, (x*w + 1)/y, w, z]
            sage: S.mutate('(x*w+1)/y')
            sage: S.cluster()
            [x, y, w, z]

            sage: S = ClusterSeed(['A',4], user_labels=[[1,2],[2,3],[4,5],[5,6]]);
            sage: S.cluster()
            [x_1_2, x_2_3, x_4_5, x_5_6]
            sage: S.mutate('[1,2]');
            sage: S.cluster()
            [(x_2_3 + 1)/x_1_2, x_2_3, x_4_5, x_5_6]

            sage: S = ClusterSeed(['A',4], user_labels=[[1,2],[2,3],[4,5],[5,6]],user_labels_prefix='P');
            sage: S.cluster()
            [P_1_2, P_2_3, P_4_5, P_5_6]
            sage: S.mutate('[1,2]')
            sage: S.cluster()
            [(P_2_3 + 1)/P_1_2, P_2_3, P_4_5, P_5_6]
            sage: S.mutate('P_4_5')
            sage: S.cluster()
            [(P_2_3 + 1)/P_1_2, P_2_3, (P_2_3*P_5_6 + 1)/P_4_5, P_5_6]

            sage: S = ClusterSeed(['A',4])
            sage: S.mutate([0,1,0,1,0,2,1])
            sage: T = ClusterSeed(S)
            sage: S.use_fpolys(False)
            sage: S.use_g_vectors(False)
            sage: S.use_c_vectors(False)
            sage: S._C
            sage: S._G
            sage: S._F
            sage: S.g_matrix()
            [ 0 -1  0  0]
            [ 1  1  1  0]
            [ 0  0 -1  0]
            [ 0  0  1  1]
            sage: S.c_matrix()
            [ 1 -1  0  0]
            [ 1  0  0  0]
            [ 1  0 -1  1]
            [ 0  0  0  1]
            sage: S.f_polynomials() == T.f_polynomials()
            True

            sage: S.cluster() == T.cluster()
            True

            sage: S._mut_path
            [0, 1, 0, 1, 0, 2, 1]

        """

        # check for sanitizable data
        if not isinstance(inplace, bool):
            raise ValueError('The second parameter must be boolean.  To mutate at a sequence of length 2, input it as a list.')

        if inplace:
            seed = self
        else:
            seed = ClusterSeed( self)

        # If we get a string, execute as a function
        if isinstance(sequence, str) and len(sequence) > 1 and sequence[0] is not '_':
            if sequence is 'green':
                sequence = self.first_green_vertex()
            elif sequence is 'red':
                sequence = self.first_red_vertex()
            elif sequence is 'urban' or sequence is 'urban_renewal':
                sequence = self.first_urban_renewal()
            elif sequence is 'all_urbans' or sequence is 'all_urban_renewals':
                sequence = self.urban_renewals()
            elif hasattr(self, sequence):
                sequence = getattr(self, sequence)()
            elif hasattr(self.quiver(), sequence):
                sequence = getattr(self.quiver(), sequence)()
            # If we are given a list in string format
            elif sequence[0] == '[' and sequence[-1] == ']':
                # convert to list
                from ast import literal_eval
                temp_list = literal_eval(sequence)

                sequence = self._user_labels_prefix
                for j in temp_list:
                    if isinstance(j, Integer):
                        sequence = sequence+"_"+j.str()
                    elif isinstance(j, int):
                        sequence = sequence+"_"+repr(j)
                    else:
                        sequence = sequence+"_"+j

        # If we get a function, execute it
        if hasattr(sequence, '__call__'):
            # function should return either integer or sequence
            sequence = sequence(seed)


        if sequence is None:
            raise ValueError('Not mutating: No vertices given.')

        if seed._use_fpolys:
            IE = seed._init_exch.values()
        else:
            IE = []

        n, m = seed.n(), seed.m()
        V = list(range(n)) + IE

        if seed._use_fpolys and isinstance(sequence, str):
            sequence = seed.cluster_index(sequence)
            if sequence is None:
                raise ValueError("Variable provided is not in our cluster")

        if (sequence in list(range(n))) or (sequence in IE):
            seqq = [sequence]
        else:
            seqq = sequence



        if isinstance(seqq, tuple):
            seqq = list( seqq )
        if not isinstance(seqq, list):
            raise ValueError('The quiver can only be mutated at a vertex or at a sequence of vertices')

        # remove ineligible vertices
        #if any( v not in V for v in seqq ):
            #v = filter( lambda v: v not in V, seqq )[0]
            #raise ValueError('The quiver cannot be mutated at the vertex ' + str( v ))

        seq = iter(seqq)

        for k in seq:

            if k in range(n):
                pass
            elif seed._use_fpolys:
                k = seed.cluster_index(k)
                if k is None:
                    raise ValueError("Variable provided is not in our cluster")
            else:
                raise ValueError('Why wasnt this caught earlier? Cannot mutate in direction ' + str(k) + '.')

            if seed._use_fpolys:
                seed._f_mutate(k)

            if seed._use_g_vec:
                seed._g_mutate(k)

            if seed._use_d_vec:
                seed._d_mutate(k)

            seed._BC.mutate(k)
            seed._M = copy(seed._BC[:n+m,:n])
            self._M.set_immutable()

            if seed._use_c_vec:
                seed._C = seed._BC[n+m:2*n+m,:n+m]

            if seed._track_mut:
                # delete involutive mutations
                if len(seed._mut_path) == 0 or seed._mut_path[len(self._mut_path)-1] != k:
                    seed._mut_path.append(k)
                else:
                    seed._mut_path.pop()

            # a mutation invalidates the cluster although it can be recomputed by F-polys and g-vectors
            # moving this into the for loop in case it does some mutations in 'seq' before finding a ValueError
            seed._cluster = None

            seed._quiver = None

        if not inplace:
            return seed

    def cluster_index(self, cluster_str):
        r"""
        Returns the index of a cluster if use_fpolys is on
        
        INPUT:

        - ``cluster_str`` -- The string to look for in the cluster

        OUTPUT:

        Returns an integer or None if the string is not a cluster variable
 
        EXAMPLES::

            sage: S = ClusterSeed(['A',4],user_labels=['x','y','z','w']); S.mutate('x')
            sage: S.cluster_index('x')
            sage: S.cluster_index('(y+1)/x')
            0

        """
        if self._use_fpolys and isinstance(cluster_str, str):
            c = FractionField(self._R)(cluster_str)
            cluster_str = ClusterVariable( FractionField(self._R), c.numerator(), c.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable',xdim=self._n )
            if cluster_str in self.cluster():
                return self.cluster().index(cluster_str)

        return None

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
            ....:     print(T.b_matrix())
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
        seed = ClusterSeed(self)

        new_clust_var = []
        seed_sequence = []

        for v in sequence:
            seed = seed.mutate(v, inplace=False)
            new_clust_var.append(seed.cluster()[v])
            seed_sequence.append(seed)

        if show_sequence:
            self.quiver().mutation_sequence2(sequence=sequence, show_sequence=True, fig_size=fig_size )

        if return_output=='seed':
            return seed_sequence
        elif return_output=='matrix':
            return [ seed._M for seed in seed_sequence ]
        elif return_output=='var':
            return new_clust_var
        else:
            raise ValueError('The parameter `return_output` can only be `seed`, `matrix`, or `var`.')

    def mutation_analysis(self, options=['all'], filter=None):
        r"""
        Runs an analysis of all potential mutation options. Note that this might take a long time on large seeds.

        Notes: Edges are only returned if we have a non-valued quiver. Green and red vertices are only returned if the cluster is principal.

        INPUT:

        - ``options`` -- (default: ['all']) a list of mutation options.
        - ``filter`` -- (default: None) A vertex or interval of vertices to limit our search to

        Possible options are:

        - ``"all"`` - All options below
        - ``"edges"`` - Number of edges (works with skew-symmetric quivers)
        - ``"edge_diff"`` - Edges added/deleted (works with skew-symmetric quivers)
        - ``"green_vertices"`` - List of green vertices (works with principals)
        - ``"green_vertices_diff"`` - Green vertices added/removed (works with principals)
        - ``"red_vertices"`` - List of red vertices (works with principals)
        - ``"red_vertices_diff"`` - Red vertices added/removed (works with principals)
        - ``"urban_renewals"`` - List of urban renewal vertices
        - ``"urban_renewals_diff"`` - Urban renewal vertices added/removed
        - ``"sources"`` - List of source vertices
        - ``"sources_diff"`` - Source vertices added/removed
        - ``"sinks"`` - List of sink vertices
        - ``"sinks_diff"`` - Sink vertices added/removed
        - ``"denominators"`` - List of all denominators of the cluster variables

        OUTPUT:

        Outputs a dictionary indexed by the vertex numbers. Each vertex will itself also be a
        dictionary with each desired option included as a key in the dictionary. As an example
        you would get something similar to: {0: {'edges': 1}, 1: {'edges': 2}}. This represents
        that if you were to do a mutation at the current seed then mutating at vertex 0 would result
        in a quiver with 1 edge and mutating at vertex 0 would result in a quiver with 2 edges.

        EXAMPLES::

            sage: B = [[0, 4, 0, -1],[-4,0, 3, 0],[0, -3, 0, 1],[1, 0, -1, 0]]
            sage: S = ClusterSeed(matrix(B)); S.mutate([2,3,1,2,1,3,0,2])
            sage: S.mutation_analysis()
            {0: {'d_matrix': [ 0  0  1  0]
              [ 0 -1  0  0]
              [ 0  0  0 -1]
              [-1  0  0  0],
              'denominators': [1, 1, x0, 1],
              'edge_diff': 6,
              'edges': 13,
              'green_vertices': [0, 1, 3],
              'green_vertices_diff': {'added': [0], 'removed': []},
              'red_vertices': [2],
              'red_vertices_diff': {'added': [], 'removed': [0]},
              'sinks': [],
              'sinks_diff': {'added': [], 'removed': [2]},
              'sources': [],
              'sources_diff': {'added': [], 'removed': []},
              'urban_renewals': [],
              'urban_renewals_diff': {'added': [], 'removed': []}},
             1: {'d_matrix': [ 1  4  1  0]
              [ 0  1  0  0]
              [ 0  0  0 -1]
              [ 1  4  0  0],
              'denominators': [x0*x3, x0^4*x1*x3^4, x0, 1],
              'edge_diff': 2,
              'edges': 9,
              'green_vertices': [0, 3],
              'green_vertices_diff': {'added': [0], 'removed': [1]},
              'red_vertices': [1, 2],
              'red_vertices_diff': {'added': [1], 'removed': [0]},
              'sinks': [2],
              'sinks_diff': {'added': [], 'removed': []},
              'sources': [],
              'sources_diff': {'added': [], 'removed': []},
              'urban_renewals': [],
              'urban_renewals_diff': {'added': [], 'removed': []}},
             2: {'d_matrix': [ 1  0  0  0]
              [ 0 -1  0  0]
              [ 0  0  0 -1]
              [ 1  0  1  0],
              'denominators': [x0*x3, 1, x3, 1],
              'edge_diff': 0,
              'edges': 7,
              'green_vertices': [1, 2, 3],
              'green_vertices_diff': {'added': [2], 'removed': []},
              'red_vertices': [0],
              'red_vertices_diff': {'added': [], 'removed': [2]},
              'sinks': [],
              'sinks_diff': {'added': [], 'removed': [2]},
              'sources': [2],
              'sources_diff': {'added': [2], 'removed': []},
              'urban_renewals': [],
              'urban_renewals_diff': {'added': [], 'removed': []}},
             3: {'d_matrix': [ 1  0  1  1]
              [ 0 -1  0  0]
              [ 0  0  0  1]
              [ 1  0  0  1],
              'denominators': [x0*x3, 1, x0, x0*x2*x3],
              'edge_diff': -1,
              'edges': 6,
              'green_vertices': [1],
              'green_vertices_diff': {'added': [], 'removed': [3]},
              'red_vertices': [0, 2, 3],
              'red_vertices_diff': {'added': [3], 'removed': []},
              'sinks': [2],
              'sinks_diff': {'added': [], 'removed': []},
              'sources': [1],
              'sources_diff': {'added': [1], 'removed': []},
              'urban_renewals': [],
              'urban_renewals_diff': {'added': [], 'removed': []}}}

            sage: S = ClusterSeed(['A',3]).principal_extension()
            sage: S.mutation_analysis()
            {0: {'d_matrix': [ 1  0  0]
              [ 0 -1  0]
              [ 0  0 -1],
              'denominators': [x0, 1, 1],
              'green_vertices': [1, 2],
              'green_vertices_diff': {'added': [], 'removed': [0]},
              'red_vertices': [0],
              'red_vertices_diff': {'added': [0], 'removed': []},
              'sinks': [],
              'sinks_diff': {'added': [], 'removed': [1]},
              'sources': [4, 5],
              'sources_diff': {'added': [], 'removed': [3]},
              'urban_renewals': [],
              'urban_renewals_diff': {'added': [], 'removed': []}},
             1: {'d_matrix': [-1  0  0]
              [ 0  1  0]
              [ 0  0 -1],
              'denominators': [1, x1, 1],
              'green_vertices': [0, 2],
              'green_vertices_diff': {'added': [], 'removed': [1]},
              'red_vertices': [1],
              'red_vertices_diff': {'added': [1], 'removed': []},
              'sinks': [0, 2, 4],
              'sinks_diff': {'added': [0, 2, 4], 'removed': [1]},
              'sources': [1, 3, 5],
              'sources_diff': {'added': [1], 'removed': [4]},
              'urban_renewals': [],
              'urban_renewals_diff': {'added': [], 'removed': []}},
             2: {'d_matrix': [-1  0  0]
              [ 0 -1  0]
              [ 0  0  1],
              'denominators': [1, 1, x2],
              'green_vertices': [0, 1],
              'green_vertices_diff': {'added': [], 'removed': [2]},
              'red_vertices': [2],
              'red_vertices_diff': {'added': [2], 'removed': []},
              'sinks': [],
              'sinks_diff': {'added': [], 'removed': [1]},
              'sources': [3, 4],
              'sources_diff': {'added': [], 'removed': [5]},
              'urban_renewals': [],
              'urban_renewals_diff': {'added': [], 'removed': []}}}

        """

        V = list(range(self._n))

        if filter is None:
            filter = V
        if filter in V:
            filter = [filter]

        # setup our initial information for differences later on
        if 'edge_diff' in options or ('all' in options and self._M.is_skew_symmetric()):
            initial_edges = self.quiver().number_of_edges()
        if 'green_vertices_diff' in options or ('all' in options and self._use_c_vec):
            initial_green_vertices = self.green_vertices()
        if 'red_vertices_diff' in options or ('all' in options and self._use_c_vec):
            initial_red_vertices = self.red_vertices()
        if 'urban_renewals_diff' in options or 'all' in options:
            initial_urban_renewals= self.urban_renewals()
        if 'sources_diff' in options or 'all' in options:
            initial_sources = self.quiver().sources()
        if 'sinks_diff' in options or 'all' in options:
            initial_sinks = self.quiver().sinks()

        #instantiate our dictionary
        analysis = {}
        for i in filter:
            #instantiate our dictionary
            analysis[i] = {}

            #run mutations not in place as we just want an analysis
            current_mutation = self.mutate(i,inplace=False)

            if ('edges' in options or 'all' in options) and self._M.is_skew_symmetric():
                analysis[i]['edges'] = current_mutation.quiver().number_of_edges()
            if ('edge_diff' in options or 'all' in options) and self._M.is_skew_symmetric():
                analysis[i]['edge_diff'] = current_mutation.quiver().number_of_edges() - initial_edges

            if ('green_vertices' in options or 'all' in options) and self._use_c_vec:
                analysis[i]['green_vertices'] = current_mutation.green_vertices()
            if ('green_vertices_diff' in options or 'all' in options) and self._use_c_vec:
                analysis[i]['green_vertices_diff'] = {}
                new_green_vertices = current_mutation.green_vertices()
                analysis[i]['green_vertices_diff']['added'] = list(set(new_green_vertices) - set(initial_green_vertices))
                analysis[i]['green_vertices_diff']['removed'] = list(set(initial_green_vertices) - set(new_green_vertices))

            if ('red_vertices' in options or 'all' in options) and self._use_c_vec:
                analysis[i]['red_vertices'] = current_mutation.red_vertices()
            if ('red_vertices_diff' in options or 'all' in options) and self._use_c_vec:
                analysis[i]['red_vertices_diff'] = {}
                new_red_vertices = current_mutation.red_vertices()
                analysis[i]['red_vertices_diff']['added'] = list(set(new_red_vertices) - set(initial_red_vertices))
                analysis[i]['red_vertices_diff']['removed'] = list(set(initial_red_vertices) - set(new_red_vertices))

            if 'urban_renewals' in options or 'all' in options:
                analysis[i]['urban_renewals'] = current_mutation.urban_renewals()
            if 'urban_renewals_diff' in options or 'all' in options:
                analysis[i]['urban_renewals_diff'] = {}
                new_urban_renewals = current_mutation.urban_renewals()
                analysis[i]['urban_renewals_diff']['added'] = list(set(new_urban_renewals) - set(initial_urban_renewals))
                analysis[i]['urban_renewals_diff']['removed'] = list(set(initial_urban_renewals) - set(new_urban_renewals))

            if 'sources' in options or 'all' in options:
                analysis[i]['sources'] = current_mutation.quiver().sources()
            if 'sources_diff' in options or 'all' in options:
                analysis[i]['sources_diff'] = {}
                new_sources = current_mutation.quiver().sources()
                analysis[i]['sources_diff']['added'] = list(set(new_sources) - set(initial_sources))
                analysis[i]['sources_diff']['removed'] = list(set(initial_sources) - set(new_sources))

            if 'sinks' in options or 'all' in options:
                analysis[i]['sinks'] = current_mutation.quiver().sinks()
            if 'sinks_diff' in options or 'all' in options:
                analysis[i]['sinks_diff'] = {}
                new_sinks = current_mutation.quiver().sinks()
                analysis[i]['sinks_diff']['added'] = list(set(new_sinks) - set(initial_sinks))
                analysis[i]['sinks_diff']['removed'] = list(set(initial_sinks) - set(new_sinks))

            if ('denominators' in options or 'all' in options) and self._use_fpolys:
                analysis[i]['denominators'] = []
                for vari in current_mutation.cluster():
                    analysis[i]['denominators'].append(vari.denominator())

            if ('d_matrix' in options or 'all' in options) and (self._use_d_vec or self._use_fpolys):
                analysis[i]['d_matrix'] = current_mutation.d_matrix()

        return analysis

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

        """
        from sage.combinat.cluster_algebra_quiver.mutation_class import _principal_part
        eval_dict = dict( [ ( self.y(i), 1 ) for i in range(self._m) ] )

        seed = ClusterSeed( _principal_part( self._M ), is_principal = True, user_labels=self._user_labels, user_labels_prefix=self._user_labels_prefix, frozen=None) 
        seed.use_c_vectors(self._use_c_vec)
        seed.use_fpolys(self._use_fpolys)
        seed.use_g_vectors(self._use_g_vec)
        seed.use_d_vectors(self._use_d_vec)
        seed.track_mutations(self._track_mut)
        if self._use_fpolys:
            self.cluster()
            seed._cluster = [self._cluster[k].subs(eval_dict)
                             for k in range(self._n)]
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

        # We give the indexing set of the Cartan matrix to be [1, 2, ..., n]
        rs = CartanMatrix(A, index_set=list(range(1,A.ncols()+1))).root_space()
        almost_positive_coroots = rs.almost_positive_roots()

        sign = [-1 if all(x <= 0 for x in self.b_matrix()[i]) else 1
                for i in range(self._n)]
        C = matrix([[sign[j] * alpha[j + 1] for j in range(self._n)]
                    for alpha in almost_positive_coroots])

        M = self._M.stack(C)
        seed = ClusterSeed(M, is_principal = False, user_labels=self._user_labels, user_labels_prefix=self._user_labels_prefix, frozen=None) 
        seed.use_c_vectors(self._use_c_vec)
        seed.use_fpolys(self._use_fpolys)
        seed.use_g_vectors(self._use_g_vec)
        seed.use_d_vectors(self._use_d_vec)
        seed.track_mutations(self._track_mut)

        seed._mutation_type = self._mutation_type
        return seed

    def principal_extension(self):
        r"""
        Returns the principal extension of self, yielding a 2n-by-n matrix.  Raises an error if the input seed has a non-square exchange matrix.  In this case, the method instead adds n frozen variables to any previously frozen variables.
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
 
            sage: S = ClusterSeed(['A',4],user_labels=['a','b','c','d'])
            sage: T= S.principal_extension()
            sage: T.cluster()
            [a, b, c, d]
            sage: T.coefficients()
            [y0, y1, y2, y3]
            sage: S2 = ClusterSeed(['A',4],user_labels={0:'a',1:'b',2:'c',3:'d'})
            sage: S2 == S
            True
            sage: T2 = S2.principal_extension()
            sage: T2 == T
            True
            """
        from sage.matrix.all import identity_matrix
        if self._m != 0:
            raise ValueError("The b-matrix is not square.")
        M = self._M.stack(identity_matrix(self._n))
        is_principal = (self._m == 0)
        if self._user_labels:
            if isinstance(self._user_labels,list):
                self._user_labels = self._user_labels + ['y%s'%i for i in range(self._n)]
            elif isinstance(self._user_labels,dict):
                self._user_labels.update( {(i+self._n):'y%s'%i for i in range(self._n)}  )
        seed = ClusterSeed(M, is_principal = is_principal, user_labels=self._user_labels, user_labels_prefix=self._user_labels_prefix, frozen=None) 
        seed.use_c_vectors(self._use_c_vec)
        seed.use_fpolys(self._use_fpolys)
        seed.use_g_vectors(self._use_g_vec)
        seed.use_d_vectors(self._use_d_vec)
        seed.track_mutations(self._track_mut)

        #### This should fix principal_extension resetting boolean flags.  Might need to update user labels to include new principals with y's.    -G
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

    def set_cluster( self, cluster, force=False ):
        r"""
        Sets the cluster for ``self`` to ``cluster``.

        Warning: Initialization may lead to inconsistent data.

        INPUT:

        - ``cluster`` -- an iterable defining a cluster for ``self``.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: cluster = S.cluster()
            sage: S.mutate([1,2,1])
            sage: S.cluster()
            [x0, (x1 + 1)/x2, (x0*x2 + x1 + 1)/(x1*x2)]
            sage: cluster2 = S.cluster()

            sage: S.set_cluster(cluster)
            Warning: using set_cluster at this point could lead to inconsistent seed data.

            sage: S.set_cluster(cluster, force=True)
            sage: S.cluster()
            [x0, x1, x2]
            sage: S.set_cluster(cluster2, force=True)
            sage: S.cluster()
            [x0, (x1 + 1)/x2, (x0*x2 + x1 + 1)/(x1*x2)]

            sage: S = ClusterSeed(['A',3]); S.use_fpolys(False)
            sage: S.set_cluster([1,1,1])
            Warning: clusters not being tracked so this command is ignored.
        """

        if len(cluster) < self._n+self._m:
            raise ValueError('The number of given cluster variables is wrong')
        if self._use_fpolys:        
            if any(c not in FractionField(self._R) for c in cluster):
                raise ValueError('The cluster variables are not all contained in %s'%FractionField(self._R))
            if not force:  # if already have f_polynomials, using set_cluster might yield data inconsistent with them.
                print("Warning: using set_cluster at this point could lead to inconsistent seed data.")
            else:
                self._cluster = [FractionField(self._R)(x)
                                 for x in cluster][0:self._n]
                self._is_principal = None
        else:
             print("Warning: clusters not being tracked so this command is ignored.") 

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
 
            sage: S = ClusterSeed(['B',3],user_labels=[[1,2],[2,3],[3,4]],user_labels_prefix='p')
            sage: S.mutate([0,1])
            sage: S.cluster()
            [(p_2_3 + 1)/p_1_2, (p_1_2*p_3_4^2 + p_2_3 + 1)/(p_1_2*p_2_3), p_3_4]
 
            sage: S.reset_cluster()
            sage: S.cluster()
            [p_1_2, p_2_3, p_3_4]
            sage: S.g_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: S.f_polynomials()
            [1, 1, 1]
        """
        if self._use_g_vec:
            self._G = matrix.identity(self._n)
        if self._use_fpolys:
            self._F = dict([(i,self._U(1)) for i in self._init_exch.values()])
        if self._use_fpolys:
            self.set_cluster(self._R.gens(), force=True)
 
    def reset_coefficients( self ):
        r"""
        Resets the coefficients of ``self`` to the frozen variables but keeps the current cluster.
        Raises an error if the number of frozen variables is different than the number of exchangeable variables.

        WARNING: This command to be phased out since 'use_c_vectors() does this more effectively.

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
        n, m = self._n, self._m
        if not n == m:
            raise ValueError("The numbers of cluster variables "
                             "and of frozen variables do not coincide.")
        newM = copy(self._M)
        for i in range(m):
            for j in range(n):
                if i == j:
                    newM[i + n, j] = 1
                else:
                    newM[i + n, j] = 0
        self._M = newM
        self._M.set_immutable()
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [1])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [0])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [1, 0])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [0, 1])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [1, 0, 1])
            (A seed for a cluster algebra of rank 2 of type ['A', [1, 1], 1], [0, 1, 0])
        """

        # Variable to track the depth
        depth_counter = 0
        n = self._n
        timer = time.time()

        # set up our initial cluster and grab variables
        if up_to_equivalence:
            cl = Set( self.cluster() )
        else:
            cl = tuple( self.cluster() )

        # If we are tracking return paths
        if return_paths:
            yield (self,[])
        else:
            yield self


        # instantiate the variables
        clusters = {}
        clusters[ cl ] = [ self, list(range(n)), [] ]

        # we get bigger the first time
        gets_bigger = True

        # If we are showing depth, show some statistics
        if show_depth:
            timer2 = time.time()
            dc = str(depth_counter)
            dc += ' ' * (5-len(dc))
            nr = str(len(clusters))
            nr += ' ' * (10-len(nr))
            print("Depth: %s found: %s Time: %.2f s" % (dc, nr, timer2-timer))

        # Each time we get bigger and we haven't hit the full depth
        while gets_bigger and depth_counter < depth:
            gets_bigger = False

            # set the keys
            keys = clusters.keys()

            # Our keys are cluster variables, so for each cluster:
            for key in keys:
                # sd is the cluster data
                sd = clusters[key]

                # another way to do a for loop for each item
                while sd[1]:
                    i = sd[1].pop()

                    # If we aren't only sinking the source
                    if not only_sink_source or all( entry >= 0 for entry in sd[0]._M.row( i ) ) or all( entry <= 0 for entry in sd[0]._M.row( i ) ):
                        # do an inplace mutation on our cluster (sd[0])
                        sd2  = sd[0].mutate( i, inplace=False )

                        # set up our new cluster variables
                        if up_to_equivalence:
                            cl2 = Set(sd2.cluster())
                        else:
                            cl2 = tuple(sd2.cluster())
                        if cl2 in clusters:
                            if not up_to_equivalence and i in clusters[cl2][1]:
                                clusters[cl2][1].remove(i)
                        else:
                            gets_bigger = True
                            if only_sink_source:
                                orbits = list(range(n))
                            else:
                                orbits = [ index for index in range(n) if index > i or sd2._M[index,i] != 0 ]

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
                print("Depth: %s found: %s Time: %.2f s" % (dc,nr,timer2-timer))

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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
            [x0, x1, x2]
            [x0, x1, (x1 + 1)/x2]
            [x0, (x0*x2 + 1)/x1, x2]
            [(x1 + 1)/x0, x1, x2]

        A finite type example where the depth is returned while computing::

            sage: it = S.cluster_class_iter(show_depth=True)
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
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
            sage: for T in it: print(T)
            x0
            x1
            x2
            (x1 + 1)/x2
            (x0*x2 + 1)/x1
            (x1 + 1)/x0

            sage: S.mutate([0,1])
            sage: it2 = S.variable_class_iter(depth=1)
            sage: for T in it2: print(T)
            (x1 + 1)/x0
            (x0*x2 + x1 + 1)/(x0*x1)
            x2
            (x1^2 + x0*x2 + 2*x1 + 1)/(x0*x1*x2)
            x1
            (x0*x2 + 1)/x1

        Infinite type examples::

            sage: S = ClusterSeed(['A',[1,1],1])
            sage: it = S.variable_class_iter(depth=2)
            sage: for T in it: print(T)
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
                    print("Found a bipartite seed - restarting the depth counter at zero and constructing the variable class using its bipartite belt.")
                depth_counter = 0
                end = False
                seed2 = ClusterSeed(seed)
                for c in seed.cluster():
                    if c not in var_class:
                        yield ClusterVariable( FractionField(seed._R), c.numerator(), c.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable',xdim=seed._n )
                var_class = var_class.union( seed.cluster())

                init_cluster = set(seed.cluster())
                while not end and depth_counter < depth:
                    depth_counter += 1
                    seed.mutate(bipartition[0])
                    seed.mutate(bipartition[1])
                    if set(seed.cluster()) in [set(seed2.cluster()),init_cluster]:
                        end = True
                    if not end:
                        for c in seed.cluster():
                            if c not in var_class:
                                yield ClusterVariable( FractionField(seed._R), c.numerator(), c.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable',xdim=seed._n )
                        var_class = var_class.union( seed.cluster() )
                        seed2.mutate(bipartition[1])
                        seed2.mutate(bipartition[0])
                        if set(seed2.cluster()) in [set(seed.cluster()),init_cluster]:
                            end = True
                        if not end:
                            for c in seed2.cluster():
                                if c not in var_class:
                                    yield ClusterVariable(FractionField(seed._R), c.numerator(), c.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable',xdim=seed._n )
                            var_class = var_class.union(seed2.cluster())
                return
            else:
                for c in seed.cluster():
                    if c not in var_class:
                        yield ClusterVariable( FractionField(seed._R), c.numerator(), c.denominator(), mutation_type=self._mutation_type, variable_type='cluster variable',xdim=seed._n)
                var_class = var_class.union(seed.cluster())

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

    def is_finite(self):
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

        - Might fail to work if it is used within different Sage processes simultaneously (that happened in the doctesting).

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

    @rename_keyword(deprecation=19572, method='algorithm')
    def greedy(self, a1, a2, algorithm='by_recursion'):
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
            if algorithm == 'by_recursion':
                ans = self.x(0)**(-a1)*self.x(1)**(-a2)
                for p in range(max(a2, 0)+1):
                    for q in range(max(a1, 0)+1):
                        if p != 0 or q != 0:
                            ans += self._R(coeff_recurs(p, q, a1, a2, b, c))*self.x(0)**(b*p-a1)*self.x(1)**(c*q-a2)
                return(ans)
            elif algorithm == 'by_combinatorics':
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
            elif algorithm == 'just_numbers':
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

    def oriented_exchange_graph(self):
        """
        Return the oriented exchange graph of ``self`` as a directed
        graph.

        The seed must be a cluster seed for a cluster algebra of
        finite type with principal coefficients (the corresponding
        quiver must have mutable vertices 0,1,...,n-1).

        EXAMPLES::

            sage: S = ClusterSeed(['A', 2]).principal_extension()
            sage: G = S.oriented_exchange_graph(); G
            Digraph on 5 vertices
            sage: G.out_degree_sequence()
            [2, 1, 1, 1, 0]

            sage: S = ClusterSeed(['B', 2]).principal_extension()
            sage: G = S.oriented_exchange_graph(); G
            Digraph on 6 vertices
            sage: G.out_degree_sequence()
            [2, 1, 1, 1, 1, 0]

        TESTS::

            sage: S = ClusterSeed(['A',[2,2],1])
            sage: S.oriented_exchange_graph()
            Traceback (most recent call last):
            ...
            TypeError: only works for finite mutation type

            sage: S = ClusterSeed(['A', 2])
            sage: S.oriented_exchange_graph()
            Traceback (most recent call last):
            ...
            TypeError: only works for principal coefficients
        """
        if not self._mutation_type.is_finite():
            raise TypeError('only works for finite mutation type')

        if not self._is_principal:
            raise TypeError('only works for principal coefficients')

        covers = []
        n = self.n()
        stack = [self]
        known_clusters = []
        while stack:
            i = stack.pop()
            Vari = tuple(sorted(i.cluster()))
            B = i.b_matrix()
            for k in range(n):
                # check if green
                if all(B[i2][k] >= 0 for i2 in range(n, 2 * n)):
                    j = i.mutate(k, inplace=False)
                    Varj = tuple(sorted(j.cluster()))
                    covers.append((Vari, Varj))
                    if not(Varj in known_clusters):
                        known_clusters += [Varj]
                        stack.append(j)

        return DiGraph(covers)

    def find_upper_bound(self, verbose=False):
        r"""
        Return the upper bound of the given cluster algebra as a quotient_ring.

        The upper bound is the intersection of the Laurent polynomial
        rings of the initial cluster and its neighboring clusters.  As
        such, it always contains both the cluster algebra and the
        upper cluster algebra.  This function uses the algorithm from
        [MaMu]_.

        When the initial seed is totally coprime (for example, when
        the unfrozen part of the exchange matrix has full rank), the
        upper bound is equal to the upper cluster algebra by
        [BeFoZe]_.

        .. WARNING::

            The computation time grows rapidly with the size
            of the seed and the number of steps.  For most seeds
            larger than four vertices, the algorithm may take an
            infeasible amount of time.  Additionally, it will run
            forever without terminating whenever the upper bound is
            infinitely-generated (such as the example in [Sp]_).

        INPUT:

        - ``verbose`` -- (default: ``False``) if ``True``, prints output
          during the computation.

        EXAMPLES:

        - finite type::

            sage: S = ClusterSeed(['A',3])
            sage: S.find_upper_bound()
            Quotient of Multivariate Polynomial Ring in x0, x1, x2, x0p, x1p, x2p, z0 over Rational Field by the ideal (x0*x0p - x1 - 1, x1*x1p - x0*x2 - 1, x2*x2p - x1 - 1, x0*z0 - x2p, x1*z0 + z0 - x0p*x2p, x2*z0 - x0p, x1p*z0 + z0 - x0p*x1p*x2p + x1 + 1)

        - Markov::

            sage: B = matrix([[0,2,-2],[-2,0,2],[2,-2,0]])
            sage: S = ClusterSeed(B)
            sage: S.find_upper_bound()
            Quotient of Multivariate Polynomial Ring in x0, x1, x2, x0p, x1p, x2p, z0 over Rational Field by the ideal (x0*x0p - x2^2 - x1^2, x1*x1p - x2^2 - x0^2, x2*x2p - x1^2 - x0^2, x0p*x1p*x2p - x0*x1*x2p - x0*x2*x1p - x1*x2*x0p - 2*x0*x1*x2, x0^3*z0 - x1p*x2p + x1*x2, x0*x1*z0 - x2p - x2, x1^3*z0 - x0p*x2p + x0*x2, x0*x2*z0 - x1p - x1, x1*x2*z0 - x0p - x0, x2^3*z0 - x0p*x1p + x0*x1)

        REFERENCES:

        .. [BeFoZe] Berenstein-Fomin-Zelevinsky, Cluster algebras III: Upper bounds and double Bruhat cells, Duke Math. J., 126(1):1-52, 2005.
        .. [MaMu] Matherne-Muller, Computing upper cluster algebras, Int. Math. Res. Not. IMRN, 2015, 3121-3149.
        .. [Sp] Speyer, An infinitely generated upper cluster algebra, :arxiv:`1305.6867`.
        """
        rank = self.n()

        xvars = ['x{}'.format(t) for t in range(rank)]
        xpvars = ['x{}p'.format(t) for t in range(rank)]
        gens = xvars + xpvars
        initial_product = '*'.join(g for g in xvars)

        lower_var = self.cluster()
        for t in range(self.b_matrix().nrows()):
            lower_var += [(self.mutate(t, inplace=False)).cluster()[t]]

        deep_gens = [initial_product]
        for t in range(rank):
            neighbor_product = '*'.join(xpvars[s] if s == t else xvars[s]
                                        for s in range(rank))
            deep_gens += [neighbor_product]

        rels = ["-{}*{}+{}".format(gens[t], gens[t + rank],
                                   lower_var[t + rank].numerator())
                for t in range(rank)]

        while True:
            R = PolynomialRing(QQ, gens, order='invlex')
            I = R.ideal(rels)
            J = R.ideal(initial_product)
            if verbose:
                msg = 'Computing relations among {} generators'
                print(msg.format(len(gens)))
            start = time.time()
            ISat = I.saturation(J)[0]
            spend = time.time() - start
            if verbose:
                msg = 'Computed {} relations in {} seconds'
                print(msg.format(len(ISat.gens()), spend))
            deep_ideal = R.ideal(deep_gens) + ISat
            initial_product_ideal = R.ideal(initial_product) + ISat
            if verbose:
                print('Attempting to find a new element of the upper bound')
            start = time.time()
            M = initial_product_ideal.saturation(deep_ideal)[0]
            spend = time.time() - start
            if M == initial_product_ideal:
                if verbose:
                    print('Verified that there are no new elements in', spend, 'seconds')
                    print('Returning a presentation for the upper bound')
                return R.quotient_ring(ISat)
            else:
                gens.append('z' + str(len(gens) - 2 * rank))
                new_gen_found = False
                i = 0
                M_gens = M.gens()
                while (not new_gen_found) and i < len(M_gens):
                    f = initial_product_ideal.reduce(M_gens[i])
                    if f != 0:
                        rels.append('z' + str(len(gens) - 2 * rank - 1) + '*' + initial_product + '-(' + str(f) + ')')
                        new_gen_found = True
                        if verbose:
                            print('Found a new element in', spend, 'seconds!')
                            print('')
                    i += 1


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
        from sage.arith.all import binomial
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
                if not uv_okay:
                        return False
        return True


def get_green_vertices(C):
    r"""
    Get the green vertices from a matrix. Will go through each clumn and return
    the ones where no entry is greater than 0.

    INPUT:

    - ``C`` -- The C matrix to check
 
    EXAMPLES::
 
        sage: from sage.combinat.cluster_algebra_quiver.cluster_seed import get_green_vertices
        sage: S = ClusterSeed(['A',4]); S.mutate([1,2,3,2,0,1,2,0,3])
        sage: get_green_vertices(S.c_matrix())
        [0, 3]

    """
    return [ i for (i,v) in enumerate(C.columns()) if any(x > 0 for x in v) ]
    ## old code commented out
    #import numpy as np
    #max_entries = [ np.max(np.array(C.column(i))) for i in range(C.ncols()) ]
    #return [i for i in range(C.ncols()) if max_entries[i] > 0]

def get_red_vertices(C):
    r"""
    Get the red vertices from a matrix. Will go through each clumn and return
    the ones where no entry is less than 0.

    INPUT:

    - ``C`` -- The C matrix to check

    EXAMPLES::
        sage: from sage.combinat.cluster_algebra_quiver.cluster_seed import get_red_vertices
        sage: S = ClusterSeed(['A',4]); S.mutate([1,2,3,2,0,1,2,0,3])
        sage: get_red_vertices(S.c_matrix())
        [1, 2]

    """
    return [ i for (i,v) in enumerate(C.columns()) if any(x < 0 for x in v) ]
    ## old code commented out
    #import numpy as np
    #min_entries = [ np.min(np.array(C.column(i))) for i in range(C.ncols()) ]
    #return [i for i in range(C.ncols()) if min_entries[i] < 0]

class ClusterVariable(FractionFieldElement):
    r"""
    This class is a thin wrapper for cluster variables in cluster seeds.

    It provides the extra feature to store if a variable is frozen or not.

    - the associated positive root::

        sage: S = ClusterSeed(['A',3])
        sage: for T in S.variable_class_iter():
        ....:     print("{} {}".format(T, T.almost_positive_root()))
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
    def __init__( self, parent, numerator, denominator, coerce=True, reduce=True, mutation_type=None, variable_type=None, xdim=0 ):
        r"""
        Initializes a cluster variable in the same way that elements in the field of rational functions are initialized.

        .. see also:: :class:`Fraction Field of Multivariate Polynomial Ring`

        TESTS::

            sage: S = ClusterSeed(['A',2])
            sage: for f in S.cluster():
            ....:     print(type(f))
            <class 'sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterVariable'>
            <class 'sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterVariable'>

            sage: S.variable_class()
            [(x0 + x1 + 1)/(x0*x1), (x1 + 1)/x0, (x0 + 1)/x1, x1, x0]
        """
        FractionFieldElement.__init__( self, parent, numerator, denominator, coerce=coerce, reduce=reduce )
        self._n = xdim;
        self._mutation_type = mutation_type
        self._variable_type = variable_type

    def almost_positive_root( self ):
        r"""
        Returns the *almost positive root* associated to ``self`` if ``self`` is of finite type.

        EXAMPLES::

            sage: S = ClusterSeed(['A',3])
            sage: for T in S.variable_class_iter():
            ....:     print("{} {}".format(T, T.almost_positive_root()))
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
                    return sum( [ root[i]*Phiplus[ i+1 ] for i in range(self._n) ] )
            else:
                raise ValueError('The cluster algebra for %s is not of finite type.'%self._repr_())
