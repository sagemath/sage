r"""
Degree sequences

The present module implements the ``DegreeSequences`` class, whose instances
represent the integer sequences of length `n`::

    sage: DegreeSequences(6)
    Degree sequences on 6 elements

With the object ``DegreeSequences(n)``, one can :

    * Check whether a sequence is indeed a degree sequence ::

        sage: DS = DegreeSequences(5)
        sage: [4, 3, 3, 3, 3] in DS
        True
        sage: [4, 4, 0, 0, 0] in DS
        False

    * List all the possible degree sequences of length `n`::

        sage: for seq in DegreeSequences(4):
        ...       print seq
        [0, 0, 0, 0]
        [1, 1, 0, 0]
        [2, 1, 1, 0]
        [3, 1, 1, 1]
        [1, 1, 1, 1]
        [2, 2, 1, 1]
        [2, 2, 2, 0]
        [3, 2, 2, 1]
        [2, 2, 2, 2]
        [3, 3, 2, 2]
        [3, 3, 3, 3]

.. NOTE::

    Given a degree sequence, one can obtain a graph realizing it by using
    :meth:`sage.graphs.graph_generators.graphs.DegreeSequence`. For instance ::

        sage: ds = [3, 3, 2, 2, 2, 2, 2, 1, 1, 0]
        sage: g = graphs.DegreeSequence(ds)
        sage: g.degree_sequence()
        [3, 3, 2, 2, 2, 2, 2, 1, 1, 0]

Definitions
~~~~~~~~~~~

A sequence of integers `d_1,...,d_n` is said to be a *degree sequence* (or
*graphic* sequence) if there exists a graph in which vertex `i` is of degree
`d_i`. It is often required to be *non-increasing*, i.e. that `d_1 \geq ... \geq
d_n`.

An integer sequence need not necessarily be a degree sequence. Indeed, in a
degree sequence of length `n` no integer can be larger than `n-1` -- the degree
of a vertex is at most `n-1` -- and the sum of them is at most `n(n-1)`.

Degree sequences are completely characterized by a result from Erdos and Gallai:

**Erdos and Gallai:** *The sequence of integers* `d_1\geq ... \geq d_n` *is a
degree sequence if and only if* `\sum_i d_i` is even and `\forall i`

.. MATH::
    \sum_{j\leq i}d_j \leq j(j-1) + \sum_{j>i}\min(d_j,i)

Alternatively, a degree sequence can be defined recursively :

**Havel and Hakimi:** *The sequence of integers* `d_1\geq ... \geq d_n` *is a
degree sequence if and only if* `d_2-1,...,d_{d_1+1}-1, d_{d_1+2}, ...,d_n` *is
also a degree sequence.*

Or equivalently :

**Havel and Hakimi (bis):** *If there is a realization of an integer sequence as
a graph (i.e. if the sequence is a degree sequence), then it can be realized in
such a way that the vertex of maximum degree* `\Delta` *is adjacent to the*
`\Delta` *vertices of highest degree (except itself, of course).*


Algorithms
~~~~~~~~~~

**Checking whether a given sequence is a degree sequence**

This is tested using Erdos and Gallai's criterion. It is also checked that the
given sequence is non-increasing and has length `n`.

**Iterating through the sequences of length** `n`

From Havel and Hakimi's recursive definition of a degree sequence, one can build
an enumeration algorithm as done in [RCES]_. It consists in trying to **extend**
a current degree sequence on `n` elements into a degree sequence on `n+1`
elements by adding a vertex of degree larger than those already present in the
sequence. This can be seen as **reversing** the reduction operation described in
Havel and Hakimi's characterization. This operation can appear in several
different ways :

    * Extensions of a degree sequence that do **not** change the value of the
      maximum element

        * If the maximum element of a given degree sequence is `0`, then one can
          remove it to reduce the sequence, following Havel and Hakimi's
          rule. Conversely, if the maximum element of the (current) sequence is
          `0`, then one can always extend it by adding a new element of degree
          `0` to the sequence.

          .. MATH::
              0, 0, 0 \xrightarrow{Extension} {\bf 0}, 0, 0, 0 \xrightarrow{Extension} {\bf 0}, 0, 0, ..., 0, 0, 0 \xrightarrow{Reduction} 0, 0, 0, 0 \xrightarrow{Reduction} 0, 0, 0

        * If there are at least `\Delta+1` elements of (maximum) degree `\Delta`
          in a given degree sequence, then one can reduce it by removing a
          vertex of degree `\Delta` and decreasing the values of `\Delta`
          elements of value `\Delta` to `\Delta-1`. Conversely, if the maximum
          element of the (current) sequence is `d>0`, then one can add a new
          element of degree `d` to the sequence if it can be linked to `d`
          elements of (current) degree `d-1`. Those `d` vertices of degree `d-1`
          hence become vertices of degree `d`, and so `d` elements of degree
          `d-1` are removed from the sequence while `d+1` elements of degree `d`
          are added to it.

          .. MATH::
              3, 2, 2, 2, 1 \xrightarrow{Extension} {\bf 3}, 3, (2+1), (2+1), (2+1), 1 =  {\bf 3}, 3, 3, 3, 3, 1 \xrightarrow{Reduction} 3, 2, 2, 2, 1

    * Extension of a degree sequence that changes the value of the maximum
      element :

        * In the general case, i.e. when the number of elements of value
          `\Delta,\Delta-1` is small compared to `\Delta` (i.e. the maximum
          element of a given degree sequence), reducing a sequence strictly
          decreases the value of the maximum element. According to Havel and
          Hakimi's characterization there is only **one** way to reduce a
          sequence, but reversing this operation is more complicated than in the
          previous cases. Indeed, the following extensions are perfectly valid
          according to the reduction rule.

          .. MATH::
              2,1,1,0,0\xrightarrow{Extension} {\bf 3}, (2+1), (1+1), (1+1), 0, 0 = 3, 3, 2, 2, 0, 0 \xrightarrow{Reduction} 2, 1, 1, 0, 0\\
              2,1,1,0,0\xrightarrow{Extension} {\bf 3}, (2+1), (1+1), 1, (0+1), 0 = 3, 3, 2, 1, 1, 0 \xrightarrow{Reduction} 2, 1, 1, 0, 0\\
              2,1,1,0,0\xrightarrow{Extension} {\bf 3}, (2+1), 1, 1, (0+1), (0+1) = 3, 3, 1, 1, 1, 1 \xrightarrow{Reduction} 2, 1, 1, 0, 0\\
              ...

          In order to extend a current degree sequence while strictly increasing
          its maximum degree, it is equivalent to pick a set `I` of elements of
          the degree sequence with `|I|>\Delta` in such a way that the
          `(d_i+1)_{i\in I}` are the `|I|` maximum elements of the sequence
          `(d_i+\genfrac{}{}{0pt}{}{1\text{ if }i\in I}{0\text{ if }i\not \in
          I})_{1\leq i \leq n}`, and to add to this new sequence an element of
          value `|I|`. The non-increasing sequence containing the elements `|I|`
          and `(d_i+\genfrac{}{}{0pt}{}{1\text{ if }i\in I}{0\text{ if }i\not
          \in I})_{1\leq i \leq n}` can be reduced to `(d_i)_{1\leq i \leq n}`
          by Havel and Hakimi's rule.

          .. MATH::
              ... 1, 1, 2, {\bf 2}, {\bf 2}, 2, 2, 3, 3, \underline{3}, {\bf 3}, {\bf 3}, {\bf 4}, {\bf 6}, ... \xrightarrow{Extension} ... 1, 1, 2, 2, 2, 3, 3, \underline{3}, {\bf 3}, {\bf 3}, {\bf 4}, {\bf 4}, {\bf 5}, {\bf 7}, ...

          The number of possible sets `I` having this property (i.e. the number
          of possible extensions of a sequence) is smaller than it
          seems. Indeed, by definition, if `j\not \in I` then for all `i\in I`
          the inequality `d_j\leq d_i+1` holds. Hence, each set `I` is entirely
          determined by the largest element `d_k` of the sequence that it does
          **not** contain (hence `I` contains `\{1,...,k-1\}`), and by the
          cardinalities of `\{i\in I:d_i= d_k\}` and `\{i\in I:d_i= d_k-1\}`.

          .. MATH::
              I = \{i \in I : d_i= d_k \} \cup \{i \in I : d_i= d_k-1 \} \cup \{i : d_i> d_k \}

          The number of possible extensions is hence at most cubic, and is
          easily enumerated.

About the implementation
~~~~~~~~~~~~~~~~~~~~~~~~

In the actual implementation of the enumeration algorithm, the degree sequence
is stored differently for reasons of efficiency.

Indeed, when enumerating all the degree sequences of length `n`, Sage first
allocates an array ``seq`` of `n+1` integers where ``seq[i]`` is the number of
elements of value ``i`` in the current sequence. Obviously, ``seq[n]=0`` holds
in permanence : it is useful to allocate a larger array than necessary to
simplify the code. The ``seq`` array is a global variable.

The recursive function ``enum(depth, maximum)`` is the one building the list of
sequences. It builds the list of degree sequences of length `n` which *extend*
the sequence currently stored in ``seq[0]...seq[depth-1]``. When it is called,
``maximum`` must be set to the maximum value of an element in the partial
sequence ``seq[0]...seq[depth-1]``.

If during its run the function ``enum`` heavily works on the content of the
``seq`` array, the value of ``seq`` is the **same** before and after the run of
``enum``.

**Extending the current partial sequence**

The two cases for which the maximum degree of the partial sequence does not
change are easy to detect. It is (sligthly) harder to enumerate all the sets `I`
corresponding to possible extensions of the partial sequence. As said
previously, to each set `I` one can associate an integer ``current_box`` such
that `I` contains all the `i` satisfying `d_i>current\_box`. The variable
``taken`` represents the number of all such elements `i`, so that when
enumerating all possible sets `I` in the algorithm we have the equality

.. MATH::
    I = \text{taken }+\text{ number of elements of value }current\_box+ \text{ number of elements of value }current\_box-1

References
~~~~~~~~~~

  .. [RCES] Alley CATs in search of good homes
    Ruskey, R. Cohen, P. Eades, A. Scott
    Congressus numerantium, 1994
    Pages 97--110


Author
~~~~~~

Nathann Cohen

Tests
~~~~~

The sequences produced by random graphs *are* degree sequences::

    sage: n = 30
    sage: DS = DegreeSequences(30)
    sage: for i in range(10):
    ...      g = graphs.RandomGNP(n,.2)
    ...      if not g.degree_sequence() in DS:
    ...          print "Something is very wrong !"

Checking that we indeed enumerate *all* the degree sequences for `n=5`::

    sage: ds1 = Set([tuple(g.degree_sequence()) for g in graphs(5)])
    sage: ds2 = Set(map(tuple,list(DegreeSequences(5))))
    sage: ds1 == ds2
    True

Checking the consistency of enumeration and test::

    sage: DS = DegreeSequences(6)
    sage: all(seq in DS for seq in DS)
    True

.. WARNING::

    For the moment, iterating over all degree sequences involves building the
    list of them first, then iterate on this list.  This is obviously bad, as it
    requires uselessly a **lot** of memory for large values of `n`.

    As soon as the ``yield`` keyword is available in Cython this should be
    changed. Updating the code does not require more than a couple of minutes.

"""

##############################################################################
#       Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.libs.gmp.all cimport mpz_t
from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer
include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'
include "sage/ext/interrupt.pxi"


cdef unsigned char * seq
cdef list sequences

class DegreeSequences:

    def __init__(self, n):
        r"""
        Degree Sequences

        An instance of this class represents the degree sequences of graphs on a
        given number `n` of vertices. It can be used to list and count them, as
        well as to test whether a sequence is a degree sequence. For more
        information, please refer to the documentation of the
        :mod:`DegreeSequence<sage.combinat.degree_sequences>` module.

        EXAMPLE::

            sage: DegreeSequences(8)
            Degree sequences on 8 elements
            sage: [3,3,2,2,2,2,2,2] in DegreeSequences(8)
            True

        """
        self._n = n

    def __contains__(self, seq):
        """
        Checks whether a given integer sequence is the degree sequence
        of a graph on `n` elements

        EXAMPLE::

            sage: [3,3,2,2,2,2,2,2] in DegreeSequences(8)
            True

        TESTS:

        :trac:`15503`::

            sage: [2,2,2,2,1,1,1] in DegreeSequences(7)
            False
        """
        cdef int n = self._n
        if len(seq)!=n:
            return False

        # Is the sum even ?
        if sum(seq)%2 == 1:
            return False

        # Partial represents the left side of Erdos and Gallai's inequality,
        # i.e. the sum of the i first integers.
        cdef int partial = 0
        cdef int i,d,dd, right

        # Temporary variable to ensure that the sequence is indeed
        # non-increasing
        cdef int prev = n-1

        for i,d in enumerate(seq):

            # Non-increasing ?
            if d > prev:
                return False
            else:
                prev = d

            # Updating the partial sum
            partial += d

            # Evaluating the right hand side
            right = i*(i+1)
            for dd in seq[i+1:]:
                right += min(dd,i+1)

            # Comparing the two
            if partial > right:
                return False

        return True

    def __repr__(self):
        """
        Representing the element

        TEST::

            sage: DegreeSequences(6)
            Degree sequences on 6 elements
        """
        return "Degree sequences on "+str(self._n)+" elements"

    def __iter__(self):
        """
        Iterate over all the degree sequences.

        TODO: THIS SHOULD BE UPDATED AS SOON AS THE YIELD KEYWORD APPEARS IN
        CYTHON. See comment in the class' documentation.

        EXAMPLE::

            sage: DS = DegreeSequences(6)
            sage: all(seq in DS for seq in DS)
            True
        """

        init(self._n)
        return iter(sequences)

    def __dealloc__():
        """
        Freeing the memory
        """
        if seq != NULL:
            sage_free(seq)

cdef init(int n):
    """
    Initializes the memory and starts the enumeration algorithm.
    """
    global seq
    global N
    global sequences

    if n == 0:
        return [[]]
    elif n == 1:
        return [[0]]

    sig_on()
    seq = <unsigned char *> sage_malloc((n+1)*sizeof(unsigned char))
    memset(seq,0,(n+1)*sizeof(unsigned char))
    sig_off()

    # We begin with one vertex of degree 0
    seq[0] = 1

    N = n
    sequences = []
    enum(1,0)
    sage_free(seq)
    return sequences

cdef inline add_seq():
     """
     This function is called whenever a sequence is found.

     Build the degree sequence corresponding to the current state of the
     algorithm and adds it to the sequences list.
     """
     global sequences
     global N
     global seq

     cdef list s = []
     cdef int i, j

     for N > i >= 0:
         for 0<= j < seq[i]:
             s.append(i)

     sequences.append(s)

cdef void enum(int k, int M):
    """
    Main function. For an explanation of the algorithm please refer to the
    class' documentation.

    INPUT:

    * ``k`` -- depth of the partial degree sequence
    * ``M`` -- value of a maximum element in the partial degree sequence
    """
    cdef int i,j
    global seq
    cdef int taken = 0
    cdef int current_box
    cdef int n_current_box
    cdef int n_previous_box
    cdef int new_vertex

    # Have we found a new degree sequence ? End of recursion !
    if k == N:
        add_seq()
        return

    sig_on()

    #############################################
    # Creating vertices of Vertices of degree M #
    #############################################

    # If 0 is the current maximum degree, we can always extend the degree
    # sequence with another 0
    if M == 0:

        seq[0] += 1
        enum(k+1, M)
        seq[0] -= 1

    # We need not automatically increase the degree at each step. In this case,
    # we have no other choice but to link the new vertex of degree M to vertices
    # of degree M-1, which will become vertices of degree M too.
    elif seq[M-1] >= M:

        seq[M]   += M+1
        seq[M-1] -= M

        enum(k+1, M)

        seq[M]   -= M+1
        seq[M-1] += M

    ###############################################
    # Creating vertices of Vertices of degree > M #
    ###############################################

    for M >= current_box > 0:

        # If there is not enough vertices in the boxes available
        if taken + (seq[current_box] - 1) + seq[current_box-1] <= M:
            taken += seq[current_box]
            seq[current_box+1] += seq[current_box]
            seq[current_box] = 0
            continue

        # The degree of the new vertex will be taken + i + j where :
        #
        # * i is the number of vertices taken in the *current* box
        # * j the number of vertices taken in the *previous* one

        n_current_box = seq[current_box]
        n_previous_box = seq[current_box-1]

        # Note to self, and others :
        #
        # In the following lines, there are many incrementation/decrementation
        # that *may* be replaced by only +1 and -1 and save some
        # instructions. This would involve adding several "if", and I feared it
        # would make the code even uglier. If you are willing to give it a try,
        # **please check the results** ! It is trickier that it seems ! Even
        # changing the lower bounds in the for loops would require tests
        # afterwards.

        for max(0,((M+1)-n_previous_box-taken)) <= i < n_current_box:
            seq[current_box] -= i
            seq[current_box+1] += i

            for max(0,((M+1)-taken-i)) <= j <= n_previous_box:
                seq[current_box-1] -= j
                seq[current_box] += j

                new_vertex = taken + i + j
                seq[new_vertex] += 1
                enum(k+1,new_vertex)
                seq[new_vertex] -= 1

                seq[current_box-1] += j
                seq[current_box] -= j

            seq[current_box] += i
            seq[current_box+1] -= i

        taken += n_current_box
        seq[current_box] = 0
        seq[current_box+1] += n_current_box

    # Corner case
    #
    # Now current_box = 0. All the vertices of nonzero degree are taken, we just
    # want to know how many vertices of degree 0 will be neighbors of the new
    # vertex.
    for max(0,((M+1)-taken)) <= i <= seq[0]:

        seq[1] += i
        seq[0] -= i
        seq[taken+i] += 1

        enum(k+1, taken+i)

        seq[taken+i] -= 1
        seq[1] -= i
        seq[0] += i

    # Shift everything back to normal ! ( cell N is always equal to 0)
    for 1 <= i < N:
        seq[i] = seq[i+1]

    sig_off()
