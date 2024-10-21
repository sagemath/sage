# sage_setup: distribution = sagemath-modules
# sage.doctest: optional - numpy
"""
Hidden Markov Models -- Utility functions

AUTHOR:

   - William Stein, 2010-03
"""

#############################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL) v2+.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


from sage.structure.element import Matrix
from sage.misc.flatten  import flatten

cdef class HMM_Util:
    """
    A class used in order to share cdef's methods between different files.
    """
    cpdef normalize_probability_TimeSeries(self, TimeSeries T, Py_ssize_t i, Py_ssize_t j):
        """
        This function is used internally by the Hidden Markov Models code.

        Replace entries of T[i:j] in place so that they are all
        nonnegative and sum to 1.  Negative entries are replaced by 0 and
        T[i:j] is then rescaled to ensure that the sum of the entries in
        each row is equal to 1.  If all entries are 0, replace them
        by 1/(j-i).

        INPUT:

        - ``T`` -- a :class:`TimeSeries`
        - ``i`` -- nonnegative integer
        - ``j`` -- nonnegative integer

        OUTPUT: ``T`` is modified

        EXAMPLES::

            sage: import sage.stats.hmm.util
            sage: T = stats.TimeSeries([.1, .3, .7, .5])
            sage: u = sage.stats.hmm.util.HMM_Util()
            sage: u.normalize_probability_TimeSeries(T,0,3)
            sage: T
            [0.0909, 0.2727, 0.6364, 0.5000]
            sage: u.normalize_probability_TimeSeries(T,0,4)
            sage: T
            [0.0606, 0.1818, 0.4242, 0.3333]
            sage: abs(T.sum()-1) < 1e-8    # might not exactly equal 1 due to rounding
            True
        """
        # One single bounds check only
        if i < 0 or j < 0 or i > T._length or j > T._length:
            raise IndexError

        if j-i <= 0:
            # Nothing to do
            return

        cdef Py_ssize_t k

        # Replace negative entries by 0, summing entries
        cdef double s = 0, t
        for k in range(i,j):
            if T._values[k] < 0:
                T._values[k] = 0
            else:
                s += T._values[k]

        if s == 0:
            # If sum is 0, make all entries 1/(j-i).
            t = 1.0/(j-i)
            for k in range(i,j):
                T._values[k] = t
        else:
            # Normalise so sum is 1.
            for k in range(i, j):
                T._values[k] /= s

    cpdef TimeSeries initial_probs_to_TimeSeries(self, pi, bint normalize):
        """
        This function is used internally by the __init__ methods of
        various Hidden Markov Models.

        INPUT:

        - ``pi`` -- vector, list, or :class:`TimeSeries`
        - ``normalize`` -- if ``True``, replace negative entries by 0 and
          rescale to ensure that the sum of the entries in each row is
          equal to 1.  If the sum of the entries in a row is 0, replace them
          all by `1/N`.

        OUTPUT: a :class:`TimeSeries` of length `N`

        EXAMPLES::

            sage: import sage.stats.hmm.util
            sage: u = sage.stats.hmm.util.HMM_Util()
            sage: u.initial_probs_to_TimeSeries([0.1,0.2,0.9], True)
            [0.0833, 0.1667, 0.7500]
            sage: u.initial_probs_to_TimeSeries([0.1,0.2,0.9], False)
            [0.1000, 0.2000, 0.9000]
        """
        cdef TimeSeries T
        if isinstance(pi, TimeSeries):
            T = pi
        else:
            if not isinstance(pi, list):
                pi = list(pi)
            T = TimeSeries(pi)
        if normalize:
            # Now normalize
            self.normalize_probability_TimeSeries(T, 0, T._length)
        return T

    cpdef TimeSeries state_matrix_to_TimeSeries(self, A, int N, bint normalize):
        """
        This function is used internally by the ``__init__`` methods of
        Hidden Markov Models to make a transition matrix from ``A``.


        INPUT:

        - ``A`` -- matrix, list, list of lists, or :class:`TimeSeries`
        - ``N`` -- number of states
        - ``normalize`` -- if ``True``, replace negative entries by 0 and
          rescale to ensure that the sum of the entries in each row is
          equal to 1.  If the sum of the entries in a row is 0, replace them
          all by `1/N`.

        OUTPUT: a :class:`TimeSeries`

        EXAMPLES::

            sage: import sage.stats.hmm.util
            sage: u = sage.stats.hmm.util.HMM_Util()
            sage: u.state_matrix_to_TimeSeries([[.1,.7],[3/7,4/7]], 2, True)
            [0.1250, 0.8750, 0.4286, 0.5714]
            sage: u.state_matrix_to_TimeSeries([[.1,.7],[3/7,4/7]], 2, False)
            [0.1000, 0.7000, 0.4286, 0.5714]
        """
        cdef TimeSeries T
        if isinstance(A, TimeSeries):
            T = A
        elif isinstance(A, Matrix):
            T = TimeSeries(A.list())
        elif isinstance(A, list):
            T = TimeSeries(flatten(A))
        else:
            T = TimeSeries(A)
        cdef Py_ssize_t i
        if normalize:
            # Set to 0 negative rows and make sure sum of entries in each
            # row is 1.
            if len(T) != N*N:
                raise ValueError("number of entries of transition matrix A must be the square of the number of entries of pi")
            for i in range(N):
                self.normalize_probability_TimeSeries(T, i*N, (i+1)*N)
        return T
