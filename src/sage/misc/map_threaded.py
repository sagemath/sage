"""
Threaded map function
"""

def map_threaded(function, sequence):
    """
    Apply the function to the elements in the sequence by threading
    recursively through all sub-sequences in the sequence.

    EXAMPLES::

        sage: map_threaded(log, [[1,2], [3,e]])
        [[0, log(2)], [log(3), 1]]
        sage: map_threaded(log, [(1,2), (3,e)])
        [[0, log(2)], [log(3), 1]]
        sage: map_threaded(N, [[1,2], [3,e]])
        [[1.00000000000000, 2.00000000000000], [3.00000000000000, 2.71828182845905]]
        sage: map_threaded((x^2).function(x), [[1,2,3,5], [2,10]])
        [[1, 4, 9, 25], [4, 100]]

    map_threaded also works on any object with an apply_map method, e.g.,
    on matrices::

        sage: map_threaded(lambda x: x^2, matrix([[1,2], [3,4]]))
        [ 1  4]
        [ 9 16]

    AUTHORS:

    - William Stein (2007-12); based on feedback from Peter Doyle.
    """
    if hasattr(sequence, 'apply_map'):
        return sequence.apply_map(function)
    return [map_threaded(function, x) if isinstance(x, (list, tuple)) \
             else function(x) for x in sequence]
