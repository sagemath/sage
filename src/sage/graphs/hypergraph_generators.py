r"""
Hypergraph generators

At the moment this module only implement one method, which calls Brendan McKay's
Nauty (`<http://cs.anu.edu.au/~bdm/nauty/>`_) to enumerate hypergraphs up to
isomorphism.
"""
from __future__ import print_function

class HypergraphGenerators():
    r"""
    A class consisting of constructors for common hypergraphs.
    """

    def nauty(self, number_of_sets, number_of_vertices,
              multiple_sets = False,
              vertex_min_degree = None, vertex_max_degree = None,
              set_max_size = None, set_min_size = None,
              regular = False, uniform = False,
              max_intersection = None,
              connected = False,
              options="", debug=False):
        r"""
        Enumerates hypergraphs up to isomorphism using Nauty.

        INPUT:

        - ``number_of_sets``, ``number_of_vertices`` (integers)

        - ``multiple_sets`` (boolean) -- whether to allow several sets
          of the hypergraph to be equal (set to ``False`` by default).

        - ``vertex_min_degree``, ``vertex_max_degree`` (integers) -- define the
          maximum and minimum degree of an element from the ground set (i.e. the
          number of sets which contain it). Set to ``None`` by default.

        - ``set_min_size``, ``set_max_size`` (integers) -- define the maximum
          and minimum size of a set. Set to ``None`` by default.

        - ``regular`` (integer) -- if set to an integer value `k`, requires the
          hypergraphs to be `k`-regular. It is actually a shortcut for the
          corresponing min/max values.

        - ``uniform`` (integer) -- if set to an integer value `k`, requires the
          hypergraphs to be `k`-uniform. It is actually a shortcut for the
          corresponing min/max values.

        - ``max_intersection`` (integer) -- constraints the maximum cardinality
          of the intersection of two sets fro the hypergraphs. Set to ``None``
          by default.

        - ``connected`` (boolean) -- whether to require the hypergraphs to be
          connected. Set to ``False`` by default.

        - ``debug`` (boolean) -- if ``True`` the first line of genbg's output to
          standard error is captured and the first call to the generator's
          ``next()`` function will return this line as a string.  A line leading
          with ">A" indicates a successful initiation of the program with some
          information on the arguments, while a line beginning with ">E"
          indicates an error with the input.

        - ``options`` (string) -- anything else that should be forwarded as
          input to Nauty's genbg. See its documentation for more information :
          `<http://cs.anu.edu.au/~bdm/nauty/>`_.

          .. NOTE::

              For genbg the *first class* elements are vertices, and *second
              class* elements are the hypergraph's sets.

        OUTPUT:

        A tuple of tuples.

        EXAMPLES:

        Small hypergraphs::

            sage: list(hypergraphs.nauty(4,2))
            [((), (0,), (1,), (0, 1))]

        Only connected ones::

            sage: list(hypergraphs.nauty(2,2, connected = True))
            [((0,), (0, 1))]

        Non-empty sets only::

            sage: list(hypergraphs.nauty(3,2, set_min_size = 1))
            [((0,), (1,), (0, 1))]

        The Fano Plane, as the only 3-uniform hypergraph with 7 sets and 7
        vertices::

            sage: fano = next(hypergraphs.nauty(7, 7, uniform=3, max_intersection=1))
            sage: print(fano)
            ((0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (2, 4, 5), (2, 3, 6), (1, 4, 6))

        The Fano Plane, as the only 3-regular hypergraph with 7 sets and 7
        vertices::

            sage: fano = next(hypergraphs.nauty(7, 7, regular=3, max_intersection=1))
            sage: print(fano)
            ((0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (2, 4, 5), (2, 3, 6), (1, 4, 6))
        """
        import subprocess

        nauty_input = options

        if connected:
            nauty_input += " -c"

        if not multiple_sets:
            nauty_input += " -z"

        if not max_intersection is None:
            nauty_input += " -Z"+str(max_intersection)

        # degrees and sizes
        if not regular is False:
            vertex_max_degree = vertex_min_degree = regular
        if vertex_max_degree is None:
            vertex_max_degree = number_of_sets
        if vertex_min_degree is None:
            vertex_min_degree = 0

        if not uniform is False:
            set_max_size = set_min_size = uniform
        if set_max_size is None:
            set_max_size = number_of_vertices
        if set_min_size is None:
            set_min_size = 0

        nauty_input += " -d"+str(vertex_min_degree)+":"+str(set_min_size)
        nauty_input += " -D"+str(vertex_max_degree)+":"+str(set_max_size)


        nauty_input +=  " "+str(number_of_vertices) +" "+str(number_of_sets)+" "

        sp = subprocess.Popen("genbg {0}".format(nauty_input), shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, close_fds=True)

        if debug:
            yield sp.stderr.readline()

        gen = sp.stdout
        total = number_of_sets + number_of_vertices
        while True:
            try:
                s = next(gen)
            except StopIteration:
                raise StopIteration("Exhausted list of graphs from nauty geng")

            from sage.graphs.graph import Graph
            G = Graph(s[:-1], format='graph6')

            yield tuple( tuple( x for x in G.neighbors(v)) for v in range(number_of_vertices, total))

    def CompleteUniform(self, n, k):
        r"""
        Return the complete `k`-uniform hypergraph on `n` points.

        INPUT:

        - ``k,n`` -- nonnegative integers with `k\leq n`

        EXAMPLE::

            sage: h = hypergraphs.CompleteUniform(5,2); h
            Incidence structure with 5 points and 10 blocks
            sage: len(h.packing())
            2
        """
        from sage.combinat.designs.incidence_structures import IncidenceStructure
        from itertools import combinations
        return IncidenceStructure(list(combinations(range(n),k)))

hypergraphs = HypergraphGenerators()
