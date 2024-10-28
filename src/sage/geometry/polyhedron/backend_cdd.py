r"""
The cdd backend for polyhedral computations
"""
# ****************************************************************************
#       Copyright (C) 2011-2014 Volker Braun <vbraun.name@gmail.com>
#                     2018      Timo Kaufmann <timokau@zoho.com>
#                     2018      Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from subprocess import Popen, PIPE
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.features.cddlib import CddExecutable

from .base import Polyhedron_base
from .base_QQ import Polyhedron_QQ


class Polyhedron_cdd(Polyhedron_base):
    r"""
    Base class for the cdd backend.
    """
    def _init_from_Vrepresentation(self, vertices, rays, lines, verbose=False):
        """
        Construct polyhedron from V-representation data.

        INPUT:

        - ``vertices`` -- list of point. Each point can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``rays`` -- list of rays. Each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``lines`` -- list of lines. Each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: Polyhedron(vertices=[(0,0)], rays=[(1,1)],               # indirect doctest
            ....:            lines=[(1,-1)], backend='cdd', base_ring=QQ)
            A 2-dimensional polyhedron in QQ^2 defined as the
            convex hull of 1 vertex, 1 ray, 1 line
        """
        from .cdd_file_format import cdd_Vrepresentation
        s = cdd_Vrepresentation(self._cdd_type, vertices, rays, lines)
        s = self._run_cdd(s, '--redcheck', verbose=verbose)
        s = self._run_cdd(s, '--repall', verbose=verbose)
        self._init_from_cdd_output(s)
        if not self.base_ring().is_exact():
            # cdd's parser cannot handle the full output of --repall, so we
            # need to extract the first block before we feed it back into cdd
            s = s.splitlines()
            s = s[:s.index('end')+1]
            s = '\n'.join(s)
            t = self._run_cdd(s, '--rep', verbose=verbose)

            def parse(intro, data):
                count = int(data[0][0])
                if count != len(self._cdd_V_to_sage_V):
                    # Upstream claims that nothing can be done about these
                    # cases/that they are features not bugs. Imho, cddlib is
                    # not really suitable for automatic parsing of its output,
                    # the implementation backed by doubles has not really been
                    # optimized for numerical stability, and makes some
                    # somewhat random numerical choices. (But I am not an
                    # expert in that field by any means.) See also
                    # https://github.com/cddlib/cddlib/pull/7.
                    from warnings import warn
                    warn("This polyhedron data is numerically complicated; cdd could not convert between the inexact V and H representation without loss of data. The resulting object might show inconsistencies.")
            Polyhedron_cdd._parse_block(t.splitlines(), 'V-representation', parse)

    def _init_from_Hrepresentation(self, ieqs, eqns, verbose=False):
        """
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ieqs`` -- list of inequalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``eqns`` -- list of equalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: Polyhedron(ieqs=[(0,1,1)], eqns=[(0,1,-1)],  # indirect doctest
            ....:            backend='cdd', base_ring=QQ)
            A 1-dimensional polyhedron in QQ^2 defined as the
            convex hull of 1 vertex and 1 ray

        TESTS:

        The polyhedron with zero inequalities can be initialized from Hrepresentation;
        see :issue:`29899`::

            sage: Polyhedron(ieqs=[], ambient_dim=5, backend='cdd')
            A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 1 vertex and 5 lines
        """
        from .cdd_file_format import cdd_Hrepresentation
        # We have to add a trivial inequality, in case the polyhedron is the universe.
        ieqs = tuple(ieqs) + ((1,) + tuple(0 for _ in range(self.ambient_dim())),)
        s = cdd_Hrepresentation(self._cdd_type, ieqs, eqns)
        s = self._run_cdd(s, '--redcheck', verbose=verbose)
        s = self._run_cdd(s, '--repall', verbose=verbose)
        self._init_from_cdd_output(s)
        if not self.base_ring().is_exact():
            if len(self._Vrepresentation) == 0:
                # cdd (reasonably) refuses to handle empty polyhedra, so we
                # skip this check
                return
            # cdd's parser cannot handle the full output of --repall, so we
            # need to extract the first block before we feed it back into cdd
            s = s.splitlines()
            s = s[:s.index('end')+1]
            s = '\n'.join(s)
            t = self._run_cdd(s, '--rep', verbose=verbose)

            def parse(intro, data):
                count = int(data[0][0])
                infinite_count = len([d for d in data[1:] if d[0] == '1' and all(c == '0' for c in d[1:])])
                if count - infinite_count != len(self._Hrepresentation):
                    # Upstream claims that nothing can be done about these
                    # cases/that they are features not bugs. Imho, cddlib is
                    # not really suitable for automatic parsing of its output,
                    # the implementation backed by doubles has not really been
                    # optimized for numerical stability, and makes some
                    # somewhat random numerical choices. (But I am not an
                    # expert in that field by any means.)
                    from warnings import warn
                    warn("This polyhedron data is numerically complicated; cdd could not convert between the inexact V and H representation without loss of data. The resulting object might show inconsistencies.")
            Polyhedron_cdd._parse_block(t.splitlines(), 'H-representation', parse)

    def _run_cdd(self, cdd_input_string, cmdline_arg, verbose=False):
        if verbose:
            print('---- CDD input -----')
            print(cdd_input_string)

        cdd_proc = Popen([CddExecutable(self._cdd_executable).absolute_filename(), cmdline_arg],
                         stdin=PIPE, stdout=PIPE, stderr=PIPE,
                         encoding='latin-1')
        ans, err = cdd_proc.communicate(input=cdd_input_string)

        if verbose:
            print('---- CDD output -----')
            print(ans)
            print(err)
        if 'Error:' in ans + err:
            # cdd reports errors on stdout and misc information on stderr
            raise ValueError(ans.strip())
        return ans

    @classmethod
    def _parse_block(cls, cddout, header, parser):
        r"""
        Parse a block of cdd data identified by ``header`` by invoking
        ``parser`` on it.

        EXAMPLES::

            sage: cddout = r'''
            ....: unrelated
            ....: HEADER
            ....: intro 0 1 2
            ....: begin
            ....: data 0 1 2
            ....: data 3 4 5
            ....: end
            ....: unrelated
            ....: '''.splitlines()
            sage: from sage.geometry.polyhedron.backend_cdd import Polyhedron_cdd
            sage: def parser(intro, data):
            ....:     print("INTRO:", intro)
            ....:     print("DATA:", data)
            sage: Polyhedron_cdd._parse_block(cddout, 'HEADER', parser)
            INTRO: [['intro', '0', '1', '2']]
            DATA: [['data', '0', '1', '2'], ['data', '3', '4', '5']]
        """
        try:
            block = cddout[cddout.index(header)+1:]
        except ValueError:
            # section is missing in the cdd output
            return

        intro = block[:block.index('begin')]
        intro = [i.strip().split() for i in intro]
        data = block[block.index('begin')+1:block.index('end')]
        data = [d.strip().split() for d in data]
        parser(intro, data)

    def _init_from_cdd_output(self, cddout):
        """
        Initialize ourselves with the output from cdd.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,1],[1,1]], backend='cdd', base_ring=QQ)  # indirect doctest
            sage: p.vertices()
            (A vertex at (0, 0), A vertex at (1, 0), A vertex at (0, 1), A vertex at (1, 1))

        Check that :issue:`29176` is fixed::

            sage: e = [[11582947.657000002, 5374.38, 4177.06, 1.0], [11562795.9322, 5373.62, 4168.38, 1.0]]
            sage: p = Polyhedron(ieqs=e); p
            A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 1 vertex, 2 rays, 1 line
            sage: p.incidence_matrix()
            [1 1]
            [1 0]
            [0 1]
            [1 1]

            sage: P = [[-2687.19, -2088.53], [-2686.81, -2084.19]]
            sage: V = VoronoiDiagram(P)
            sage: R = V.regions()
            sage: V.points()[0], R[V.points()[0]]
            (P(-2687.19000000000, -2088.53000000000),
             A 2-dimensional polyhedron in RDF^2 defined as the convex hull of 1 vertex, 1 ray, 1 line)
            sage: V.points()[1], R[V.points()[1]]
            (P(-2686.81000000000, -2084.19000000000),
             A 2-dimensional polyhedron in RDF^2 defined as the convex hull of 1 vertex, 1 ray, 1 line)

        Check that :issue:`31253` is fixed::

            sage: P = polytopes.permutahedron(2, backend='cdd')
            sage: P.Hrepresentation()
            (An inequality (0, 1) x - 1 >= 0,
             An inequality (1, 0) x - 1 >= 0,
             An equation (1, 1) x - 3 == 0)
            sage: Q = Polyhedron(P.vertices(), backend='cdd')
            sage: Q.Hrepresentation()
            (An inequality (-1, 0) x + 2 >= 0,
             An inequality (1, 0) x - 1 >= 0,
             An equation (1, 1) x - 3 == 0)
            sage: [x.ambient_Hrepresentation() for x in P.facets()]
            [(An inequality (1, 0) x - 1 >= 0, An equation (1, 1) x - 3 == 0),
             (An inequality (0, 1) x - 1 >= 0, An equation (1, 1) x - 3 == 0)]
        """
        cddout = cddout.splitlines()

        def parse_indices(count, cdd_indices, cdd_indices_to_sage_indices=None):
            cdd_indices = [int(x) for x in cdd_indices]
            if cdd_indices_to_sage_indices is None:
                cdd_indices_to_sage_indices = {i: i-1 for i in cdd_indices}
            if count < 0:
                assert cdd_indices_to_sage_indices is not None, "Did not expect negative counts here"
                count = -count
                cdd_indices = list(set(cdd_indices_to_sage_indices.keys()) - set(cdd_indices))
                assert count in [len(cdd_indices), len(cdd_indices) - 1]
            assert count == len(cdd_indices)
            return [cdd_indices_to_sage_indices[i] for i in cdd_indices if cdd_indices_to_sage_indices[i] is not None]

        def parse_linearities(intro):
            for entries in intro:
                if entries and entries.pop(0) == 'linearity':
                    return parse_indices(int(entries.pop(0)), entries)
            return []

        def parse_H_representation(intro, data):
            if '_Hrepresentation' in self.__dict__:
                raise NotImplementedError("cannot replace internal representation as this breaks caching")
            self._Hrepresentation = []
            # we drop some entries in cdd's output and this changes the numbering; this dict keeps track of that
            self._cdd_H_to_sage_H = {}
            equations = parse_linearities(intro)
            data[0].pop(2)  # ignore data type, we know the base ring already
            count, dimension = map(int, data.pop(0))
            assert self.ambient_dim() == dimension - 1, "Unexpected ambient dimension"
            assert len(data) == count, "Unexpected number of lines"
            R = self.base_ring()
            from itertools import chain
            # We add equations to the end of the Hrepresentation.
            for i in chain(
                    (j for j in range(len(data)) if j not in equations),
                    equations):
                line = data[i]
                coefficients = [R(x) for x in line]
                if coefficients[0] != 0 and all(e == 0 for e in coefficients[1:]):
                    # cddlib sometimes includes an implicit plane at infinity: 1 0 0 ... 0
                    # We do not care about this entry.
                    self._cdd_H_to_sage_H[i+1] = None
                    continue

                self._cdd_H_to_sage_H[i+1] = len(self._Hrepresentation)
                if i in equations:
                    self.parent()._make_Equation(self, coefficients)
                else:
                    self.parent()._make_Inequality(self, coefficients)

            self._Hrepresentation = tuple(self._Hrepresentation)

        def parse_V_representation(intro, data):
            if '_Vrepresentation' in self.__dict__:
                raise NotImplementedError("cannot replace internal representation as this breaks caching")
            self._Vrepresentation = []
            # we drop some entries in cdd's output and this changes the numbering; this dict keeps track of that
            self._cdd_V_to_sage_V = {}
            lines = parse_linearities(intro)
            data[0].pop(2)  # ignore data type, we know the base ring already
            count, dimension = map(int, data.pop(0))
            assert self.ambient_dim() == dimension - 1, "Unexpected ambient dimension"
            assert len(data) == count, "Unexpected number of lines"
            has_vertex = False
            for i, line in enumerate(data):
                kind = line.pop(0)
                coefficients = map(self.base_ring(), line)
                self._cdd_V_to_sage_V[i+1] = len(self._Vrepresentation)
                if i in lines:
                    self.parent()._make_Line(self, coefficients)
                elif kind == '0':
                    self.parent()._make_Ray(self, coefficients)
                else:
                    self.parent()._make_Vertex(self, coefficients)
                    has_vertex = True
            if self._Vrepresentation and not has_vertex:
                # when the Polyhedron consists only of lines/rays from the
                # origin, cddlib does not output the single vertex at the
                # origin so we have to add it here as the Polyhedron class
                # expects it to be there.
                self.parent()._make_Vertex(self, [self.base_ring().zero()] * self.ambient_dim())
            self._Vrepresentation = tuple(self._Vrepresentation)

        def parse_adjacency(intro, data, M, N, cdd_indices_to_sage_indices, cdd_indices_to_sage_indices2=None):
            # This function is also used to parse the incidence matrix.
            if cdd_indices_to_sage_indices2 is None:
                cdd_indices_to_sage_indices2 = cdd_indices_to_sage_indices
            ret = matrix(ZZ, M, N, 0)
            data.pop(0)
            data.reverse()
            for adjacencies in data:
                assert adjacencies[2] == ':', "Not a line of adjacency data"
                cdd_vertex = int(adjacencies[0])
                count = int(adjacencies[1])

                # cdd sometimes prints implicit adjacencies for the plane at
                # infinity at the end of the output (even though it's not part
                # of the V/H representation) so we ignore indices that we do
                # not know about.
                if cdd_vertex not in cdd_indices_to_sage_indices:
                    cdd_indices_to_sage_indices[cdd_vertex] = None
                v = cdd_indices_to_sage_indices[cdd_vertex]
                if v is None:
                    continue
                for w in parse_indices(count, adjacencies[3:], cdd_indices_to_sage_indices2):
                    if w is None:
                        continue
                    ret[v, w] = 1
            return ret

        def parse_vertex_adjacency(intro, data):
            if '_V_adjacency_matrix' in self.__dict__:
                raise NotImplementedError("cannot replace internal representation as this breaks caching")
            N = len(self._Vrepresentation)
            self._V_adjacency_matrix = parse_adjacency(intro, data, N, N, self._cdd_V_to_sage_V)
            for i, v in enumerate(self._Vrepresentation):
                # cdd reports that lines are never adjacent to anything.
                # we disagree, they are adjacent to everything.
                if v.is_line():
                    for j in range(len(self._Vrepresentation)):
                        self._V_adjacency_matrix[i ,j] = 1
                        self._V_adjacency_matrix[j, i] = 1
                self._V_adjacency_matrix[i,i] = 0
            self._V_adjacency_matrix.set_immutable()
            self.vertex_adjacency_matrix.set_cache(self._V_adjacency_matrix)

        def parse_facet_adjacency(intro, data):
            if '_H_adjacency_matrix' in self.__dict__:
                raise NotImplementedError("cannot replace internal representation as this breaks caching")
            N = len(self._Hrepresentation)
            self._H_adjacency_matrix = parse_adjacency(intro, data, N, N, self._cdd_H_to_sage_H)
            self._H_adjacency_matrix.set_immutable()
            self.facet_adjacency_matrix.set_cache(self._H_adjacency_matrix)

        def parse_incidence_matrix(intro, data):
            if 'incidence_matrix' in self.__dict__:
                raise NotImplementedError("cannot replace internal representation as this breaks caching")
            N = len(self._Hrepresentation)
            M = len(self._Vrepresentation)
            inc_mat = parse_adjacency(intro, data, M, N, self._cdd_V_to_sage_V, self._cdd_H_to_sage_H)
            inc_mat.set_immutable()
            self.incidence_matrix.set_cache(inc_mat)

        Polyhedron_cdd._parse_block(cddout, 'H-representation', parse_H_representation)
        Polyhedron_cdd._parse_block(cddout, 'V-representation', parse_V_representation)
        Polyhedron_cdd._parse_block(cddout, 'Facet adjacency', parse_facet_adjacency)
        Polyhedron_cdd._parse_block(cddout, 'Vertex adjacency', parse_vertex_adjacency)
        Polyhedron_cdd._parse_block(cddout, 'Vertex incidence', parse_incidence_matrix)


class Polyhedron_QQ_cdd(Polyhedron_cdd, Polyhedron_QQ):
    """
    Polyhedra over QQ with cdd.

    INPUT:

    - ``parent`` -- the parent, an instance of
      :class:`~sage.geometry.polyhedron.parent.Polyhedra`

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``

    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: parent = Polyhedra(QQ, 2, backend='cdd')
        sage: from sage.geometry.polyhedron.backend_cdd import Polyhedron_QQ_cdd
        sage: Polyhedron_QQ_cdd(parent, [ [(1,0),(0,1),(0,0)], [], []], None, verbose=False)
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices

    TESTS:

    Check that :issue:`19803` is fixed::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: P_cdd = Polyhedra(QQ, 3, 'cdd')
        sage: P_cdd([[],[],[]], None)
        The empty polyhedron in QQ^3
        sage: Polyhedron(vertices=[], backend='cdd', base_ring=QQ)
        The empty polyhedron in QQ^0
    """

    _cdd_type = 'rational'

    _cdd_executable = 'cddexec_gmp'

    def __init__(self, parent, Vrep, Hrep, **kwds):
        """
        The Python constructor.

        See :class:`Polyhedron_base` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron(backend='cdd', base_ring=QQ)
            sage: type(p)
            <class 'sage.geometry.polyhedron.parent.Polyhedra_QQ_cdd_with_category.element_class'>
            sage: TestSuite(p).run()
        """
        Polyhedron_cdd.__init__(self, parent, Vrep, Hrep, **kwds)
