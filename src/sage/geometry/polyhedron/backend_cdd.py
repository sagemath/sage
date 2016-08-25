"""
The cdd backend for polyhedral computations
"""
from __future__ import print_function
from __future__ import absolute_import

from subprocess import Popen, PIPE
from sage.rings.all import ZZ, QQ, RDF
from sage.misc.all import SAGE_TMP, tmp_filename, union, cached_method, prod
from sage.matrix.constructor import matrix

from .base import Polyhedron_base
from .base_QQ import Polyhedron_QQ
from .base_RDF import Polyhedron_RDF



#########################################################################
class Polyhedron_cdd(Polyhedron_base):
    """
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

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        EXAMPLES::

            sage: Polyhedron(vertices=[(0,0)], rays=[(1,1)],
            ...              lines=[(1,-1)], backend='cdd', base_ring=QQ)  # indirect doctest
            A 2-dimensional polyhedron in QQ^2 defined as the
            convex hull of 1 vertex, 1 ray, 1 line
        """
        from .cdd_file_format import cdd_Vrepresentation
        s = cdd_Vrepresentation(self._cdd_type, vertices, rays, lines)
        self._init_from_cdd_input(s, '--reps', verbose)


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

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        EXAMPLES::

            sage: Polyhedron(ieqs=[(0,1,1)], eqns=[(0,1,-1)],
            ...              backend='cdd', base_ring=QQ)  # indirect doctest
            A 1-dimensional polyhedron in QQ^2 defined as the
            convex hull of 1 vertex and 1 ray
        """
        from .cdd_file_format import cdd_Hrepresentation
        s = cdd_Hrepresentation(self._cdd_type, ieqs, eqns)
        self._init_from_cdd_input(s, '--reps', verbose)


    def _init_facet_adjacency_matrix(self, verbose=False):
        """
        Compute the facet adjacency matrix in case it has not been
        computed during initialization.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], backend='cdd', base_ring=QQ)
            sage: '_H_adjacency_matrix' in p.__dict__
            False
            sage: p._init_facet_adjacency_matrix()
            sage: p._H_adjacency_matrix
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        self._init_from_cdd_input(self.cdd_Hrepresentation(),
                                  '--adjacency', verbose)


    def _init_vertex_adjacency_matrix(self, verbose=False):
        """
        Compute the vertex adjacency matrix in case it has not been
        computed during initialization.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], backend='cdd', base_ring=QQ)
            sage: '_V_adjacency_matrix' in p.__dict__
            False
            sage: p._init_vertex_adjacency_matrix()
            sage: p._V_adjacency_matrix
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        self._init_from_cdd_input(self.cdd_Vrepresentation(),
                                  '--adjacency', verbose)


    def _init_from_cdd_input(self, cdd_input_string, cmdline_arg='--all', verbose=False):
        """
        Internal method: run cdd on a cdd H- or V-representation
        and initialize ourselves with the output.

        TESTS::

            sage: p = Polyhedron(vertices=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],
            ...                  backend='cdd', base_ring=QQ)
            sage: from sage.geometry.polyhedron.cdd_file_format import cdd_Vrepresentation
            sage: s = cdd_Vrepresentation('rational', [[0,0,1],[0,1,0],[1,0,0]], [], [])
            sage: p._init_from_cdd_input(s)
            sage: p.dim()
            2

            sage: point_list = [[0.132, -1.028, 0.028],[0.5, 0.5, -1.5],
            ....:    [-0.5, 1.5, -0.5],[0.5, 0.5, 0.5],[1.5, -0.5, -0.5],
            ....:    [-0.332, -0.332, -0.668],[-1.332, 0.668, 0.332],
            ....:    [-0.932, 0.068, 0.932],[-0.38, -0.38, 1.38],
            ....:    [-0.744, -0.12, 1.12],[-0.7781818182, -0.12, 0.9490909091],
            ....:    [0.62, -1.38, 0.38],[0.144, -1.04, 0.04],
            ....:    [0.1309090909, -1.0290909091, 0.04]]
            sage: Polyhedron(point_list)
            A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 14 vertices
            sage: Polyhedron(point_list, base_ring=QQ)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 14 vertices
        """
        if verbose:
            print('---- CDD input -----')
            print(cdd_input_string)

        cdd_proc = Popen([self._cdd_executable, cmdline_arg],
                         stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = cdd_proc.communicate(input=cdd_input_string)

        if verbose:
            print('---- CDD output -----')
            print(ans)
            print(err)
        if 'Error:' in ans + err:
            # cdd reports errors on stdout and misc information on stderr
            raise ValueError(ans.strip())
        self._init_from_cdd_output(ans)


    def _init_from_cdd_output(self, cdd_output_string):
        """
        Initialize ourselves with the output from cdd.

        TESTS::

            sage: from sage.geometry.polyhedron.cdd_file_format import cdd_Vrepresentation
            sage: s = cdd_Vrepresentation('rational',[[0,0],[1,0],[0,1],[1,1]], [], [])
            sage: from subprocess import Popen, PIPE
            sage: cdd_proc = Popen(['cdd_both_reps_gmp', '--all'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
            sage: ans, err = cdd_proc.communicate(input=s)
            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,1],[1,1]], backend='cdd', base_ring=QQ)
            sage: p._init_from_cdd_output(ans)
            sage: p.vertices()
            (A vertex at (0, 0), A vertex at (1, 0), A vertex at (0, 1), A vertex at (1, 1))
        """
        cddout=cdd_output_string.splitlines()
        suppressed_vertex = False   # whether cdd suppressed the vertex in output
        parent = self.parent()

        # nested function
        def expect_in_cddout(expected_string):
            l = cddout.pop(0).strip()
            if l!=expected_string:
                raise ValueError('Error while parsing cdd output: expected "'
                                   +expected_string+'" but got "'+l+'".\n')
        # nested function
        def cdd_linearities():
            l = cddout[0].split()
            if l[0] != "linearity":
                return []
            cddout.pop(0)
            assert len(l) == int(l[1])+2, "Not enough linearities given"
            return [int(i)-1 for i in l[2:]]  # make indices pythonic

        # nested function
        def cdd_convert(string, field=self.field()):
            """
            Converts the cdd output string to a numerical value.
            """
            return [field(x) for x in string.split()]

        # nested function
        def find_in_cddout(expected_string):
            """
            Find the expected string in a list of strings, and
            truncates ``cddout`` to start at that point. Returns
            ``False`` if search fails.
            """
            for pos in range(0,len(cddout)):
                l = cddout[pos].strip();
                if l==expected_string:
                    # must not assign to cddout in nested function
                    for i in range(0,pos+1):
                        cddout.pop(0)
                    return True
            return False

        if find_in_cddout('V-representation'):
            self._Vrepresentation = []
            lines = cdd_linearities()
            expect_in_cddout('begin')
            l = cddout.pop(0).split()
            assert self.ambient_dim() == int(l[1])-1,  "Different ambient dimension?"
            suppressed_vertex = True
            for i in range(int(l[0])):
                l = cddout.pop(0).strip()
                l_type = l[0]
                l = l[1:]
                if i in lines:
                    parent._make_Line(self, cdd_convert(l));
                elif l_type == '0':
                    parent._make_Ray(self, cdd_convert(l));
                else:
                    parent._make_Vertex(self, cdd_convert(l));
                    suppressed_vertex = False
            if suppressed_vertex and self.n_Vrepresentation()>0:
                # cdd does not output the vertex if it is only the origin
                parent._make_Vertex(self, [0] * self.ambient_dim())
            self._Vrepresentation = tuple(self._Vrepresentation)
            expect_in_cddout('end')

        if find_in_cddout('H-representation'):
            self._Hrepresentation = []
            equations = cdd_linearities()
            expect_in_cddout('begin')
            l = cddout.pop(0).split()
            assert self.ambient_dim() == int(l[1])-1, "Different ambient dimension?"
            for i in range(int(l[0])):
                l = cddout.pop(0)
                if i in equations:
                    parent._make_Equation(self, cdd_convert(l));
                else:
                    parent._make_Inequality(self, cdd_convert(l));
            self._Hrepresentation = tuple(self._Hrepresentation)
            expect_in_cddout('end')

        # nested function
        def cdd_adjacencies():
            l = cddout.pop(0).split()
            assert l[2] == ':', "Not a line of the adjacency data?"
            return [int(i)-1 for i in l[3:]]

        if find_in_cddout('Vertex graph'):
            n = len(self._Vrepresentation);
            if suppressed_vertex:
                n_cdd=n-1;
            else:
                n_cdd=n;
            self._V_adjacency_matrix = matrix(ZZ, n, n, 0)
            expect_in_cddout('begin')
            l = cddout.pop(0).split()
            assert int(l[0]) == n_cdd, "Not enough V-adjacencies in cdd output?"
            for i in range(n_cdd):
                for a in cdd_adjacencies():
                    self._V_adjacency_matrix[i,a] = 1
                # cdd reports that lines are never adjacent to anything.
                # I disagree, they are adjacent to everything!
                if self._Vrepresentation[i].is_line():
                    for j in range(n):
                        self._V_adjacency_matrix[i,j] = 1
                        self._V_adjacency_matrix[j,i] = 1
                    self._V_adjacency_matrix[i,i] = 0
            if suppressed_vertex: # cdd implied that there is only one vertex
                for i in range(n-1):
                    self._V_adjacency_matrix[i,n-1] = 1
                    self._V_adjacency_matrix[n-1,i] = 1
            self._V_adjacency_matrix.set_immutable()
            expect_in_cddout('end')

        if find_in_cddout('Facet graph'):
            n = len(self._Hrepresentation);
            self._H_adjacency_matrix = matrix(ZZ, n, n, 0)
            expect_in_cddout('begin')
            l = cddout.pop(0).split()
            assert int(l[0]) == n, "Not enough H-adjacencies in cdd output?"
            for i in range(n):
                for a in cdd_adjacencies():
                    self._H_adjacency_matrix[i,a] = 1
            self._H_adjacency_matrix.set_immutable()
            expect_in_cddout('end')


#########################################################################
class Polyhedron_QQ_cdd(Polyhedron_cdd, Polyhedron_QQ):
    """
    Polyhedra over QQ with cdd

    INPUT:

    - ``parent`` -- the parent, an instance of
      :class:`~sage.geometry.polyhedron.parent.Polyhedra`.

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: parent = Polyhedra(QQ, 2, backend='cdd')
        sage: from sage.geometry.polyhedron.backend_cdd import Polyhedron_QQ_cdd
        sage: Polyhedron_QQ_cdd(parent, [ [(1,0),(0,1),(0,0)], [], []], None, verbose=False)
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
    """

    _cdd_type = 'rational'

    _cdd_executable = 'cdd_both_reps_gmp'

    def __init__(self, parent, Vrep, Hrep, **kwds):
        """
        The Python constructor.

        See :class:`Polyhedron_base` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron(backend='cdd', base_ring=QQ)
            sage: type(p)
            <class 'sage.geometry.polyhedron.backend_cdd.Polyhedra_QQ_cdd_with_category.element_class'>
            sage: TestSuite(p).run()
        """
        Polyhedron_cdd.__init__(self, parent, Vrep, Hrep, **kwds)


#########################################################################
class Polyhedron_RDF_cdd(Polyhedron_cdd, Polyhedron_RDF):
    """
    Polyhedra over RDF with cdd

    INPUT:

    - ``ambient_dim`` -- integer. The dimension of the ambient space.

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: parent = Polyhedra(RDF, 2, backend='cdd')
        sage: from sage.geometry.polyhedron.backend_cdd import Polyhedron_RDF_cdd
        sage: Polyhedron_RDF_cdd(parent, [ [(1,0),(0,1),(0,0)], [], []], None, verbose=False)
        A 2-dimensional polyhedron in RDF^2 defined as the convex hull of 3 vertices
    """
    _cdd_type = 'real'

    _cdd_executable = 'cdd_both_reps'

    def __init__(self, parent, Vrep, Hrep, **kwds):
        """
        The Python constructor.

        See :class:`Polyhedron_base` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron(backend='cdd', base_ring=RDF)
            sage: type(p)
            <class 'sage.geometry.polyhedron.backend_cdd.Polyhedra_RDF_cdd_with_category.element_class'>
            sage: TestSuite(p).run()
        """
        Polyhedron_cdd.__init__(self, parent, Vrep, Hrep, **kwds)

