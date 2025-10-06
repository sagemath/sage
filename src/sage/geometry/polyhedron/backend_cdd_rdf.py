# sage.doctest: needs sage.rings.real_double

r"""
The cdd backend for polyhedral computations, floating point version
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


from .backend_cdd import Polyhedron_cdd
from .base_RDF import Polyhedron_RDF


class Polyhedron_RDF_cdd(Polyhedron_cdd, Polyhedron_RDF):
    """
    Polyhedra over RDF with cdd.

    INPUT:

    - ``ambient_dim`` -- integer; the dimension of the ambient space

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``

    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: parent = Polyhedra(RDF, 2, backend='cdd')
        sage: from sage.geometry.polyhedron.backend_cdd_rdf import Polyhedron_RDF_cdd
        sage: Polyhedron_RDF_cdd(parent, [ [(1,0),(0,1),(0,0)], [], []], None, verbose=False)
        A 2-dimensional polyhedron in RDF^2 defined as the convex hull of 3 vertices

    TESTS:

    Checks that :issue:`24877` is fixed::

        sage: n1 = 1045602428815736513789288687833080060779
        sage: n2 = 76591188009721216624438400001815308369088648782156930777145
        sage: n3 = 141046287872967162025203834781636948939209065735662536571684677443277621519222367249160281646288602157866548267640061850035
        sage: n4 = 169296796161110084211548448622149955145002732358082778064645608216077666698460018565094060494217
        sage: verts = [[159852/261157, 227425/261157],
        ....:  [9/10, 7/10],
        ....:  [132/179, 143/179],
        ....:  [8/11, -59/33],
        ....:  [174/167, 95/167],
        ....:  [3/2, -1/2],
        ....:  [-1162016360399650274197433414376009691155/n1,
        ....:    1626522696050475596930360993440360903664/n1],
        ....:  [-112565666321600055047037445519656973805313121630809713051718/n2,
        ....:    -15318574020578896781701071673537253327221557273483622682053/n2],
        ....:  [-222823992658914823798345935660863293259608796350232624336738123149601409997996952470726909468671437285720616325991022633438/n3,
        ....:   (-20857694835570598502487921801401627779907095024585170129381924208334510445828894861553290291713792691651471189597832832973*5)/n3],
        ....:  [-100432602675156818915933977983765863676402454634873648118147187022041830166292457614016362515164/n4,
        ....:   -429364759737031049317769174492863890735634068814210512342503744054527903830844433491149538512537/n4]]
        sage: P = Polyhedron(verts, base_ring=RDF)
        sage: len(P.faces(1))
        10
        sage: P.n_vertices()
        10
        sage: P.n_facets()
        10

    Check that :issue:`19803` is fixed::

        sage: from sage.geometry.polyhedron.parent import Polyhedra
        sage: P_cdd = Polyhedra(RDF, 3, 'cdd')
        sage: P_cdd([[],[],[]], None)
        The empty polyhedron in RDF^3
        sage: Polyhedron(vertices=[], backend='cdd', base_ring=RDF)
        The empty polyhedron in RDF^0
    """
    _cdd_type = 'real'

    _cdd_executable = 'cddexec'

    def __init__(self, parent, Vrep, Hrep, **kwds):
        """
        The Python constructor.

        See :class:`Polyhedron_base` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron(backend='cdd', base_ring=RDF)
            sage: type(p)
            <class 'sage.geometry.polyhedron.parent.Polyhedra_RDF_cdd_with_category.element_class'>
            sage: TestSuite(p).run()
        """
        Polyhedron_cdd.__init__(self, parent, Vrep, Hrep, **kwds)

    def _init_from_Vrepresentation_and_Hrepresentation(self, Vrep, Hrep, verbose=False):
        """
        Construct polyhedron from Vrepresentation and Hrepresentation data.

        See :class:`Polyhedron_base` for a description of ``Vrep`` and ``Hrep``.

        .. NOTE::

            The representation is assumed to be correct.

            As long as cdd can obtain a consistent object with Vrepresentation
            or Hrepresentation no warning is raised. Consistency is checked by
            comparing the output length of Vrepresentation and Hrepresentation
            with the input.

            In comparison, the "normal" initialization from Vrepresentation over RDF
            expects the output length to be consistent with the computed length
            when re-feeding cdd the outputted Hrepresentation.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra_RDF_cdd
            sage: from sage.geometry.polyhedron.backend_cdd_rdf import Polyhedron_RDF_cdd
            sage: parent = Polyhedra_RDF_cdd(RDF, 1, 'cdd')
            sage: Vrep = [[[0.0], [1.0]], [], []]
            sage: Hrep = [[[0.0, 1.0], [1.0, -1.0]], []]
            sage: p = Polyhedron_RDF_cdd(parent, Vrep, Hrep,                    # indirect doctest
            ....:                        Vrep_minimal=True, Hrep_minimal=True)
            sage: p
            A 1-dimensional polyhedron in RDF^1 defined as the convex hull of 2 vertices

        TESTS:

        Test that :issue:`29568` is fixed::

            sage: # needs sage.groups
            sage: P = polytopes.buckyball(exact=False)
            sage: Q = P + P.center()
            sage: P.is_combinatorially_isomorphic(Q)
            True
            sage: R = 2*P
            sage: P.is_combinatorially_isomorphic(R)
            True

        The polyhedron with zero inequalities works correctly; see :issue:`29899`::

            sage: Vrep = [[], [], [[1.0]]]
            sage: Hrep = [[], []]
            sage: p = Polyhedron_RDF_cdd(parent, Vrep, Hrep,                    # indirect doctest
            ....:                        Vrep_minimal=True, Hrep_minimal=True)
            sage: p
            A 1-dimensional polyhedron in RDF^1 defined as the convex hull of 1 vertex and 1 line

        Test that :issue:`30330` is fixed::

            sage: P1 = polytopes.regular_polygon(5, exact=False)
            sage: P2 = Polyhedron()
            sage: P1*P2
            The empty polyhedron in RDF^2
        """
        def parse_Vrep(intro, data):
            count = int(data[0][0])
            if count != len(vertices) + len(rays) + len(lines):
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

        def parse_Hrep(intro, data):
            count = int(data[0][0])
            infinite_count = len([d for d in data[1:] if d[0] == '1' and all(c == '0' for c in d[1:])])
            if count - infinite_count != len(ieqs) + len(eqns):
                # Upstream claims that nothing can be done about these
                # cases/that they are features not bugs. Imho, cddlib is
                # not really suitable for automatic parsing of its output,
                # the implementation backed by doubles has not really been
                # optimized for numerical stability, and makes some
                # somewhat random numerical choices. (But I am not an
                # expert in that field by any means.)
                from warnings import warn
                warn("This polyhedron data is numerically complicated; cdd could not convert between the inexact V and H representation without loss of data. The resulting object might show inconsistencies.")

        def try_init(rep):
            if rep == "Vrep":
                from .cdd_file_format import cdd_Vrepresentation
                s = cdd_Vrepresentation(self._cdd_type, vertices, rays, lines)
            else:
                # We have to add a trivial inequality, in case the polyhedron is the universe.
                new_ieqs = ieqs + ((1,) + tuple(0 for _ in range(self.ambient_dim())),)

                from .cdd_file_format import cdd_Hrepresentation
                s = cdd_Hrepresentation(self._cdd_type, new_ieqs, eqns)

            s = self._run_cdd(s, '--redcheck', verbose=verbose)
            s = self._run_cdd(s, '--repall', verbose=verbose)
            Polyhedron_cdd._parse_block(s.splitlines(), 'V-representation', parse_Vrep)
            Polyhedron_cdd._parse_block(s.splitlines(), 'H-representation', parse_Hrep)
            self._init_from_cdd_output(s)

        from warnings import catch_warnings, simplefilter

        vertices, rays, lines = (tuple(x) for x in Vrep)
        ieqs, eqns = (tuple(x) for x in Hrep)

        if not (vertices or rays or lines):
            # cdd refuses to handle empty polyhedra.
            self._init_empty_polyhedron()
            return

        # We prefer the shorter representation.
        # Note that for the empty polyhedron we prefer Hrepresentation.
        prim = "Hrep" if len(ieqs) <= len(vertices) + len(rays) else "Vrep"
        sec = "Vrep" if len(ieqs) <= len(vertices) + len(rays) else "Hrep"

        with catch_warnings():
            # Raise an error and try the other representation in case of
            # numerical inconsistency.
            simplefilter("error")
            try:
                try_init(prim)
            except UserWarning:
                simplefilter("once")  # Only print the first warning.
                try_init(sec)
