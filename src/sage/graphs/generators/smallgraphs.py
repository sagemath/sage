# -*- coding: utf-8 -*-
r"""
Small graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""
#*****************************************************************************
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


# import from Sage library
from sage.graphs.graph import Graph
from math import sin, cos, pi
from sage.graphs.graph_plot import _circle_embedding, _line_embedding

#######################################################################
#   Named Graphs
#######################################################################

def HarriesGraph(embedding=1):
    r"""
    Returns the Harries Graph.

    The Harries graph is a Hamiltonian 3-regular graph on 70
    vertices. See the :wikipedia:`Wikipedia page on the Harries
    graph <Harries_graph>`.

    The default embedding here is to emphasize the graph's 4 orbits.
    This graph actually has a funny construction. The following
    procedure gives an idea of it, though not all the adjacencies
    are being properly defined.

    #. Take two disjoint copies of a :meth:`Petersen graph
       <PetersenGraph>`. Their vertices will form an orbit of the
       final graph.

    #. Subdivide all the edges once, to create 15+15=30 new
       vertices, which together form another orbit.

    #. Create 15 vertices, each of them linked to 2 corresponding
       vertices of the previous orbit, one in each of the two
       subdivided Petersen graphs. At the end of this step all
       vertices from the previous orbit have degree 3, and the only
       vertices of degree 2 in the graph are those that were just
       created.

    #. Create 5 vertices connected only to the ones from the
       previous orbit so that the graph becomes 3-regular.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.HarriesGraph()
        sage: g.order()
        70
        sage: g.size()
        105
        sage: g.girth()
        10
        sage: g.diameter()
        6
        sage: g.show(figsize=[10, 10])   # long time
        sage: graphs.HarriesGraph(embedding=2).show(figsize=[10, 10])   # long time

    TESTS::

        sage: graphs.HarriesGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.

    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(70, [-29, -19, -13, 13, 21, -27, 27, 33, -13, 13,
                             19, -21, -33, 29], 5)
    g.name("Harries Graph")

    if embedding == 1:
        gpos = g.get_pos()
        ppos = PetersenGraph().get_pos()

        # The graph's four orbits
        o = [None]*4
        o[0] = [0, 2, 6, 8, 14, 16, 20, 22, 28, 30, 34, 36, 42, 44, 48, 50,
                56, 58, 62, 64]
        o[1] = [1, 3, 5, 7, 9, 13, 15, 17, 19, 21, 23, 27, 29, 31, 33, 35,
                37, 41, 43, 45, 47, 49, 51, 55, 57, 59, 61, 63, 65, 69]
        o[2] = [60, 10, 12, 4, 24, 26, 18, 38, 40, 32, 52, 54, 46, 66, 68]
        o[3] = [11, 25, 39, 53, 67]

        # Correspondence between the vertices of one of the two Petersen
        # graphs on o[0] and the vertices of a standard Petersen graph
        # object
        g_to_p = {0: 0, 2: 1, 42: 5, 44: 8, 14: 7, 16: 2, 56: 9, 58: 6,
                  28: 4, 30: 3}

        # Correspondence between the vertices of the other Petersen graph
        # on o[0] and the vertices of the first one
        g_to_g = {64: 44, 34: 0, 36: 28, 6: 2, 8: 58, 48: 16, 50: 30,
                  20: 14, 22: 56, 62: 42}

        # Position for the vertices from the first copy
        for v, i in g_to_p.iteritems():
            gpos[v] = ppos[i]

        # Position for the vertices in the second copy. Moves the first,
        # too.
        offset = 3.5
        for v, i in g_to_g.iteritems():
            x, y = gpos[i]
            gpos[v] = (x + offset*0.5, y)
            gpos[i] = (x - offset*0.5, y)

        # Vertices from o[1]. These are actually the "edges" of the
        # copies of Petersen.
        for v in o[1]:
            p1, p2 = [gpos[x] for x in g.neighbors(v) if x in o[0]]
            gpos[v] = ((p1[0] + p2[0])/2, (p1[1] + p2[1])/2)

        # 15 vertices from o[2]
        for i, v in enumerate(o[2]):
            gpos[v] = (-1.75 + i*.25, 2)

        # 5 vertices from o[3]
        for i, v in enumerate(o[3]):
            gpos[v] = (-1 + i*.5, 2.5)

        return g

    elif embedding == 2:
        return g
    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def HarriesWongGraph(embedding=1):
    r"""
    Returns the Harries-Wong Graph.

    See the :wikipedia:`Wikipedia page on the Harries-Wong graph
    <Harries-Wong_graph>`.

    *About the default embedding:*

    The default embedding is an attempt to emphasize the graph's
    8 (!!!) different orbits. In order to understand this better,
    one can picture the graph as being built in the following way:

        #. One first creates a 3-dimensional cube (8 vertices, 12
           edges), whose vertices define the first orbit of the
           final graph.

        #. The edges of this graph are subdivided once, to create 12
           new vertices which define a second orbit.

        #. The edges of the graph are subdivided once more, to
           create 24 new vertices giving a third orbit.

        #. 4 vertices are created and made adjacent to the vertices
           of the second orbit so that they have degree
           3. These 4 vertices also define a new orbit.

        #. In order to make the vertices from the third orbit
           3-regular (they all miss one edge), one creates a binary
           tree on 1 + 3 + 6 + 12 vertices. The leaves of this new
           tree are made adjacent to the 12 vertices of the third
           orbit, and the graph is now 3-regular. This binary tree
           contributes 4 new orbits to the Harries-Wong graph.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.HarriesWongGraph()
        sage: g.order()
        70
        sage: g.size()
        105
        sage: g.girth()
        10
        sage: g.diameter()
        6
        sage: orbits = g.automorphism_group(orbits=True)[-1]
        sage: g.show(figsize=[15, 15], partition=orbits)   # long time

    Alternative embedding::

        sage: graphs.HarriesWongGraph(embedding=2).show()

    TESTS::

        sage: graphs.HarriesWongGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    L = [9, 25, 31, -17, 17, 33, 9, -29, -15, -9, 9, 25, -25, 29, 17, -9,
         9, -27, 35, -9, 9, -17, 21, 27, -29, -9, -25, 13, 19, -9, -33,
         -17, 19, -31, 27, 11, -25, 29, -33, 13, -13, 21, -29, -21, 25,
         9, -11, -19, 29, 9, -27, -19, -13, -35, -9, 9, 17, 25, -9, 9, 27,
         -27, -21, 15, -9, 29, -29, 33, -9, -25]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(70, L, 1)
    g.name("Harries-Wong graph")

    if embedding == 1:
        d = g.get_pos()

        # Binary tree (left side)
        d[66] = (-9.5, 0)
        _line_embedding(g, [37, 65, 67], first=(-8, 2.25),
                last=(-8, -2.25))
        _line_embedding(g, [36, 38, 64, 24, 68, 30], first=(-7, 3),
                last=(-7, -3))
        _line_embedding(g, [35, 39, 63, 25, 59, 29, 11, 5, 55, 23, 69, 31],
                first=(-6, 3.5), last=(-6, -3.5))

        # Cube, corners: [9, 15, 21, 27, 45, 51, 57, 61]
        _circle_embedding(g, [61, 9], center=(0, -1.5), shift=.2,
                radius=4)
        _circle_embedding(g, [27, 15], center=(0, -1.5), shift=.7,
                radius=4*.707)
        _circle_embedding(g, [51, 21], center=(0, 2.5), shift=.2,
                radius=4)
        _circle_embedding(g, [45, 57], center=(0, 2.5), shift=.7,
                radius=4*.707)

        # Cube, subdivision
        _line_embedding(g, [21, 22, 43, 44, 45], first=d[21], last=d[45])
        _line_embedding(g, [21, 4, 3, 56, 57], first=d[21], last=d[57])
        _line_embedding(g, [57, 12, 13, 14, 15], first=d[57], last=d[15])
        _line_embedding(g, [15, 6, 7, 8, 9], first=d[15], last=d[9])
        _line_embedding(g, [9, 10, 19, 20, 21], first=d[9], last=d[21])
        _line_embedding(g, [45, 54, 53, 52, 51], first=d[45], last=d[51])
        _line_embedding(g, [51, 50, 49, 58, 57], first=d[51], last=d[57])
        _line_embedding(g, [51, 32, 33, 34, 61], first=d[51], last=d[61])
        _line_embedding(g, [61, 62, 41, 40, 27], first=d[61], last=d[27])
        _line_embedding(g, [9, 0, 1, 26, 27], first=d[9], last=d[27])
        _line_embedding(g, [27, 28, 47, 46, 45], first=d[27], last=d[45])
        _line_embedding(g, [15, 16, 17, 60, 61], first=d[15], last=d[61])

        # Top vertices
        _line_embedding(g, [2, 18, 42, 48], first=(-1, 7), last=(3, 7))

        return g

    elif embedding == 2:
        return g
    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def WellsGraph():
    r"""
    Returns the Wells graph.

    For more information on the Wells graph (also called Armanios-Wells graph),
    see `this page <http://www.win.tue.nl/~aeb/graphs/Wells.html>`_.

    The implementation follows the construction given on page 266 of
    [BCN89]_. This requires to create intermediate graphs and run a small
    isomorphism test, while everything could be replaced by a pre-computed list
    of edges : I believe that it is better to keep "the recipe" in the code,
    however, as it is quite unlikely that this could become the most
    time-consuming operation in any sensible algorithm, and .... "preserves
    knowledge", which is what open-source software is meant to do.

    EXAMPLES::

        sage: g = graphs.WellsGraph(); g
        Wells graph: Graph on 32 vertices
        sage: g.order()
        32
        sage: g.size()
        80
        sage: g.girth()
        5
        sage: g.diameter()
        4
        sage: g.chromatic_number()
        4
        sage: g.is_regular(k=5)
        True

    REFERENCES:

    .. [BCN89] A. E. Brouwer, A. M. Cohen, A. Neumaier,
      Distance-Regular Graphs,
      Springer, 1989.
    """
    from platonic_solids import DodecahedralGraph
    from basic import CompleteBipartiteGraph

    # Following the construction from the book "Distance-regular graphs"
    dodecahedron = DodecahedralGraph()

    # Vertices at distance 3 in the Dodecahedron
    distance3 = dodecahedron.distance_graph([3])

    # Building the graph whose line graph is the dodecahedron.
    b = CompleteBipartiteGraph(5,5)
    b.delete_edges([(0,5), (1,6), (2,7), (3,8), (4,9)])

    # Computing the isomorphism between the two
    b = b.line_graph(labels = False)
    _, labels = distance3.is_isomorphic(b, certify = True)

    # The relabeling that the books claims to exist.
    for v,new_name in labels.items():
        x,y = new_name
        labels[v] = (x%5,y%5)

    dodecahedron.relabel(labels)

    # Checking that the above computations indeed produces a good labeling.
    for u in dodecahedron:
        for v in dodecahedron:
            if u == v:
                continue

            if (u[0] != v[0]) and (u[1] != v[1]):
                continue

            if dodecahedron.distance(u,v) != 3:
                raise ValueError("There is something wrong going on !")

    # The graph we will return, starting from the dodecahedron
    g = dodecahedron

    # Good ! Now adding 12 new vertices
    for i in range(5):
        g.add_edge((i,'+'),('inf','+'))
        g.add_edge((i,'-'),('inf','-'))
        for k in range(5):
            if k == i:
                continue
            g.add_edge((i,'+'),(i,k))
            g.add_edge((i,'-'),(k,i))

    g.name("Wells graph")

    # Giving our graph a "not-so-bad" layout
    g.relabel({
            (1, 3): 8, (3, 0): 18, (3, '+'): 22, (2, 1): 13,
            (1, '+'): 10, (0, 3): 2, (2, '+'): 16, ('inf', '-'): 31,
            (4, 0): 24, (1, 2): 7, (4, '+'): 28, (0, '-'): 5,
            (0, 4): 3, (4, 1): 25, (2, '-'): 17, (3, 2): 20,
            (3, '-'): 23, (1, '-'): 11, (1, 4): 9, (2, 3): 14,
            ('inf', '+'): 30, (4, 2): 26, (1, 0): 6, (0, 1): 0,
            (3, 1): 19, (0, 2): 1, (2, 0): 12, (4, '-'): 29,
            (0, '+'): 4, (4, 3): 27, (3, 4): 21, (2, 4): 15})

    p = [(1, 29, 20, 13, 12, 28, 14, 7),
         (2, 5, 30, 23, 18, 4, 31, 22),
         (3, 17, 21, 9, 24, 16, 27, 25),
         (6, 10, 8, 15, 0, 11, 19, 26)]

    from sage.graphs.graph_plot import _circle_embedding
    _circle_embedding(g, p[0], radius = 1)
    _circle_embedding(g, p[1], radius = .9)
    _circle_embedding(g, p[2], radius = .8)
    _circle_embedding(g, p[3], radius = .7)

    return g


def HallJankoGraph(from_string=True):
    r"""
    Returns the Hall-Janko graph.

    For more information on the Hall-Janko graph, see its
    :wikipedia:`Wikipedia page <Hall-Janko_graph>`.

    The construction used to generate this graph in Sage is by
    a 100-point permutation representation of the Janko group `J_2`,
    as described in version 3 of the ATLAS of Finite Group
    representations, in particular on the page `ATLAS: J2
    -- Permutation representation on 100 points
    <http://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/J2G1-p100B0>`_.

    INPUT:

    - ``from_string`` (boolean) -- whether to build the graph from
      its sparse6 string or through GAP. The two methods return the
      same graph though doing it through GAP takes more time. It is
      set to ``True`` by default.

    EXAMPLES::

        sage: g = graphs.HallJankoGraph()
        sage: g.is_regular(36)
        True
        sage: g.is_vertex_transitive()
        True

    Is it really strongly regular with parameters 14, 12? ::

        sage: nu = set(g.neighbors(0))
        sage: for v in range(1, 100):
        ....:     if v in nu:
        ....:         expected = 14
        ....:     else:
        ....:         expected = 12
        ....:     nv = set(g.neighbors(v))
        ....:     nv.discard(0)
        ....:     if len(nu & nv) != expected:
        ....:         print "Something is wrong here!!!"
        ....:         break

    Some other properties that we know how to check::

        sage: g.diameter()
        2
        sage: g.girth()
        3
        sage: factor(g.characteristic_polynomial())
        (x - 36) * (x - 6)^36 * (x + 4)^63

    TESTS::

        sage: gg = graphs.HallJankoGraph(from_string=False) # long time
        sage: g == gg # long time
        True
    """

    string = (":~?@c__E@?g?A?w?A@GCA_?CA`OWF`W?EAW?@?_OD@_[GAgcIaGGB@OcIA"
              "wCE@o_K_?GB@?WGAouC@OsN_?GB@O[GB`A@@_e?@OgLB_{Q_?GC@O[GAOs"
              "OCWGBA?kKBPA@?_[KB_{OCPKT`o_RD`]A?o[HBOwODW?DA?cIB?wRDP[X`"
              "ogKB_{QD@]B@o_KBPWXE`mC@o_JB?{PDPq@?oWGA_{OCPKTDp_YEwCA@_c"
              "IBOwOC`OX_OGB@?WPDPcYFg?C@_gKBp?SE@cYF`{_`?SGAOoOC`_\\FwCE"
              "A?gKBO{QD@k[FqI??_OFA_oQE@k\\Fq?`GgCB@pGRD@_XFP{a_?SE@ocIA"
              "ooNCPOUEqU@?oODA?cJB_{UEqYC@_kLC@CREPk]GAGbHgCA@?SMBpCSD`["
              "YFq?`Ga]BA?gPC`KSD`_\\Fa?cHWGB@?[IAooPD`[WF@s^HASeIg?@@OcP"
              "C`KYF@w^GQ[h`O[HAooMC@CQCpSVEPk\\GaSeIG?FA?kLB_{OC`OVE@cYG"
              "QUA@?WLBp?PC`KVEqKgJg?DA?sMBpCSDP[WEQKfIay@?_KD@_[GC`SUE@k"
              "[FaKdHa[k_?OLC@CRD@WVEpo^HAWfIAciIqoo_?CB@?kMCpOUE`o\\GAKg"
              "IQgq_?GD@_[GB?{OCpWVE@cYFACaHAWhJR?q_?CC@_kKBpC\\GACdHa[kJ"
              "a{o_?CA?oOFBpGRD@o\\GaKdIQonKrOt_?WHA`?PC`KTD`k]FqSeIaolJr"
              "CqLWCA@OkKCPGRDpcYGAKdIAgjJAsmJr?t__OE@ogJB_{XEps`HA[gIQwn"
              "KWKGAOoMBpGUE`k[Fa?aHqckJbSuLw?@?_SHA_kLC@OTFPw^GaOkLg?B@?"
              "[HA_{PDP_XFaCbHa[gIqooKRWx_?CFBpOTE@cZFPw^GACcHQgoKrSvMwWG"
              "BOwQCp_YFP{`HASfJAwnKRSx_OSSDP[WEq?aGqSfIQsoKR_zNWCE@o_HA_"
              "sREPg^GAGcHQWfIAciKbOxNg?A@__IAooMC`KTD`g\\GAKcIasoKrOtLb["
              "wMbyCA?cKBp?TD`[WE`s^GQGbHqcjJrK{NRw~_oODA?sNC@CQCpOZF@s]G"
              "QOfIaolJrGsLbk}_?OFA_sRD@SVE`k[HQcjJa{qLb[xMb|?_OOFA?cIAos"
              "RDP_ZFa?aGqOfIAsuMbk{Ns@@OsQAA_sPDPWXE`o\\FqKdIQkkJrCuLr_x"
              "Mro}NsDAPG?@@OWFApKUE@o`IQolKRKsLrc|NsQC@OWGAOgJCpOWE`o_GQ"
              "KiIqwnKr_~OcLCPS]A?oWHA_oMBpKSDP[\\FagjKBWxMbk{OSQ@@O_IAoo"
              "LBpCSD`g\\FaGbHQWgIQgmKRKwMRl?PgGC@OWHB@KSE@c[FqCaGqSeIAkk"
              "KBCqLBSuMBpGQWCA@?cKBOwRDPWVE@k^GqOfJr?pKbKtLrs}OSHDQwKIBO"
              "wPD@WWEQ?`HQWfIQglKBOtLbo}Ns@@OsTE_?kLCpWWHA[gIqomKBGwMRgz"
              "NBw~OSPDPc\\H_?CFAOoLCPSVE`o\\GAOeJAwpKbKtMrx?Qcq??OKFA?gJ"
              "B`?QDpcYEpo]FqKfIAgjJB?qKr_{NS@A__SE@o_HBO{PC`OTD`{_HaciIq"
              "{vMbt?OcPFQCeB@?SKBOwRD@SXE`k[FPw`HQ_lKRKxNRxBPC\\HQclK_?K"
              "EB?sOC`OTDa?`GqWgJRCrNBw~OSHFQStMRtDQ_?KC@OoQE`k_GaOdHa[gI"
              "q{tMBg|Nb|?OcPMSDDQSwCB@_cJB_{OCpOVFP{dHa[jJQwqKrk}NsHBQCd"
              "MRtMA?oSEA_wPDp_YEpo]GAOeIq{pLBk}NsLEQCtNTDU??OKEA_oLC@[[G"
              "aKnKBOtLbk~OCPFQStNSDLSTgGKC@GSD`[WEpw_GQGcIAciJAwpKb_xMbk"
              "~QShJRc|R`_wNCPcZF@s^GAGbHA_hJR?qKrOvMRg|NsDEPsxTTgCB@?gJB"
              "?sMC@CUDp_]FqCaHQcjJQwtLrhCPS\\IRCtQTw?B@?SHA_wPC`_aGqOiJa"
              "{oKRKvMRpFQChKRtXVUTi??ocNC@KUE@cYFaGdHa_mJrKsLb[yMro|OcXI"
              "RdPTTddZaOgJB@?UEPk[FQCfIaolJrSvMBczNR|AOsXFQCtOTtaB@?WGAP"
              "?TEPo\\GAGdHqgmKBCqLR[xMb|?PC`HQs|TTt`XUtu@?o[HB?sNCPGXF@{"
              "_GQKcIqolJb_yNCLDPs`MRtDRTTdYUwSEA?kLB`CWF@s]FqGgIqooLRgzN"
              "RxFQSlMSDDQTDXVUTi@?_KDAOoLBpKUEQOfIa{oLB_xMrt?Os\\HQcpMST"
              "HSTtl[VT}A@ocJBOwSD`_XEpo_Ha_mJrKtLbgzNSTGQspLRtDUUDp\\WG["
              "HB`CQCp[WFQGgIQgkJQ{rLbc{Nc@APsdLRt@PSt\\WUtt_Wn")

    if from_string:
        g = Graph(string, loops = False, multiedges = False)
    else:

        # The following construction is due to version 3 of the ATLAS of
        # Finite Group Representations, specifically the page at
        # http://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/J2G1-p100B0 .

        from sage.interfaces.gap import gap
        gap.eval("g1 := (1,84)(2,20)(3,48)(4,56)(5,82)(6,67)(7,55)(8,41)"
                 "(9,35)(10,40)(11,78)(12,100)(13,49)(14,37)(15,94)(16,76)"
                 "(17,19)(18,44)(21,34)(22,85)(23,92)(24,57)(25,75)(26,28)"
                 "(27,64)(29,90)(30,97)(31,38)(32,68)(33,69)(36,53)(39,61)"
                 "(42,73)(43,91)(45,86)(46,81)(47,89)(50,93)(51,96)(52,72)"
                 "(54,74)(58,99)(59,95)(60,63)(62,83)(65,70)(66,88)(71,87)"
                 "(77,98)(79,80);")

        gap.eval("g2 := (1,80,22)(2,9,11)(3,53,87)(4,23,78)(5,51,18)"
                 "(6,37,24)(8,27,60)(10,62,47)(12,65,31)(13,64,19)"
                 "(14,61,52)(15,98,25)(16,73,32)(17,39,33)(20,97,58)"
                 "(21,96,67)(26,93,99)(28,57,35)(29,71,55)(30,69,45)"
                 "(34,86,82)(38,59,94)(40,43,91)(42,68,44)(46,85,89)"
                 "(48,76,90)(49,92,77)(50,66,88)(54,95,56)(63,74,72)"
                 "(70,81,75)(79,100,83);")

        gap.eval("G := Group([g1,g2]);")
        edges = gap('Orbit(G,[1,5],OnSets)').sage()
        g = Graph([(int(u), int(v)) for u,v in edges])
        g.relabel()

    _circle_embedding(g, range(100))
    g.name("Hall-Janko graph")
    return g

def Balaban10Cage(embedding=1):
    r"""
    Returns the Balaban 10-cage.

    The Balaban 10-cage is a 3-regular graph with 70 vertices and
    105 edges. See its :wikipedia:`Wikipedia page
    <Balaban_10-cage>`.

    The default embedding gives a deeper understanding of the
    graph's automorphism group. It is divided into 4 layers (each
    layer being a set of points at equal distance from the drawing's
    center). From outside to inside:

    - L1: The outer layer (vertices which are the furthest from the
      origin) is actually the disjoint union of two cycles of length
      10.

    - L2: The second layer is an independent set of 20 vertices.

    - L3: The third layer is a matching on 10 vertices.

    - L4: The inner layer (vertices which are the closest from the
      origin) is also the disjoint union of two cycles of length 10.

    This graph is not vertex-transitive, and its vertices are
    partitioned into 3 orbits: L2, L3, and the union of L1 of L4
    whose elements are equivalent.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to be either 1 or 2.

    EXAMPLES::

        sage: g = graphs.Balaban10Cage()
        sage: g.girth()
        10
        sage: g.chromatic_number()
        2
        sage: g.diameter()
        6
        sage: g.is_hamiltonian()
        True
        sage: g.show(figsize=[10,10])   # long time

    TESTS::

        sage: graphs.Balaban10Cage(embedding='foo')
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    L = [-9, -25, -19, 29, 13, 35, -13, -29, 19, 25, 9, -29, 29, 17, 33,
          21, 9,-13, -31, -9, 25, 17, 9, -31, 27, -9, 17, -19, -29, 27,
          -17, -9, -29, 33, -25,25, -21, 17, -17, 29, 35, -29, 17, -17,
          21, -25, 25, -33, 29, 9, 17, -27, 29, 19, -17, 9, -27, 31, -9,
          -17, -25, 9, 31, 13, -9, -21, -33, -17, -29, 29]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(70, L, 1)
    g.name("Balaban 10-cage")

    if embedding == 2:
        return g
    elif embedding != 1:
        raise ValueError("The value of embedding must be 1 or 2.")

    L3 = [5, 24, 35, 46, 29, 40, 51, 34, 45, 56]
    _circle_embedding(g, L3, center=(0,0), radius = 4.3)

    L2  = [6, 4, 23, 25, 60, 36, 1, 47, 28, 30, 39, 41, 50, 52, 33, 9, 44,
            20, 55, 57]
    _circle_embedding(g, L2, center=(0,0), radius = 5, shift=-.5)


    L1a = [69, 68, 67, 66, 65, 64, 63, 62, 61, 0]
    L1b = [19, 18, 17, 16, 15, 14, 13, 12, 11, 10]
    _circle_embedding(g, L1a, center=(0,0), radius = 6, shift = 3.25)
    _circle_embedding(g, L1b, center=(0,0), radius = 6, shift = -1.25)

    L4a = [37, 2, 31, 38, 53, 32, 21, 54, 3, 22]
    _circle_embedding(g, L4a, center=(0,0), radius = 3, shift = 1.9)

    L4b = [26, 59, 48, 27, 42, 49, 8, 43, 58, 7]
    _circle_embedding(g, L4b, center=(0,0), radius = 3, shift = 1.1)

    return g

def Balaban11Cage(embedding = 1):
    r"""
    Returns the Balaban 11-cage.

    For more information, see this :wikipedia:`Wikipedia article on
    the Balaban 11-cage <Balaban_11-cage>`.

    INPUT:

    - ``embedding`` -- three embeddings are available, and can be
      selected by setting ``embedding`` to be 1, 2, or 3.

      - The first embedding is the one appearing on page 9 of the
        Fifth Annual Graph Drawing Contest report [FAGDC]_. It
        separates vertices based on their eccentricity (see
        :meth:`eccentricity()
        <sage.graphs.generic_graph.GenericGraph.eccentricity>`).

      - The second embedding has been produced just for Sage and is
        meant to emphasize the automorphism group's 6 orbits.

      - The last embedding is the default one produced by the
        :meth:`LCFGraph` constructor.

    .. NOTE::

        The vertex labeling changes according to the value of
        ``embedding=1``.

    EXAMPLES:

    Basic properties::

        sage: g = graphs.Balaban11Cage()
        sage: g.order()
        112
        sage: g.size()
        168
        sage: g.girth()
        11
        sage: g.diameter()
        8
        sage: g.automorphism_group().cardinality()
        64

    Our many embeddings::

        sage: g1 = graphs.Balaban11Cage(embedding=1)
        sage: g2 = graphs.Balaban11Cage(embedding=2)
        sage: g3 = graphs.Balaban11Cage(embedding=3)
        sage: g1.show(figsize=[10,10])   # long time
        sage: g2.show(figsize=[10,10])   # long time
        sage: g3.show(figsize=[10,10])   # long time

    Proof that the embeddings are the same graph::

        sage: g1.is_isomorphic(g2) # g2 and g3 are obviously isomorphic
        True

    TESTS::

        sage: graphs.Balaban11Cage(embedding='xyzzy')
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1, 2, or 3.

    REFERENCES:

    .. [FAGDC] Fifth Annual Graph Drawing Contest
       P. Eaded, J. Marks, P.Mutzel, S. North
       http://www.merl.com/papers/docs/TR98-16.pdf
    """
    if embedding == 1:
        pos_dict = {}
        for j in range(8):
            for i in range(8):
                pos_dict[str(j) + str(i)]= [
                        0.8 * float(cos(2*((8*j + i)*pi/64 + pi/128))),
                        0.8 * float(sin(2*((8*j + i)*pi/64 + pi/128)))
                ]
            for i in range(4):
                pos_dict['1' + str(j) + str(i)] = [
                        1.1 * float(cos(2*((4*j + i)*pi/32 + pi/64))),
                        1.1 * float(sin(2*((4*j + i)*pi/32 + pi/64)))
                ]
            for i in range(2):
                pos_dict['1' + str(j) + str(i + 4)] = [
                        1.4 * float(cos(2*((2*j + i)*pi/16 + pi/32))),
                        1.4 * float(sin(2*((2*j + i)*pi/16 + pi/32)))
                ]

        edge_dict = {
            "00": ["11"], "01": ["10"],   "02": ["53"], "03": ["52"],
            "11": ["20"], "10": ["21"],   "53": ["22"], "52": ["23"],
            "20": ["31"], "21": ["30"],   "22": ["33"], "23": ["32"],
            "31": ["40"], "30": ["41"],   "33": ["43"], "32": ["42"],
            "40": ["50"], "41": ["51"],   "43": ["12"], "42": ["13"],
            "50": ["61"], "51": ["60"],   "12": ["63"], "13": ["62"],
            "61": ["70"], "60": ["71"],   "63": ["72"], "62": ["73"],
            "70": ["01"], "71": ["00"],   "72": ["03"], "73": ["02"],

            "04": ["35"], "05": ["34"],   "06": ["37"], "07": ["36"],
            "35": ["64"], "34": ["65"],   "37": ["66"], "36": ["67"],
            "64": ["55"], "65": ["54"],   "66": ["17"], "67": ["16"],
            "55": ["45"], "54": ["44"],   "17": ["46"], "16": ["47"],
            "45": ["74"], "44": ["75"],   "46": ["76"], "47": ["77"],
            "74": ["25"], "75": ["24"],   "76": ["27"], "77": ["26"],
            "25": ["14"], "24": ["15"],   "27": ["56"], "26": ["57"],
            "14": ["05"], "15": ["04"],   "56": ["07"], "57": ["06"],

            "100": ["03", "04"],   "110": ["10", "12"],
            "101": ["01", "06"],   "111": ["11", "13"],
            "102": ["00", "07"],   "112": ["14", "16"],
            "103": ["02", "05"],   "113": ["15", "17"],

            "120": ["22", "24"],   "130": ["33", "36"],
            "121": ["20", "26"],   "131": ["32", "37"],
            "122": ["21", "27"],   "132": ["31", "34"],
            "123": ["23", "25"],   "133": ["30", "35"],

            "140": ["43", "45"],   "150": ["50", "52"],
            "141": ["40", "46"],   "151": ["51", "53"],
            "142": ["41", "47"],   "152": ["54", "56"],
            "143": ["42", "44"],   "153": ["55", "57"],

            "160": ["60", "66"],   "170": ["73", "76"],
            "161": ["63", "65"],   "171": ["72", "77"],
            "162": ["62", "64"],   "172": ["71", "74"],
            "163": ["61", "67"],   "173": ["70", "75"],

            "104": ["100", "102", "105"],   "114": ["110", "111", "115"],
            "105": ["101", "103", "104"],   "115": ["112", "113", "114"],

            "124": ["120", "121", "125"],   "134": ["130", "131", "135"],
            "125": ["122", "123", "124"],   "135": ["132", "133", "134"],

            "144": ["140", "141", "145"],   "154": ["150", "151", "155"],
            "145": ["142", "143", "144"],   "155": ["152", "153", "154"],

            "164": ["160", "161", "165"],   "174": ["170", "171", "175"],
            "165": ["162", "163", "164"],   "175": ["172", "173", "174"]
        }

        return Graph(edge_dict, pos=pos_dict, name="Balaban 11-cage")

    elif embedding == 2 or embedding == 3:
        L = [44, 26, -47, -15, 35, -39, 11, -27, 38, -37, 43, 14, 28, 51,
             -29, -16, 41, -11, -26, 15, 22, -51, -35, 36, 52, -14, -33,
             -26, -46, 52, 26, 16, 43, 33, -15, 17, -53, 23, -42, -35, -28,
             30, -22, 45, -44, 16, -38, -16, 50, -55, 20, 28, -17, -43,
             47, 34, -26, -41, 11, -36, -23, -16, 41, 17, -51, 26, -33,
             47, 17, -11, -20, -30, 21, 29, 36, -43, -52, 10, 39, -28, -17,
             -52, 51, 26, 37, -17, 10, -10, -45, -34, 17, -26, 27, -21,
             46, 53, -10, 29, -50, 35, 15, -47, -29, -41, 26, 33, 55, -17,
             42, -26, -36, 16]

        from sage.graphs.generators.families import LCFGraph
        g = LCFGraph(112, L, 1)
        g.name("Balaban 11-cage")

        if embedding == 3:
            return g

        v1 = [34, 2, 54, 43, 66, 20, 89, 100, 72, 76, 6, 58, 16, 78, 74,
              70, 36, 94, 27, 25, 10, 8, 45, 60, 14, 64, 80, 82, 109, 107,
              49, 98]
        v2 = [88, 3, 19, 55, 67, 42, 101, 33, 77, 5, 17, 57, 69, 71, 73,
              75, 11, 61, 28, 9, 37, 26, 46, 95, 13, 63, 81, 83, 108, 106,
              48, 97]
        l1 = [35, 93, 1, 24, 53, 7, 44, 59, 15, 65, 79, 21, 110, 90, 50,
              99]
        l2 = [87, 4, 18, 56, 68, 41, 102, 32, 12, 62, 29, 84, 38, 105, 47,
              96]

        d = g.get_pos()
        for i,v in enumerate(v1):
            d[v] = (-2, 16.5-i)

        for i,v in enumerate(l1):
            d[v] = (-10, 8-i)

        for i,v in enumerate(l2):
            d[v] = (10, 8.5-i)

        for i,v in enumerate(v2):
            d[v] = (2, 16.5-i)

        for i,v in enumerate([0, 111, 92, 91, 52, 51, 23, 22]):
            d[v] = (-20, 14.5-4*i)

        for i,v in enumerate([104, 103, 86, 85, 40, 39, 31, 30]):
            d[v] = (20, 14.5-4*i)

        return g

    else:
        raise ValueError("The value of embedding must be 1, 2, or 3.")

def BidiakisCube():
    r"""
    Returns the Bidiakis cube.

    For more information, see this
    `Wikipedia article on the Bidiakis cube <http://en.wikipedia.org/wiki/Bidiakis_cube>`_.

    EXAMPLES:

    The Bidiakis cube is a 3-regular graph having 12 vertices and 18
    edges. This means that each vertex has a degree of 3. ::

        sage: g = graphs.BidiakisCube(); g
        Bidiakis cube: Graph on 12 vertices
        sage: g.show()  # long time
        sage: g.order()
        12
        sage: g.size()
        18
        sage: g.is_regular(3)
        True

    It is a Hamiltonian graph with diameter 3 and girth 4::

        sage: g.is_hamiltonian()
        True
        sage: g.diameter()
        3
        sage: g.girth()
        4

    It is a planar graph with characteristic polynomial
    `(x - 3) (x - 2) (x^4) (x + 1) (x + 2) (x^2 + x - 4)^2` and
    chromatic number 3::

        sage: g.is_planar()
        True
        sage: bool(g.characteristic_polynomial() == expand((x - 3) * (x - 2) * (x^4) * (x + 1) * (x + 2) * (x^2 + x - 4)^2))
        True
        sage: g.chromatic_number()
        3
    """
    edge_dict = {
        0:[1,6,11], 1:[2,5], 2:[3,10], 3:[4,9], 4:[5,8],
        5:[6], 6:[7], 7:[8,11], 8:[9], 9:[10], 10:[11]}
    pos_dict = {
        0: [0, 1],
        1: [0.5, 0.866025403784439],
        2: [0.866025403784439, 0.500000000000000],
        3: [1, 0],
        4: [0.866025403784439, -0.5],
        5: [0.5, -0.866025403784439],
        6: [0, -1],
        7: [-0.5, -0.866025403784439],
        8: [-0.866025403784439, -0.5],
        9: [-1, 0],
        10: [-0.866025403784439, 0.5],
        11: [-0.5, 0.866025403784439]}
    return Graph(edge_dict, pos=pos_dict, name="Bidiakis cube")

def BiggsSmithGraph(embedding=1):
    r"""
    Returns the Biggs-Smith graph.

    For more information, see this :wikipedia:`Wikipedia article on
    the Biggs-Smith graph <Biggs-Smith_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to be 1 or 2.

    EXAMPLES:

    Basic properties::

        sage: g = graphs.BiggsSmithGraph()
        sage: g.order()
        102
        sage: g.size()
        153
        sage: g.girth()
        9
        sage: g.diameter()
        7
        sage: g.automorphism_group().cardinality()
        2448
        sage: g.show(figsize=[10, 10])   # long time

    The other embedding::

        sage: graphs.BiggsSmithGraph(embedding=2).show()

    TESTS::

        sage: graphs.BiggsSmithGraph(embedding='xyzzy')
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.

    """
    L = [16, 24, -38, 17, 34, 48, -19, 41, -35, 47, -20, 34, -36,
         21, 14, 48, -16, -36, -43, 28, -17, 21, 29, -43, 46, -24,
         28, -38, -14, -50, -45, 21, 8, 27, -21, 20, -37, 39, -34,
         -44, -8, 38, -21, 25, 15, -34, 18, -28, -41, 36, 8, -29,
         -21, -48, -28, -20, -47, 14, -8, -15, -27, 38, 24, -48, -18,
         25, 38, 31, -25, 24, -46, -14, 28, 11, 21, 35, -39, 43, 36,
         -38, 14, 50, 43, 36, -11, -36, -24, 45, 8, 19, -25, 38, 20,
         -24, -14, -21, -8, 44, -31, -38, -28, 37]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(102, L, 1)
    g.name("Biggs-Smith graph")

    if embedding == 1:

        orbs = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 0],
                [17, 101, 25, 66, 20, 38, 53, 89, 48, 75, 56, 92, 45, 78,
                 34, 28, 63],
                [18, 36, 26, 65, 19, 37, 54, 90, 47, 76, 55, 91, 46, 77,
                 35, 27, 64],
                [21, 39, 52, 88, 49, 74, 57, 93, 44, 79, 33, 29, 62, 83,
                 100, 24, 67],
                [22, 97, 51, 96, 50, 95, 58, 94, 59, 80, 60, 81, 61, 82,
                 99, 23, 98],
                [30, 86, 84, 72, 70, 68, 42, 40, 31, 87, 85, 73, 71, 69,
                 43, 41, 32]]

        # central orbits
        _circle_embedding(g, orbs[1], center=(-.4, 0), radius=.2)
        _circle_embedding(g, orbs[3], center=(.4, 0), radius=.2, shift=4)

        # lower orbits
        _circle_embedding(g, orbs[0], center=(-.9, -.5), radius=.3,
                shift=2)
        _circle_embedding(g, orbs[2], center=(-.9, .5), radius=.3)

        # upper orbits
        _circle_embedding(g, orbs[4], center=(.9, -.5), radius=.3, shift=4)
        _circle_embedding(g, orbs[5], center=(.9, .5), radius=.3, shift=-2)

    elif embedding == 2:
        pass
    else:
        raise ValueError("The value of embedding must be 1 or 2.")

    return g

def BlanusaFirstSnarkGraph():
    r"""
    Returns the first Blanusa Snark Graph.

    The Blanusa graphs are two snarks on 18 vertices and 27 edges. For more
    information on them, see the :wikipedia:`Blanusa_snarks`.

    .. SEEALSO::

        * :meth:`~sage.graphs.graph_generators.GraphGenerators.BlanusaSecondSnarkGraph`.

    EXAMPLES::

        sage: g = graphs.BlanusaFirstSnarkGraph()
        sage: g.order()
        18
        sage: g.size()
        27
        sage: g.diameter()
        4
        sage: g.girth()
        5
        sage: g.automorphism_group().cardinality()
        8
    """
    g = Graph({17:[4,7,1],0:[5],
               3:[8],13:[9],12:[16],
               10:[15],11:[6],14:[2]},
              name="Blanusa First Snark Graph")

    g.add_cycle(range(17))
    _circle_embedding(g, range(17), shift=0.25)
    g.get_pos()[17] = (0,0)
    return g

def BlanusaSecondSnarkGraph():
    r"""
    Returns the second Blanusa Snark Graph.

    The Blanusa graphs are two snarks on 18 vertices and 27 edges. For more
    information on them, see the :wikipedia:`Blanusa_snarks`.

    .. SEEALSO::

        * :meth:`~sage.graphs.graph_generators.GraphGenerators.BlanusaFirstSnarkGraph`.

    EXAMPLES::

        sage: g = graphs.BlanusaSecondSnarkGraph()
        sage: g.order()
        18
        sage: g.size()
        27
        sage: g.diameter()
        4
        sage: g.girth()
        5
        sage: g.automorphism_group().cardinality()
        4
    """
    g = Graph({0:[(0,0),(1,4),1],1:[(0,3),(1,1)],(0,2):[(0,5)],
               (0,6):[(0,4)],(0,7):[(0,1)],(1,7):[(1,2)],
               (1,0):[(1,6)],(1,3):[(1,5)]},
              name="Blanusa Second Snark Graph")

    g.add_cycle([(0,i) for i in range(5)])
    g.add_cycle([(1,i) for i in range(5)])
    g.add_cycle([(0,5),(0,6),(0,7),(1,5),(1,6),(1,7)])

    _circle_embedding(g,
                      [(0,(2*i)%5) for i in range(5)],
                      center = (-1.5,0),
                      shift = .5)
    _circle_embedding(g,
                      [(1,(2*i)%5) for i in range(5)],
                      center = (1.5,0))

    _circle_embedding(g,
                      [(0,i) for i in range(5,8)]+[0]*4,
                      center = (-1.2,0),
                      shift = 2.5,
                      radius = 2.2)
    _circle_embedding(g,
                      [(1,i) for i in range(5,8)]+[0]*4,
                      center = (1.2,0),
                      shift = -1,
                      radius = 2.2)

    _circle_embedding(g,[0,1], shift=.5)
    g.relabel()
    return g

def BrinkmannGraph():
    r"""
    Returns the Brinkmann graph.

    For more information, see the
    `Wikipedia article on the Brinkmann graph <http://en.wikipedia.org/wiki/Brinkmann_graph>`_.

    EXAMPLES:

    The Brinkmann graph is a 4-regular graph having 21 vertices and 42
    edges. This means that each vertex has degree 4. ::

        sage: G = graphs.BrinkmannGraph(); G
        Brinkmann graph: Graph on 21 vertices
        sage: G.show()  # long time
        sage: G.order()
        21
        sage: G.size()
        42
        sage: G.is_regular(4)
        True

    It is an Eulerian graph with radius 3, diameter 3, and girth 5. ::

        sage: G.is_eulerian()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        3
        sage: G.girth()
        5

    The Brinkmann graph is also Hamiltonian with chromatic number 4::

        sage: G.is_hamiltonian()
        True
        sage: G.chromatic_number()
        4

    Its automorphism group is isomorphic to `D_7`::

        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(7))
        True
    """
    edge_dict = {
        0: [2,5,7,13],
        1: [3,6,7,8],
        2: [4,8,9],
        3: [5,9,10],
        4: [6,10,11],
        5: [11,12],
        6: [12,13],
        7: [15,20],
        8: [14,16],
        9: [15,17],
        10: [16,18],
        11: [17,19],
        12: [18,20],
        13: [14,19],
        14: [17,18],
        15: [18,19],
        16: [19,20],
        17: [20]}
    pos_dict = {
        0: [0, 4],
        1: [3.12732592987212, 2.49395920743493],
        2: [3.89971164872729, -0.890083735825258],
        3: [1.73553495647023, -3.60387547160968],
        4: [-1.73553495647023, -3.60387547160968],
        5: [-3.89971164872729, -0.890083735825258],
        6: [-3.12732592987212, 2.49395920743493],
        7: [0.867767478235116, 1.80193773580484],
        8: [1.94985582436365, 0.445041867912629],
        9: [1.56366296493606, -1.24697960371747],
        10: [0, -2],
        11: [-1.56366296493606, -1.24697960371747],
        12: [-1.94985582436365, 0.445041867912629],
        13: [-0.867767478235116, 1.80193773580484],
        14: [0.433883739117558, 0.900968867902419],
        15: [0.974927912181824, 0.222520933956314],
        16: [0.781831482468030, -0.623489801858733],
        17: [0, -1],
        18: [-0.781831482468030, -0.623489801858733],
        19: [-0.974927912181824, 0.222520933956315],
        20: [-0.433883739117558, 0.900968867902419]}
    return Graph(edge_dict, pos=pos_dict, name="Brinkmann graph")

def BrouwerHaemersGraph():
    r"""
    Returns the Brouwer-Haemers Graph.

    The Brouwer-Haemers is the only strongly regular graph of parameters
    `(81,20,1,6)`. It is build in Sage as the Affine Orthogonal graph
    `VO^-(6,3)`. For more information on this graph, see its `corresponding page
    on Andries Brouwer's website
    <http://www.win.tue.nl/~aeb/graphs/Brouwer-Haemers.html>`_.

    EXAMPLE::

        sage: g = graphs.BrouwerHaemersGraph()
        sage: g
        Brouwer-Haemers: Graph on 81 vertices

    It is indeed strongly regular with parameters `(81,20,1,6)`::

        sage: g.is_strongly_regular(parameters = True) # long time
        (81, 20, 1, 6)

    Its has as eigenvalues `20,2` and `-7`::

        sage: set(g.spectrum()) == {20,2,-7}
        True
    """
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.modules.free_module import VectorSpace
    from sage.matrix.constructor import Matrix
    from sage.matrix.constructor import identity_matrix

    d = 4
    q = 3
    F = FiniteField(q,"x")
    V = VectorSpace(F,d)
    M = Matrix(F,identity_matrix(d))
    M[1,1]=-1
    G = Graph([map(tuple,V), lambda x,y:(V(x)-V(y))*(M*(V(x)-V(y))) == 0], loops = False)
    G.relabel()
    ordering = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                18, 19, 20, 21, 22, 23, 24, 25, 26, 48, 49, 50, 51, 52, 53,
                45, 46, 47, 30, 31, 32, 33, 34, 35, 27, 28, 29, 39, 40, 41,
                42, 43, 44, 36, 37, 38, 69, 70, 71, 63, 64, 65, 66, 67, 68,
                78, 79, 80, 72, 73, 74, 75, 76, 77, 60, 61, 62, 54, 55, 56,
                57, 58, 59]
    _circle_embedding(G, ordering)
    G.name("Brouwer-Haemers")
    return G

def BuckyBall():
    r"""
    Create the Bucky Ball graph.

    This graph is a 3-regular 60-vertex planar graph. Its vertices
    and edges correspond precisely to the carbon atoms and bonds
    in buckminsterfullerene.  When embedded on a sphere, its 12
    pentagon and 20 hexagon faces are arranged exactly as the
    sections of a soccer ball.

    EXAMPLES:

    The Bucky Ball is planar. ::

        sage: g = graphs.BuckyBall()
        sage: g.is_planar()
        True

    The Bucky Ball can also be created by extracting the 1-skeleton
    of the Bucky Ball polyhedron, but this is much slower. ::

        sage: g = polytopes.buckyball().vertex_graph()
        sage: g.remove_loops()
        sage: h = graphs.BuckyBall()
        sage: g.is_isomorphic(h)
        True

    The graph is returned along with an attractive embedding. ::

        sage: g = graphs.BuckyBall()
        sage: g.plot(vertex_labels=False, vertex_size=10).show() # long time
    """
    edges = [(0, 2), (0, 48), (0, 59), (1, 3), (1, 9), (1, 58),
             (2, 3), (2, 36), (3, 17), (4, 6), (4, 8), (4, 12),
             (5, 7), (5, 9), (5, 16), (6, 7), (6, 20), (7, 21),
             (8, 9), (8, 56), (10, 11), (10, 12), (10, 20), (11, 27),
             (11, 47), (12, 13), (13, 46), (13, 54), (14, 15), (14, 16),
             (14, 21), (15, 25), (15, 41), (16, 17), (17, 40), (18, 19),
             (18, 20), (18, 26), (19, 21), (19, 24), (22, 23), (22, 31),
             (22, 34), (23, 25), (23, 38), (24, 25), (24, 30), (26, 27),
             (26, 30), (27, 29), (28, 29), (28, 31), (28, 35), (29, 44),
             (30, 31), (32, 34), (32, 39), (32, 50), (33, 35), (33, 45),
             (33, 51), (34, 35), (36, 37), (36, 40), (37, 39), (37, 52),
             (38, 39), (38, 41), (40, 41), (42, 43), (42, 46), (42, 55),
             (43, 45), (43, 53), (44, 45), (44, 47), (46, 47), (48, 49),
             (48, 52), (49, 53), (49, 57), (50, 51), (50, 52), (51, 53),
             (54, 55), (54, 56), (55, 57), (56, 58), (57, 59), (58, 59)
             ]
    g = Graph()
    g.add_edges(edges)
    g.name("Bucky Ball")

    pos = {
        0 :  (1.00000000000000, 0.000000000000000),
        1 :  (-1.00000000000000, 0.000000000000000),
        2 :  (0.500000000000000, 0.866025403784439),
        3 :  (-0.500000000000000, 0.866025403784439),
        4 :  (-0.252886764483159, -0.146004241548845),
        5 :  (-0.368953972399043, 0.0928336233191176),
        6 :  (-0.217853192651371, -0.0480798425451855),
        7 :  (-0.255589950938772, 0.0495517623332213),
        8 :  (-0.390242139418333, -0.225306404242310),
        9 :  (-0.586398703939125, -0.0441575936410641),
        10:  (-0.113926229169631, -0.101751920396670),
        11:  (-0.0461308635969359, -0.0928422349110366),
        12:  (-0.150564961379772, -0.164626477859040),
        13:  (-0.0848818904865275, -0.246123271631605),
        14:  (-0.170708060452244, 0.196571509298384),
        15:  (-0.0672882312715990, 0.212706320404226),
        16:  (-0.264873262319233, 0.273106701265196),
        17:  (-0.254957754106411, 0.529914971178085),
        18:  (-0.103469165775548, 0.00647061768205703),
        19:  (-0.113590051906687, 0.0655812470455896),
        20:  (-0.145082862532183, -0.0477870484199328),
        21:  (-0.179962687765901, 0.103901506225732),
        22:  (0.0573383021786124, 0.0863716172289798),
        23:  (0.0311566333625530, 0.149538968816603),
        24:  (-0.0573383021786121, 0.0863716172289799),
        25:  (-0.0311566333625527, 0.149538968816603),
        26:  (-0.0517345828877740, 0.00161765442051429),
        27:  (-0.0244663616211774, -0.0456122902452611),
        28:  (0.0517345828877743, 0.00161765442051431),
        29:  (0.0244663616211777, -0.0456122902452611),
        30:  (-0.0272682212665964, 0.0439946358247470),
        31:  (0.0272682212665968, 0.0439946358247470),
        32:  (0.179962687765901, 0.103901506225732),
        33:  (0.145082862532184, -0.0477870484199329),
        34:  (0.113590051906687, 0.0655812470455895),
        35:  (0.103469165775548, 0.00647061768205698),
        36:  (0.254957754106411, 0.529914971178085),
        37:  (0.264873262319233, 0.273106701265196),
        38:  (0.0672882312715993, 0.212706320404226),
        39:  (0.170708060452245, 0.196571509298384),
        40:  (1.59594559789866e-16, 0.450612808484620),
        41:  (2.01227923213310e-16, 0.292008483097691),
        42:  (0.0848818904865278, -0.246123271631605),
        43:  (0.150564961379773, -0.164626477859040),
        44:  (0.0461308635969362, -0.0928422349110366),
        45:  (0.113926229169631, -0.101751920396670),
        46:  (1.66533453693773e-16, -0.207803012451463),
        47:  (1.80411241501588e-16, -0.131162494091179),
        48:  (0.586398703939126, -0.0441575936410641),
        49:  (0.390242139418333, -0.225306404242310),
        50:  (0.255589950938772, 0.0495517623332212),
        51:  (0.217853192651372, -0.0480798425451855),
        52:  (0.368953972399044, 0.0928336233191175),
        53:  (0.252886764483159, -0.146004241548845),
        54:  (-0.104080710079810, -0.365940324584313),
        55:  (0.104080710079811, -0.365940324584313),
        56:  (-0.331440949832714, -0.485757377537020),
        57:  (0.331440949832715, -0.485757377537021),
        58:  (-0.500000000000000, -0.866025403784438),
        59:  (0.500000000000000, -0.866025403784439)
    }

    g.set_pos(pos)

    return g

def DoubleStarSnark():
    r"""
    Returns the double star snark.

    The double star snark is a 3-regular graph on 30 vertices. See
    the :wikipedia:`Wikipedia page on the double star snark
    <Double-star_snark>`.

    EXAMPLES::

        sage: g = graphs.DoubleStarSnark()
        sage: g.order()
        30
        sage: g.size()
        45
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian()
        False
        sage: g.automorphism_group().cardinality()
        80
        sage: g.show()
    """

    d = { 0: [1, 14, 15]
        , 1: [0, 2, 11]
        , 2: [1, 3, 7]
        , 3: [2, 4, 18]
        , 4: [3, 5, 14]
        , 5: [10, 4, 6]
        , 6: [5, 21, 7]
        , 7: [8, 2, 6]
        , 8: [9, 13, 7]
        , 9: [24, 8, 10]
        , 10: [9, 11, 5]
        , 11: [1, 10, 12]
        , 12: [11, 27, 13]
        , 13: [8, 12, 14]
        , 14: [0, 4, 13]
        , 15: [0, 16, 29]
        , 16: [15, 20, 23]
        , 17: [25, 18, 28]
        , 18: [3, 17, 19]
        , 19: [18, 26, 23]
        , 20: [16, 28, 21]
        , 21: [20, 6, 22]
        , 22: [26, 21, 29]
        , 23: [16, 24, 19]
        , 24: [25, 9, 23]
        , 25: [24, 17, 29]
        , 26: [27, 19, 22]
        , 27: [12, 26, 28]
        , 28: [17, 27, 20]
        , 29: [25, 22, 15]
        }

    g = Graph(d, pos={}, name="Double star snark")
    _circle_embedding(g, range(15), radius=2)
    _circle_embedding(g, range(15, 30), radius=1.4)

    return g

def MeredithGraph():
    r"""
    Returns the Meredith Graph

    The Meredith Graph is a 4-regular 4-connected non-hamiltonian graph. For
    more information on the Meredith Graph, see the :wikipedia:`Meredith_graph`.

    EXAMPLES::

        sage: g = graphs.MeredithGraph()
        sage: g.is_regular(4)
        True
        sage: g.order()
        70
        sage: g.size()
        140
        sage: g.radius()
        7
        sage: g.diameter()
        8
        sage: g.girth()
        4
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian() # long time
        False
    """
    g = Graph(name="Meredith Graph")
    g.add_vertex(0)

    # Edges between copies of K_{4,3}
    for i in range(5):
        g.add_edge(('outer',i,3),('outer',(i+1)%5,0))
        g.add_edge(('inner',i,3),('inner',(i+2)%5,0))
        g.add_edge(('outer',i,1),('inner',i      ,1))
        g.add_edge(('outer',i,2),('inner',i      ,2))

    # Edges inside of the K_{4,3}s.
    for i in range(5):
        for j in range(4):
            for k in range(3):
                g.add_edge(('inner',i,j),('inner',i,k+4))
                g.add_edge(('outer',i,j),('outer',i,k+4))

    _circle_embedding(g, sum([[('outer',i,j) for j in range(4)]+10*[0] for i in range(5)],[]), radius = 1, shift = 2)
    _circle_embedding(g, sum([[('outer',i,j) for j in range(4,7)]+10*[0] for i in range(5)],[]), radius = 1.2, shift = 2.2)
    _circle_embedding(g, sum([[('inner',i,j) for j in range(4)]+7*[0] for i in range(5)],[]), radius = .6, shift = 1.24)
    _circle_embedding(g, sum([[('inner',i,j) for j in range(4,7)]+5*[0] for i in range(5)],[]), radius = .4, shift = 1.05)

    g.delete_vertex(0)
    g.relabel()
    return g

def KittellGraph():
    r"""
    Returns the Kittell Graph.

    For more information on the Kittell Graph, see the `corresponding Wolfram
    page <http://mathworld.wolfram.com/KittellGraph.html>`_.

    EXAMPLES::

        sage: g = graphs.KittellGraph()
        sage: g.order()
        23
        sage: g.size()
        63
        sage: g.radius()
        3
        sage: g.diameter()
        4
        sage: g.girth()
        3
        sage: g.chromatic_number()
        4
    """
    g = Graph({0: [1, 2, 4, 5, 6, 7], 1: [0, 2, 7, 10, 11, 13],
               2: [0, 1, 11, 4, 14], 3: [16, 12, 4, 5, 14], 4: [0, 2, 3, 5, 14],
               5: [0, 16, 3, 4, 6], 6: [0, 5, 7, 15, 16, 17, 18],
               7: [0, 1, 6, 8, 13, 18], 8: [9, 18, 19, 13, 7],
               9: [8, 10, 19, 20, 13], 10: [1, 9, 11, 13, 20, 21],
               11: [1, 2, 10, 12, 14, 15, 21], 12: [11, 16, 3, 14, 15],
               13: [8, 1, 10, 9, 7], 14: [11, 12, 2, 3, 4],
               15: [6, 11, 12, 16, 17, 21, 22],
               16: [3, 12, 5, 6, 15], 17: [18, 19, 22, 6, 15],
               18: [8, 17, 19, 6, 7], 19: [8, 9, 17, 18, 20, 22],
               20: [9, 10, 19, 21, 22], 21: [10, 11, 20, 22, 15],
               22: [17, 19, 20, 21, 15]},
              name = "Kittell Graph")

    _circle_embedding(g, range(3), shift=.75)
    _circle_embedding(g, range(3,13), radius = .4)
    _circle_embedding(g, range(15,22), radius = .2, shift=-.15)
    pos = g.get_pos()
    pos[13] = (-.65,-.35)
    pos[14] = (.65,-.35)
    pos[22] = (0,0)

    return g

def CameronGraph():
    r"""
    Returns the Cameron graph.

    The Cameron graph is strongly regular with parameters `v = 231, k = 30,
    \lambda = 9, \mu = 3`.

    For more information on the Cameron graph, see
    `<http://www.win.tue.nl/~aeb/graphs/Cameron.html>`_.

    EXAMPLES::

        sage: g = graphs.CameronGraph()
        sage: g.order()
        231
        sage: g.size()
        3465
        sage: g.is_strongly_regular(parameters = True) # long time
        (231, 30, 9, 3)

    """
    from sage.groups.perm_gps.permgroup_named import MathieuGroup
    from itertools import combinations
    g = Graph(name="Cameron Graph")
    sets = MathieuGroup(22).orbit((1,2,3,7,10,20), action = "OnSets")
    for s in sets:
        for a,b,c,d in combinations(set(s),4):
            g.add_edges([((a,b),(c,d)),((a,c),(b,d)), ((a,d),(b,c))])

    g.relabel()
    ordering = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 18, 19, 20,
                21, 24, 25, 26, 27, 29, 31, 34, 35, 38, 39, 96, 97, 101, 105,
                51, 117, 198, 32, 196, 201, 131, 167, 199, 197, 86, 102, 195,
                200, 186, 144, 202, 177, 44, 53, 58, 45, 48, 54, 43, 57, 50,
                46, 59, 133, 169, 104, 188, 118, 208, 157, 52, 207, 209, 132,
                204, 13, 187, 33, 203, 70, 145, 103, 168, 178, 87, 124, 123,
                125, 111, 120, 116, 119, 112, 95, 114, 115, 137, 218, 213, 108,
                76, 77, 74, 62, 64, 67, 63, 68, 69, 61, 41, 75, 73, 66, 71, 72,
                60, 22, 230, 151, 184, 138, 193, 109, 228, 174, 214, 219, 93,
                126, 143, 150, 146, 224, 181, 16, 223, 171, 90, 135, 106, 205,
                211, 121, 148, 160, 216, 222, 190, 36, 55, 185, 175, 94, 139,
                110, 215, 152, 220, 229, 194, 40, 128, 99, 141, 173, 154, 82,
                156, 164, 159, 28, 127, 158, 65, 162, 163, 153, 161, 155, 140,
                98, 47, 113, 84, 180, 30, 129, 179, 183, 165, 176, 142, 100,
                49, 134, 210, 170, 147, 91, 37, 206, 182, 191, 56, 136, 225,
                221, 149, 227, 217, 17, 107, 172, 212, 122, 226, 23, 85, 42,
                80, 92, 81, 89, 78, 83, 88, 79, 130, 192, 189, 166]

    _circle_embedding(g, ordering)
    return g

def ChvatalGraph():
    r"""
    Returns the Chvatal graph.

    Chvatal graph is one of the few known graphs to satisfy Grunbaum's
    conjecture that for every m, n, there is an m-regular,
    m-chromatic graph of girth at least n. For more information, see this
    `Wikipedia article on the Chvatal graph <http://en.wikipedia.org/wiki/Chv%C3%A1tal_graph>`_.

    EXAMPLES:

    The Chvatal graph has 12 vertices and 24 edges. It is a 4-regular,
    4-chromatic graph with radius 2, diameter 2, and girth 4. ::

        sage: G = graphs.ChvatalGraph(); G
        Chvatal graph: Graph on 12 vertices
        sage: G.order(); G.size()
        12
        24
        sage: G.degree()
        [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
        sage: G.chromatic_number()
        4
        sage: G.radius(); G.diameter(); G.girth()
        2
        2
        4
    """
    import networkx
    pos_dict = {}
    for i in range(5, 10):
        x = float(cos((pi / 2) + ((2 * pi) / 5) * i))
        y = float(sin((pi / 2) + ((2 * pi) / 5) * i))
        pos_dict[i] = (x, y)
    for i in range(5):
        x = float(2 * (cos((pi / 2) + ((2 * pi) / 5) * (i - 5))))
        y = float(2 * (sin((pi / 2) + ((2 * pi) / 5) * (i - 5))))
        pos_dict[i] = (x, y)
    pos_dict[10] = (0.5, 0)
    pos_dict[11] = (-0.5, 0)

    return Graph(networkx.chvatal_graph(), pos=pos_dict, name="Chvatal graph")

def ClebschGraph():
    r"""
    Return the Clebsch graph.

    EXAMPLES::

        sage: g = graphs.ClebschGraph()
        sage: g.automorphism_group().cardinality()
        1920
        sage: g.girth()
        4
        sage: g.chromatic_number()
        4
        sage: g.diameter()
        2
        sage: g.show(figsize=[10, 10]) # long time
    """
    g = Graph(pos={})
    x = 0
    for i in range(8):
        g.add_edge(x % 16, (x + 1) % 16)
        g.add_edge(x % 16, (x + 6) % 16)
        g.add_edge(x % 16, (x + 8) % 16)
        x += 1
        g.add_edge(x % 16, (x + 3) % 16)
        g.add_edge(x % 16, (x + 2) % 16)
        g.add_edge(x % 16, (x + 8) % 16)
        x += 1

    _circle_embedding(g, range(16), shift=.5)
    g.name("Clebsch graph")

    return g

def CoxeterGraph():
    r"""
    Return the Coxeter graph.

    See the :wikipedia:`Wikipedia page on the Coxeter graph
    <Coxeter_graph>`.

    EXAMPLES::

        sage: g = graphs.CoxeterGraph()
        sage: g.automorphism_group().cardinality()
        336
        sage: g.girth()
        7
        sage: g.chromatic_number()
        3
        sage: g.diameter()
        4
        sage: g.show(figsize=[10, 10]) # long time
    """
    g = Graph({
            27: [6, 22, 14],
            24: [0, 7, 18],
            25: [8, 15, 2],
            26: [10, 16, 23],
            }, pos={})

    g.add_cycle(range(24))
    g.add_edges([(5, 11), (9, 20), (12, 1), (13, 19), (17, 4), (3, 21)])

    _circle_embedding(g, range(24))
    _circle_embedding(g, [24, 25, 26], radius=.5)
    g.get_pos()[27] = (0, 0)

    g.name("Coxeter Graph")

    return g

def DesarguesGraph():
    """
    Returns the Desargues graph.

    PLOTTING: The layout chosen is the same as on the cover of [1].

    EXAMPLE::

        sage: D = graphs.DesarguesGraph()
        sage: L = graphs.LCFGraph(20,[5,-5,9,-9],5)
        sage: D.is_isomorphic(L)
        True
        sage: D.show()  # long time

    REFERENCE:

    - [1] Harary, F. Graph Theory. Reading, MA: Addison-Wesley,
      1994.
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    G = GeneralizedPetersenGraph(10,3)
    G.name("Desargues Graph")
    return G

def DurerGraph():
    r"""
    Returns the Dürer graph.

    For more information, see this
    `Wikipedia article on the Dürer graph <http://en.wikipedia.org/wiki/D%C3%BCrer_graph>`_.

    EXAMPLES:

    The Dürer graph is named after Albrecht Dürer. It is a planar graph
    with 12 vertices and 18 edges. ::

        sage: G = graphs.DurerGraph(); G
        Durer graph: Graph on 12 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        12
        sage: G.size()
        18

    The Dürer graph has chromatic number 3, diameter 4, and girth 3. ::

        sage: G.chromatic_number()
        3
        sage: G.diameter()
        4
        sage: G.girth()
        3

    Its automorphism group is isomorphic to `D_6`. ::

        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(6))
        True
    """
    edge_dict = {
        0: [1,5,6],
        1: [2,7],
        2: [3,8],
        3: [4,9],
        4: [5,10],
        5: [11],
        6: [8,10],
        7: [9,11],
        8: [10],
        9: [11]}
    pos_dict = {
        0: [2, 0],
        1: [1, 1.73205080756888],
        2: [-1, 1.73205080756888],
        3: [-2, 0],
        4: [-1, -1.73205080756888],
        5: [1, -1.73205080756888],
        6: [1, 0],
        7: [0.5, 0.866025403784439],
        8: [-0.5, 0.866025403784439],
        9: [-1, 0],
        10: [-0.5, -0.866025403784439],
        11: [0.5, -0.866025403784439]}
    return Graph(edge_dict, pos=pos_dict, name="Durer graph")

def DyckGraph():
    """
    Returns the Dyck graph.

    For more information, see the `MathWorld article on the Dyck graph
    <http://mathworld.wolfram.com/DyckGraph.html>`_ or the `Wikipedia
    article on the Dyck graph <http://en.wikipedia.org/wiki/Dyck_graph>`_.

    EXAMPLES:

    The Dyck graph was defined by Walther von Dyck in 1881. It has `32`
    vertices and `48` edges, and is a cubic graph (regular of degree `3`)::

        sage: G = graphs.DyckGraph(); G
        Dyck graph: Graph on 32 vertices
        sage: G.order()
        32
        sage: G.size()
        48
        sage: G.is_regular()
        True
        sage: G.is_regular(3)
        True

    It is non-planar and Hamiltonian, as well as bipartite (making it a
    bicubic graph)::

        sage: G.is_planar()
        False
        sage: G.is_hamiltonian()
        True
        sage: G.is_bipartite()
        True

    It has radius `5`, diameter `5`, and girth `6`::

        sage: G.radius()
        5
        sage: G.diameter()
        5
        sage: G.girth()
        6

    Its chromatic number is `2` and its automorphism group is of order
    `192`::

        sage: G.chromatic_number()
        2
        sage: G.automorphism_group().cardinality()
        192

    It is a non-integral graph as it has irrational eigenvalues::

        sage: G.characteristic_polynomial().factor()
        (x - 3) * (x + 3) * (x - 1)^9 * (x + 1)^9 * (x^2 - 5)^6

    It is a toroidal graph, and its embedding on a torus is dual to an
    embedding of the Shrikhande graph (:meth:`ShrikhandeGraph
    <GraphGenerators.ShrikhandeGraph>`).
    """
    pos_dict = {}
    for i in range(8):
        pos_dict[i] = [float(cos((2*i) * pi/8)),
                       float(sin((2*i) * pi/8))]
        pos_dict[8 + i]  = [0.75 * pos_dict[i][0],
                            0.75 * pos_dict[i][1]]
        pos_dict[16 + i] = [0.50 * pos_dict[i][0],
                            0.50 * pos_dict[i][1]]
        pos_dict[24 + i] = [0.25 * pos_dict[i][0],
                            0.25 * pos_dict[i][1]]

    edge_dict = {
        0O00: [0O07, 0O01,   0O10], 0O10: [0O00,   0O27, 0O21],
        0O01: [0O00, 0O02,   0O11], 0O11: [0O01,   0O20, 0O22],
        0O02: [0O01, 0O03,   0O12], 0O12: [0O02,   0O21, 0O23],
        0O03: [0O02, 0O04,   0O13], 0O13: [0O03,   0O22, 0O24],
        0O04: [0O03, 0O05,   0O14], 0O14: [0O04,   0O23, 0O25],
        0O05: [0O04, 0O06,   0O15], 0O15: [0O05,   0O24, 0O26],
        0O06: [0O05, 0O07,   0O16], 0O16: [0O06,   0O25, 0O27],
        0O07: [0O06, 0O00,   0O17], 0O17: [0O07,   0O26, 0O20],

        0O20: [0O17, 0O11,   0O30], 0O30: [0O20,   0O35, 0O33],
        0O21: [0O10, 0O12,   0O31], 0O31: [0O21,   0O36, 0O34],
        0O22: [0O11, 0O13,   0O32], 0O32: [0O22,   0O37, 0O35],
        0O23: [0O12, 0O14,   0O33], 0O33: [0O23,   0O30, 0O36],
        0O24: [0O13, 0O15,   0O34], 0O34: [0O24,   0O31, 0O37],
        0O25: [0O14, 0O16,   0O35], 0O35: [0O25,   0O32, 0O30],
        0O26: [0O15, 0O17,   0O36], 0O36: [0O26,   0O33, 0O31],
        0O27: [0O16, 0O10,   0O37], 0O37: [0O27,   0O34, 0O32],
    }

    return Graph(edge_dict, pos=pos_dict, name="Dyck graph")

def HortonGraph():
    r"""
    Returns the Horton Graph.

    The Horton graph is a cubic 3-connected non-hamiltonian graph. For more
    information, see the :wikipedia:`Horton_graph`.

    EXAMPLES::

        sage: g = graphs.HortonGraph()
        sage: g.order()
        96
        sage: g.size()
        144
        sage: g.radius()
        10
        sage: g.diameter()
        10
        sage: g.girth()
        6
        sage: g.automorphism_group().cardinality()
        96
        sage: g.chromatic_number()
        2
        sage: g.is_hamiltonian() # not tested -- veeeery long
        False
    """
    g = Graph(name = "Horton Graph")

    # Each group of the 6 groups of vertices is based on the same 3-regular
    # graph.
    from sage.graphs.generators.families import LCFGraph
    lcf = LCFGraph(16,[5,-5],8)
    lcf.delete_edge(15,0)
    lcf.delete_edge(7,8)

    for i in range(6):
        for u,v in lcf.edges(labels=False):
            g.add_edge((i,u),(i,v))

    # Modifying the groups and linking them together
    for i in range(3):
        g.add_edge((2*i,0),(2*i+1,7))
        g.add_edge((2*i+1,8),(2*i,7))
        g.add_edge((2*i,15),(2*i+1,0))
        g.add_edge((2*i,8),1)
        g.add_edge((2*i+1,14),2)
        g.add_edge((2*i+1,10),0)

    # Embedding
    for i in range(6):
        _circle_embedding(g, [(i,j) for j in range(16)], center=(cos(2*i*pi/6),sin(2*i*pi/6)), radius=.3)

    for i in range(3):
        g.delete_vertex((2*i+1,15))

    _circle_embedding(g, range(3), radius=.2, shift=-0.75)

    g.relabel()

    return g

def EllinghamHorton54Graph():
    r"""
    Returns the Ellingham-Horton 54-graph.

    For more information, see the :wikipedia:`Wikipedia page on the
    Ellingham-Horton graphs <Ellingham-Horton_graph>`

    EXAMPLE:

    This graph is 3-regular::

        sage: g = graphs.EllinghamHorton54Graph()
        sage: g.is_regular(k=3)
        True

    It is 3-connected and bipartite::

        sage: g.vertex_connectivity() # not tested - too long
        3
        sage: g.is_bipartite()
        True

    It is not Hamiltonian::

        sage: g.is_hamiltonian() # not tested - too long
        False

    ... and it has a nice drawing ::

        sage: g.show(figsize=[10, 10]) # not tested - too long

    TESTS::

        sage: g.show() # long time
    """
    from sage.graphs.generators.basic import CycleGraph
    up = CycleGraph(16)
    low = 2*CycleGraph(6)

    for v in range(6):
        low.add_edge(v, v + 12)
        low.add_edge(v + 6, v + 12)
    low.add_edge(12, 15)
    low.delete_edge(1, 2)
    low.delete_edge(8, 7)
    low.add_edge(1, 8)
    low.add_edge(7, 2)


    # The set of vertices on top is 0..15
    # Bottom left is 16..33
    # Bottom right is 34..52
    # The two other vertices are 53, 54
    g = up + 2*low
    g.name("Ellingham-Horton 54-graph")
    g.set_pos({})

    g.add_edges([(15, 4), (3, 8), (7, 12), (11, 0), (2, 13), (5, 10)])
    g.add_edges([(30, 6), (29, 9), (48, 14), (47, 1)])
    g.add_edge(32, 52)
    g.add_edge(50, 52)
    g.add_edge(33, 53)
    g.add_edge(51, 53)
    g.add_edge(52, 53)

    # Top
    _circle_embedding(g, range(16), center=(0, .5), shift=.5, radius=.5)

    # Bottom-left
    _circle_embedding(g, range(16, 22), center=(-1.5, -1))
    _circle_embedding(g, range(22, 28), center=(-1.5, -1), radius=.5)
    _circle_embedding(g, range(28, 34), center=(-1.5, -1), radius=.7)

    # Bottom right
    _circle_embedding(g, range(34, 40), center=(1.5, -1))
    _circle_embedding(g, range(40, 46), center=(1.5, -1), radius=.5)
    _circle_embedding(g, range(46, 52), center=(1.5, -1), radius=.7)

    d = g.get_pos()
    d[52] = (-.3, -2.5)
    d[53] = (.3, -2.5)
    d[31] = (-2.2, -.9)
    d[28] = (-.8, -.9)
    d[46] = (2.2, -.9)
    d[49] = (.8, -.9)


    return g

def EllinghamHorton78Graph():
    r"""
    Returns the Ellingham-Horton 78-graph.

    For more information, see the :wikipedia:`Wikipedia page on the
    Ellingham-Horton graphs
    <http://en.wikipedia.org/wiki/Ellingham%E2%80%93Horton_graph>`

    EXAMPLE:

    This graph is 3-regular::

        sage: g = graphs.EllinghamHorton78Graph()
        sage: g.is_regular(k=3)
        True

    It is 3-connected and bipartite::

        sage: g.vertex_connectivity() # not tested - too long
        3
        sage: g.is_bipartite()
        True

    It is not Hamiltonian::

        sage: g.is_hamiltonian() # not tested - too long
        False

    ... and it has a nice drawing ::

        sage: g.show(figsize=[10,10]) # not tested - too long

    TESTS::

        sage: g.show(figsize=[10, 10]) # not tested - too long
    """
    g = Graph({
            0: [1, 5, 60], 1: [2, 12], 2: [3, 7], 3: [4, 14], 4: [5, 9],
            5: [6], 6: [7, 11], 7: [15], 8: [9, 13, 22], 9: [10],
            10: [11, 72], 11: [12], 12: [13], 13: [14], 14: [72],
            15: [16, 20], 16: [17, 27], 17: [18, 22], 18: [19, 29],
            19: [20, 24], 20: [21], 21: [22, 26], 23: [24, 28, 72],
            24: [25], 25: [26, 71], 26: [27], 27: [28], 28: [29],
            29: [69], 30: [31, 35, 52], 31: [32, 42], 32: [33, 37],
            33: [34, 43], 34: [35, 39], 35: [36], 36: [41, 63],
            37: [65, 66], 38: [39, 59, 74], 39: [40], 40: [41, 44],
            41: [42], 42: [74], 43: [44, 74], 44: [45], 45: [46, 50],
            46: [47, 57], 47: [48, 52], 48: [49, 75], 49: [50, 54],
            50: [51], 51: [52, 56], 53: [54, 58, 73], 54: [55],
            55: [56, 59], 56: [57], 57: [58], 58: [75], 59: [75],
            60: [61, 64], 61: [62, 71], 62: [63, 77], 63: [67],
            64: [65, 69], 65: [77], 66: [70, 73], 67: [68, 73],
            68: [69, 76], 70: [71, 76], 76: [77]}, pos={})

    _circle_embedding(g, range(15), center=(-2.5, 1.5))
    _circle_embedding(g, range(15, 30), center=(-2.5, -1.5))
    _circle_embedding(g, [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
        42, 74, 43, 44], center=(2.5, 1.5))
    _circle_embedding(g, [45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
        57, 58, 75, 59], center=(2.5, -1.5))

    d = g.get_pos()

    d[76] = (-.2, -.1)
    d[77] = (.2, .1)
    d[38] = (2.2, .1)
    d[52] = (2.3, -.1)
    d[15] = (-2.1, -.1)
    d[72] = (-2.1, .1)

    _line_embedding(g, [60, 61, 62, 63], first=(-1, 2), last=(1, 2))
    _line_embedding(g, [64, 65, 37], first=(-.5, 1.5), last=(1.2, 1.5))
    _line_embedding(g, [66, 73, 67, 68, 69], first=(1.2, -2),
            last=(-.8, -2))
    _line_embedding(g, [66, 70, 71], first=(.7, -1.5), last=(-1, -1.5))

    g.name("Ellingham-Horton 78-graph")

    return g

def ErreraGraph():
    r"""
    Returns the Errera graph.

    For more information, see this
    `Wikipedia article on the Errera graph <http://en.wikipedia.org/wiki/Errera_graph>`_.

    EXAMPLES:

    The Errera graph is named after Alfred Errera. It is a planar graph
    on 17 vertices and having 45 edges. ::

        sage: G = graphs.ErreraGraph(); G
        Errera graph: Graph on 17 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        17
        sage: G.size()
        45

    The Errera graph is Hamiltonian with radius 3, diameter 4, girth 3,
    and chromatic number 4. ::

        sage: G.is_hamiltonian()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        4
        sage: G.girth()
        3
        sage: G.chromatic_number()
        4

    Each vertex degree is either 5 or 6. That is, if `f` counts the
    number of vertices of degree 5 and `s` counts the number of vertices
    of degree 6, then `f + s` is equal to the order of the Errera
    graph. ::

        sage: D = G.degree_sequence()
        sage: D.count(5) + D.count(6) == G.order()
        True

    The automorphism group of the Errera graph is isomorphic to the
    dihedral group of order 20. ::

        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(10))
        True
    """
    edge_dict = {
        0: [1,7,14,15,16],
        1: [2,9,14,15],
        2: [3,8,9,10,14],
        3: [4,9,10,11],
        4: [5,10,11,12],
        5: [6,11,12,13],
        6: [7,8,12,13,16],
        7: [13,15,16],
        8: [10,12,14,16],
        9: [11,13,15],
        10: [12],
        11: [13],
        13: [15],
        14: [16]}
    return Graph(edge_dict, name="Errera graph")

def FlowerSnark():
    """
    Returns a Flower Snark.

    A flower snark has 20 vertices. It is part of the class of
    biconnected cubic graphs with edge chromatic number = 4, known as
    snarks. (i.e.: the Petersen graph). All snarks are not Hamiltonian,
    non-planar and have Petersen graph graph minors.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the nodes are
    drawn 0-14 on the outer circle, and 15-19 in an inner pentagon.

    REFERENCES:

    - [1] Weisstein, E. (1999). "Flower Snark - from Wolfram
      MathWorld". [Online] Available:
      http://mathworld.wolfram.com/FlowerSnark.html [2007, February 17]

    EXAMPLES: Inspect a flower snark::

        sage: F = graphs.FlowerSnark()
        sage: F
        Flower Snark: Graph on 20 vertices
        sage: F.graph6_string()
        'ShCGHC@?GGg@?@?Gp?K??C?CA?G?_G?Cc'

    Now show it::

        sage: F.show() # long time
    """
    pos_dict = {}
    for i in range(15):
        x = float(2.5*(cos((pi/2) + ((2*pi)/15)*i)))
        y = float(2.5*(sin((pi/2) + ((2*pi)/15)*i)))
        pos_dict[i] = (x,y)
    for i in range(15,20):
        x = float(cos((pi/2) + ((2*pi)/5)*i))
        y = float(sin((pi/2) + ((2*pi)/5)*i))
        pos_dict[i] = (x,y)
    return Graph({0:[1,14,15],1:[2,11],2:[3,7],3:[2,4,16],4:[5,14], \
                        5:[6,10],6:[5,7,17],8:[7,9,13],9:[10,18],11:[10,12], \
                        12:[13,19],13:[14],15:[19],16:[15,17],18:[17,19]}, \
                        pos=pos_dict, name="Flower Snark")


def FolkmanGraph():
    """
    Returns the Folkman graph.

    See the :wikipedia:`Wikipedia page on the Folkman Graph
    <Folkman_graph>`.

    EXAMPLE::

        sage: g = graphs.FolkmanGraph()
        sage: g.order()
        20
        sage: g.size()
        40
        sage: g.diameter()
        4
        sage: g.girth()
        4
        sage: g.charpoly().factor()
        (x - 4) * (x + 4) * x^10 * (x^2 - 6)^4
        sage: g.chromatic_number()
        2
        sage: g.is_eulerian()
        True
        sage: g.is_hamiltonian()
        True
        sage: g.is_vertex_transitive()
        False
        sage: g.is_bipartite()
        True
    """
    from sage.graphs.generators.families import LCFGraph
    g= LCFGraph(20, [5, -7, -7, 5], 5)
    g.name("Folkman Graph")
    return g


def FosterGraph():
    """
    Returns the Foster graph.

    See the :wikipedia:`Wikipedia page on the Foster Graph
    <Foster_graph>`.

    EXAMPLE::

        sage: g = graphs.FosterGraph()
        sage: g.order()
        90
        sage: g.size()
        135
        sage: g.diameter()
        8
        sage: g.girth()
        10
        sage: g.automorphism_group().cardinality()
        4320
        sage: g.is_hamiltonian()
        True
    """
    from sage.graphs.generators.families import LCFGraph
    g= LCFGraph(90, [17, -9, 37, -37, 9, -17], 15)
    g.name("Foster Graph")
    return g


def FranklinGraph():
    r"""
    Returns the Franklin graph.

    For more information, see this
    `Wikipedia article on the Franklin graph <http://en.wikipedia.org/wiki/Franklin_graph>`_.

    EXAMPLES:

    The Franklin graph is named after Philip Franklin. It is a
    3-regular graph on 12 vertices and having 18 edges. ::

        sage: G = graphs.FranklinGraph(); G
        Franklin graph: Graph on 12 vertices
        sage: G.is_regular(3)
        True
        sage: G.order()
        12
        sage: G.size()
        18

    The Franklin graph is a Hamiltonian, bipartite graph with radius 3,
    diameter 3, and girth 4. ::

        sage: G.is_hamiltonian()
        True
        sage: G.is_bipartite()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        3
        sage: G.girth()
        4

    It is a perfect, triangle-free graph having chromatic number 2. ::

        sage: G.is_perfect()
        True
        sage: G.is_triangle_free()
        True
        sage: G.chromatic_number()
        2
    """
    edge_dict = {
        0: [1,5,6],
        1: [2,7],
        2: [3,8],
        3: [4,9],
        4: [5,10],
        5: [11],
        6: [7,9],
        7: [10],
        8: [9,11],
        10: [11]}
    pos_dict = {
        0: [2, 0],
        1: [1, 1.73205080756888],
        2: [-1, 1.73205080756888],
        3: [-2, 0],
        4: [-1, -1.73205080756888],
        5: [1, -1.73205080756888],
        6: [1, 0],
        7: [0.5, 0.866025403784439],
        8: [-0.5, 0.866025403784439],
        9: [-1, 0],
        10: [-0.5, -0.866025403784439],
        11: [0.5, -0.866025403784439]}
    return Graph(edge_dict, pos=pos_dict, name="Franklin graph")

def FruchtGraph():
    """
    Returns a Frucht Graph.

    A Frucht graph has 12 nodes and 18 edges. It is the smallest cubic
    identity graph. It is planar and it is Hamiltonian.

    This constructor is dependent on NetworkX's numeric labeling.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the first
    seven nodes are on the outer circle, with the next four on an inner
    circle and the last in the center.

    REFERENCES:

    - [1] Weisstein, E. (1999). "Frucht Graph - from Wolfram
      MathWorld". [Online] Available:
      http://mathworld.wolfram.com/FruchtGraph.html [2007, February 17]

    EXAMPLES::

        sage: FRUCHT = graphs.FruchtGraph()
        sage: FRUCHT
        Frucht graph: Graph on 12 vertices
        sage: FRUCHT.graph6_string()
        'KhCKM?_EGK?L'
        sage: (graphs.FruchtGraph()).show() # long time
    """
    pos_dict = {}
    for i in range(7):
        x = float(2*(cos((pi/2) + ((2*pi)/7)*i)))
        y = float(2*(sin((pi/2) + ((2*pi)/7)*i)))
        pos_dict[i] = (x,y)
    pos_dict[7] = (0,1)
    pos_dict[8] = (-1,0)
    pos_dict[9] = (0,-1)
    pos_dict[10] = (1,0)
    pos_dict[11] = (0,0)
    import networkx
    G = networkx.frucht_graph()
    return Graph(G, pos=pos_dict, name="Frucht graph")

def GoldnerHararyGraph():
    r"""
    Return the Goldner-Harary graph.

    For more information, see this
    `Wikipedia article on the Goldner-Harary graph <http://en.wikipedia.org/wiki/Goldner%E2%80%93Harary_graph>`_.

    EXAMPLES:

    The Goldner-Harary graph is named after A. Goldner and Frank Harary.
    It is a planar graph having 11 vertices and 27 edges. ::

        sage: G = graphs.GoldnerHararyGraph(); G
        Goldner-Harary graph: Graph on 11 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        11
        sage: G.size()
        27

    The Goldner-Harary graph is chordal with radius 2, diameter 2, and
    girth 3. ::

        sage: G.is_chordal()
        True
        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        3

    Its chromatic number is 4 and its automorphism group is isomorphic to
    the dihedral group `D_6`. ::

        sage: G.chromatic_number()
        4
        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(6))
        True
    """
    edge_dict = {
        0: [1,3,4],
        1: [2,3,4,5,6,7,10],
        2: [3,7],
        3: [7,8,9,10],
        4: [3,5,9,10],
        5: [10],
        6: [7,10],
        7: [8,10],
        8: [10],
        9: [10]}

    pos = {
        0: (-2, 0),
        1: (0, 1.5),
        2: (2, 0),
        3: (0, -1.5),
        4: (-1.5, 0),
        5: (-0.5, 0.5),
        6: (0.5, 0.5),
        7: (1.5, 0),
        8: (0.5, -0.5),
        9: (-0.5, -0.5),
        10: (0, 0)}

    return Graph(edge_dict, pos = pos, name="Goldner-Harary graph")

def GrayGraph(embedding=1):
    r"""
    Returns the Gray graph.

    See the :wikipedia:`Wikipedia page on the Gray Graph
    <Gray_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.GrayGraph()
        sage: g.order()
        54
        sage: g.size()
        81
        sage: g.girth()
        8
        sage: g.diameter()
        6
        sage: g.show(figsize=[10, 10])   # long time
        sage: graphs.GrayGraph(embedding = 2).show(figsize=[10, 10])   # long time

    TESTS::

        sage: graphs.GrayGraph(embedding = 3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1, 2, or 3.
    """

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(54, [-25,7,-7,13,-13,25], 9)
    g.name("Gray graph")

    if embedding == 1:
        o = g.automorphism_group(orbits=True)[-1]
        _circle_embedding(g, o[0], center=(0, 0), radius=1)
        _circle_embedding(g, o[1], center=(0, 0), radius=.6, shift=-.5)

    elif embedding != 2:
        raise ValueError("The value of embedding must be 1, 2, or 3.")

    return g

def GrotzschGraph():
    r"""
    Returns the Grötzsch graph.

    The Grötzsch graph is an example of a triangle-free graph with
    chromatic number equal to 4. For more information, see this
    `Wikipedia article on Grötzsch graph <http://en.wikipedia.org/wiki/Gr%C3%B6tzsch_graph>`_.

    REFERENCE:

    - [1] Weisstein, Eric W. "Grotzsch Graph."
      From MathWorld--A Wolfram Web Resource.
      http://mathworld.wolfram.com/GroetzschGraph.html

    EXAMPLES:

    The Grötzsch graph is named after Herbert Grötzsch. It is a
    Hamiltonian graph with 11 vertices and 20 edges. ::

        sage: G = graphs.GrotzschGraph(); G
        Grotzsch graph: Graph on 11 vertices
        sage: G.is_hamiltonian()
        True
        sage: G.order()
        11
        sage: G.size()
        20

    The Grötzsch graph is triangle-free and having radius 2, diameter 2,
    and girth 4. ::

        sage: G.is_triangle_free()
        True
        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        4

    Its chromatic number is 4 and its automorphism group is isomorphic
    to the dihedral group `D_5`. ::

        sage: G.chromatic_number()
        4
        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(5))
        True
    """
    g = Graph()
    g.add_vertices(range(11))

    edges = [];
    for u in range(1,6):
        edges.append( (0,u) )

    edges.append( (10,6) )

    for u in range(6,10):
        edges.append( (u,u+1) )
        edges.append( (u,u-4) )

    edges.append( (10,1) )

    for u in range(7,11):
        edges.append( (u,u-6) )

    edges.append((6,5))

    g.add_edges(edges)

    pos = {}
    pos[0] = (0,0)
    for u in range(1,6):
        theta = (u-1)*2*pi/5
        pos[u] = (float(5*sin(theta)),float(5*cos(theta)))
        pos[u+5] = (2*pos[u][0], 2*pos[u][1])

    g.set_pos(pos)
    g.name("Grotzsch graph")
    return g

def HeawoodGraph():
    """
    Returns a Heawood graph.

    The Heawood graph is a cage graph that has 14 nodes. It is a cubic
    symmetric graph. (See also the Moebius-Kantor graph). It is
    nonplanar and Hamiltonian. It has diameter = 3, radius = 3, girth =
    6, chromatic number = 2. It is 4-transitive but not 5-transitive.

    This constructor is dependent on NetworkX's numeric labeling.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the nodes are
    positioned in a circular layout with the first node appearing at
    the top, and then continuing counterclockwise.

    REFERENCES:

    - [1] Weisstein, E. (1999). "Heawood Graph - from Wolfram
      MathWorld". [Online] Available:
      http://mathworld.wolfram.com/HeawoodGraph.html [2007, February 17]

    EXAMPLES::

        sage: H = graphs.HeawoodGraph()
        sage: H
        Heawood graph: Graph on 14 vertices
        sage: H.graph6_string()
        'MhEGHC@AI?_PC@_G_'
        sage: (graphs.HeawoodGraph()).show() # long time
    """
    pos_dict = {}
    for i in range(14):
        x = float(cos((pi/2) + (pi/7)*i))
        y = float(sin((pi/2) + (pi/7)*i))
        pos_dict[i] = (x,y)
    import networkx
    G = networkx.heawood_graph()
    return Graph(G, pos=pos_dict, name="Heawood graph")

def HerschelGraph():
    r"""
    Returns the Herschel graph.

    For more information, see this
    `Wikipedia article on the Herschel graph <http://en.wikipedia.org/wiki/Herschel_graph>`_.

    EXAMPLES:

    The Herschel graph is named after Alexander Stewart Herschel. It is
    a planar, bipartite graph with 11 vertices and 18 edges. ::

        sage: G = graphs.HerschelGraph(); G
        Herschel graph: Graph on 11 vertices
        sage: G.is_planar()
        True
        sage: G.is_bipartite()
        True
        sage: G.order()
        11
        sage: G.size()
        18

    The Herschel graph is a perfect graph with radius 3, diameter 4, and
    girth 4. ::

        sage: G.is_perfect()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        4
        sage: G.girth()
        4

    Its chromatic number is 2 and its automorphism group is
    isomorphic to the dihedral group `D_6`. ::

        sage: G.chromatic_number()
        2
        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(6))
        True
    """
    edge_dict = {
        0: [1,3,4],
        1: [2,5,6],
        2: [3,7],
        3: [8,9],
        4: [5,9],
        5: [10],
        6: [7,10],
        7: [8],
        8: [10],
        9: [10]}
    pos_dict = {
        0: [2, 0],
        1: [0, 2],
        2: [-2, 0],
        3: [0, -2],
        4: [1, 0],
        5: [0.5, 0.866025403784439],
        6: [-0.5, 0.866025403784439],
        7: [-1, 0],
        8: [-0.5, -0.866025403784439],
        9: [0.5, -0.866025403784439],
        10: [0, 0]}
    return Graph(edge_dict, pos=pos_dict, name="Herschel graph")

def HigmanSimsGraph(relabel=True):
    r"""
    The Higman-Sims graph is a remarkable strongly regular
    graph of degree 22 on 100 vertices.  For example, it can
    be split into two sets of 50 vertices each, so that each
    half induces a subgraph isomorphic to the
    Hoffman-Singleton graph
    (:meth:`~HoffmanSingletonGraph`).
    This can be done in 352 ways (see [BROUWER-HS-2009]_).

    Its most famous property is that the automorphism
    group has an index 2 subgroup which is one of the
    26 sporadic groups. [HIGMAN1968]_

    The construction used here follows [HAFNER2004]_.

    INPUT:

    - ``relabel`` - default: ``True``.  If ``True`` the
      vertices will be labeled with consecutive integers.
      If ``False`` the labels are strings that are three
      digits long. "xyz" means the vertex is in group
      x (zero through three), pentagon or pentagram y
      (zero through four), and is vertex z (zero
      through four) of that pentagon or pentagram.
      See [HAFNER2004]_ for more.

    OUTPUT:

    The Higman-Sims graph.

    EXAMPLES:

    A split into the first 50 and last 50 vertices
    will induce two copies of the Hoffman-Singleton graph,
    and we illustrate another such split, which is obvious
    based on the construction used. ::

        sage: H = graphs.HigmanSimsGraph()
        sage: A = H.subgraph(range(0,50))
        sage: B = H.subgraph(range(50,100))
        sage: K = graphs.HoffmanSingletonGraph()
        sage: K.is_isomorphic(A) and K.is_isomorphic(B)
        True
        sage: C = H.subgraph(range(25,75))
        sage: D = H.subgraph(range(0,25)+range(75,100))
        sage: K.is_isomorphic(C) and K.is_isomorphic(D)
        True

    The automorphism group contains only one nontrivial
    proper normal subgroup, which is of index 2 and is
    simple.  It is known as the Higman-Sims group.  ::

        sage: H = graphs.HigmanSimsGraph()
        sage: G = H.automorphism_group()
        sage: g=G.order(); g
        88704000
        sage: K = G.normal_subgroups()[1]
        sage: K.is_simple()
        True
        sage: g//K.order()
        2

    REFERENCES:

        .. [BROUWER-HS-2009] `Higman-Sims graph
           <http://www.win.tue.nl/~aeb/graphs/Higman-Sims.html>`_.
           Andries E. Brouwer, accessed 24 October 2009.
        .. [HIGMAN1968] A simple group of order 44,352,000,
           Math.Z. 105 (1968) 110-113. D.G. Higman & C. Sims.
        .. [HAFNER2004] `On the graphs of Hoffman-Singleton and
           Higman-Sims
           <http://www.combinatorics.org/Volume_11/PDF/v11i1r77.pdf>`_.
           The Electronic Journal of Combinatorics 11 (2004), #R77,
           Paul R. Hafner, accessed 24 October 2009.

    AUTHOR:

        - Rob Beezer (2009-10-24)
    """
    HS = Graph()
    HS.name('Higman-Sims graph')

    # Four groups of either five pentagons, or five pentagrams
    # 4 x 5 x 5 = 100 vertices
    # First digit is "group", second is "penta{gon|gram}", third is "vertex"
    vlist = ['%d%d%d'%(g,p,v)
                    for g in range(4) for p in range(5) for v in range(5)]
    for avertex in vlist:
        HS.add_vertex(avertex)

    # Edges: Within groups 0 and 2, joined as pentagons
    # Edges: Within groups 1 and 3, joined as pentagrams
    for g in range(4):
        shift = 1
        if g in [1,3]:
            shift += 1
        for p in range(5):
            for v in range(5):
                HS.add_edge(('%d%d%d'%(g,p,v), '%d%d%d'%(g,p,(v+shift)%5)))

    # Edges: group 0 to group 1
    for x in range(5):
        for m in range(5):
            for c in range(5):
                y = (m*x+c)%5
                HS.add_edge(('0%d%d'%(x,y), '1%d%d'%(m,c)))

    # Edges: group 1 to group 2
    for m in range(5):
        for A in range(5):
            for B in range(5):
                c = (2*(m-A)*(m-A)+B)%5
                HS.add_edge(('1%d%d'%(m,c), '2%d%d'%(A,B)))

    # Edges: group 2 to group 3
    for A in range(5):
        for a in range(5):
            for b in range(5):
                B = (2*A*A+3*a*A-a*a+b)%5
                HS.add_edge(('2%d%d'%(A,B), '3%d%d'%(a,b)))

    # Edges: group 3 to group 0
    for a in range(5):
        for b in range(5):
            for x in range(5):
                y = ((x-a)*(x-a)+b)%5
                HS.add_edge(('3%d%d'%(a,b), '0%d%d'%(x,y)))

    # Edges: group 0 to group 2
    for x in range(5):
        for A in range(5):
            for B in range(5):
                y = (3*x*x+A*x+B+1)%5
                HS.add_edge(('0%d%d'%(x,y), '2%d%d'%(A,B)))
                y = (3*x*x+A*x+B-1)%5
                HS.add_edge(('0%d%d'%(x,y), '2%d%d'%(A,B)))

    # Edges: group 1 to group 3
    for m in range(5):
        for a in range(5):
            for b in range(5):
                c = (m*(m-a)+b+2)%5
                HS.add_edge(('1%d%d'%(m,c), '3%d%d'%(a,b)))
                c = (m*(m-a)+b-2)%5
                HS.add_edge(('1%d%d'%(m,c), '3%d%d'%(a,b)))

    # Rename to integer vertex labels, creating dictionary
    # Or not, and create identity mapping
    if relabel:
        vmap = HS.relabel(return_map=True)
    else:
        vmap={}
        for v in vlist:
            vmap[v] = v
    # Layout vertices in a circle
    # In the order given in vlist
    # Using labels from vmap
    pos_dict = {}
    for i in range(100):
        x = float(cos((pi/2) + ((2*pi)/100)*i))
        y = float(sin((pi/2) + ((2*pi)/100)*i))
        pos_dict[vmap[vlist[i]]] = (x,y)
    HS.set_pos(pos_dict)
    return HS

def HoffmanSingletonGraph():
    r"""
    Returns the Hoffman-Singleton graph.

    The Hoffman-Singleton graph is the Moore graph of degree 7,
    diameter 2 and girth 5. The Hoffman-Singleton theorem states that
    any Moore graph with girth 5 must have degree 2, 3, 7 or 57. The
    first three respectively are the pentagon, the Petersen graph, and
    the Hoffman-Singleton graph. The existence of a Moore graph with
    girth 5 and degree 57 is still open.

    A Moore graph is a graph with diameter `d` and girth
    `2d + 1`. This implies that the graph is regular, and
    distance regular.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. A novel algorithm written by
    Tom Boothby gives a random layout which is pleasing to the eye.

    REFERENCES:

    .. [GodsilRoyle] Godsil, C. and Royle, G. Algebraic Graph Theory.
      Springer, 2001.

    EXAMPLES::

        sage: HS = graphs.HoffmanSingletonGraph()
        sage: Set(HS.degree())
        {7}
        sage: HS.girth()
        5
        sage: HS.diameter()
        2
        sage: HS.num_verts()
        50

    Note that you get a different layout each time you create the graph.
    ::

        sage: HS.layout()[1]
        (-0.844..., 0.535...)
        sage: graphs.HoffmanSingletonGraph().layout()[1]
        (-0.904..., 0.425...)

    """
    H = Graph({ \
    'q00':['q01'], 'q01':['q02'], 'q02':['q03'], 'q03':['q04'], 'q04':['q00'], \
    'q10':['q11'], 'q11':['q12'], 'q12':['q13'], 'q13':['q14'], 'q14':['q10'], \
    'q20':['q21'], 'q21':['q22'], 'q22':['q23'], 'q23':['q24'], 'q24':['q20'], \
    'q30':['q31'], 'q31':['q32'], 'q32':['q33'], 'q33':['q34'], 'q34':['q30'], \
    'q40':['q41'], 'q41':['q42'], 'q42':['q43'], 'q43':['q44'], 'q44':['q40'], \
    'p00':['p02'], 'p02':['p04'], 'p04':['p01'], 'p01':['p03'], 'p03':['p00'], \
    'p10':['p12'], 'p12':['p14'], 'p14':['p11'], 'p11':['p13'], 'p13':['p10'], \
    'p20':['p22'], 'p22':['p24'], 'p24':['p21'], 'p21':['p23'], 'p23':['p20'], \
    'p30':['p32'], 'p32':['p34'], 'p34':['p31'], 'p31':['p33'], 'p33':['p30'], \
    'p40':['p42'], 'p42':['p44'], 'p44':['p41'], 'p41':['p43'], 'p43':['p40']})
    for j in range(5):
        for i in range(5):
            for k in range(5):
                con = (i+j*k)%5
                H.add_edge(('q%d%d'%(k,con),'p%d%d'%(j,i)))
    H.name('Hoffman-Singleton graph')
    from sage.combinat.permutation import Permutations
    from sage.misc.prandom import randint
    P = Permutations([1,2,3,4])
    qpp = [0] + list(P[randint(0,23)])
    ppp = [0] + list(P[randint(0,23)])
    qcycle = lambda i,s : ['q%s%s'%(i,(j+s)%5) for j in qpp]
    pcycle = lambda i,s : ['p%s%s'%(i,(j+s)%5) for j in ppp]
    l = 0
    s = 0
    D = []
    while l < 5:
        for q in qcycle(l,s):
            D.append(q)
        vv = 'p%s'%q[1]
        s = int([v[-1] for v in H.neighbors(q) if v[:2] == vv][0])
        for p in pcycle(l,s):
            D.append(p)
        vv = 'q%s'%(int(p[1])+1)
        v = [v[-1] for v in H.neighbors(p) if v[:2] == vv]
        if len(v):
            s = int(v[0])
        l+=1
    map = H.relabel(return_map=True)
    pos_dict = {}
    for i in range(50):
        x = float(cos((pi/2) + ((2*pi)/50)*i))
        y = float(sin((pi/2) + ((2*pi)/50)*i))
        pos_dict[map[D[i]]] = (x,y)
    H.set_pos(pos_dict)
    return H

def HoffmanGraph():
    r"""
    Returns the Hoffman Graph.

    See the :wikipedia:`Wikipedia page on the Hoffman graph
    <Hoffman_graph>`.

    EXAMPLES::

        sage: g = graphs.HoffmanGraph()
        sage: g.is_bipartite()
        True
        sage: g.is_hamiltonian() # long time
        True
        sage: g.radius()
        3
        sage: g.diameter()
        4
        sage: g.automorphism_group().cardinality()
        48
    """
    g = Graph({
            0: [1, 7, 8, 13],
            1: [2, 9, 14],
            2: [3, 8, 10],
            3: [4, 9, 15],
            4: [5, 10, 11],
            5: [6, 12, 14],
            6: [7, 11, 13],
            7: [12, 15],
            8: [12, 14],
            9: [11, 13],
            10: [12, 15],
            11: [14],
            13: [15]})
    g.set_pos({})
    _circle_embedding(g, range(8))
    _circle_embedding(g, range(8, 14), radius=.7, shift=.5)
    _circle_embedding(g, [14, 15], radius=.1)

    g.name("Hoffman Graph")

    return g

def HoltGraph():
    r"""
    Returns the Holt graph (also called the Doyle graph)

    See the :wikipedia:`Wikipedia page on the Holt graph
    <Holt_graph>`.

    EXAMPLES::

        sage: g = graphs.HoltGraph();g
        Holt graph: Graph on 27 vertices
        sage: g.is_regular()
        True
        sage: g.is_vertex_transitive()
        True
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian() # long time
        True
        sage: g.radius()
        3
        sage: g.diameter()
        3
        sage: g.girth()
        5
        sage: g.automorphism_group().cardinality()
        54
    """
    g = Graph(loops=False, name = "Holt graph", pos={})
    for x in range(9):
        for y in range(3):
            g.add_edge((x,y),((4*x+1)%9,(y-1)%3))
            g.add_edge((x,y),((4*x-1)%9,(y-1)%3))
            g.add_edge((x,y),((7*x+7)%9,(y+1)%3))
            g.add_edge((x,y),((7*x-7)%9,(y+1)%3))

    for j in range(0,6,2):
        _line_embedding(g, [(x,j/2) for x in range(9)],
                        first=(cos(2*j*pi/6),sin(2*j*pi/6)),
                        last=(cos(2*(j+1)*pi/6),sin(2*(j+1)*pi/6)))

    return g

def KrackhardtKiteGraph():
    """
    Returns a Krackhardt kite graph with 10 nodes.

    The Krackhardt kite graph was originally developed by David
    Krackhardt for the purpose of studying social networks. It is used
    to show the distinction between: degree centrality, betweeness
    centrality, and closeness centrality. For more information read the
    plotting section below in conjunction with the example.

    REFERENCES:

    - [1] Kreps, V. (2002). "Social Network Analysis".  [Online] Available:
      http://www.orgnet.com/sna.html

    This constructor depends on NetworkX numeric labeling.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the graph is
    drawn left to right, in top to bottom row sequence of [2, 3, 2, 1,
    1, 1] nodes on each row. This places the fourth node (3) in the
    center of the kite, with the highest degree. But the fourth node
    only connects nodes that are otherwise connected, or those in its
    clique (i.e.: Degree Centrality). The eighth (7) node is where the
    kite meets the tail. It has degree = 3, less than the average, but
    is the only connection between the kite and tail (i.e.: Betweenness
    Centrality). The sixth and seventh nodes (5 and 6) are drawn in the
    third row and have degree = 5. These nodes have the shortest path
    to all other nodes in the graph (i.e.: Closeness Centrality).
    Please execute the example for visualization.

    EXAMPLE: Construct and show a Krackhardt kite graph

    ::

        sage: g = graphs.KrackhardtKiteGraph()
        sage: g.show() # long time
    """
    pos_dict = {0:(-1,4),1:(1,4),2:(-2,3),3:(0,3),4:(2,3),5:(-1,2),6:(1,2),7:(0,1),8:(0,0),9:(0,-1)}
    import networkx
    G = networkx.krackhardt_kite_graph()
    return Graph(G, pos=pos_dict, name="Krackhardt Kite Graph")

def LjubljanaGraph(embedding=1):
    r"""
    Returns the Ljubljana Graph.

    The Ljubljana graph is a bipartite 3-regular graph on 112
    vertices and 168 edges. It is not vertex-transitive as it has
    two orbits which are also independent sets of size 56. See the
    :wikipedia:`Wikipedia page on the Ljubljana Graph
    <Ljubljana_graph>`.

    The default embedding is obtained from the Heawood graph.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.LjubljanaGraph()
        sage: g.order()
        112
        sage: g.size()
        168
        sage: g.girth()
        10
        sage: g.diameter()
        8
        sage: g.show(figsize=[10, 10])   # long time
        sage: graphs.LjubljanaGraph(embedding=2).show(figsize=[10, 10])   # long time

    TESTS::

        sage: graphs.LjubljanaGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    L = [47, -23, -31, 39, 25, -21, -31, -41, 25, 15, 29, -41, -19, 15,
         -49, 33, 39, -35, -21, 17, -33, 49, 41, 31, -15, -29, 41, 31,
         -15, -25, 21, 31, -51, -25, 23, 9, -17, 51, 35, -29, 21, -51,
         -39, 33, -9, -51, 51, -47, -33, 19, 51, -21, 29, 21, -31, -39]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(112, L, 2)
    g.name("Ljubljana graph")

    if embedding == 1:
        dh = HeawoodGraph().get_pos()

        # Correspondence between the vertices of the Heawood Graph and
        # 8-sets of the Ljubljana Graph.

        d = {
            0: [1, 21, 39, 57, 51, 77, 95, 107],
            1: [2, 22, 38, 58, 50, 78, 94, 106],
            2: [3, 23, 37, 59, 49, 79, 93, 105],
            3: [4, 24, 36, 60, 48, 80, 92, 104],
            4: [5, 25, 35, 61, 15, 81, 91, 71],
            9: [6, 26, 44, 62, 16, 82, 100, 72],
            10: [7, 27, 45, 63, 17, 83, 101, 73],
            11: [8, 28, 46, 64, 18, 84, 102, 74],
            12: [9, 29, 47, 65, 19, 85, 103, 75],
            13: [10, 30, 0, 66, 20, 86, 56, 76],
            8: [11, 31, 111, 67, 99, 87, 55, 43],
            7: [12, 32, 110, 68, 98, 88, 54, 42],
            6: [13, 33, 109, 69, 97, 89, 53, 41],
            5: [14, 34, 108, 70, 96, 90, 52, 40]
            }

        # The vertices of each 8-set are plotted on a circle, and the
        # circles are slowly shifted to obtain a symmetric drawing.

        for i, (u, vertices) in enumerate(d.iteritems()):
            _circle_embedding(g, vertices, center=dh[u], radius=.1,
                    shift=8.*i/14)

        return g

    elif embedding == 2:
        return g

    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def M22Graph():
    r"""
    Returns the M22 graph.

    The `M_{22}` graph is the unique strongly regular graph with parameters
    `v = 77, k = 16, \lambda = 0, \mu = 4`.

    For more information on the `M_{22}` graph, see
    `<http://www.win.tue.nl/~aeb/graphs/M22.html>`_.

    EXAMPLES::

        sage: g = graphs.M22Graph()
        sage: g.order()
        77
        sage: g.size()
        616
        sage: g.is_strongly_regular(parameters = True)
        (77, 16, 0, 4)
    """
    from sage.groups.perm_gps.permgroup_named import MathieuGroup
    sets = map(tuple,MathieuGroup(22).orbit((1,2,3,7,10,20), action = "OnSets"))
    g = Graph([sets, lambda x,y : not any(xx in y for xx in x)], name="M22 Graph")
    g.relabel()
    ordering = [0, 1, 3, 4, 5, 6, 7, 10, 12, 19, 20, 31, 2, 24, 35, 34, 22, 32,
                36, 23, 27, 25, 40, 26, 16, 71, 61, 63, 50, 68, 39, 52, 48, 44,
                69, 28, 9, 64, 60, 17, 38, 49, 45, 65, 14, 70, 72, 21, 43, 56,
                33, 73, 58, 55, 41, 29, 66, 54, 76, 46, 67, 11, 51, 47, 62, 53,
                15, 8, 18, 13, 59, 37, 30, 57, 75, 74, 42]

    _circle_embedding(g, ordering)

    return g

def MarkstroemGraph():
    r"""
    Returns the Markström Graph.

    The Markström Graph is a cubic planar graph with no cycles of length 4 nor
    8, but containing cycles of length 16. For more information on the Markström
    Graph, see the `corresponding Wolfram page
    <http://mathworld.wolfram.com/MarkstroemGraph.html>`_.

    EXAMPLES::

        sage: g = graphs.MarkstroemGraph()
        sage: g.order()
        24
        sage: g.size()
        36
        sage: g.is_planar()
        True
        sage: g.is_regular(3)
        True
        sage: g.subgraph_search(graphs.CycleGraph(4)) is None
        True
        sage: g.subgraph_search(graphs.CycleGraph(8)) is None
        True
        sage: g.subgraph_search(graphs.CycleGraph(16))
        Subgraph of (Markstroem Graph): Graph on 16 vertices
    """
    g = Graph(name="Markstroem Graph")

    g.add_cycle(range(9))
    g.add_path([0,9,10,11,2,1,11])
    g.add_path([3,12,13,14,5,4,14])
    g.add_path([6,15,16,17,8,7,17])
    g.add_cycle([10,9,18])
    g.add_cycle([12,13,19])
    g.add_cycle([15,16,20])
    g.add_cycle([21,22,23])
    g.add_edges([(19,22),(18,21),(20,23)])

    _circle_embedding(g, sum([[9+3*i+j for j in range(3)]+[0]*2 for i in range(3)],[]), radius=.6, shift=.7)
    _circle_embedding(g, [18,19,20], radius=.35, shift=.25)
    _circle_embedding(g, [21,22,23], radius=.15, shift=.25)
    _circle_embedding(g, range(9))

    return g

def McGeeGraph(embedding=2):
    r"""
    Returns the McGee Graph.

    See the :wikipedia:`Wikipedia page on the McGee Graph
    <McGee_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.McGeeGraph()
        sage: g.order()
        24
        sage: g.size()
        36
        sage: g.girth()
        7
        sage: g.diameter()
        4
        sage: g.show()
        sage: graphs.McGeeGraph(embedding=1).show()

    TESTS::

        sage: graphs.McGeeGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    L = [47, -23, -31, 39, 25, -21, -31, -41, 25, 15, 29, -41, -19, 15,
         -49, 33, 39, -35, -21, 17, -33, 49, 41, 31, -15, -29, 41, 31,
         -15, -25, 21, 31, -51, -25, 23, 9, -17, 51, 35, -29, 21, -51,
         -39, 33, -9, -51, 51, -47, -33, 19, 51, -21, 29, 21, -31, -39]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(24, [12, 7, -7], 8)
    g.name('McGee graph')

    if embedding == 1:
        return g

    elif embedding == 2:

        o = [[7, 2, 13, 8, 19, 14, 1, 20],
             [5, 4, 11, 10, 17, 16, 23, 22],
             [3, 12, 9, 18, 15, 0, 21, 6]]

        _circle_embedding(g, o[0], radius=1.5)
        _circle_embedding(g, o[1], radius=3, shift=-.5)
        _circle_embedding(g, o[2], radius=2.25, shift=.5)

        return g

    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def McLaughlinGraph():
    r"""
    Returns the McLaughlin Graph.

    The McLaughlin Graph is the unique strongly regular graph of parameters
    `(275, 112, 30, 56)`.

    For more information on the McLaughlin Graph, see its web page on `Andries
    Brouwer's website <http://www.win.tue.nl/~aeb/graphs/McL.html>`_ which gives
    the definition that this method implements.

    .. NOTE::

        To create this graph you must have the gap_packages spkg installed.

    EXAMPLES::

        sage: g = graphs.McLaughlinGraph()           # optional gap_packages
        sage: g.is_strongly_regular(parameters=True) # optional gap_packages
        (275, 112, 30, 56)
        sage: set(g.spectrum()) == {112, 2, -28}     # optional gap_packages
        True
    """
    from sage.combinat.designs.block_design import WittDesign
    from itertools import combinations
    from sage.sets.set import Set

    blocks = WittDesign(23).blocks()
    blocks = map(Set, blocks)
    B = [b for b in blocks if 0 in b]
    C = [b for b in blocks if not 0 in b]
    g = graph.Graph()
    for b in B:
        for x in range(23):
            if not x in b:
                g.add_edge(b, x)

    for b in C:
        for x in b:
            g.add_edge(b, x)

    for b, bb in combinations(B, 2):
        if len(b & bb) == 1:
            g.add_edge(b, bb)

    for c, cc in combinations(C, 2):
        if len(c & cc) == 1:
            g.add_edge(c, cc)

    for b in B:
        for c in C:
            if len(b & c) == 3:
                g.add_edge(b, c)

    g.relabel()
    g.name("McLaughlin")
    return g

def MoebiusKantorGraph():
    """
    Returns a Moebius-Kantor Graph.

    A Moebius-Kantor graph is a cubic symmetric graph. (See also the
    Heawood graph). It has 16 nodes and 24 edges. It is nonplanar and
    Hamiltonian. It has diameter = 4, girth = 6, and chromatic number =
    2. It is identical to the Generalized Petersen graph, P[8,3].

    PLOTTING: See the plotting section for the generalized Petersen graphs.

    REFERENCES:

    - [1] Weisstein, E. (1999). "Moebius-Kantor Graph - from
      Wolfram MathWorld". [Online] Available:
      http://mathworld.wolfram.com/Moebius-KantorGraph.html [2007,
      February 17]

    EXAMPLES::

        sage: MK = graphs.MoebiusKantorGraph()
        sage: MK
        Moebius-Kantor Graph: Graph on 16 vertices
        sage: MK.graph6_string()
        'OhCGKE?O@?ACAC@I?Q_AS'
        sage: (graphs.MoebiusKantorGraph()).show() # long time
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    G=GeneralizedPetersenGraph(8,3)
    G.name("Moebius-Kantor Graph")
    return G

def MoserSpindle():
    r"""
    Returns the Moser spindle.

    For more information, see this
    `MathWorld article on the Moser spindle <http://mathworld.wolfram.com/MoserSpindle.html>`_.

    EXAMPLES:

    The Moser spindle is a planar graph having 7 vertices and 11 edges. ::

        sage: G = graphs.MoserSpindle(); G
        Moser spindle: Graph on 7 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        7
        sage: G.size()
        11

    It is a Hamiltonian graph with radius 2, diameter 2, and girth 3. ::

        sage: G.is_hamiltonian()
        True
        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        3

    The Moser spindle has chromatic number 4 and its automorphism
    group is isomorphic to the dihedral group `D_4`. ::

        sage: G.chromatic_number()
        4
        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(4))
        True
    """
    edge_dict = {
        0: [1,4,5,6],
        1: [2,5],
        2: [3,5],
        3: [4,6],
        4: [6]}
    pos_dict = {
        0: [0, 2],
        1: [-1.90211303259031, 0.618033988749895],
        2: [-1.17557050458495, -1.61803398874989],
        3: [1.17557050458495, -1.61803398874989],
        4: [1.90211303259031, 0.618033988749895],
        5: [1, 0],
        6: [-1, 0]}
    return Graph(edge_dict, pos=pos_dict, name="Moser spindle")


def NauruGraph(embedding=2):
    """
    Returns the Nauru Graph.

    See the :wikipedia:`Wikipedia page on the Nauru Graph
    <Nauru_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.NauruGraph()
        sage: g.order()
        24
        sage: g.size()
        36
        sage: g.girth()
        6
        sage: g.diameter()
        4
        sage: g.show()
        sage: graphs.NauruGraph(embedding=1).show()

    TESTS::

        sage: graphs.NauruGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
        sage: graphs.NauruGraph(embedding=1).is_isomorphic(g)
        True
    """

    if embedding == 1:
        from sage.graphs.generators.families import LCFGraph
        g = LCFGraph(24, [5, -9, 7, -7, 9, -5], 4)
        g.name('Nauru Graph')
        return g
    elif embedding == 2:
        from sage.graphs.generators.families import GeneralizedPetersenGraph
        g = GeneralizedPetersenGraph(12, 5)
        g.name("Nauru Graph")
        return g
    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def PappusGraph():
    """
    Returns the Pappus graph, a graph on 18 vertices.

    The Pappus graph is cubic, symmetric, and distance-regular.

    EXAMPLES::

        sage: G = graphs.PappusGraph()
        sage: G.show()  # long time
        sage: L = graphs.LCFGraph(18, [5,7,-7,7,-7,-5], 3)
        sage: L.show()  # long time
        sage: G.is_isomorphic(L)
        True
    """
    pos_dict = {}
    for i in range(6):
        pos_dict[i] = [float(cos(pi/2 + ((2*pi)/6)*i)),\
                       float(sin(pi/2 + ((2*pi)/6)*i))]
        pos_dict[6 + i] = [(2/3.0)*float(cos(pi/2 + ((2*pi)/6)*i)),\
                           (2/3.0)*float(sin(pi/2 + ((2*pi)/6)*i))]
        pos_dict[12 + i] = [(1/3.0)*float(cos(pi/2 + ((2*pi)/6)*i)),\
                            (1/3.0)*float(sin(pi/2 + ((2*pi)/6)*i))]
    return Graph({0:[1,5,6],1:[2,7],2:[3,8],3:[4,9],4:[5,10],\
                        5:[11],6:[13,17],7:[12,14],8:[13,15],9:[14,16],\
                        10:[15,17],11:[12,16],12:[15],13:[16],14:[17]},\
                       pos=pos_dict, name="Pappus Graph")

def PoussinGraph():
    r"""
    Returns the Poussin Graph.

    For more information on the Poussin Graph, see its corresponding `Wolfram
    page <http://mathworld.wolfram.com/PoussinGraph.html>`_.

    EXAMPLES::

        sage: g = graphs.PoussinGraph()
        sage: g.order()
        15
        sage: g.is_planar()
        True
    """
    g = Graph({2:[7,8,3,4],1:[7,6],0:[6,5,4],3:[5]},name="Poussin Graph")

    g.add_cycle(range(3))
    g.add_cycle(range(3,9))
    g.add_cycle(range(9,14))
    g.add_path([8,12,7,11,6,10,5,9,3,13,8,12])
    g.add_edges([(14,i) for i in range(9,14)])
    _circle_embedding(g, range(3), shift=.75)
    _circle_embedding(g, range(3,9), radius=.4, shift=0)
    _circle_embedding(g, range(9,14), radius=.2, shift=.4)
    g.get_pos()[14] = (0,0)

    return g

def PetersenGraph():
    """
    The Petersen Graph is a named graph that consists of 10 vertices
    and 15 edges, usually drawn as a five-point star embedded in a
    pentagon.

    The Petersen Graph is a common counterexample. For example, it is
    not Hamiltonian.

    PLOTTING: See the plotting section for the generalized Petersen graphs.

    EXAMPLES: We compare below the Petersen graph with the default
    spring-layout versus a planned position dictionary of [x,y]
    tuples::

        sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
        sage: petersen_spring.show() # long time
        sage: petersen_database = graphs.PetersenGraph()
        sage: petersen_database.show() # long time
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    P=GeneralizedPetersenGraph(5,2)
    P.name("Petersen graph")
    return P


def RobertsonGraph():
    """
    Returns the Robertson graph.

    See the :wikipedia:`Wikipedia page on the Robertson Graph
    <Robertson_graph>`.

    EXAMPLE::

        sage: g = graphs.RobertsonGraph()
        sage: g.order()
        19
        sage: g.size()
        38
        sage: g.diameter()
        3
        sage: g.girth()
        5
        sage: g.charpoly().factor()
        (x - 4) * (x - 1)^2 * (x^2 + x - 5) * (x^2 + x - 1) * (x^2 - 3)^2 * (x^2 + x - 4)^2 * (x^2 + x - 3)^2
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian()
        True
        sage: g.is_vertex_transitive()
        False
    """
    from sage.graphs.generators.families import LCFGraph
    lcf = [8, 4, 7, 4, 8, 5, 7, 4, 7, 8, 4, 5, 7, 8, 4, 8, 4, 8, 4]
    g = LCFGraph(19, lcf, 1)
    g.name("Robertson Graph")
    return g


def SchlaefliGraph():
    r"""
    Returns the Schläfli graph.

    The Schläfli graph is the only strongly regular graphs of parameters
    `(27,16,10,8)` (see [GodsilRoyle]_).

    For more information, see the :wikipedia:`Wikipedia article on the
    Schläfli graph <Schläfli_graph>`.

    .. SEEALSO::

        :meth:`Graph.is_strongly_regular` -- tests whether a graph is strongly
        regular and/or returns its parameters.

    .. TODO::

        Find a beautiful layout for this beautiful graph.

    EXAMPLE:

    Checking that the method actually returns the Schläfli graph::

        sage: S = graphs.SchlaefliGraph()
        sage: S.is_strongly_regular(parameters = True)
        (27, 16, 10, 8)

    The graph is vertex-transitive::

        sage: S.is_vertex_transitive()
        True

    The neighborhood of each vertex is isomorphic to the complement of the
    Clebsch graph::

        sage: neighborhood = S.subgraph(vertices = S.neighbors(0))
        sage: graphs.ClebschGraph().complement().is_isomorphic(neighborhood)
        True
    """
    from sage.graphs.graph import Graph
    G = Graph('ZBXzr|}^z~TTitjLth|dmkrmsl|if}TmbJMhrJX]YfFyTbmsseztKTvyhDvw')
    order = [1,8,5,10,2,6,11,15,17,13,18,12,9,24,25,3,26,7,16,20,23,0,21,14,22,4,19]
    _circle_embedding(G, order)
    G.name("Schläfli graph")
    return G

def ShrikhandeGraph():
    """
    Returns the Shrikhande graph.

    For more information, see the `MathWorld article on the Shrikhande graph
    <http://mathworld.wolfram.com/ShrikhandeGraph.html>`_ or the
    :wikipedia:`Wikipedia article on the Shrikhande graph <Shrikhande_graph>`.

    .. SEEALSO::

        :meth:`Graph.is_strongly_regular` -- tests whether a graph is strongly
        regular and/or returns its parameters.

    EXAMPLES:

    The Shrikhande graph was defined by S. S. Shrikhande in 1959. It has
    `16` vertices and `48` edges, and is strongly regular of degree `6` with
    parameters `(2,2)`::

        sage: G = graphs.ShrikhandeGraph(); G
        Shrikhande graph: Graph on 16 vertices
        sage: G.order()
        16
        sage: G.size()
        48
        sage: G.is_regular(6)
        True
        sage: set([ len([x for x in G.neighbors(i) if x in G.neighbors(j)])
        ....:     for i in range(G.order())
        ....:     for j in range(i) ])
        set([2])

    It is non-planar, and both Hamiltonian and Eulerian::

        sage: G.is_planar()
        False
        sage: G.is_hamiltonian()
        True
        sage: G.is_eulerian()
        True

    It has radius `2`, diameter `2`, and girth `3`::

        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        3

    Its chromatic number is `4` and its automorphism group is of order
    `192`::

        sage: G.chromatic_number()
        4
        sage: G.automorphism_group().cardinality()
        192

    It is an integral graph since it has only integral eigenvalues::

        sage: G.characteristic_polynomial().factor()
        (x - 6) * (x - 2)^6 * (x + 2)^9

    It is a toroidal graph, and its embedding on a torus is dual to an
    embedding of the Dyck graph (:meth:`DyckGraph <GraphGenerators.DyckGraph>`).
    """
    pos_dict = {}
    for i in range(8):
        pos_dict[i] = [float(cos((2*i) * pi/8)),
                       float(sin((2*i) * pi/8))]
        pos_dict[8 + i] = [0.5 * pos_dict[i][0],
                           0.5 * pos_dict[i][1]]
    edge_dict = {
        0O00: [0O06, 0O07, 0O01, 0O02,   0O11, 0O17],
        0O01: [0O07, 0O00, 0O02, 0O03,   0O12, 0O10],
        0O02: [0O00, 0O01, 0O03, 0O04,   0O13, 0O11],
        0O03: [0O01, 0O02, 0O04, 0O05,   0O14, 0O12],
        0O04: [0O02, 0O03, 0O05, 0O06,   0O15, 0O13],
        0O05: [0O03, 0O04, 0O06, 0O07,   0O16, 0O14],
        0O06: [0O04, 0O05, 0O07, 0O00,   0O17, 0O15],
        0O07: [0O05, 0O06, 0O00, 0O01,   0O10, 0O16],

        0O10: [0O12, 0O13, 0O15, 0O16,   0O07, 0O01],
        0O11: [0O13, 0O14, 0O16, 0O17,   0O00, 0O02],
        0O12: [0O14, 0O15, 0O17, 0O10,   0O01, 0O03],
        0O13: [0O15, 0O16, 0O10, 0O11,   0O02, 0O04],
        0O14: [0O16, 0O17, 0O11, 0O12,   0O03, 0O05],
        0O15: [0O17, 0O10, 0O12, 0O13,   0O04, 0O06],
        0O16: [0O10, 0O11, 0O13, 0O14,   0O05, 0O07],
        0O17: [0O11, 0O12, 0O14, 0O15,   0O06, 0O00]
    }

    return Graph(edge_dict, pos=pos_dict, name="Shrikhande graph")

def SylvesterGraph():
    """
    Returns the Sylvester Graph.

    This graph is obtained from the Hoffman Singleton graph by considering the
    graph induced by the vertices at distance two from the vertices of an (any)
    edge.

    For more information on the Sylvester graph, see
    `<http://www.win.tue.nl/~aeb/graphs/Sylvester.html>`_.

    .. SEEALSO::

        * :meth:`~sage.graphs.graph_generators.GraphGenerators.HoffmanSingletonGraph`.

    EXAMPLE::

        sage: g = graphs.SylvesterGraph(); g
        Sylvester Graph: Graph on 36 vertices
        sage: g.order()
        36
        sage: g.size()
        90
        sage: g.is_regular(k=5)
        True
    """
    g = HoffmanSingletonGraph()
    e = g.edge_iterator(labels = False).next()
    g.delete_vertices(g.neighbors(e[0]) + g.neighbors(e[1]))
    g.relabel()
    ordering = [0, 1, 2, 4, 5, 9, 16, 35, 15, 18, 20, 30, 22, 6, 33, 32, 14,
                10, 28, 29, 7, 24, 23, 26, 19, 12, 13, 21, 11, 31, 3, 27, 25,
                17, 8, 34]
    _circle_embedding(g,ordering, shift=.5)
    g.name("Sylvester Graph")
    return g

def SimsGewirtzGraph():
    """
    Returns the Sims-Gewirtz Graph.

    This graph is obtained from the Higman Sims graph by considering the graph
    induced by the vertices at distance two from the vertices of an (any)
    edge. It is the only strongly regular graph with parameters `v = 56, k = 10,
    \lambda = 0, \mu = 2`

    For more information on the Sylvester graph, see
    `<http://www.win.tue.nl/~aeb/graphs/Sims-Gewirtz.html>`_ or its
    :wikipedia:`Wikipedia page <Gewirtz graph>`.

    .. SEEALSO::

        * :meth:`~sage.graphs.graph_generators.GraphGenerators.HigmanSimsGraph`.

    EXAMPLE::

        sage: g = graphs.SimsGewirtzGraph(); g
        Sims-Gewirtz Graph: Graph on 56 vertices
        sage: g.order()
        56
        sage: g.size()
        280
        sage: g.is_strongly_regular(parameters = True)
        (56, 10, 0, 2)

    """
    g = HigmanSimsGraph()
    e = g.edge_iterator(labels = False).next()
    g.delete_vertices(g.neighbors(e[0]) + g.neighbors(e[1]))
    g.relabel()
    ordering = [0, 2, 3, 4, 6, 7, 8, 17, 1, 41, 49, 5, 22, 26, 11, 27, 15, 47,
                53, 52, 38, 43, 44, 18, 20, 32, 19, 42, 54, 36, 51, 30, 33, 35,
                37, 28, 34, 12, 29, 23, 55, 25, 40, 24, 9, 14, 48, 39, 45, 16,
                13, 21, 31, 50, 10, 46]
    _circle_embedding(g,ordering)
    g.name("Sims-Gewirtz Graph")
    return g

def SousselierGraph():
    r"""
    Returns the Sousselier Graph.

    The Sousselier graph is a hypohamiltonian graph on 16 vertices and 27
    edges. For more information, see the corresponding `Wikipedia page (in
    French) <http://fr.wikipedia.org/wiki/Graphe_de_Sousselier>`_.

    EXAMPLES::

        sage: g = graphs.SousselierGraph()
        sage: g.order()
        16
        sage: g.size()
        27
        sage: g.radius()
        2
        sage: g.diameter()
        3
        sage: g.automorphism_group().cardinality()
        2
        sage: g.is_hamiltonian()
        False
        sage: g.delete_vertex(g.random_vertex())
        sage: g.is_hamiltonian()
        True
    """
    g = Graph(name="Sousselier Graph")

    g.add_cycle(range(15))
    g.add_path([12,8,3,14])
    g.add_path([9,5,0,11])
    g.add_edge(6,2)
    g.add_edges([(15,i) for i in range(15) if i%3==1])

    _circle_embedding(g, range(15), shift=-.25)
    g.get_pos()[15] = (0,0)

    return g

def SzekeresSnarkGraph():
    r"""
    Returns the Szekeres Snark Graph.

    The Szekeres graph is a snark with 50 vertices and 75 edges. For more
    information on this graph, see the :wikipedia:`Szekeres_snark`.

    EXAMPLES::

        sage: g = graphs.SzekeresSnarkGraph()
        sage: g.order()
        50
        sage: g.size()
        75
        sage: g.chromatic_number()
        3
    """
    g = Graph(name="Szekeres Snark Graph")

    for i in range(5):
        g.add_cycle([(i,j) for j in range(9)])
        g.delete_edge((i,0),(i,8))
        g.add_edge((i,1),i)
        g.add_edge((i,4),i)
        g.add_edge((i,7),i)
        g.add_edge((i,0),(i,5))
        g.add_edge((i,8),(i,3))

        g.add_edge((i,0),((i+1)%5,8))
        g.add_edge((i,6),((i+2)%5,2))
        _circle_embedding(g, [(i,j) for j in range(9)],
                          radius=.3,
                          center=(cos(2*(i+.25)*pi/5),sin(2*(i+.25)*pi/5)),
                          shift=5.45+1.8*i)

    _circle_embedding(g, range(5), radius=1, shift=.25)

    g.relabel()
    return g


def ThomsenGraph():
    """
    Returns the Thomsen Graph.

    The Thomsen Graph is actually a complete bipartite graph with `(n1, n2) =
    (3, 3)`. It is also called the Utility graph.

    PLOTTING: See CompleteBipartiteGraph.

    EXAMPLES::

        sage: T = graphs.ThomsenGraph()
        sage: T
        Thomsen graph: Graph on 6 vertices
        sage: T.graph6_string()
        'EFz_'
        sage: (graphs.ThomsenGraph()).show() # long time
    """
    pos_dict = {0:(-1,1),1:(0,1),2:(1,1),3:(-1,0),4:(0,0),5:(1,0)}
    import networkx
    G = networkx.complete_bipartite_graph(3,3)
    return Graph(G, pos=pos_dict, name="Thomsen graph")

def TietzeGraph():
    r"""
    Returns the Tietze Graph.

    For more information on the Tietze Graph, see the
    :wikipedia:`Tietze's_graph`.

    EXAMPLES::

        sage: g = graphs.TietzeGraph()
        sage: g.order()
        12
        sage: g.size()
        18
        sage: g.diameter()
        3
        sage: g.girth()
        3
        sage: g.automorphism_group().cardinality()
        12
        sage: g.automorphism_group().is_isomorphic(groups.permutation.Dihedral(6))
        True
    """
    g = Graph([(0,9),(3,10),(6,11),(1,5),(2,7),(4,8),(7,2)], name="Tietze Graph")
    g.add_cycle(range(9))
    g.add_cycle([9,10,11])
    _circle_embedding(g,range(9))
    _circle_embedding(g,[9,10,11],radius=.5)

    return g

def Tutte12Cage():
    r"""
    Returns Tutte's 12-Cage.

    See the :wikipedia:`Wikipedia page on the Tutte 12-Cage
    <Tutte_12-cage>`.

    EXAMPLES::

        sage: g = graphs.Tutte12Cage()
        sage: g.order()
        126
        sage: g.size()
        189
        sage: g.girth()
        12
        sage: g.diameter()
        6
        sage: g.show()
    """
    L = [17, 27, -13, -59, -35, 35, -11, 13, -53, 53, -27, 21, 57, 11,
         -21, -57, 59, -17]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(126, L, 7)
    g.name("Tutte 12-Cage")
    return g

def TutteCoxeterGraph(embedding=2):
    r"""
    Returns the Tutte-Coxeter graph.

    See the :wikipedia:`Wikipedia page on the Tutte-Coxeter Graph
    <Tutte-Coxeter_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.TutteCoxeterGraph()
        sage: g.order()
        30
        sage: g.size()
        45
        sage: g.girth()
        8
        sage: g.diameter()
        4
        sage: g.show()
        sage: graphs.TutteCoxeterGraph(embedding=1).show()

    TESTS::

        sage: graphs.TutteCoxeterGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(30, [-13, -9, 7, -7, 9, 13], 5)
    g.name("Tutte-Coxeter graph")

    if embedding == 1:
        d = {
            0: [1, 3, 5, 7, 29],
            1: [2, 4, 6, 28, 0],
            2: [8, 18, 26, 22, 12],
            3: [9, 13, 23, 27, 17],
            4: [11, 15, 21, 25, 19],
            5: [10, 14, 24, 20, 16]
            }

        _circle_embedding(g, d[0], center=(-1, 1), radius=.25)
        _circle_embedding(g, d[1], center=(1, 1), radius=.25)
        _circle_embedding(g, d[2], center=(-.8, 0), radius=.25, shift=2.5)
        _circle_embedding(g, d[3], center=(1.2, 0), radius=.25)
        _circle_embedding(g, d[4], center=(-1, -1), radius=.25, shift=2)
        _circle_embedding(g, d[5], center=(1, -1), radius=.25)

        return g

    elif embedding == 2:
        return g

    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def TutteGraph():
    r"""
    Returns the Tutte Graph.

    The Tutte graph is a 3-regular, 3-connected, and planar non-hamiltonian
    graph. For more information on the Tutte Graph, see the
    :wikipedia:`Tutte_graph`.

    EXAMPLES::

        sage: g = graphs.TutteGraph()
        sage: g.order()
        46
        sage: g.size()
        69
        sage: g.is_planar()
        True
        sage: g.vertex_connectivity() # long
        3
        sage: g.girth()
        4
        sage: g.automorphism_group().cardinality()
        3
        sage: g.is_hamiltonian()
        False
    """
    g = Graph(name="Tutte Graph")
    from sage.graphs.graph_plot import _circle_embedding

    g.add_cycle([(i,j) for i in range(3) for j in range(3) ])
    for i in range(3):
        g.add_cycle([(i,j) for j in range(9)])
        g.add_cycle([(i,j) for j in range(9,14)])
        g.add_edge((i,5),0)
        g.add_edge((i,13),(i,3))
        g.add_edge((i,12),(i,1))
        g.add_edge((i,11),(i,8))
        g.add_edge((i,10),(i,7))
        g.add_edge((i,6),(i,14))
        g.add_edge((i,4),(i,14))
        g.add_edge((i,9),(i,14))

    _circle_embedding(g, [(i,j) for i in range(3)  for j in range(6)], shift=.5)
    _circle_embedding(g, [(i,14) for i in range(3) ], radius=.3,shift=.25)

    for i in range(3):
        _circle_embedding(g, [(i,j) for j in range(3,9)]+[0]*5,
                          shift=3.7*(i-2)+.75,
                          radius=.4,
                          center=(.6*cos(2*(i+.25)*pi/3),.6*sin(2*(i+.25)*pi/3)))
        _circle_embedding(g, [(i,j) for j in range(9,14)],
                          shift=1.7*(i-2)+1,
                          radius=.2,
                          center=(.6*cos(2*(i+.25)*pi/3),.6*sin(2*(i+.25)*pi/3)))

    g.get_pos()[0] = (0,0)

    return g

def WagnerGraph():
    """
    Returns the Wagner Graph.

    See the :wikipedia:`Wikipedia page on the Wagner Graph
    <Wagner_graph>`.

    EXAMPLES::

        sage: g = graphs.WagnerGraph()
        sage: g.order()
        8
        sage: g.size()
        12
        sage: g.girth()
        4
        sage: g.diameter()
        2
        sage: g.show()
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(8, [4], 8)
    g.name("Wagner Graph")
    return g

def WatkinsSnarkGraph():
    r"""
    Returns the Watkins Snark Graph.

    The Watkins Graph is a snark with 50 vertices and 75 edges. For more
    information, see the :wikipedia:`Watkins_snark`.

    EXAMPLES::

        sage: g = graphs.WatkinsSnarkGraph()
        sage: g.order()
        50
        sage: g.size()
        75
        sage: g.chromatic_number()
        3
    """
    g = Graph(name="Watkins Snark Graph")

    for i in range(5):
        g.add_cycle([(i,j) for j in range(9)])
        _circle_embedding(g,
                          [(i,j) for j in range(4)]+[0]*2+[(i,4)]+[0]*2+[(i,j) for j in range(5,9)],
                          radius=.3,
                          center=(cos(2*(i+.25)*pi/5),sin(2*(i+.25)*pi/5)),
                          shift=2.7*i+7.55)
        g.add_edge((i,5),((i+1)%5,0))
        g.add_edge((i,8),((i+2)%5,3))
        g.add_edge((i,1),i)
        g.add_edge((i,7),i)
        g.add_edge((i,4),i)
        g.add_edge((i,6),(i,2))

    _circle_embedding(g, range(5), shift=.25, radius=1.1)
    return g

def WienerArayaGraph():
    r"""
    Returns the Wiener-Araya Graph.

    The Wiener-Araya Graph is a planar hypohamiltonian graph on 42 vertices and
    67 edges. For more information on the Wiener-Araya Graph, see its
    corresponding `Wolfram Page
    <http://mathworld.wolfram.com/Wiener-ArayaGraph.html>`_ or its `(french)
    Wikipedia page <http://fr.wikipedia.org/wiki/Graphe_de_Wiener-Araya>`_.

    EXAMPLES::

        sage: g = graphs.WienerArayaGraph()
        sage: g.order()
        42
        sage: g.size()
        67
        sage: g.girth()
        4
        sage: g.is_planar()
        True
        sage: g.is_hamiltonian() # not tested -- around 30s long
        False
        sage: g.delete_vertex(g.random_vertex())
        sage: g.is_hamiltonian()
        True
    """
    g = Graph(name="Wiener-Araya Graph")
    from sage.graphs.graph_plot import _circle_embedding

    g.add_cycle([(0,i) for i in range(4)])
    g.add_cycle([(1,i) for i in range(12)])
    g.add_cycle([(2,i) for i in range(20)])
    g.add_cycle([(3,i) for i in range(6)])
    _circle_embedding(g, [(0,i) for i in range(4)], shift=.5)
    _circle_embedding(g,
                      sum([[(1,3*i),(1,3*i+1)]+[0]*3+[(1,3*i+2)]+[0]*3 for i in range(4)],[]),
                      shift=4,
                      radius=.65)
    _circle_embedding(g, [(2,i) for i in range(20)], radius=.5)
    _circle_embedding(g, [(3,i) for i in range(6)], radius=.3, shift=.5)

    for i in range(4):
        g.delete_edge((1,3*i),(1,3*i+1))
        g.add_edge((1,3*i),(0,i))
        g.add_edge((1,3*i+1),(0,i))
        g.add_edge((2,5*i+2),(1,3*i))
        g.add_edge((2,5*i+3),(1,3*i+1))
        g.add_edge((2,(5*i+5)%20),(1,3*i+2))
        g.add_edge((2,(5*i+1)%20),(3,i+(i>=1)+(i>=3)))
        g.add_edge((2,(5*i+4)%20),(3,i+(i>=1)+(i>=3)))

    g.delete_edge((3,1),(3,0))
    g.add_edge((3,1),(2,4))
    g.delete_edge((3,4),(3,3))
    g.add_edge((3,4),(2,14))
    g.add_edge((3,1),(3,4))

    g.get_pos().pop(0)
    g.relabel()
    return g
