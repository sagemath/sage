r"""
Projective planarity

This module contains the 35 minimal forbidden minors for projective plane
and a function for checking if a graph `G` has one of them as a minor.

EXAMPLES::

The Peterson graph is a known projective planar graph so it doesn't have a
forbidden minor::

    sage: P = graphs.PetersenGraph()
    sage: _ = get_p2_forbidden_minor(P); type(_)        # long
    <class 'NoneType'>

K_{4,4} has a projective plane crossing number of 2. One of the minimal
forbidden minors is K_{4,4} - e, so we get a 1-to-1 dictionary from
:meth:`~Graph.minor`::

    sage: K44 = graphs.CompleteBipartiteGraph(4, 4)
    sage: minor_map = get_p2_forbidden_minor(K44); minor_map
    {0: [0], 1: [1], 2: [2], 3: [4], 4: [8], 5: [6], 6: [5], 7: [7]}

AUTHORS:

- Juan M. Lazaro Ruiz (2025-01-27): initial version
- Steve Schluchter (2025-04-06). gently edited

"""

# ****************************************************************************
#       Copyright (C) 2024 Juan M. Lazaro Ruiz <juan.m.lazaro.ruiz@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

P2_FORBIDDEN_MINORS = [
    'KFz_????wF?[',
    'J~{???F@oM?',
    'I~{?GKF@w',
    'JFz_?AB_sE?',
    'I~{?CME`_',
    'H~}CKMF',
    'G^~EMK',
    'H^|ACME',
    'Himp`cr',
    'Iimp_CpKO',
    'IFz@GCdHO',
    'IBz__aB_o',
    'FQ~~w',
    'GlvJ`k',
    'HilKH`J',
    'GjlKJs',
    'HhI]ECZ',
    'HiMIKSp',
    'HFwO]Kf',
    'I]q?a?n@o',
    'IHIWuFGo_',
    'IXJWMC`Eg',
    'GFzfF?',
    'I]o__OF@o',
    'G?^vf_',
    'H?]ufBo',
    'GlrHhs',
    'HhIWuRB',
    'IXCO]FGb?',
    'Fvz~o',
    'GlfH]{',
    'Hl`HGvV',
    'HhcIHmv',
    'IhEGICRiw',
    'JhEIDSD?ga_'
]