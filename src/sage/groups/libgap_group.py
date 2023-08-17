# sage_setup: distribution = sagemath-gap
"""
Generic LibGAP-based Group

This is useful if you need to use a GAP group implementation in Sage
that does not have a dedicated Sage interface.

If you want to implement your own group class, you should not derive
from this but directly from
:class:`~sage.groups.libgap_wrapper.ParentLibGAP`.

EXAMPLES::

    sage: F.<a,b> = FreeGroup()
    sage: G_gap = libgap.Group([ (a*b^2).gap() ])
    sage: from sage.groups.libgap_group import GroupLibGAP
    sage: G = GroupLibGAP(G_gap);  G
    Group([ a*b^2 ])
    sage: type(G)
    <class 'sage.groups.libgap_group.GroupLibGAP_with_category'>
    sage: G.gens()
    (a*b^2,)
"""

##############################################################################
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################


from sage.groups.group import Group
from sage.groups.libgap_wrapper import ParentLibGAP, ElementLibGAP
from sage.groups.libgap_mixin import GroupMixinLibGAP

class GroupLibGAP(GroupMixinLibGAP, Group, ParentLibGAP):

    Element = ElementLibGAP

    def __init__(self, *args, **kwds):
        """
        Group interface for LibGAP-based groups.

        INPUT:

        Same as :class:`~sage.groups.libgap_wrapper.ParentLibGAP`.

        TESTS::

            sage: F.<a,b> = FreeGroup()
            sage: G_gap = libgap.Group([ (a*b^2).gap() ])
            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: G = GroupLibGAP(G_gap);  G
            Group([ a*b^2 ])
            sage: g = G.gen(0);  g
            a*b^2
            sage: TestSuite(G).run(skip=['_test_pickling', '_test_elements'])
            sage: TestSuite(g).run(skip=['_test_pickling'])
        """
        ParentLibGAP.__init__(self, *args, **kwds)
        Group.__init__(self)
