r"""
Group Algebras
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def GroupAlgebras(base_ring):
    r"""
    The category of group algebras over ``base_ring``.

    EXAMPLES::

        sage: C = GroupAlgebras(QQ); C
        Category of group algebras over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of hopf algebras with basis over Rational Field,
         Category of monoid algebras over Rational Field]

    This is just an alias for::

        sage: C is Groups().Algebras(QQ)
        True

    TESTS::

        sage: TestSuite(GroupAlgebras(ZZ)).run()
    """
    from sage.categories.all import Groups
    return Groups().Algebras(base_ring)

