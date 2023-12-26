# sage_setup: distribution = sagemath-objects
"""
Pure python code for abstract base class for objects with generators
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


def multiplicative_iterator(M):
    from sage.rings.infinity import infinity
    G = M.gens()
    if len(G) == 0:
        yield M(1)
        return

    stop = [g.multiplicative_order() for g in G]
    for i in range(len(stop)):
        if stop[i] is infinity:
            raise ArithmeticError("%s is not finite." % M)
        stop[i] = stop[i] - 1
    z = M(1)
    yield z
    cnt = [0] * len(G)
    while cnt != stop:
        z = z * G[0]
        cnt[0] = cnt[0] + 1
        i = 0
        while i < len(cnt)-1 and cnt[i] > stop[i]:
            cnt[i] = 0
            cnt[i+1] = cnt[i+1] + 1
            z = z * G[i+1]
            i += 1
        yield z


def abelian_iterator(M):
    from sage.rings.infinity import infinity
    G = M.gens()
    if len(G) == 0:
        yield M(0)
        return

    stop = [g.additive_order() for g in G]
    for i in range(len(stop)):
        if stop[i] is infinity:
            raise ArithmeticError("%s is not finite." % M)
        stop[i] = stop[i] - 1
    z = M(0)
    yield z
    cnt = [0] * len(G)
    while cnt != stop:
        cnt[0] = cnt[0] + 1
        z = z + G[0]
        i = 0
        while i < len(cnt)-1 and cnt[i] > stop[i]:
            cnt[i] = 0
            cnt[i+1] = cnt[i+1] + 1
            z = z + G[i+1]
            i += 1
        yield z
