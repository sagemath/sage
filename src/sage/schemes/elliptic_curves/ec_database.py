r"""
Tables of elliptic curves of given rank

The default database of curves contains the following data:

+------+------------------+--------------------+
| Rank | Number of curves | Maximal conductor  |
+======+==================+====================+
| 0    | 30427            |            9999    |
+------+------------------+--------------------+
| 1    | 31871            |            9999    |
+------+------------------+--------------------+
| 2    |  2388            |            9999    |
+------+------------------+--------------------+
| 3    |   836            |          119888    |
+------+------------------+--------------------+
| 4    |     1            |          234446    |
+------+------------------+--------------------+
| 5    |     1            |        19047851    |
+------+------------------+--------------------+
| 6    |     1            |      5187563742    |
+------+------------------+--------------------+
| 7    |     1            |    382623908456    |
+------+------------------+--------------------+
| 8    |     1            | 457532830151317    |
+------+------------------+--------------------+

AUTHOR:
- William Stein (2007-10-07): initial version

See also the functions cremona_curves() and cremona_optimal_curves()
which enable easy looping through the Cremona elliptic curve database.

"""

import os

from constructor import EllipticCurve

class EllipticCurves:
    def rank(self, rank, tors=0, n=10, labels=False):
        r"""
        Return a list of at most `n` non-isogenous curves with given
        rank and torsion order.

        INPUT:

        - ``rank`` (int) -- the desired rank

        - ``tors`` (int, default 0) -- the desired torsion order (ignored if 0)

        - ``n`` (int, default 10) -- the maximum number of curves returned.

        - ``labels`` (bool, default False) -- if True, return Cremona
          labels instead of curves.

        OUTPUT:

        (list) A list at most `n` of elliptic curves of required rank.

        EXAMPLES::

            sage: elliptic_curves.rank(n=5, rank=3, tors=2, labels=True)
            ['59450i1', '59450i2', '61376c1', '61376c2', '65481c1']

        ::

            sage: elliptic_curves.rank(n=5, rank=0, tors=5, labels=True)
            ['11a1', '11a3', '38b1', '50b1', '50b2']

        ::

            sage: elliptic_curves.rank(n=5, rank=1, tors=7, labels=True)
            ['574i1', '4730k1', '6378c1']

        ::

            sage: e = elliptic_curves.rank(6)[0]; e.ainvs(), e.conductor()
            ((1, 1, 0, -2582, 48720), 5187563742)
            sage: e = elliptic_curves.rank(7)[0]; e.ainvs(), e.conductor()
            ((0, 0, 0, -10012, 346900), 382623908456)
            sage: e = elliptic_curves.rank(8)[0]; e.ainvs(), e.conductor()
            ((0, 0, 1, -23737, 960366), 457532830151317)

        """
        from sage.env import SAGE_SHARE
        db = os.path.join(SAGE_SHARE,'ellcurves')
        data = os.path.join(db,'rank%s'%rank)
        if not os.path.exists(data):
            return []
        v = []
        tors = int(tors)
        for w in open(data).readlines():
            N,iso,num,ainvs,r,t = w.split()
            if tors and tors != int(t):
                continue
            label = '%s%s%s'%(N,iso,num)
            if labels:
                v.append(label)
            else:
                E = EllipticCurve(eval(ainvs))
                E._set_rank(r)
                E._set_torsion_order(t)
                E._set_conductor(N)
                E._set_cremona_label(label)
                v.append(E)
            if len(v) >= n:
                break
        return v

elliptic_curves = EllipticCurves()


