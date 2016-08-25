r"""
Tables of elliptic curves of given rank

The default database of curves contains the following data:

+------+------------------+--------------------+
| Rank | Number of curves | Maximal conductor  |
+======+==================+====================+
| 0    | 30427            |               9999 |
+------+------------------+--------------------+
| 1    | 31871            |               9999 |
+------+------------------+--------------------+
| 2    |  2388            |               9999 |
+------+------------------+--------------------+
| 3    |   836            |             119888 |
+------+------------------+--------------------+
| 4    |    10            |            1175648 |
+------+------------------+--------------------+
| 5    |     5            |           37396136 |
+------+------------------+--------------------+
| 6    |     5            |         6663562874 |
+------+------------------+--------------------+
| 7    |     5            |       896913586322 |
+------+------------------+--------------------+
| 8    |     6            |    457532830151317 |
+------+------------------+--------------------+
| 9    |     7            |      ~9.612839e+21 |
+------+------------------+--------------------+
| 10   |     6            |      ~1.971057e+21 |
+------+------------------+--------------------+
| 11   |     6            |      ~1.803406e+24 |
+------+------------------+--------------------+
| 12   |     1            |      ~2.696017e+29 |
+------+------------------+--------------------+
| 14   |     1            |      ~3.627533e+37 |
+------+------------------+--------------------+
| 15   |     1            |      ~1.640078e+56 |
+------+------------------+--------------------+
| 17   |     1            |      ~2.750021e+56 |
+------+------------------+--------------------+
| 19   |     1            |      ~1.373776e+65 |
+------+------------------+--------------------+
| 20   |     1            |      ~7.381324e+73 |
+------+------------------+--------------------+
| 21   |     1            |      ~2.611208e+85 |
+------+------------------+--------------------+
| 22   |     1            |      ~2.272064e+79 |
+------+------------------+--------------------+
| 23   |     1            |      ~1.139647e+89 |
+------+------------------+--------------------+
| 24   |     1            |      ~3.257638e+95 |
+------+------------------+--------------------+
| 28   |     1            |     ~3.455601e+141 |
+------+------------------+--------------------+

Note that lists for r>=4 are not exhaustive; there may well be curves
of the given rank with conductor less than the listed maximal conductor,
which are not included in the tables.

AUTHORS:
- William Stein (2007-10-07): initial version
- Simon Spicer (2014-10-24): Added examples of more high-rank curves

See also the functions cremona_curves() and cremona_optimal_curves()
which enable easy looping through the Cremona elliptic curve database.

"""
from __future__ import absolute_import

import os

from .constructor import EllipticCurve

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
            ((1, -1, 0, -106384, 13075804), 249649566346838)

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


