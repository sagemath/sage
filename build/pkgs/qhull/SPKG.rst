qhull: Compute convex hulls, Delaunay triangulations, Voronoi diagrams
======================================================================

Description
-----------

From the README.txt of Qhull:

Qhull computes convex hulls, Delaunay triangulations, Voronoi diagrams,
furthest-site Voronoi diagrams, and halfspace intersections about a
point. It runs in 2-d, 3-d, 4-d, or higher. It implements the Quickhull
algorithm for computing convex hulls. Qhull handles round-off errors
from floating point arithmetic. It can approximate a convex hull.

The program includes options for hull volume, facet area, partial hulls,
input transformations, randomization, tracing, multiple output formats,
and execution statistics.

Further notes:

The qhull library is already shipped with the Python library scipy (from
version 1.4), see

-  http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
-  http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html
-  http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Voronoi.html

There is also the Python interface Pyhull available on PyPI
https://pypi.python.org/pypi/pyhull (see also documentation at
http://pythonhosted.org/pyhull/).


Upstream Contact
----------------

http://www.qhull.org/html

C. Bradford Barber bradb@shore.net or qhull@qhull.org

Dependencies
------------

Can be compiled with Qt support, but the Sage version currently doesn't
try to do this.

License
-------

Not a standard license, but Sage compatible. See the COPYING.txt file in
the source directory for details.
