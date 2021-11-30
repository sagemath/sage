ratpoints: Find rational points on hyperelliptic curves
=======================================================

Description
-----------

Michael Stoll's program which searches for rational points on
hyperelliptic curves.

NOTE: the ratpoints package has been assimilated by PARI/GP. Therefore,
this package (as Sage package) is deprecated. In the future, it will be
removed from Sage.


Upstream Contact
----------------

-  Author: Michael Stoll
-  Email: Michael.Stoll@uni-bayreuth.de
-  Website: http://www.mathe2.uni-bayreuth.de/stoll/programs/


Note on SSE2 instructions
~~~~~~~~~~~~~~~~~~~~~~~~~

-  On several architectures, the SSE2 instructions used by ratpoints
   cause
   compiler errors. In the case that ratpoints fails to build with SSE2
   instructions enabled, the build is repeated with SSE2 disabled.
