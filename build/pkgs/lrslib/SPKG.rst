lrslib: Reverse search algorithm for vertex enumeration and convex hull problems
================================================================================

Description
-----------

lrslib implements the linear reverse search algorithm of Avis and
Fukuda.

See the homepage (http://cgm.cs.mcgill.ca/~avis/C/lrs.html) for details.

We use an autotoolized version from
https://github.com/mkoeppe/lrslib/tree/autoconfiscation

License
-------

lrslib is released under a GPL v2+ license.


Upstream Contact
----------------

David Avis, avis at cs dot mcgill dot edu.

Dependencies
------------

To build and install the "plrs" binary, a multi-thread version of lrs,
need to first install the full Boost package ("sage -i boost").

If the package finds an MPI C++ compiler script (mpic++), it also builds
and installs the "mplrs" binary, a distributed version of lrs using MPI.

(Sage currently does not make use of plrs and mplrs.)


Special Update/Build Instructions
---------------------------------
