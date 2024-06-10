tachyon: A ray tracing system
=============================

Description
-----------

Tachyon is a raytracer developed by John E. Stone. Tachyon supports the
typical ray tracer features, most of the common geometric primitives,
shading and texturing modes, etc. It also supports less common features
such as HDR image output, ambient occlusion lighting, and support for
various triangle mesh and volumetric texture formats beneficial for
molecular visualization (e.g. rendering VMD scenes).

Currently not all of Tachyon's functionality is exported by the Sage
interface.

License
-------

Copyright (c) 1994-2010 John E. Stone
All rights reserved.


Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

3. The name of the author may not be used to endorse or promote
   products
   derived from this software without specific prior written permission.


Upstream Contact
----------------

- http://jedi.ks.uiuc.edu/~johns/raytracer/
- http://www.photonlimited.com/~johns/raytracer/
- John Stone <johns@ks.uiuc.edu>

Dependencies
------------

This spkg depends on:

-  libpng


Special Update/Build Instructions
---------------------------------

-  Delete the scenes directory, which has lots of cool examples.
-  Delete the msvc directory, which is also large and not used within
   Sage.
-  The CVS subdirectories are currently (almost) empty, but should
   otherwise be deleted.

-  The upstream files had strange permissions, i.e. some source files
   were executable, while almost all files weren't world-readable.

-  There's seems to be some crap like ``tachyon.html.tar.gz`` and a few
   ``.#*`` files I haven't [yet] deleted, since they're not that large.

-  TODO: Check whether building multi-threaded versions on MacOS X
   meanwhile works. (This was said to fail with an old beta.)

-  TODO: Use ``patch`` instead of copying over pre-patched files.
-  TODO: [Optionally] also install some of the documentation.
-  TODO: I doubt the CFLAGS set for AIX and HP-UX won't get overridden
   by the created Makefile, but that's a minor issue. -leif
