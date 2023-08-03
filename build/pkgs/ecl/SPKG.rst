ecl: An implementation of the Common Lisp language
==================================================

Description
-----------

ECL is an implementation of the Common Lisp language as defined by the
ANSI X3J13 specification. The most relevant features:

-  A bytecodes compiler and interpreter.
-  Compiles Lisp also with any C/C++ compiler.
-  It can build standalone executables and libraries.
-  ASDF, Sockets, Gray streams, MOP, and other useful components.
-  Extremely portable.
-  A reasonable license.

ECL supports the operating systems Linux, FreeBSD, NetBSD, OpenBSD,
Solaris and Windows, running on top of the Intel, Sparc, Alpha and
PowerPC processors. Porting to other architectures should be rather
easy.

Website: https://common-lisp.net/project/ecl/

License
-------

-  LGPL V2+ or compatible - for details see

   https://common-lisp.net/project/ecl/static/manual/Copyrights.html#Copyright-of-ECL


Upstream Contact
----------------

-  the ECL mailing list - see https://mailman.common-lisp.net/listinfo/ecl-devel

Special Update/Build Instructions
---------------------------------

-  Note: for the time being, ECL is built single threaded library as it
   seems to interact badly with the pexpect interface and Sage's signal
   handling when built multithreaded.

-  Do NOT quote SAGE_LOCAL when setting CPPFLAGS and/or LDFLAGS,
   in spkg-install as this caused the build to break. See
   :issue:`10187#comment:117`

-  TODO: Add the ECL test suite, and an spkg-check file to run it.
-  TODO: Make ECL use Sage's Boehm GC on MacOS X as well (but perhaps
   put some changes from ECL's into Sage's Boehm GC), then remove
   the src/src/gc directory, too.
