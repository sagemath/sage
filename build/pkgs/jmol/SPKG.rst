
Jmol for Sage
=============

Description
-----------

This provides files necessary for Jmol(java) and JSmol (javascript) to
operate from the command line and the Notebook. It does not contain the
Notebook javascript library jmol_lib.js or changes to Notebook or Sage
code.

License
-------

GPLv2+


Upstream Contact
----------------

-  Bob Hanson
-  e-mail: hansonr@stolaf.edu
-  Homepage: https://www.stolaf.edu/people/hansonr/
-  Development page: https://github.com/BobHanson/Jmol-SwingJS
-  Download page: https://sourceforge.net/projects/jmol/files/Jmol/

Dependencies
------------

No build-time dependencies.

The commandline jmol requires java at runtime.


Special Build Instructions
--------------------------

To avoid depending on ``unzip`` at build time, we have to repack the
tarball, see ``spkg-src``. We take the opportunity to remove some
unnecessary subdirectories, see
http://wiki.jmol.org/index.php/Jmol_JavaScript_Object#In_detail
