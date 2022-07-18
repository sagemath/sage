readline: Command line editing library
======================================

Description
-----------

The GNU Readline library provides a set of functions for use by
applications that allow users to edit command lines as they are typed
in. Both Emacs and vi editing modes are available. The Readline library
includes additional functions to maintain a list of previously-entered
command lines, to recall and perhaps reedit those lines, and perform
csh-like history expansion on previous commands.

Website: http://tiswww.case.edu/php/chet/readline/rltop.html

License
-------

-  GPL V3+


Upstream Contact
----------------

-  Chet Ramey at http://cnswww.cns.cwru.edu/~chet

Special Update/Build Instructions
---------------------------------

We build readline using ncurses. Readline needs to be told to link with
libtinfo (part of ncurses), this is what the patch 0002-ltinfo.patch
does.

Patches
-------

-  0001-macports.patch: Changes to shobj.conf for OS/X, from macports:

   https://trac.macports.org/browser/trunk/dports/devel/readline/files/patch-shobj-conf.diff

-  0002-ltinfo.patch: We build readline using ncurses, and for that it
   needs to be told to link with libtinfo (part of ncurses).
