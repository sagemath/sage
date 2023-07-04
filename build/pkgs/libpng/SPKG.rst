libpng: Bitmap image support
============================

Description
-----------

libpng is the official PNG reference library. It supports almost all PNG
features, is extensible, and has been extensively tested for over 13
years. The home site for development versions (i.e., may be buggy or
subject to change or include experimental features) is
http://libpng.sourceforge.net/, and the place to go for questions about
the library is the png-mng-implement mailing list.

Website: http://www.libpng.org/pub/png/libpng.html

License
-------

The libpng license - see
http://www.libpng.org/pub/png/src/libpng-LICENSE.txt


Upstream Contact
----------------

https://libpng.sourceforge.io

The png mailing lists - see
http://www.libpng.org/pub/png/pngmisc.html#lists

Special Update/Build Instructions
---------------------------------

-  On old versions of Darwin, the symbolic links libpng.\* created by
   libpng16 may
   interfere with a system-wide libPng.dylib.

   -- the following is very likely to be obsolete in 2014 ---

   This system-wide library is likely to be a different version and on
   top of that, the symbols exported there are prefixed with "_cg"
   (for "Core Graphics"). So even if by chance the functionalities of
   the two libraries were interchangeable, libraries or applications
   looking for one and being presented the other won't find the symbols
   they expect. Note the uppercase "P" which could prevent this
   conflict; unfortunately, the default filesystem used by Apple is
   case-insensitive.

   Note there would be no problem if the system-wide library was not
   looked for when Sage is being built or run, but that's not the case
   either; it is at least looked for by the "ImageIO" framework:

   -  when Python is built with Mac OS extensions, fixed in #4008;
   -  when Mercurial is built because it uses $EDITOR, cf. #4678;
   -  when R is built and it finds ``-lpng``, cf. #4409 and #11696.

   -- this is no longer done, as of #27186 ---

   As not all of these problems are easily dealt with and new ones may
   arise, we chose to delete the $SAGE_LOCAL/lib/libpng.\* symlinks.
   Therefore, some packages like Tachyon, which by default look for
   ``-lpng`` are patched to look for ``-lpng16`` instead.
