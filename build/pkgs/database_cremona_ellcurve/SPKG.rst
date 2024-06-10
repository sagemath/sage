database_cremona_ellcurve: Database of elliptic curves
======================================================

Description
-----------

John Cremona's database of elliptic curves

See https://github.com/JohnCremona/ecdata

This is an optional package, not included by default.

License
-------

Public Domain

Upstream Contact
----------------

-  Author: John Cremona
-  Email: john.cremona@gmail.com
-  Website: http://homepages.warwick.ac.uk/staff/J.E.Cremona/


Update Instructions
-------------------

Get an up-to-date copy of the git repository ecdata from
https://github.com/JohnCremona/ecdata.

If the cremona database has already been installed, remove
``SAGE_DATA/cremona/cremona.db``. Then run

The build script expects to find the files in subfolders allcurves,
allgens, degphi and allbsd of the ecdata folder. It extracts them and
builds the new cremona.db file from the contents.

Finally, copy ``SAGE_DATA/cremona/cremona.db`` to the src directory of
the spkg.
