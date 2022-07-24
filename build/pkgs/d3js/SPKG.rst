d3js: JavaScript library for manipulating documents based on data
=================================================================

Description
-----------

D3.js is a JavaScript library for manipulating documents based on data.
The file d3.min.js will be placed into the ${SAGE_SHARE}/d3js/
directory.

License
-------

BSD 3-Clause License


Upstream Contact
----------------

- Author: Mike Bostock (http://bost.ocks.org/mike/)
- Home page: http://d3js.org/

Special Update/Build Instructions
---------------------------------

Two kind of archives can be downloaded from d3.js website: one with all
source code and tests that weights 2,9M (both in zip and tar.gz formats)
and one with the final javascript scripts which weights 121K (zip format
only). Since testing requires node.js that is not shipped with Sage, we
currenlty ship the final js only. Hence we have to transform it from zip
to tar.gz format. Running sage-src should do all the repackaging job.
