bzip2: High-quality data compressor
===================================

Description
-----------

bzip2 is a freely available, patent free, high-quality data compressor.

It typically compresses files to within 10% to 15% of the best available
techniques (the PPM family of statistical compressors), whilst being
around twice as fast at compression and six times faster at
decompression.

License
-------

BSD-style


Upstream Contact
----------------

-  Website http://bzip.org/
-  Author: Julian Seward <julian@bzip.org>

Special Update/Build Instructions
---------------------------------

This package must not be bzip2 compressed, so create it using ::

    tar c bzip2-1.0.6 | gzip --best >bzip2-1.0.6.spkg

The build system has been autotoolized based on a patch by the Suse folk
at
http://ftp.uni-kl.de/pub/linux/suse/people/sbrabec/bzip2/for_downstream/bzip2-1.0.6-autoconfiscated.patch

See patches/autotools and spkg-src for details.
