#!/usr/bin/env bash
# The file python3.inv contains the database of Sphinx hyperlink targets used by
# the intersphinx extension. See
#
#    http://sphinx-doc.org/ext/intersphinx.html
#
# To be able to compile Sage without accessing the net, we use a local copy of
# this database. Here is how to update it by downloading the file
# for the latest stable Python version.
#
# Likewise for other intersphinx targets.

set -x

if command -v wget > /dev/null 2>&1 ; then
    DOWNLOAD="wget -O -"
elif command -v curl > /dev/null 2>&1 ; then
    # On OS X, curl is installed by default, but not wget.
    DOWNLOAD=curl
else
    echo "Error: neither wget nor curl is installed."
    return 1
fi

rm -f python.inv python2.inv python3.inv
$DOWNLOAD https://docs.python.org/3/objects.inv > _vendor/python.inv
$DOWNLOAD https://docs.scipy.org/doc/scipy/reference/objects.inv > _vendor/scipy.inv
$DOWNLOAD https://flintlib.org/doc/objects.inv > _vendor/flint.inv
