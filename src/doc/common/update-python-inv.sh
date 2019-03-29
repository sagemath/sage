#!/usr/bin/env bash
# The files python2.inv and python3.inv contain the database of Sphinx hyperlink targets
# used by the intersphinx extension. See
#    http://sphinx-doc.org/ext/intersphinx.html
# To be able to compile Sage without accessing the net, we use a local copy of
# this database. Here is how to update it:

if command -v wget > /dev/null 2>&1 ; then
    rm -f python.inv python2.inv python3.inv
    wget https://docs.python.org/release/`sage --python2 --version 2>&1 | cut -d " " -f 2`/objects.inv -O - > python2.inv
    wget https://docs.python.org/release/`sage --python3 --version 2>&1 | cut -d " " -f 2`/objects.inv -O - > python3.inv
elif command -v curl > /dev/null 2>&1 ; then
    # On OS X, curl is installed by default, but not wget.
    rm -f python.inv python2.inv python3.inv
    curl https://docs.python.org/release/`sage --python2 --version 2>&1 | cut -d " " -f 2`/objects.inv > python2.inv
    curl https://docs.python.org/release/`sage --python3 --version 2>&1 | cut -d " " -f 2`/objects.inv > python3.inv
else
    echo "Error: neither wget nor curl is installed."
    return 1
fi
