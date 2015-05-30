"""
Saving Sage objects to a file (deprecated)
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import cPickle
import os
from sage.env import SAGE_ROOT

PATH = os.path.join(SAGE_ROOT, "db")

USE_DB = False

def path():
    from sage.misc.superseded import deprecation
    deprecation(17653, 'The sage.misc.db module is deprecated, use the load/save functions from sage.structure.sage_object instead')
    from sage.misc.misc import sage_makedirs
    sage_makedirs(PATH)

def save(x, filename, bzip2=False, gzip=False):
    """
    save(x, filename):

    Saves x to a file.  Pretty much the only constraint on x is that
    it have no circular references (it must be Python pickle-able).
    This uses the pickle module, so data you save is *guaranteed*
    to be readable by future versions of Python.

    INPUT:
       x -- almost arbitrary object
       filename -- a string

    OUTPUT:
       Creates a file named filename, from which the object x
       can be reconstructed.
    """
    from sage.misc.superseded import deprecation
    deprecation(17653, 'The sage.misc.db module is deprecated, use the load/save functions from sage.structure.sage_object instead')

    o=open(filename,"w")
    # Note: don't use protocol 2 here (use 1), since loading doesn't work
    # on my extension types.
    cPickle.dump(x,o,1)
    o.close()
    if bzip2:
        os.system("bzip2 -f %s"%filename)
    if gzip:
        os.system("gzip -f %s"%filename)


def load(filename, bzip2=False, gzip=False):
    """
    load(filename):

    Loads an object from filename and returns it.

    INPUT:
       filename -- a string that defines a valid file.  If the
          file doesn't exist then an IOError exception is raised.

    OUTPUT:
       An almost arbitrary object.
    """
    from sage.misc.superseded import deprecation
    deprecation(17653, 'The sage.misc.db module is deprecated, use the load/save functions from sage.structure.sage_object instead')

    if bzip2:
        os.system("bunzip2 -f -k %s"%(filename + ".bz2"))
    if gzip:
        os.system("cat %s.gz | gunzip -f > %s"%(filename,filename))
    assert os.path.exists(filename)
    o = open(filename,"r")
    X = cPickle.load(o)
    if bzip2 or gzip:
        os.remove(filename)
    return X


def save_db(x):
    """
    Save x to the database.  x must define a filename method.
    """
    from sage.misc.superseded import deprecation
    deprecation(17653, 'The sage.misc.db module is deprecated, use the load/save functions from sage.structure.sage_object instead')

    path()
    fn = PATH + x.filename()
    save(x,fn)
    os.system("bzip2 -f %s"%fn)


def load_db(x):
    """
    Load x from the database.  x must define a filename method.
    """
    from sage.misc.superseded import deprecation
    deprecation(17653, 'The sage.misc.db module is deprecated, use the load/save functions from sage.structure.sage_object instead')

    fn = PATH + x.filename()
    if os.path.exists(fn + ".bz2"):
        print("Loading {} from {}.".format(x, x.filename()))
        os.system("bunzip2 -f -k {}".format(fn + ".bz2"))
        o = open(fn, "r")
        x = cPickle.load(o)
        os.remove(fn)
        return x
    else:
        return None
