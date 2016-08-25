r"""
Object persistence

You can load and save most Sage object to disk using the load and
save member functions and commands.

.. note::

   It is impossible to save certain Sage objects to disk. For example,
   if `x` is a MAGMA object, i.e., a wrapper around an object
   that is defined in MAGMA, there is no way to save `x` it to
   disk, since MAGMA doesn't support saving of individual objects to
   disk.


-  Versions: Loading and saving of objects is guaranteed to work
   even if the version of Python changes. Saved objects can be loaded
   in future versions of Python. However, if the data structure that
   defines the object, e.g., in Sage code, changes drastically (or
   changes name or disappears), then the object might not load
   correctly or work correctly.

-  Objects are zlib compressed for space efficiency.
"""
from __future__ import absolute_import

from sage.structure.sage_object import save, load, \
     loads, dumps, SageObject

from .misc import SAGE_DB

def load_sage_object(cls, dic):   # not used
    X = cls.__new__(cls)
    try:
        X.__setstate__(dic)
    except AttributeError:
        X.__dict__ = dic
    return X

from six.moves import cPickle
def load_sage_element(cls, parent, dic_pic):
    X = cls.__new__(cls)
    X._set_parent(parent)
    X.__dict__ = cPickle.loads(dic_pic)
    return X

def db(name):
    r"""
    Load object with given name from the Sage database. Use x.db(name)
    or db_save(x, name) to save objects to the database.

    The database directory is ``$HOME/.sage/db``.
    """
    return load('%s/%s'%(SAGE_DB,name))

def db_save(x, name=None):
    r"""
    Save x to the Sage database.

    The database directory is ``$HOME/.sage/db``.
    """
    try:
        x.db(name)
    except AttributeError:
        save(x, '%s/%s'%(SAGE_DB,name))

