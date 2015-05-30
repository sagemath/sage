"""
Local copy of Sloane On-Line Encyclopedia of Integer Sequences

The SloaneEncyclopedia object provides access to a local copy of the database
containing only the sequences and their names. To use this you must download
and install the database using ``SloaneEncyclopedia.install()``, or
``SloaneEncyclopedia.install_from_gz()`` if you have already downloaded the
database manually.

To look up a sequence, type

::

    sage: SloaneEncyclopedia[60843]               # optional - sloane_database
    [1, 6, 21, 107]

To get the name of a sequence, type

::

    sage: SloaneEncyclopedia.sequence_name(1)     # optional - sloane_database
    'Number of groups of order n.'

To search locally for a particular subsequence, type

::

    sage: SloaneEncyclopedia.find([1,2,3,4,5], 1)    # optional - sloane_database
    [(15, [1, 2, 3, 4, 5, 7, 7, 8, 9, 11, 11, 13, 13, 16, 16, 16, 17, 19, 19, 23, 23, 23, 23, 25, 25, 27, 27, 29, 29, 31, 31, 32, 37, 37, 37, 37, 37, 41, 41, 41, 41, 43, 43, 47, 47, 47, 47, 49, 49, 53, 53, 53, 53, 59, 59, 59, 59, 59, 59, 61, 61, 64, 64, 64, 67, 67, 67, 71, 71, 71, 71, 73])]

The default maximum number of results is 30, but to return up to
100, type

::

    sage: SloaneEncyclopedia.find([1,2,3,4,5], 100)    # optional - sloane_database
    [(15, [1, 2, 3, 4, 5, 7, 7, 8, 9, 11, 11, ...

Results in either case are of the form [ (number, list) ].


.. SEEALSO::

    - If you want to get more informations relative to a sequence (references,
      links, examples, programs, ...), you can use the On-Line Encyclopedia of
      Integer Sequences provided by the :mod:`OEIS <sage.databases.oeis>`
      module.
    - Some infinite OEIS sequences are implemented in Sage, via the
      :mod:`sloane_functions <sage.combinat.sloane_functions>` module.


AUTHORS:

- Steven Sivek (2005-12-22): first version

- Steven Sivek (2006-02-07): updated to correctly handle the new
  search form on the Sloane website, and it's now also smarter about
  loading the local database in that it doesn't convert a sequence
  from string form to a list of integers until absolutely necessary.
  This seems to cut the loading time roughly in half.

- Steven Sivek (2009-12-22): added the SloaneEncyclopedia functions
  install() and install_from_gz() so users can get the latest versions
  of the OEIS without having to get an updated spkg; added
  sequence_name() to return the description of a sequence; and changed
  the data type for elements of each sequence from int to Integer.

- Thierry Monteil (2012-02-10): deprecate dead code and update related doc and
  tests.

Classes and methods
-------------------
"""

#*****************************************************************************
#
#      Sage: Copyright (C) 2005-2006 William Stein <wstein@gmail.com>
#                               and  Steven Sivek  <ssivek@mit.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import bz2, os, re, urllib

from sage.misc.all import verbose
from sage.env import SAGE_SHARE
import sage.rings.integer_ring
ZZ = sage.rings.integer_ring.IntegerRing()
from sage.misc.superseded import deprecation

class SloaneEncyclopediaClass:
    """
    A local copy of the Sloane Online Encyclopedia of Integer Sequences
    that contains only the sequence numbers and the sequences
    themselves.
    """
    def __init__(self):
        """
        Initialize the database but do not load any of the data.
        """
        self.__path__ = os.path.join(SAGE_SHARE, 'sloane')
        self.__file__ = os.path.join(self.__path__, 'sloane-oeis.bz2')
        self.__file_names__ = os.path.join(self.__path__, 'sloane-names.bz2')
        self.__loaded__ = False
        self.__loaded_names__ = False

    def __repr__(self):
        """
        String representation of this database. OUTPUT: str
        """
        return "Local copy of Sloane Online Encyclopedia of Integer Sequences"

    def __iter__(self):
        """
        Returns an iterator through the encyclopedia. Elements are of the
        form [number, sequence].
        """
        for i in self.__data__:
            yield [i, self[i]]

    def __getitem__(self, N):
        """
        Return sequence N in the encyclopedia. If sequence N does not
        exist, return [].

        INPUT:


        -  ``N`` - int


        OUTPUT: list
        """
        self.load()
        if not N in self.__data__: # sequence N does not exist
            return []
        if self.__data__[N][1] is None: # list N has not been created yet
            list = self.__data__[N][2].strip(',').split(',')
            self.__data__[N][1] = [ZZ(n) for n in list]
        return self.__data__[N][1]

    def __len__(self):
        """
        Return the number of sequences in the encyclopedia.
        """
        self.load()
        return len(self.__data__)


    def find(self, seq, maxresults=30):
        """
        Return a list of all sequences which have seq as a subsequence, up
        to maxresults results. Sequences are returned in the form (number,
        list).

        INPUT:


        -  ``seq`` - list

        -  ``maxresults`` - int


        OUTPUT: list of 2-tuples (i, v), where v is a sequence with seq as
        a subsequence.
        """
        self.load()

        answer, nanswer = [], 0
        pattern = re.sub(r'[\[\]]', ',', str(seq).replace(' ',''))
        for i in self.__data__:
            if self.__data__[i][2].find(pattern) != -1:
                answer.append((i, self[i]))
                nanswer = nanswer + 1
                if nanswer == maxresults:
                    return answer

        return answer

    def install(self, oeis_url="http://oeis.org/stripped.gz", names_url="http://oeis.org/names.gz", overwrite=False):
        """
        Download and install the online encyclopedia, raising an IOError if
        either step fails.

        INPUT:

        - ``oeis_url`` - string (default: "http://oeis.org...")
          The URL of the stripped.gz encyclopedia file.

        - ``names_url`` - string (default: "http://oeis.org...")
          The URL of the names.gz encyclopedia file.  If you do not want to
          download this file, set names_url=None.

        - ``overwrite`` - boolean (default: False) If the encyclopedia is
          already installed and overwrite=True, download and install the latest
          version over the installed one.
        """
        ### See if the encyclopedia already exists
        if not overwrite and os.path.exists(self.__file__):
            raise IOError("Sloane encyclopedia is already installed")

        tm = verbose("Downloading stripped version of Sloane encyclopedia")
        try:
            fname, _ = urllib.urlretrieve(oeis_url);
        except IOError as msg:
            raise IOError("%s\nError fetching the following website:\n    %s\nTry checking your internet connection."%(msg, oeis_url))

        if not names_url is None:
            try:
                nname, _ = urllib.urlretrieve(names_url);
            except IOError as msg:
                raise IOError("%s\nError fetching the following website:\n    %s\nTry checking your internet connection."%(msg, names_url))
        else:
            nname = None
        verbose("Finished downloading", tm)

        self.install_from_gz(fname, nname, overwrite)
        # Delete the temporary downloaded files
        os.remove(fname)
        if not nname is None:
            os.remove(nname)

    def install_from_gz(self, stripped_file, names_file, overwrite=False):
        """
        Install the online encyclopedia from a local stripped.gz file.

        INPUT:

        - ``stripped_file`` - string. The name of the stripped.gz OEIS file.

        - ``names_file`` - string.  The name of the names.gz OEIS file, or
          None if the user does not want it installed.

        - ``overwrite`` - boolean (default: False) If the encyclopedia is
          already installed and overwrite=True, install 'filename' over the
          old encyclopedia.
        """
        if not overwrite and os.path.exists(self.__file__):
            raise IOError("Sloane encyclopedia is already installed")

        copy_gz_file(stripped_file, self.__file__)

        if not names_file is None:
            copy_gz_file(names_file, self.__file_names__)
        else:
            # Delete old copies of names.gz since their sequence numbers
            # probably won't match the newly installed stripped.gz
            if os.path.exists(self.__file_names__):
                os.remove(self.__file_names__)

        # Remove the old database from memory so the new one will be
        # automatically loaded next time the user tries to access it
        self.unload()


    def load(self):
        """
        Load the entire encyclopedia into memory from a file. This is done
        automatically if the user tries to perform a lookup or a search.
        """
        if self.__loaded__ == True:
            return
        try:
            file_seq = bz2.BZ2File(self.__file__, 'r')
        except IOError:
            raise IOError("The Sloane Encyclopedia database must be installed.  Use e.g. 'SloaneEncyclopedia.install()' to download and install it.")

        self.__data__ = {}

        tm = verbose("Loading Sloane encyclopedia from disk")
        entry = re.compile(r'A(?P<num>\d{6}) ,(?P<body>.*),$');
        for L in file_seq:
            if len(L) == 0:
                continue
            m = entry.search(L)
            if m:
                seqnum = int(m.group('num'))
                msg = m.group('body').strip()
                self.__data__[seqnum] = [seqnum, None, ','+msg+',', None]
        file_seq.close()

        try:
            file_names = bz2.BZ2File(self.__file_names__, 'r')
            entry = re.compile(r'A(?P<num>\d{6}) (?P<body>.*)$');
            for L in file_names:
                if len(L) == 0: continue
                m = entry.search(L)
                if m:
                    seqnum = int(m.group('num'))
                    self.__data__[seqnum][3] = m.group('body').strip()
            file_names.close()
            self.__loaded_names__ = True
        except KeyError:
            ### Some sequence in the names file isn't in the database
            raise KeyError("Sloane OEIS sequence and name files do not match.  Try reinstalling, e.g. SloaneEncyclopedia.install(overwrite=True).")
        except IOError as msg:
            ### The names database is not installed
            self.__loaded_names__ = False

        verbose("Finished loading", tm)
        self.__loaded__ = True

    def sequence_name(self, N):
        """
        Return the name of sequence N in the encyclopedia. If sequence N
        does not exist, return ''.  If the names database is not installed,
        raise an IOError.

        INPUT:

        -  ``N`` - int

        OUTPUT: string

        EXAMPLES:

        sage: SloaneEncyclopedia.sequence_name(1) # optional - sloane_database
        'Number of groups of order n.'
        """
        self.load()
        if not self.__loaded_names__:
            raise IOError("The Sloane OEIS names file is not installed.  Try reinstalling, e.g. SloaneEncyclopedia.install(overwrite=True).")

        if not N in self.__data__: # sequence N does not exist
            return ''
        return self.__data__[N][3]

    def unload(self):
        """
        Remove the database from memory.
        """
        if self.__loaded__ == False:
            return
        del self.__data__
        self.__loaded__ = False
        self.__loaded_names__ = False

SloaneEncyclopedia = SloaneEncyclopediaClass()

def copy_gz_file(gz_source, bz_destination):
    """
    Decompress a gzipped file and install the bzipped verson.  This is
    used by SloaneEncyclopedia.install_from_gz to install several
    gzipped OEIS database files.

    INPUT:

    - ``gz_source`` - string. The name of the gzipped file.

    - ``bz_destination`` - string.  The name of the newly compressed file.
    """
    import gzip
    from sage.misc.misc import sage_makedirs

    # Read the gzipped input
    try:
        gz_input = gzip.open(gz_source, 'r')
        db_text = gz_input.read()
        gz_input.close()
    except IOError as msg:
        raise IOError("Error reading gzipped input file:\n%s"%msg)

    # Write the bzipped output
    try:
        sage_makedirs(os.path.dirname(bz_destination))
        bz2_output = bz2.BZ2File(bz_destination, 'w')
        bz2_output.write(db_text)
        bz2_output.close()
    except IOError as msg:
        raise IOError("Error writing bzipped output file:\n%s"%msg)

def parse_sequence(text=''):
    r"""
    This internal function was only used by the sloane_find function,
    which is now deprecated.

    TESTS::
        sage: from sage.databases.sloane import parse_sequence
        sage: parse_sequence()
        doctest:...: DeprecationWarning: The function parse_sequence is not used anymore (2012-01-01).
        See http://trac.sagemath.org/10358 for details.
    """
    deprecation(10358, "The function parse_sequence is not used anymore (2012-01-01).")

def sloane_sequence(number=1, verbose=True):
    r"""
    This function is broken. It is replaced by the
    :mod:`OEIS <sage.databases.oeis>` module.

    Type ``oeis?`` for more information.

    TESTS::

        sage: sloane_sequence(123)
        doctest:...: DeprecationWarning: The function sloane_sequence is deprecated. Use oeis() instead (2012-01-01).
        See http://trac.sagemath.org/10358 for details.
    """
    deprecation(10358,
            "The function sloane_sequence is deprecated. "
            "Use oeis() instead (2012-01-01).")

def sloane_find(list=[], nresults=30, verbose=True):
    r"""
    This function is broken. It is replaced by the
    :mod:`OEIS <sage.databases.oeis>` module.

    Type ``oeis?`` for more information.

    TESTS::

        sage: sloane_find([1,2,3])
        doctest:...: DeprecationWarning: The function sloane_find is deprecated. Use oeis() instead (2012-01-01).
        See http://trac.sagemath.org/10358 for details.
    """
    deprecation(10358,
            "The function sloane_find is deprecated. "
            "Use oeis() instead (2012-01-01).")

