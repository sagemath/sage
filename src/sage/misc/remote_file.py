"get_remote_file"
import os
from pathlib import Path
from urllib.request import Request, urlopen
from ssl import create_default_context as default_context


def get_remote_file(filename, verbose=True) -> Path:
    """
    INPUT:

    - ``filename`` -- the URL of a file on the web, e.g.,
      ``'http://modular.math.washington.edu/myfile.txt'``

    - ``verbose`` -- whether to display download status

    OUTPUT:

    This creates a file in the temp directory and returns the absolute path
    to that file as a :class:`Path` object.

    EXAMPLES::

        sage: url = 'https://www.sagemath.org/files/loadtest.py'
        sage: g = get_remote_file(url, verbose=False)      # optional - internet
        sage: with open(g) as f: print(f.read())           # optional - internet
        print("hi from the net")
        <BLANKLINE>
        print(2 + 3)
    """
    if verbose:
        print("Attempting to load remote file: " + filename)

    from sage.misc.temporary_file import tmp_filename
    ext = os.path.splitext(filename)[1]
    temp_name = Path(tmp_filename(ext=ext))
    # IMPORTANT -- urllib takes a long time to load,
    # so do not import it in the module scope.

    req = Request(filename, headers={"User-Agent": "sage-doctest"})

    if verbose:
        print("Loading started")

    content = urlopen(req, timeout=1, context=default_context())
    with temp_name.open('wb') as f:
        f.write(content.read())

    if verbose:
        print("Loading ended")

    return temp_name
