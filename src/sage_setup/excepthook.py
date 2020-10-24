import os
import sys

def excepthook(*exc):
    """
    When an error occurs, display an error message similar to the error
    messages from ``sage-spkg``.

    In particular, ``build/make/install`` will recognize "sage" as a failed
    package, see :trac:`16774`.
    """
    stars = '*' * 72

    print(stars, file=sys.stderr)
    import traceback
    traceback.print_exception(*exc, file=sys.stderr)
    print(stars, file=sys.stderr)
    print("Error building the Sage library", file=sys.stderr)
    print(stars, file=sys.stderr)

    try:
        logfile = os.environ['SAGE_LOGFILE']  # Set by build/bin/sage-logger
    except Exception:
        pass
    else:
        print("Please email sage-devel (http://groups.google.com/group/sage-devel)", file=sys.stderr)
        print("explaining the problem and including the relevant part of the log file", file=sys.stderr)
        print("  " + logfile, file=sys.stderr)
        print("Describe your computer, operating system, etc.", file=sys.stderr)
        print(stars, file=sys.stderr)
