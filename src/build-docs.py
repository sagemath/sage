# Add the current directory to the end of sys.path
# to ensure that an installed version of sage is used first.
import sys

current_dir = sys.path.pop(0)
sys.path.append(current_dir)

from sage_docbuild.__main__ import main

if __name__ == "__main__":
    main()
