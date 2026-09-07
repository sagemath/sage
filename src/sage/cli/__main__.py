import os
import sys

# Allow importing modules from the current directory, matching python behavior
sys.path.append(os.getcwd())

from sage.cli import main

if __name__ == '__main__':
    sys.exit(main())
