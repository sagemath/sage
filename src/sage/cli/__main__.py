import os
import sys

# Allow importing modules from the current directory, matching python behavior
sys.path.append(os.getcwd())

from sage.cli import main

sys.exit(main())
