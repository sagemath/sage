#!/usr/bin/env python3
# Copy from $SAGE_DOC_SRC to $SAGE_DOC/src

from os import environ
from pathlib import Path
from shutil import copytree

from_path = environ['SAGE_DOC_SRC']
to_path = Path(environ['SAGE_DOC']) / "src"

print(f"Copying {from_path} to {to_path}")
copytree(from_path, to_path)
