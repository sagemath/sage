from __future__ import absolute_import, print_function

import argparse
import os

from . import rebuild

parser = argparse.ArgumentParser()
parser.add_argument("output_dir")
args = parser.parse_args()

output_dir = args.output_dir
if not output_dir:
    from sage.env import SAGE_SRC
    output_dir = os.path.join(SAGE_SRC, "sage", "ext", "interpreters")


rebuild(output_dir)
