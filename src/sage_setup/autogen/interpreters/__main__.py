# Usage: python -m sage_setup.autogen.interpreters <output_dir>

import argparse

from . import rebuild

parser = argparse.ArgumentParser()
parser.add_argument("output_dir", help="Output directory")
args = parser.parse_args()

rebuild(args.output_dir)
