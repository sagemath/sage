#!/usr/bin/env python3
# Usage: python tools/autogen/interpreters <output_dir>

import argparse

from internal import rebuild

parser = argparse.ArgumentParser()
parser.add_argument("output_dir", help="Output directory")
args = parser.parse_args()

rebuild(args.output_dir)
