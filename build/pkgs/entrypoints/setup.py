#!/usr/bin/env python
#
# Upstream does not have a setup.py file, so we use this instead
#

from entrypoints import __version__ as v

from setuptools import setup
setup(name="entrypoints", version=v, py_modules=["entrypoints"])
