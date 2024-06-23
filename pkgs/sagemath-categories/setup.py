#!/usr/bin/env python

from sage_setup import sage_setup

sage_setup(
    ['sagemath-categories'],
    interpreters=['Element', 'Python'])  # RDF uses gsl --> sagemath-modules
