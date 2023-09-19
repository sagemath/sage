#!/usr/bin/env python

from sage_setup import sage_setup

sage_setup(['sagemath-modules'],
           interpreters=['CDF', 'RDF', 'RR', 'CC'],
           required_modules=('gsl',))
