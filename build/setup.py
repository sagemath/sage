#!/usr/bin/env sage-bootstrap-python

from distutils.core import setup

setup(
    name='sage_bootstrap',
    description='',
    author='Volker Braun',
    author_email='vbraun.name@gmail.com',
    packages=[
        'sage_bootstrap',
        'sage_bootstrap.download',
        'sage_bootstrap.compat'
    ],
    scripts=['bin/sage-package', 'bin/sage-download-file', 'bin/sage-bootstrap-python'],
    version='1.0',
    url='https://www.sagemath.org',
)
