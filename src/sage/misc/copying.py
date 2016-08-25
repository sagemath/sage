"License"
from __future__ import absolute_import

import os
from . import pager

from sage.env import SAGE_ROOT

class License:
    def __call__(self):
        pager.pager()(str(self))

    def __repr__(self):
        return "Type license() to see the full license text."

    def __str__(self):
        return open(os.path.join(SAGE_ROOT,'COPYING.txt')).read()

license = License()
