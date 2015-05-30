"License"

import os
import pager

from sage.env import SAGE_ROOT

class License:
    def __call__(self):
        pager.pager()(self.__str__())

    def __repr__(self):
        return "Type license() to see the full license text."

    def __str__(self):
        return open(os.path.join(SAGE_ROOT,'COPYING.txt')).read()

license = License()
