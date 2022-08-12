# -*- coding: utf-8 -*-
"""
Python 2/3 compatibility utils
"""


#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************


try:
    # Python 3
    import urllib.request as urllib
    import urllib.parse as urlparse
except ImportError:
    import urllib
    import urlparse


try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
