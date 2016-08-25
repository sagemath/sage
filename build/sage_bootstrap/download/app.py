# -*- coding: utf-8 -*-
"""
Controller for the commandline actions
"""


#*****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import sys
import logging
log = logging.getLogger()


from sage_bootstrap.tarball import Tarball
from sage_bootstrap.download.mirror_list import MirrorList
from sage_bootstrap.download.transfer import Download


class Application(object):

    def __init__(self, timeout, quiet):
        import socket
        socket.setdefaulttimeout(timeout)
        self.quiet = quiet
    
    def print_fastest_mirror(self):
        print(MirrorList().fastest)

    def download_url(self, url, destination):
        Download(url, destination, progress=not self.quiet, ignore_errors=True).run()

    def download_tarball(self, tarball_filename, destination=None):
        tarball = Tarball(tarball_filename)
        tarball.download()
        if destination is not None:
            tarball.save_as(destination)
        
