# -*- coding: utf-8 -*-
"""
Set Up Input/output

Output should always be unbuffered so that it appears immediately on
the terminal.
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import sys


class UnbufferedStream(object):

    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


REAL_STDOUT = sys.stdout
REAL_STDERR = sys.stderr


def init_streams(config):
    if not config.interactive:
        sys.stdout = UnbufferedStream(REAL_STDOUT)


def flush():
    REAL_STDOUT.flush()
    REAL_STDERR.flush()
