# -*- coding: utf-8 -*-
"""
Capture output for testing purposes
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
import contextlib
import logging

from sage_bootstrap.compat import StringIO

log = logging.getLogger()


class LogCaptureHandler(logging.Handler):

    def __init__(self, log_capture):
        self.records = log_capture.records
        logging.Handler.__init__(self)

    def emit(self, record):
        self.records.append(record)


class CapturedLog(object):

    def __init__(self):
        self.records = []

    def __enter__(self):
        self.old_level = log.level
        self.old_handlers = log.handlers
        log.level = logging.INFO
        log.handlers = [LogCaptureHandler(self)]
        return self

    def __exit__(self, type, value, traceback):
        log.level = self.old_level
        log.handlers = self.old_handlers

    def messages(self):
        return tuple((rec.levelname, rec.getMessage()) for rec in self.records)


@contextlib.contextmanager
def CapturedOutput():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err
