# -*- coding: utf-8 -*-
"""
Test the printing and logging
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


import unittest
from textwrap import dedent
from .runnable import run_log_with


class LoggerTestCase(unittest.TestCase):

    def test_interactive(self):
        stdout, stderr = run_log_with('interactive:true')
        self.assertEqual(stderr.strip(), dedent("""
            WARNING [runnable|print_log:25]: This is the warning log level
            CRITICAL [runnable|print_log:26]: This is the critical log level
            ERROR [runnable|print_log:27]: This is the error log level
        """).strip())
        self.assertEqual(stdout.strip(), dedent("""
            This is the info log level
            This is printed
        """).strip())

    def test_noninteractive(self):
        stdout, stderr = run_log_with('interactive:false')
        self.assertEqual(stderr.strip(), dedent("""
            This is the info log level
            WARNING [runnable|print_log:25]: This is the warning log level
            CRITICAL [runnable|print_log:26]: This is the critical log level
            ERROR [runnable|print_log:27]: This is the error log level
        """).strip())
        self.assertEqual(stdout.strip(), dedent("""
            This is printed
        """).strip())

    def test_debug(self):
        """
        The lowest logging level
        """
        stdout, stderr = run_log_with('log:debug,interactive:true')
        self.assertEqual(stderr.strip(), dedent("""
            DEBUG [runnable|print_log:23]: This is the debug log level
            WARNING [runnable|print_log:25]: This is the warning log level
            CRITICAL [runnable|print_log:26]: This is the critical log level
            ERROR [runnable|print_log:27]: This is the error log level
        """).strip())
        self.assertEqual(stdout.strip(), dedent("""
            This is the info log level
            This is printed
        """).strip())

    def test_error(self):
        """
        The highest logging level
        """
        stdout, stderr = run_log_with('log:error,interactive:true')
        self.assertEqual(stderr.strip(), dedent("""
            CRITICAL [runnable|print_log:26]: This is the critical log level
            ERROR [runnable|print_log:27]: This is the error log level
        """).strip())
        self.assertEqual(stdout.strip(), dedent("""
            This is printed
        """).strip())
