# -*- coding: utf-8 -*-
"""
Test the is_url utility
"""

# ****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import unittest

from sage_bootstrap.util import is_url


class IsUrlTestCase(unittest.TestCase):

    def test_http(self):
        self.assertTrue(is_url('http://foo.bar/baz'))

    def test_https(self):
        self.assertTrue(is_url('https://foo.bar/baz'))

    def test_ftp(self):
        self.assertTrue(is_url('ftp://foo.bar/baz'))

    def test_nospace(self):
        self.assertFalse(is_url('http://foo. bar/baz'))

    def test_single_line(self):
        self.assertFalse(is_url('http://foo.bar/baz\nhttp://foo.bar/baz'))
