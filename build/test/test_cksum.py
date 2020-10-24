# -*- coding: utf-8 -*-

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
from sage_bootstrap.cksum import CksumAlgorithm


class CksumTestCase(unittest.TestCase):

    def test_cksum_bytes(self):
        cksum = CksumAlgorithm()
        cksum.update(b'The quick brown fox jumps over the lazy dog\n')
        self.assertEqual(cksum.hexdigest(), '2382472371')

    def test_cksum_string(self):
        cksum = CksumAlgorithm()
        cksum.update('The quick brown fox jumps over the lazy dog\n')
        self.assertEqual(cksum.hexdigest(), '2382472371')
