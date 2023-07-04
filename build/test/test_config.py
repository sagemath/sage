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
from sage_bootstrap.config import Configuration, LOG_LEVELS
from .runnable import run_config_with


class ConfigurationTestCase(unittest.TestCase):

    def test_default(self):
        """
        Test the default configuration
        """
        config = Configuration()
        self.assertEqual(config.log, 'info')
        self.assertTrue(config.interactive)

    def test_example(self):
        """
        Test all ``SAGE_BOOTSTRAP`` settings
        """
        SAGE_BOOTSTRAP = ' loG:CrItIcAl, interactive:TRUE'
        result = run_config_with(SAGE_BOOTSTRAP)
        self.assertEqual(len(result), 4)
        self.assertEqual(result['log'], u'critical')
        self.assertTrue(result['interactive'])
        self.assertEqual(result['stdout'], 'default stdout')
        self.assertEqual(result['stderr'], 'default stderr')

    def test_logging(self):
        """
        Test that the different log levels are understood
        """
        for level in LOG_LEVELS:
            self.assertEqual(
                run_config_with('LOG:{0}'.format(level.upper()))['log'],
                level)

    def test_overriding(self):
        """
        Test that overriding the isatty detection works
        """
        interactive = run_config_with('interactive:true')
        self.assertTrue(interactive['interactive'])
        self.assertEqual(interactive['stdout'], 'default stdout')
        self.assertEqual(interactive['stderr'], 'default stderr')
        in_pipe = run_config_with('interactive:false')
        self.assertFalse(in_pipe['interactive'])
        self.assertEqual(in_pipe['stdout'],
                         u"<class 'sage_bootstrap.stdio.UnbufferedStream'>")
        self.assertEqual(in_pipe['stderr'], 'default stderr')
