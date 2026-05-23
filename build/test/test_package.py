# -*- coding: utf-8 -*-
"""
Test Sage Package Handling
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
from types import SimpleNamespace
from unittest.mock import patch

from sage_bootstrap.package import Package
from sage_bootstrap.tarball import Tarball


class PackageTestCase(unittest.TestCase):

    maxDiff = None

    def test_package(self):
        pkg = Package('pari')
        self.assertTrue(pkg.name, 'pari')
        self.assertTrue(pkg.path.endswith('build/pkgs/pari'))
        self.assertEqual(pkg.tarball_pattern, 'pari-VERSION.tar.gz')
        self.assertEqual(pkg.tarball_filename, pkg.tarball.filename)
        self.assertTrue(pkg.tarball.filename.startswith('pari-') and
                        pkg.tarball.filename.endswith('.tar.gz'))
        self.assertTrue(pkg.tarball.filename.startswith('pari-') and
                        pkg.tarball.filename.endswith('.tar.gz'))
        self.assertTrue(isinstance(pkg.tarball, Tarball))

    def test_all(self):
        pari = Package('pari')
        self.assertTrue(pari in Package.all())

    @patch('subprocess.run')
    @patch('sage_bootstrap.package.os.path.exists', return_value=True)
    def test_choose_whl(self, mock_exists, mock_run):
        pkg = Package('rpds_py')
        pkg._Package__tarballs_info = [
            {'tarball': 'demo-VERSION-py3-none-any.whl', 'sha256': '1', 'sha1': None, 'upstream_url': None},
            {'tarball': 'demo-VERSION-cp313-abi3-manylinux_2_17_x86_64.whl', 'sha256': '2', 'sha1': None, 'upstream_url': None},
            {'tarball': 'demo-VERSION-cp313-cp313-manylinux_2_17_x86_64.whl', 'sha256': '3', 'sha1': None, 'upstream_url': None},
            {'tarball': 'demo-VERSION.tar.gz', 'sha256': '4', 'sha1': None, 'upstream_url': None},
        ]
        mock_run.side_effect = [
            SimpleNamespace(
                returncode=0,
                stdout='["cp313-cp313-manylinux_2_17_x86_64", "cp313-abi3-manylinux_2_17_x86_64", "py3-none-any"]\n',
                stderr='',
            ),
            SimpleNamespace(
                returncode=0,
                stdout=(
                    '{"demo-0.30.0-py3-none-any.whl": ["py3-none-any"], '
                    '"demo-0.30.0-cp313-abi3-manylinux_2_17_x86_64.whl": ["cp313-abi3-manylinux_2_17_x86_64"], '
                    '"demo-0.30.0-cp313-cp313-manylinux_2_17_x86_64.whl": ["cp313-cp313-manylinux_2_17_x86_64"]}\n'
                ),
                stderr='',
            ),
        ]

        tarball_info = pkg.find_tarball_for_platform()

        self.assertEqual(
            tarball_info['tarball'],
            'demo-VERSION-cp313-cp313-manylinux_2_17_x86_64.whl',
        )
        self.assertIn(
            'packaging.tags.sys_tags()',
            mock_run.call_args_list[0].args[0][3],
        )
        self.assertEqual(
            mock_run.call_args_list[1].args[0][-3:],
            [
                'demo-0.30.0-py3-none-any.whl',
                'demo-0.30.0-cp313-abi3-manylinux_2_17_x86_64.whl',
                'demo-0.30.0-cp313-cp313-manylinux_2_17_x86_64.whl',
            ],
        )

    @patch('subprocess.run')
    @patch('sage_bootstrap.package.os.path.exists', return_value=True)
    def test_fallback_to_tarball_when_no_wheel_matches(self, mock_exists, mock_run):
        pkg = Package('rpds_py')
        pkg._Package__tarballs_info = [
            {'tarball': 'demo-VERSION-cp313-cp313-manylinux_2_17_x86_64.whl', 'sha256': '1', 'sha1': None, 'upstream_url': None},
            {'tarball': 'demo-VERSION-cp313-cp313-win_amd64.whl', 'sha256': '2', 'sha1': None, 'upstream_url': None},
            {'tarball': 'demo-VERSION.tar.gz', 'sha256': '3', 'sha1': None, 'upstream_url': None},
        ]
        mock_run.side_effect = [
            SimpleNamespace(
                returncode=0,
                stdout='["cp313-cp313-manylinux_2_17_aarch64"]\n',
                stderr='',
            ),
            SimpleNamespace(
                returncode=0,
                stdout=(
                    '{"demo-0.30.0-cp313-cp313-manylinux_2_17_x86_64.whl": ["cp313-cp313-manylinux_2_17_x86_64"], '
                    '"demo-0.30.0-cp313-cp313-win_amd64.whl": ["cp313-cp313-win_amd64"]}\n'
                ),
                stderr='',
            ),
        ]

        tarball_info = pkg.find_tarball_for_platform()

        self.assertEqual(tarball_info['tarball'], 'demo-VERSION.tar.gz')

    @patch('subprocess.run')
    @patch('sage_bootstrap.package.os.path.exists', return_value=True)
    def test_fallback_to_tarball_when_tag_query_fails(self, mock_exists, mock_run):
        pkg = Package('rpds_py')
        pkg._Package__tarballs_info = [
            {'tarball': 'demo-VERSION-cp313-cp313-manylinux_2_17_x86_64.whl', 'sha256': '1', 'sha1': None, 'upstream_url': None},
            {'tarball': 'demo-VERSION.tar.gz', 'sha256': '2', 'sha1': None, 'upstream_url': None},
        ]
        mock_run.return_value = SimpleNamespace(
            returncode=1,
            stdout='',
            stderr='boom',
        )

        tarball_info = pkg.find_tarball_for_platform()

        self.assertEqual(tarball_info['tarball'], 'demo-VERSION.tar.gz')
