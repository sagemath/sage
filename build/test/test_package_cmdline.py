# -*- coding: utf-8 -*-
"""
Test sage-download-file commandline utility
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
import unittest
import tempfile
import subprocess
import shutil

from sage_bootstrap.env import SAGE_DISTFILES
from sage_bootstrap.download.mirror_list import MIRRORLIST_FILENAME
from sage_bootstrap.util import is_url
from sage_bootstrap.package import Package
from test.capture import CapturedLog


EXECUTABLE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    'bin',
    'sage-package',
)


class SagePackageTestCase(unittest.TestCase):

    def run_command(self, *args, **kwds):
        kwds.update(
            stdin=None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        proc = subprocess.Popen(args, **kwds)
        stdout, stderr = proc.communicate()
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        rc = proc.returncode
        return (rc, stdout, stderr)
    
    def test_config(self):
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'config')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertTrue(stdout.startswith('Configuration:\n'))
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_list(self):
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'list')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertTrue('configure' in stdout.splitlines())
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_name(self):
        pkg = Package('configure')
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'name', pkg.tarball_filename)
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertEqual(stdout.rstrip(), 'configure')
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_tarball(self):
        pkg = Package('configure')
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'tarball', pkg.name)
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertEqual(stdout.rstrip(), pkg.tarball_filename)
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_apropos(self):
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'apropos', 'python')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertTrue(stdout.startswith('Did you mean:'))
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_download(self):
        pkg = Package('configure')
        with CapturedLog() as log:
            pkg.tarball.download()
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'download', pkg.name)
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints filename to stdout
        self.assertEqual(stdout.rstrip(), pkg.tarball.upstream_fqn)
        # Prints info to stderr
        self.assertTrue(stderr.startswith('Using cached file'))

    def test_update(self):
        pkg = Package('configure')
        # The confball never has a patchlevel since we are upstream...
        self.assertEqual(pkg.patchlevel, -1)
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'update', pkg.name, pkg.version)
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints nothing to stdout
        self.assertEqual(stdout, '')
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_fix_checksum(self):
        pkg = Package('configure')
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'fix-checksum', 'configure')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertEqual(stdout.rstrip(), 'Checksum of {0} unchanged'.format(pkg.tarball_filename))
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_create(self):
        tmp = tempfile.mkdtemp()
        with open(os.path.join(tmp, 'configure.ac'), 'w+') as f:
            f.write('test')
        os.mkdir(os.path.join(tmp, 'build'))
        os.mkdir(os.path.join(tmp, 'build', 'pkgs'))
        os.mkdir(os.path.join(tmp, 'upstream'))
        with open(os.path.join(tmp, 'upstream', 'Foo-13.5.tgz'), 'w+') as f:
            f.write('tarball content')
        rc, stdout, stderr = self.run_command(
            EXECUTABLE,
            'create', 'foo',
            '--version', '13.5',
            '--tarball', 'Foo-VERSION.tgz',
            '--type', 'standard',
            env=dict(SAGE_ROOT=tmp)
        )
        self.assertEqual(rc, 0)
        with open(os.path.join(tmp, 'build', 'pkgs', 'foo', 'package-version.txt')) as f:
            self.assertEqual(f.read(), '13.5\n')
        with open(os.path.join(tmp, 'build', 'pkgs', 'foo', 'type')) as f:
            self.assertEqual(f.read(), 'standard\n')
        with open(os.path.join(tmp, 'build', 'pkgs', 'foo', 'checksums.ini')) as f:
            self.assertEqual(
                f.read(),
                'tarball=Foo-VERSION.tgz\n' + 
                'sha1=15d0e36e27c69bc758231f8e9add837f40a40cd0\n' +
                'md5=bc62fed5e35f31aeea2af95c00473d4d\n' +
                'cksum=1436769867\n'
            )
        shutil.rmtree(tmp)
        
