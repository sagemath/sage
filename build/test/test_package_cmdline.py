# -*- coding: utf-8 -*-
"""
Test sage-download-file commandline utility
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

import os
import unittest
import tempfile
import subprocess
import shutil
import logging

from sage_bootstrap.env import SAGE_DISTFILES
from sage_bootstrap.download.mirror_list import MIRRORLIST_FILENAME
from sage_bootstrap.package import Package
from test.capture import CapturedLog


log = logging.getLogger()


PATH = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    'bin'
)


EXECUTABLE = os.path.join(PATH, 'sage-package')


class SagePackageTestCase(unittest.TestCase):

    def run_command(self, *args, **kwds):
        env = dict(os.environ)
        env.update(kwds.get('env', {}))
        env['PATH'] = PATH + os.pathsep + env['PATH']
        kwds.update(
            stdin=None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
        )
        log.debug('running {}: {}'.format(args, kwds))
        proc = subprocess.Popen(args, **kwds)
        stdout, stderr = proc.communicate()
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        log.debug(u'stdout="{}", stderr="{}"'.format(stdout, stderr))
        rc = proc.returncode
        return (rc, stdout, stderr)

    def test_config(self):
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'config')
        # Prints nothing to stderr
        self.assertEqual(stderr, '')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertTrue(stdout.startswith('Configuration:\n'))

    def test_list(self):
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'list')
        # Prints nothing to stderr
        self.assertEqual(stderr, '')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertTrue('configure' in stdout.splitlines())

    def test_name(self):
        pkg = Package('configure')
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'name', pkg.tarball_filename)
        # Prints nothing to stderr
        self.assertEqual(stderr, '')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertEqual(stdout.rstrip(), 'configure')

    def test_tarball(self):
        pkg = Package('configure')
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'tarball', pkg.name)
        # Prints nothing to stderr
        self.assertEqual(stderr, '')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertEqual(stdout.rstrip(), pkg.tarball_filename)

    def test_apropos(self):
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'apropos', 'python')
        # Prints nothing to stderr
        self.assertEqual(stderr, '')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertTrue(stdout.startswith('Did you mean:'))

    def test_download(self):
        pkg = Package('configure')
        with CapturedLog() as _:
            pkg.tarball.download()
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'download', pkg.name)
        # Prints info to stderr
        self.assertTrue(stderr.startswith('Using cached file'))
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints filename to stdout
        self.assertEqual(stdout.rstrip(), pkg.tarball.upstream_fqn)

    def test_update(self):
        pkg = Package('configure')
        # The confball never has a patchlevel since we are upstream...
        self.assertEqual(pkg.patchlevel, -1)
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'update', pkg.name, pkg.version)
        # Prints nothing to stderr
        self.assertEqual(stderr, '')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints nothing to stdout
        self.assertEqual(stdout, '')

    def test_fix_checksum(self):
        pkg = Package('configure')
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'fix-checksum', 'configure')
        # Prints nothing to stderr
        self.assertEqual(stderr, '')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertEqual(
            stdout.rstrip(),
            'Checksum of {0} (tarball {1}) unchanged'.format(pkg.name, pkg.tarball_filename))

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
