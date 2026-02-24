# -*- coding: utf-8 -*-
"""
Access the List of Sage Download Mirrors
"""

#*****************************************************************************
#       Copyright (C) 2014-2016 Volker Braun <vbraun.name@gmail.com>
#                     2015      Jeroen Demeyer
#                     2023      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import contextlib
import logging
log = logging.getLogger()

from sage_bootstrap.compat import urllib, urlparse
from sage_bootstrap.env import SAGE_DISTFILES, SAGE_ROOT

from fcntl import flock, LOCK_SH, LOCK_EX
from errno import ENOLCK


def try_lock(fd, operation):
    """
    Try flock() but ignore ``ENOLCK`` errors, which could happen if the
    file system does not support locking.
    """
    try:
        flock(fd, operation)
    except IOError as e:
        if e.errno != ENOLCK:
            raise

        
class MirrorListException(RuntimeError):
    pass
        

class MirrorList(object):

    def __init__(self):
        self.sources = []
        upstream_d = os.path.join(SAGE_ROOT, '.upstream.d')
        for fname in sorted(os.listdir(upstream_d)):
            if '~' in fname or '#' in fname:
                # Ignore auto-save and backup files
                continue
            try:
                with open(os.path.join(upstream_d, fname), 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('#'):
                            continue
                        if not line:
                            continue
                        line = line.replace('${SAGE_ROOT}', SAGE_ROOT)
                        line = line.replace('${SAGE_DISTFILES}', SAGE_DISTFILES)
                        if '${SAGE_SERVER}' in line:
                            SAGE_SERVER = os.environ.get("SAGE_SERVER", "")
                            if not SAGE_SERVER:
                                continue
                            line = line.replace('${SAGE_SERVER}', SAGE_SERVER)
                        if line.endswith('mirror_list'):
                            cache_filename = os.path.join(SAGE_DISTFILES, line.rpartition('/')[2])
                            self.sources.append(MirrorList_from_url(line, cache_filename))
                        else:
                            self.sources.append([line])
            except IOError:
                # Silently ignore files that do not exist
                pass

    def __iter__(self):
        """
        Iterate through the list of mirrors.

        This is the main entry point into the mirror list. Every
        script should just use this function to try mirrors in order
        of preference. This will not just yield the official mirrors,
        but also urls for packages that are currently being tested.
        """
        for source in self.sources:
            for mirror in source:
                yield mirror


class MirrorList_from_url(object):
    
    MAXAGE = 24*60*60   # seconds

    def __init__(self, url, filename):
        self.url = url
        self.filename = filename
        self._mirrors = None

    @property
    def mirrors(self):
        if self._mirrors is not None:
            return self._mirrors

        try:
            self.mirrorfile = open(self.filename, 'r+t')
        except IOError:
            self.mirrorfile = open(self.filename, 'w+t')

        with self.mirrorfile:
            self.mirrorfd = self.mirrorfile.fileno()
            try_lock(self.mirrorfd, LOCK_SH)  # shared (read) lock
            if self._must_refresh():
                try_lock(self.mirrorfd, LOCK_EX)  # exclusive (write) lock
                # Maybe the mirror list file was updated by a different
                # process while we waited for the lock?  Check again.
                if self._must_refresh():
                    self._refresh()
            if self._mirrors is None:
                self._mirrors = self._load()

        return self._mirrors

    def _load(self, mirror_list=None):
        """
        Load and return `mirror_list` (defaults to the one on disk) as
        a list of strings
        """
        if mirror_list is None:
            try:
                self.mirrorfile.seek(0)
                mirror_list = self.mirrorfile.read()
            except IOError:
                log.critical('Failed to load the cached mirror list')
                return []
        if mirror_list == '':
            return []
        import ast
        try:
            return ast.literal_eval(mirror_list)
        except SyntaxError:
            log.critical('Downloaded mirror list has syntax error: {0}'.format(mirror_list))
            return []

    def _save(self):
        """
        Save the mirror list for (short-term) future  use.
        """
        self.mirrorfile.seek(0)
        self.mirrorfile.write(repr(self.mirrors))
        self.mirrorfile.truncate()
        self.mirrorfile.flush()

    def _port_of_mirror(self, mirror):
        if mirror.startswith('http://'):
            return 80
        if mirror.startswith('https://'):
            return 443
        if mirror.startswith('ftp://'):
            return 21
        # Sensible default (invalid mirror?)
        return 80

    def _rank_mirrors(self):
        """
        Sort the mirrors by speed, fastest being first

        This method is used by the YUM fastestmirror plugin
        """
        timed_mirrors = []
        import time
        import socket
        log.info('Searching fastest mirror')
        timeout = 1
        for mirror in self.mirrors:
            if not mirror.startswith('http'):
                log.debug('we currently can only handle http, got %s', mirror)
                continue
            port = self._port_of_mirror(mirror)
            mirror_hostname = urlparse.urlsplit(mirror).netloc
            time_before = time.time()
            try:
                sock = socket.create_connection((mirror_hostname, port), timeout)
                sock.close()
            except (IOError, socket.error, socket.timeout) as err:
                log.warning(str(err).strip() + ': ' + mirror)
                continue
            result = time.time() - time_before
            result_ms = int(1000 * result)
            log.info(str(result_ms).rjust(5) + 'ms: ' + mirror)
            timed_mirrors.append((result, mirror))
            timed_mirrors.sort()
            if len(timed_mirrors) >= 5 and timed_mirrors[4][0] < 0.3:
                # We don't need more than 5 decent mirrors
                break

        if len(timed_mirrors) == 0:
            # We cannot reach any mirror directly, most likely firewall issue
            if 'http_proxy' not in os.environ:
                log.error('Could not reach any mirror directly and no proxy set')
                raise MirrorListException('Failed to connect to any mirror, probably no internet connection')
            log.info('Cannot time mirrors via proxy, using default order')
        else:
            self._mirrors = [m[1] for m in timed_mirrors]
        log.info('Fastest mirror: ' + self.fastest)

    def _age(self):
        """
        Return the age of the cached mirror list in seconds
        """
        import time
        mtime = os.fstat(self.mirrorfd).st_mtime
        now = time.mktime(time.localtime())
        return now - mtime

    def _must_refresh(self):
        """
        Return whether we must download the mirror list.

        If and only if this method returns ``False`` is it admissible
        to use the cached mirror list.
        """
        if os.fstat(self.mirrorfd).st_size == 0:
            return True
        return self._age() > self.MAXAGE

    def _refresh(self):
        """
        Download and rank the mirror list.
        """
        log.info('Downloading the Sage mirror list')
        try:
            with contextlib.closing(urllib.urlopen(self.url)) as f:
                mirror_list = f.read().decode("ascii")
        except IOError:
            log.critical('Downloading the mirror list failed, using cached version')
        else:
            self._mirrors = self._load(mirror_list)
            self._rank_mirrors()
            self._save()

    def __iter__(self):
        """
        Iterate through the list of mirrors.

        This is the main entry point into the mirror list. Every
        script should just use this function to try mirrors in order
        of preference. This will not just yield the official mirrors,
        but also urls for packages that are currently being tested.
        """
        try:
            yield os.environ['SAGE_SERVER']
        except KeyError:
            pass
        for mirror in self.mirrors:
            if not mirror.endswith('/'):
                mirror += '/'
            yield mirror + '/'.join(['spkg', 'upstream', '${SPKG}'])

    @property
    def fastest(self):
        return next(iter(self))
