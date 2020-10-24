# -*- coding: utf-8 -*-
"""
Download files from the internet
"""


#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import sys
import logging
log = logging.getLogger()

from sage_bootstrap.stdio import flush
from sage_bootstrap.compat import urllib


class ProgressBar(object):
    """
    Progress bar as urllib reporthook
    """

    def __init__(self, stream, length=70):
        self.length = length
        self.progress = 0
        self.stream = stream

    def start(self):
        flush()    # make sure to not interleave stdout/stderr
        self.stream.write('[')
        self.stream.flush()

    def __call__(self, chunks_so_far, chunk_size, total_size):
        if total_size == -1:  # we do not know size
            n = 0 if chunks_so_far == 0 else self.length // 2
        else:
            n = chunks_so_far * chunk_size * self.length // total_size
        if n > self.length: 
            # If there is a Content-Length, this will be sent as the last progress
            return
        # n ranges from 0 to length*total (exclude), so we'll print at most length dots
        if n >= self.progress:
            self.stream.write('.' * (n-self.progress))
            self.stream.flush()
        self.progress = n
        
    def stop(self):
        missing = '.' * (self.length - self.progress)
        self.stream.write(missing + ']\n')
        self.stream.flush()

    def error_stop(self):
        missing = 'x' * (self.length - self.progress)
        self.stream.write(missing + ']\n')
        self.stream.flush()
        

class DownloadError(IOError):
    pass

        
class Download(object):
    """
    Download URL
    
    Right now, only via HTTP

    This should work for FTP as well but, in fact, hangs on python <
    3.4, see http://bugs.python.org/issue16270

    INPUT:

    - ``url`` -- string. The URL to download.

    - ``destination`` -- string or ``None`` (default). The destination
      file name to save to. If not specified, the file is written to
      stdout.

    - ``progress`` -- boolean (default: ``True``). Whether to print a
      progress bar to stderr. For testing, this can also be a stream
      to which the progress bar is being sent.

    - ``ignore_errors`` -- boolean (default: ``False``). Catch network
      errors (a message is still being logged).
    """

    def __init__(self, url, destination=None, progress=True, ignore_errors=False):
        self.url = url
        self.destination = destination or '/dev/stdout'
        self.progress = (progress is not False)
        self.progress_stream = sys.stderr if isinstance(progress, bool) else progress
        self.ignore_errors = ignore_errors

    def http_error_default(self, url, fp, errcode, errmsg, headers):
        """
        Callback for the URLopener to raise an exception on HTTP errors
        """
        fp.close()
        raise DownloadError(errcode, errmsg, url)

    def start_progress_bar(self):
        if self.progress:
            self.progress_bar = ProgressBar(self.progress_stream)
            self.progress_bar.start()

    def success_progress_bar(self):
        if self.progress:
            self.progress_bar.stop()

    def error_progress_bar(self):
        if self.progress:
            self.progress_bar.error_stop()
    
    def run(self):
        opener = urllib.FancyURLopener()
        opener.http_error_default = self.http_error_default
        self.start_progress_bar()
        try:
            if self.progress:
                filename, info = opener.retrieve(
                    self.url, self.destination, self.progress_bar)
            else:
                filename, info = opener.retrieve(
                    self.url, self.destination)
        except IOError as error:
            self.error_progress_bar()
            log.error(error)
            if not self.ignore_errors:
                raise error
        self.success_progress_bar()
