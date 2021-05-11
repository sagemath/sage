# -*- coding: utf-8 -*-
"""
Set Up Logging

Logging can be customized using the ``SAGE_BOOTSTRAP`` environment
variable. It is a comma-separated list of ``key:value`` pairs. They
are not case sensitive. Valid pairs are:

* ``log:[level]``, where ``[level]`` is one of 

    * ``debug``
    * ``info``
    * ``warning``
    * ``critical``
    * ``error``

* ``interactive:true`` or ``interactive:false``, to override isatty detection.
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
import os


LOG_LEVELS = (
    'debug',
    'info',
    'warning',
    'critical',
    'error'
)


class Configuration(object):

    _initialized = False

    log = 'info'

    interactive = os.isatty(sys.stdout.fileno())
    
    def __init__(self):
        if not Configuration._initialized:
            Configuration._init_from_environ()
        if self.log not in LOG_LEVELS:
            raise ValueError('invalid log level: {0}'.format(self.log))
        assert isinstance(self.interactive, bool)

    @classmethod
    def _init_from_environ(cls):
        env = os.environ.get('SAGE_BOOTSTRAP', '').lower()
        for pair in env.split(','):
            if not pair.strip():
                continue
            key, value = pair.split(':', 1)
            key = key.strip()
            value = value.strip()
            if key == 'log':
                cls.log = value
            elif key == 'interactive':
                if value == 'true':
                    cls.interactive = True
                elif value == 'false':
                    cls.interactive = False
                else:
                    raise ValueError('interactive value must be "true" or "false", got "{0}"'
                                     .format(value))
            else:
                raise ValueError('unknown key: "{0}"'.format(key))
        cls._initialized = True

    def __repr__(self):
        return '\n'.join([
            'Configuration:',
            '  * log = {0}'.format(self.log),
            '  * interactive = {0}'.format(self.interactive)
        ])
