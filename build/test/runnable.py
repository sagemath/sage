#!/usr/bin/env sage-bootstrap-python
# -*- coding: utf-8 -*-
"""
Utility to test running with different values for ``SAGE_BOOTSTRAP``
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


# This function's line numbers are in unit tests, try to not move it
# up or down. This is why it is up here at the beginning.
def print_log():
    import logging
    log = logging.getLogger()
    log.debug('This is the debug log level')
    log.info('This is the info log level')
    log.warning('This is the warning log level')
    log.critical('This is the critical log level')
    log.error('This is the error log level')
    print('This is printed')

# From here on the line number does not matter

import sys
import os
import json
import subprocess


def run_with(command, SAGE_BOOTSTRAP):
    env = dict(os.environ)
    env['SAGE_BOOTSTRAP'] = SAGE_BOOTSTRAP
    proc = subprocess.Popen(
        [sys.executable, __file__, command],
        env=env,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    out, err = proc.communicate()
    return out.decode('ascii'), err.decode('ascii')


def run_config_with(SAGE_BOOTSTRAP):
    out, err = run_with('print_config', SAGE_BOOTSTRAP)
    assert not err, err
    return json.loads(out)


def print_config():
    from sage_bootstrap.config import Configuration
    from sage_bootstrap.stdio import REAL_STDOUT, REAL_STDERR
    config = Configuration()
    result = dict(
        log=config.log,
        interactive=config.interactive,
        stdout='default stdout' if sys.stdout == REAL_STDOUT else str(type(sys.stdout)),
        stderr='default stderr' if sys.stderr == REAL_STDERR else str(type(sys.stderr)),
    )
    print(json.dumps(result))


def run_log_with(SAGE_BOOTSTRAP):
    return run_with('print_log', SAGE_BOOTSTRAP)


commands = dict(
    print_config=print_config,
    print_log=print_log,
)


if __name__ == '__main__':
    sys.path.insert(0,
        os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    import sage_bootstrap
    commands[sys.argv[1]]()
