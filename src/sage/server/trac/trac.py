r"""
Sage Trac Server

This module configures and launches a Trac server, if an optional
package (e.g., trac-x.y.z.spkg) is installed.
"""
from __future__ import print_function

import os
import sys

from sage.env import SAGE_LIB, SAGE_LOCAL
from sage.misc.viewer import browser

def trac_create_instance(directory = 'sage_trac', easy_setup = False):
    """
    Create a new Trac project instance if Trac is installed.

    INPUT:

    - ``directory`` - a string (default: 'sage_trac'); the name of the
      project directory

    - ``easy_setup`` - a bool (default: False); whether to use the
      project name 'Sage', the default database, enable the webadmin
      plugin, and enable the source browser for the 'sage' Mercurial
      repository.

    .. note::

        To access the webadmin panel, first create a user, give that
        user admin permissions, and log in as that user.  Create new
        accounts using htdigest (part of Apache)::

            cd <directory>/conf
            htdigest passwd <server_address> <username>

        Grant a user administrative privileges with::

            trac-admin <directory> add <username> TRAC_ADMIN
    """
    from sage.misc.sage_ostools import have_program

    if not have_program('trac-admin'):
        raise RuntimeError("trac is not installed")

    if easy_setup:
        cmd = 'trac-admin "%s" initenv "Sage" "sqlite:db/trac.db" "" ""' % directory
    else:
        cmd = 'trac-admin "%s" initenv' % directory

    e = os.system(cmd)
    if e:
        raise RuntimeError("Error creating trac environment.")

    if easy_setup:
        conf_name = os.path.abspath(os.path.join(directory, 'conf/trac.ini'))

        __import__('trac.config')
        Configuration = sys.modules['trac.config'].Configuration
        conf = Configuration(conf_name)

        conf.set('trac', 'repository_dir', SAGE_LIB)
        conf.set('trac', 'repository_type', 'hg')

        conf.set('components', 'tracext.hg.*', 'enabled')
        conf.set('components', 'webadmin.*', 'enabled')
        conf.set('hg', 'node_format', 'short')
        conf.set('hg', 'show_rev', 'yes')

        conf.save()


def trac(directory = 'sage_trac', port = 10000, address = 'localhost',
         open_viewer = False, auto_reload = False, easy_setup = False,
         options = ''):
    r"""
    Start a Trac server, creating a new project, if necessary.  An
    "optional" Sage package trac-x.y.z.spkg must already be installed.
    (Use optional_packages() to get the exact name.)

    INPUT:

    - ``directory`` - a string (default: 'sage_trac'); name of the
      project directory

    - ``port`` - an integer (default: 10000); the server's port number

    - ``address`` - a string (default: 'localhost'); the server's address

    - ``open_viewer`` - a bool (default: False); whether to open a
      browser at the server's address

    - ``auto_reload`` - a bool (default: False); whether the server
      should restart automatically when *its* sources are modified

    - ``easy_setup`` - a bool (default: False); **if** creating a new
      project, whether to enable optional plug-ins.  See
      :func:`trac_create_instance`.

    - ``options`` - a string (default: ''); command-line options to pass
      directly to tracd
    """
    from sage.misc.superseded import deprecation
    deprecation(16759, "This Sage Trac Server interface is deprecated, you can just run the Trac server outside of Sage")

    if not os.path.exists(directory):
        trac_create_instance(directory, easy_setup = easy_setup)

    url = 'http://%s:%s' % (address, port)
    if open_viewer:
        cmd = 'sleep 2 && %s ' % browser() + url + ' 1>&2 >/dev/null &'
        os.system(cmd)

    passwd = os.path.join(directory, 'conf/passwd')
    if not os.path.exists(passwd) or len(open(passwd).read()) < 2:
        print("*" * 80)
        print("Create new accounts with the htdigest command (part of Apache):")
        print("\nTo add a new user with name username:")
        print("    cd %s" % os.path.abspath(os.path.join(directory, 'conf')))
        print("    htdigest passwd %s <username>" % address)
        print("\nTo grant full admin permissions to a user:")
        print("    %s %s permission add <username> TRAC_ADMIN" % (os.path.join(SAGE_LOCAL, 'bin','trac-admin'), os.path.abspath(directory)))
        print("\nThen restart the trac server.")
        print("*" * 80)
        open(passwd, 'w').close()

    print("Open your web browser to " + url)
    print("Starting a Trac server...")

    if auto_reload:
        options += ' --auto-reload '

    cmd ='tracd %s --port %s --hostname %s --auth *,%s,%s "%s" ' % (options, port, address, passwd, address, directory)

    print(cmd)
    os.system(cmd)
