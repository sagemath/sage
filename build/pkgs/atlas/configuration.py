######################################################################
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
######################################################################

from __future__ import print_function

import platform
import os
import sys
import shutil
import glob
import subprocess


# this dictionary will hold all configuration information
conf = dict()

# Here is a list of keys and (some of) their possible values. First,
# strings:
#
# system:  Linux, SunOS, Darwin, FreeBSD, CYGWIN
# release: (Version number of the operating system, output of uname -r)
# machine: i386, x86_64,   # Linux, Darwin, *BSD
#          sun4u, i86pc    # SunOS
# processor: i386, x86_64, powerpc, sparc
# bits:    32bit, 64bit
# fortran: g95, gfortran
# ld:      GNU, Solaris, Darwin

# The following pre-defined boolean values are determined from the
# strings. The keys are distinguished by a question mark at the
# end. If possible use these keys as they guard against typos.
#
# Linux?, Solaris?, Darwin?, FreeBSD?, CYGWIN?  # system
# OS_X_Lion?                                    # release
# Intel?, PPC?, SPARC? ARM?                     # processor
# 64bit?, 32bit?                                # bits
# fortran_g95?, fortran_GNU?                    # fortran
# linker_GNU?, linker_Solaris?, linker_Darwin?  # ld


######################################################################
### Functions
######################################################################

def try_run(command, ignore=False):
    """
    Try to execute ``command`` and return its output as string.

    Return ``None`` if command not found. Return ``None`` if execution
    fails and ``ignore=False`` (default). Otherwise, return output of
    ``command`` as string. Here, output always means the concatenation
    of stdout and stderr.
    """
    f = subprocess.Popen(command, shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = f.communicate()
    rc = f.wait()
    if (not ignore) and (rc!=0):
        return None
    # concatenate stdout and stderr
    return (result[0].strip() + result[1].strip()).decode('utf8')


def cp(source_pattern, destination):
    """
    Portable implementation of "cp -p"
    """
    for filename in glob.iglob(source_pattern):
        print('Copying', filename, 'to', destination)
        shutil.copy2(filename, destination)


def ln(source, destination):
    """
    Portable implementation of "ln -sf"
    """
    if os.path.exists(destination):
        os.remove(destination)
    print('Linking', source, 'to', destination)
    os.symlink(source, destination)


def is_executable_file(fpath):
    """
    Check that `fpath` is a file and executable (ie not a directory with +x set)
    """
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def which(program):
    """
    Portable implementation of "which"
    """
    fpath, fname = os.path.split(program)
    if fpath != '':
        return program
    for path in os.environ['PATH'].split(os.pathsep):
        f = os.path.join(path, program)
        if is_executable_file(f):
            return f
    raise ValueError('The program '+program+' is not in the $PATH.')

class edit_in_place(object):
    """
    Edit a file in-place by search/replacing regular expressions.
    """
    def __init__(self, filename):
        self.filename = filename
        f = open(filename, 'r')
        self.data = f.read()
        f.close()

    def replace(self, find_regex, subs, count=0):
        self.data = glob.re.sub(find_regex, subs, self.data, count)
        return self

    def close(self):
        f = open(self.filename, 'w')
        f.write(self.data)
        f.close()


######################################################################
### uname
######################################################################

try:
    conf['system'] = os.environ['UNAME']
except KeyError:
    conf['system'] = platform.system()

try:
    conf['release'] = subprocess.check_output(['uname', '-r']).decode('utf8').strip()
except subprocess.CalledProcessError:
    conf['release'] = ""

conf['Linux?'  ] = (conf['system'] == 'Linux')
conf['Solaris?'] = (conf['system'] == 'SunOS')
conf['Darwin?' ] = (conf['system'] == 'Darwin')
conf['OS_X_Lion?'] = (conf['Darwin?'] and conf['release'].startswith('11.'))
conf['FreeBSD?'] = (conf['system'] == 'FreeBSD')
conf['CYGWIN?' ] = (conf['system'] == 'CYGWIN')

conf['machine'] = platform.machine()

conf['processor'] = platform.processor()

conf['ARM?'] = (platform.machine().startswith('arm'))
conf['Intel?'] = (platform.machine() in ('i386', 'i486', 'i586', 'i686', 'x86_64',
                                         'AMD64', 'i86pc'))
conf['IA64?']   = ((platform.processor() == 'ia64')
                   or (platform.machine() == 'ia64'))
conf['PPC?']   = (platform.processor() == 'powerpc')
conf['SPARC?'] = (platform.processor() == 'sparc')

conf['generic_binary?'] = (os.environ.get('SAGE_FAT_BINARY', 'no') == 'yes')

conf['user'] = os.environ.get('ATLAS_CONFIGURE', '')

######################################################################
### bit width
######################################################################

conf['bits'] = platform.architecture()[0]

conf['64bit?'] = (conf['bits'] == '64bit')
conf['32bit?'] = (conf['bits'] != '64bit')


######################################################################
### fortran compiler
######################################################################

fortran_version = try_run('$FC --version')
if fortran_version is None:
    print('Cannot execute fortran compiler ($FC)!')
    sys.exit(3)
if 'G95' in fortran_version:
    conf['fortran'] = 'g95'
elif 'GNU Fortran' in fortran_version:
    conf['fortran'] = 'gfortran'
else:
    print('Unknown fortran compiler version: '+fortran_version)
    conf['fortran'] = None


conf['fortran_g95?'] = (conf['fortran'] == 'g95')
conf['fortran_GNU?'] = (conf['fortran'] == 'gfortran')


if conf['fortran_g95?']:
    g95_dir = glob.glob(os.environ['SAGE_LOCAL']+'/lib/gcc-lib/*/*')
    assert len(g95_dir)==1, 'Could not find G95 dir: '+str(g95_dir)
    conf['fortran_g95_dir'] = g95_dir[0]

######################################################################
### linker
######################################################################

if conf['Solaris?']:
    # Note: Solaris linker does not accept --version
    # ignore error code as 'ld -V' likes to error on some Solaris versions
    ld_version = try_run('ld  -V', ignore=True)
else:
    ld_version = try_run('ld  -v')

if ld_version is None:
    print('Cannot execute ld!')
    sys.exit(3)
if 'GNU' in ld_version:
    conf['ld'] = 'GNU'
elif 'Solaris' in ld_version:
    conf['ld'] = 'Solaris'
elif 'Apple' in ld_version:
    conf['ld'] = 'Darwin'
else:
    print('Unknown linker: '+ld_version)
    conf['ld'] = None

conf['linker_GNU?'] = (conf['ld'] == 'GNU')
conf['linker_Solaris?'] = (conf['ld'] == 'Solaris')
conf['linker_Darwin?'] = (conf['ld'] == 'Darwin')


if conf['Solaris?'] and conf['linker_GNU?']:
    print("WARNING: You are using the GNU linker from 'binutils'")
    print("Generally it is considered better to use the Sun linker")
    print("but Sage has been built on Solaris using the GNU linker")
    print("although that was a very old version of Sage, which")
    print("never passes all the Sage test-suite.")



######################################################################
### paths, files, and environment variables
######################################################################

conf['SPKG_DIR'] = os.getcwd()
conf['SAGE_LOCAL'] = os.environ['SAGE_LOCAL']


######################################################################
### The end: print configuration
######################################################################

print("Configuration:")
for k in sorted(conf, key=str.lower):
    print('    {}: {}'.format(k, conf[k]))
