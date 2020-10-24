import os
import sys
import glob
from subprocess import check_call, CalledProcessError

try:
    SAGE_LOCAL = os.environ['SAGE_LOCAL']
    DOT_SAGE = os.environ['DOT_SAGE']
except KeyError:
    print('You need to run this script in a Sage shell ("sage -sh")')
    sys.exit(1)


TEXLIVE_INSTALL_UNIX_URL = \
    'http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz'

TEXLIVE_PROFILE_TEMPLATE = \
"""
selected_scheme scheme-minimal
TEXDIR {SAGE_LOCAL}/share/texlive
TEXMFLOCAL {SAGE_LOCAL}/share/texlive/texmf-local
TEXMFSYSCONFIG {SAGE_LOCAL}/share/texlive/texmf-config
TEXMFSYSVAR {SAGE_LOCAL}/share/texlive/texmf-var
TEXMFCONFIG {DOT_SAGE}/texlive/texmf-config
TEXMFHOME {DOT_SAGE}/texlive/texmf-home
TEXMFVAR {DOT_SAGE}/texlive/texmf-var
collection-basic 1
in_place 0
option_adjustrepo 1
option_autobackup 1
option_backupdir tlpkg/backups
option_desktop_integration 1
option_doc 1
option_file_assocs 1
option_fmt 1
option_letter 0
option_menu_integration 1
option_path 1
option_post_code 1
option_src 1
option_sys_bin {SAGE_LOCAL}/bin
option_sys_info {SAGE_LOCAL}/share/info
option_sys_man {SAGE_LOCAL}/share/man
option_w32_multi_user 1
option_write18_restricted 1
portable 0
"""


def have_texlive():
    try:
        check_call(['tlmgr', '--version'])
        return True
    except (OSError, CalledProcessError):
        return False


def download_install_script(tarball):
    try:
        from urllib import urlretrieve
    except ImportError:
        from urllib.request import urlretrieve
    urlretrieve(TEXLIVE_INSTALL_UNIX_URL, tarball)

    
def write_profile(filename):
    profile = TEXLIVE_PROFILE_TEMPLATE.format(
        SAGE_LOCAL=SAGE_LOCAL, DOT_SAGE=DOT_SAGE)
    with open(filename, 'w') as f:
        f.write(profile)
    print('TeXlive unattended install profile: {0}'.format(filename))

def first_dir(glob_pattern):
    """
    Return the first directory matching ``glob_pattern``
    """
    for dirent in glob.glob(glob_pattern):
        if os.path.isdir(dirent):
            return dirent
    raise RuntimeError('no directory found: {0}'.format(glob_pattern))
    
def install_texlive():
    import tempfile
    tmp_dir = tempfile.mkdtemp()
    tarball = os.path.join(tmp_dir, 'install-tl-unx.tar.gz')
    profile = os.path.join(tmp_dir, 'sage.profile')
    download_install_script(tarball)
    write_profile(profile)

    original_dir = os.getcwd()
    os.chdir(tmp_dir)
    check_call(['tar', '-x', '-z', '-f', tarball])
    install_dir = first_dir('install-tl-*')
    os.chdir(install_dir)

    check_call([
        './install-tl',
        '-profile',
        profile
    ])
    os.chdir(original_dir)
    

def install_packages():
    package_list_txt = os.path.join(
        os.path.dirname(__file__),
        'package-list.txt'
    )
    with open(package_list_txt) as f:
        package_list = f.read()
    packages = []
    for pkg in package_list.splitlines():
        pkg = pkg.strip()
        if len(pkg) == 0:
            continue
        if pkg.startswith('#'):
            continue
        packages.append(pkg)
    print('installing the following TeXlive packages:')
    for pkg in packages:
        print(' *  ' + pkg)
    check_call(['tlmgr', 'install'] + packages)
    check_call(['tlmgr', 'path', 'add'])
        

if __name__ == '__main__':
    if have_texlive():
        print('Using your own texlive install (see "tlmgr --version")')
    else:
        print('Performing minimal texlive install into {0}'
              .format(SAGE_LOCAL))
        install_texlive()
    install_packages()

