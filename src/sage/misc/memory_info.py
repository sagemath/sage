"""
Information about available RAM/swap

There is no portable way to figure these out, nor should you generally
have to. But GAP currently needs to allocate a cache of fixed size
upon startup, and we would like a certain fraction of the swap address
space.

EXAMPLES::

    sage: from sage.misc.memory_info import MemoryInfo, MemoryInfo_proc
    sage: mem = MemoryInfo()
    sage: mem.total_ram()          # random output
    16708194304
    sage: mem.available_ram()      # random output
    1690738688
    sage: mem.total_swap()         # random output
    15728635904
    sage: mem.available_swap()     # random output
    15340593152
"""

########################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import subprocess
from sage.structure.sage_object import SageObject

memory_info_instance = None

def MemoryInfo():
    """
    Provide information about memory

    OUTPUT:

    A class that is encapsulates memory information. If no method for
    the particular host OS is provided, reasonable guesses are given.

    EXAMPLES::

        sage: from sage.misc.memory_info import MemoryInfo, MemoryInfo_proc
        sage: mem = MemoryInfo()
        sage: mem.total_ram()       # random output
        16708194304
        sage: mem.available_ram()   # random output
        1690738688
        sage: mem.total_swap()      # random output
        15728635904
        sage: mem.available_swap()  # random output
        15340593152
    """
    global memory_info_instance
    if memory_info_instance is not None:
        return memory_info_instance
    import platform
    system = platform.system()
    if memory_info_instance is None and \
            system != 'Darwin':
        try:
            memory_info_instance = MemoryInfo_proc()
        except OSError:
            pass
    if memory_info_instance is None and \
            system == 'Darwin':
        try:
            memory_info_instance = MemoryInfo_OSX()
        except OSError:
            pass
    if memory_info_instance is None:
        memory_info_instance = MemoryInfo_guess()
    return memory_info_instance


class MemoryInfo_base(SageObject):
    """
    Base class for memory info objects.
    """

    def rlimit_address_space(self):
        """
        Return ``RLIMIT_AS``.

        OUTPUT:

        Integer. The limit in bytes or `-1` if no limit is set or cannot
        be found out.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: mem = MemoryInfo()
            sage: mem.rlimit_address_space() in ZZ
            True
        """
        import resource
        try:
            return resource.getrlimit(resource.RLIMIT_AS)[1]
        except resource.error:
            return -1

    def virtual_memory_limit(self):
        """
        Return the upper limit for virtual memory usage

        This is the value set by ``ulimit -v`` at the command line or
        a practical limit if no limit is set.

        OUTPUT:

        Integer. The virtual memory limit in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: mem = MemoryInfo()
            sage: mem.virtual_memory_limit() > 0
            True
        """
        limit = self.rlimit_address_space()
        if limit >=0:
            return limit
        else:
            avail = self.total_swap() + self.total_ram()
            import platform
            if platform.architecture()[0] == '32bit':
                # 2GB is likely the single-process address space limit
                return min(avail, 2 * 1024**3)
            else:
                return avail


class MemoryInfo_proc(MemoryInfo_base):
    """
    Provide information from ``/proc/`` pseudo-filesystem on most UNIXes

    EXAMPLES::

        sage: from sage.misc.memory_info import MemoryInfo
        sage: mem = MemoryInfo()
        sage: mem.total_ram()   # random output
        16708194304
    """

    def __init__(self):
        try:
            self._parse_proc_meminfo()
        except (IOError, ValueError):
            raise OSError('/proc/meminfo is not available')

    def _parse_proc_meminfo(self):
        """
        Parse ``/proc/meminfo``

        OUTPUT:

        A dictionary. All sizes are in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo, MemoryInfo_proc
            sage: mem = MemoryInfo()
            sage: if isinstance(mem, MemoryInfo_proc):
            ...       info = mem._parse_proc_meminfo()
            ... else:
            ...       info = None
            sage: info   # random output
            {'available_ram': 1749782528,
             'total_swap': 15728635904,
             'free_swap': 15340572672,
             'total_ram': 16708194304}
            sage: keys = set(['available_ram', 'total_swap', 'free_swap', 'total_ram'])
            sage: (info is None) or keys.issubset(info.keys())
            True
        """
        kb = 1024
        result = dict()
        meminfo = open('/proc/meminfo', 'r')
        for line in meminfo.readlines():
            line = line.split()
            if line[0].startswith('MemTotal') and line[2] == 'kB':
                result['total_ram'] = int(line[1]) * kb
            if line[0].startswith('MemFree') and line[2] == 'kB':
                result['available_ram'] = int(line[1]) * kb
            if line[0].startswith('SwapTotal') and line[2] == 'kB':
                result['total_swap'] = int(line[1]) * kb
            if line[0].startswith('SwapFree') and line[2] == 'kB':
                result['free_swap'] = int(line[1]) * kb
            if line[0].startswith('Committed_AS') and line[2] == 'kB':
                result['Committed_AS'] = int(line[1]) * kb
        meminfo.close()
        required = set(['available_ram', 'total_swap', 'free_swap', 'total_ram'])
        if not required.issubset(result.keys()):
            raise OSError('failed to parse /proc/meminfo correctly')
        return result

    def total_ram(self):
        """
        Return the total RAM size

        OUTPUT:

        Integer. The RAM size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().total_ram() > 0
            True
        """
        return self._parse_proc_meminfo()['total_ram']

    def available_ram(self):
        """
        Return the available (free) RAM size

        OUTPUT:

        Integer. The free RAM size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().available_ram() > 0
            True
        """
        return self._parse_proc_meminfo()['available_ram']

    def total_swap(self):
        """
        Return the total swap size

        OUTPUT:

        Integer. The swap size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().total_swap() >= 0
            True
        """
        return self._parse_proc_meminfo()['total_swap']

    def available_swap(self):
        """
        Return the available (free) swap size

        OUTPUT:

        Integer. The free swap size in bytes, excluding reserved swap
        space. Can be negative if the system is overcommitting memory.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().available_swap() in ZZ
            True
        """
        info = self._parse_proc_meminfo()
        try:
            return info['total_swap'] - info['Committed_AS']
        except KeyError:
            return info['free_swap']


class MemoryInfo_OSX(MemoryInfo_base):
    """
    Memory info on OSX

    TESTS::

        sage: from sage.misc.memory_info import MemoryInfo_OSX
    """

    def __init__(self):
        self._maxage = 10   # cache result for 10 seconds
        self._age = -self._maxage
        try:
            self._parse_top()
        except (IOError, ValueError, subprocess.CalledProcessError, KeyError):
            raise OSError('failed to parse OSX "top" output')

    def _parse_top_output(self, meminfo):
        """
        Pick total and available memory out of the "top" output

        INPUT:

        - ``meminfo`` -- output of "top"

        OUTPUT:

        See :meth:`_parse_top`.

        TESTS::

            sage: from sage.misc.memory_info import MemoryInfo_OSX
            sage: m = MemoryInfo_OSX.__new__(MemoryInfo_OSX)
            sage: osx_ppc = 'PhysMem:  64.7M wired, 87.3M active,  14.1M inactive,  29.3M used,  21.8M free'
            sage: m._parse_top_output(osx_ppc)
            {'available_ram': 22858956, 'total_ram': 53582232}

            sage: osx_x86 = 'PhysMem: 8861M wired, 3574M active, 678M inactive, 13G used, 19G free.'
            sage: m._parse_top_output(osx_x86)
            {'available_ram': 20401094656L, 'total_ram': 34359738368L}    # 32-bit
            {'available_ram': 20401094656, 'total_ram': 34359738368}      # 64-bit
        """
        units = { 'K': 1024, 'M':1024**2, 'G':1024**3 }
        for line in meminfo.splitlines():
            if not line.startswith('PhysMem:'):
                continue
            line = line.split()
            if not line[-1].startswith('free') or not line[-3].startswith('used'):
                raise OSError('failed to parse PhysMem: line in "top" output')
            used_ram = line[-4]
            free_ram = line[-2]
            used_ram = int( float(used_ram[:-1]) * units[used_ram[-1]])
            free_ram = int( float(free_ram[:-1]) * units[free_ram[-1]])
            return { 'total_ram': used_ram + free_ram,
                     'available_ram': free_ram }
        raise OSError('failed to parse "top" output, no PhysMem: section')

    def _parse_top(self):
        """
        Parse ``top`` output

        OUTPUT:

        A dictionary. All sizes are in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo, MemoryInfo_OSX
            sage: mem = MemoryInfo()
            sage: if isinstance(mem, MemoryInfo_OSX):
            ...       info = mem._parse_top()
            ... else:
            ...       info = None
            sage: info   # random output
            {'available_ram': 1749782528,
             'total_ram': 16708194304}
            sage: keys = set(['available_ram', 'total_ram'])
            sage: (info is None) or (set(info.keys()) == keys)
            True
        """
        import time
        if (time.time()-self._age) < self._maxage:
            return self._parse_top_cache
        meminfo = subprocess.check_output(['top', '-l', '1'],
                                          stderr=subprocess.STDOUT)
        result = self._parse_top_output(meminfo)
        self._age = time.time()
        self._parse_top_cache = result
        return result

    def total_ram(self):
        """
        Return the total RAM size

        OUTPUT:

        Integer. The RAM size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().total_ram() > 0
            True
        """
        return self._parse_top()['total_ram']

    def available_ram(self):
        """
        Return the available (free) RAM size

        OUTPUT:

        Integer. The free RAM size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().available_ram() > 0
            True
        """
        return self._parse_top()['available_ram']

    def total_swap(self):
        """
        Return the total swap size

        The OSX swap file is growing dynamically, so we just return
        twice the total ram.

        OUTPUT:

        Integer. The swap size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().total_swap() >= 0
            True
        """
        return 2*self.total_ram()

    def available_swap(self):
        """
        Return the available (free) swap size

        The OSX swap file is growing dynamically, so we just return
        twice the available ram.

        OUTPUT:

        Integer. The free swap size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().available_swap() in ZZ
            True
        """
        return 2*self.available_ram()


class MemoryInfo_guess(MemoryInfo_base):
    """
    Guess memory as a fallback.

    TESTS::

        sage: from sage.misc.memory_info import MemoryInfo_guess
        sage: mem = MemoryInfo_guess()
        sage: mem.total_ram()
        4294967296       # 64-bit
        4294967296L      # 32-bit
    """
    def total_ram(self):
        """
        Return the total RAM size

        OUTPUT:

        Integer. The RAM size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().total_ram() > 0
            True
        """
        GB = 1024 * 1024 * 1024
        return 4*GB

    def available_ram(self):
        """
        Return the available (free) RAM size

        OUTPUT:

        Integer. The free RAM size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().available_ram() > 0
            True
        """
        return self.total_ram()

    def total_swap(self):
        """
        Return the total swap size

        OUTPUT:

        Integer. The swap size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().total_swap() >= 0
            True
        """
        GB = 1024 * 1024 * 1024
        return 4*GB

    def available_swap(self):
        """
        Return the available (free) swap size

        OUTPUT:

        Integer. The free swap size in bytes.

        EXAMPLES::

            sage: from sage.misc.memory_info import MemoryInfo
            sage: MemoryInfo().available_swap() in ZZ
            True
        """
        return self.total_swap()
