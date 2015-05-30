"""
Verify that a given file probably defines a valid element class.
"""

import os

def element_verify(file, module_element=False, ring_element=False, monoid_element=False):
    try:
        r = open(file).read()
    except OSError as msg:
        print msg
        print "Invalid file!"
        return False

    sagex = (file[-4:] == '.pyx')

    basefile = os.path.split(file)[1]

    def msg(s):
        print "%s: %s"%(basefile, s)

    if not 'class' in r:
        return

    if sagex:
        if not "def __richcmp__" in r:
            msg("The following method *must* be in your file.")
            msg("def __richcmp__(left, right, int op)")

        for x in ['_richcmp(self, right, int op)']:
            if (' ' + x) in r:
                msg("The following forbidden method is in your file but must *not* be.")
                msg("             " + x)

        if not '_cmp_c_impl(left,' in r:
            msg("WARNING: You should define '_cmp_c_impl(left,'")
            msg("And be sure to also define 'def __richcmp__(left, right, int op)'")

        if module_element:
            for x in ['_add_', '_sub_', '_neg_c_impl']:
                if not (('Element ' + x) in r):
                    msg("WARNING: You should define the cdef'd method '%s'"%x)

        if monoid_element or ring_element:
            if not 'Element _mul_' in r:
                msg("WARNING: You should define the cdef'd method '_mul_'")

        if ring_element:
            if not 'Element _div_' in r:
                msg("WARNING: You should define the cdef'd method '_div_'")

    else:
        # pure python class
        if not 'def __cmp__(' in r:
            msg("WARNING: You should define 'def __cmp__(left, right)'")
            msg("which may assume the parents of left and right are identical.")

        if module_element:
            for x in ['_add_', '_sub_', '_neg_']:
                if not (('def ' + x) in r):
                    msg("WARNING: You should define the method '%s'"%x)

        if monoid_element or ring_element:
            if not 'def _mul_' in r:
                msg("WARNING: You should define the method '_mul_'")

        if ring_element:
            if not 'def _div_' in r:
                msg("WARNING: You should define the method '_div_'")






