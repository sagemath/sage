"""
Return types for PARI calls
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

class PariReturn(object):
    """
    This class represents the return value of a PARI call.
    """
    def __init__(self):
        self.name = "_ret"

    def __repr__(self):
        return self.ctype()

    def ctype(self):
        """
        Return the C type of the result of the PARI call.
        """
        raise NotImplementedError

    def assign_code(self, value):
        """
        Return code to assign the result of the PARI call in ``value``
        to the variable named ``self.name``.
        """
        s  = "        cdef {ctype} {name} = {value}\n"
        return s.format(ctype=self.ctype(), name=self.name, value=value)

    def return_code(self):
        """
        Return code to return from the Cython wrapper.
        """
        s  = "        pari_instance.clear_stack()\n"
        s += "        return {name}\n"
        return s.format(name=self.name)


class PariReturnGEN(PariReturn):
    def ctype(self):
        return "GEN"
    def return_code(self):
        s = "        return pari_instance.new_gen({name})\n"
        return s.format(name=self.name)

class PariReturnInt(PariReturn):
    def ctype(self):
        return "int"

class PariReturnLong(PariReturn):
    def ctype(self):
        return "long"

class PariReturnULong(PariReturn):
    def ctype(self):
        return "unsigned long"

class PariReturnVoid(PariReturn):
    def ctype(self):
        return "void"
    def assign_code(self, value):
        return "        {value}\n".format(value=value)
    def return_code(self):
        s = "        pari_instance.clear_stack()\n"
        return s


pari_ret_types = {
        '':  PariReturnGEN,
        'm': PariReturnGEN,
        'i': PariReturnInt,
        'l': PariReturnLong,
        'u': PariReturnULong,
        'v': PariReturnVoid,
        }
