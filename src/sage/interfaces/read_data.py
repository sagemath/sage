"""
An interface to read data files
"""

###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2010 Paul Zimmermann
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later. The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

def read_data(f,t):
   r"""
   Read data from file 'f' and class 't' (one element per line),
   and returns a list of elements.

   INPUT:

       - 'f' - a file name
       - 't' - a class (objects will be coerced to that class)

   OUTPUT:

       a list of elements of class 't'.

   EXAMPLES::

       sage: indata = tmp_filename()
       sage: f = open(indata, "w")
       sage: _ = f.write("17\n42\n")
       sage: f.close()
       sage: l = read_data(indata, ZZ); l
       [17, 42]
       sage: f = open(indata, "w")
       sage: _ = f.write("1.234\n5.678\n")
       sage: f.close()
       sage: l = read_data(indata, RealField(17)); l
       [1.234, 5.678]
   """
   fp = open(f,"r")
   l = []
   while True:
      s = fp.readline()
      s = s.rstrip()
      if s == '':
         break
      l.append(t(s))
   fp.close()
   return l
