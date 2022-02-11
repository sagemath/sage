"""
Pager for showing strings

Currently we just use the IPython pager.
If we want to use something else, we can just change this function.

Any code in sage that uses a pager should use this pager.
"""
# ---------------------------------------------------------------------------
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ---------------------------------------------------------------------------


def pager():
    import IPython.core.page
    return IPython.core.page.page
