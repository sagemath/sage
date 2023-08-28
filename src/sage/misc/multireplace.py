# sage_setup: distribution = sagemath-categories
"multi_replace"

##########################################################################
#
# multi_replace function
#
# By Xavier Defrang.
#
# From the Python cookbook:
#
#  http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/81330
#
##########################################################################


import re


#
# The simplest, lambda-based implementation
#

def multiple_replace(dic, text):
    """
    Replace in 'text' all occurrences of any key in the given
    dictionary by its corresponding value.  Returns the new string.

    EXAMPLES::

        sage: from sage.misc.multireplace import multiple_replace
        sage: txt = "This monkey really likes the bananas."
        sage: dic = {'monkey': 'penguin', 'bananas': 'fish'}
        sage: multiple_replace(dic, txt)
        'This penguin really likes the fish.'
    """
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(re.escape(k) for k in dic))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dic[mo.string[mo.start():mo.end()]], text)
