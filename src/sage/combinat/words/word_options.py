# sage_setup: distribution = sagemath-combinat
r"""
User-customizable options for words
"""
# ****************************************************************************
#       Copyright (C) 2009 Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import copy
from sage.rings.integer import Integer

word_options = {'identifier': 'word: ',
                'display': 'string',
                'truncate': True,
                'truncate_length': 40,
                'letter_separator': ',',
                'cache': True,
                'old_repr': False}


def WordOptions(**kwargs):
    """
    Set the global options for elements of the word class.
    The defaults are for words to be displayed in list notation.

    INPUT:

    - ``display`` -- 'string' (default), or 'list', words are displayed in
      string or list notation
    - ``truncate`` -- boolean (default: ``True``); whether to truncate the string
      output of long words (see truncate_length below)
    - ``truncate_length`` -- integer (default: 40); if the length of the word
      is greater than this integer, then the word is truncated
    - ``letter_separator`` -- string (default: ``','``); if the string
      representation of letters have length greater than 1, then
      the letters are separated by this string in the string
      representation of the word

    If no parameters are set, then the function returns a copy of the
    options dictionary.

    EXAMPLES::

        sage: w = Word([2,1,3,12])
        sage: u = Word("abba")
        sage: WordOptions(display='list')
        sage: w
        word: [2, 1, 3, 12]
        sage: u
        word: ['a', 'b', 'b', 'a']
        sage: WordOptions(display='string')
        sage: w
        word: 2,1,3,12
        sage: u
        word: abba
    """
    global word_options
    if kwargs == {}:
        return copy.copy(word_options)
    if 'display' in kwargs:
        if kwargs['display'] not in ['list', 'string']:
            raise ValueError("display must be either 'list' or 'string'")
        else:
            word_options['display'] = kwargs['display']
    elif 'truncate' in kwargs:
        if not isinstance(kwargs['truncate'], bool):
            raise ValueError("truncate must be ``True`` or False")
        else:
            word_options['truncate'] = kwargs['truncate']
    elif 'truncate_length' in kwargs:
        if not isinstance(kwargs['truncate_length'], (int, Integer)) or kwargs['truncate_length'] <= 0:
            raise ValueError("truncate_length must be a positive integer")
        else:
            word_options['truncate_length'] = kwargs['truncate_length']
    elif 'letter_separator' in kwargs:
        if not isinstance(kwargs['letter_separator'], str):
            raise ValueError("letter_separator must be a string")
        else:
            word_options['letter_separator'] = kwargs['letter_separator']
    elif 'identifier' in kwargs:
        if not isinstance(kwargs['identifier'], str):
            raise ValueError("identifier must be a string")
        else:
            word_options['identifier'] = kwargs['identifier']
    elif 'cache' in kwargs:
        if not isinstance(kwargs['cache'], bool):
            raise ValueError("cache must be ``True`` or False")
        else:
            word_options['cache'] = kwargs['cache']
    elif 'old_repr' in kwargs:
        if not isinstance(kwargs['old_repr'], bool):
            raise ValueError("old_repr must be ``True`` or False")
        else:
            word_options['old_repr'] = kwargs['old_repr']
