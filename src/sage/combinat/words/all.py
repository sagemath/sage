# sage_setup: distribution = sagemath-combinat
r"""
Combinatorics on words

**Main modules and their methods:**

- :ref:`sage.combinat.words.abstract_word`
- :ref:`sage.combinat.words.finite_word`
- :ref:`sage.combinat.words.infinite_word`
- :ref:`sage.combinat.words.alphabet`
- :ref:`sage.combinat.words.words`
- :ref:`sage.combinat.words.paths`
- :ref:`sage.combinat.words.morphism`
- :ref:`sage.combinat.words.shuffle_product`
- :ref:`sage.combinat.words.suffix_trees`

Main classes and functions meant to be used by the user:

    :func:`~sage.combinat.words.word.Word`,
    :class:`~sage.combinat.words.words.FiniteWords`,
    :class:`~sage.combinat.words.words.InfiniteWords`,
    :func:`~sage.combinat.words.words.Words`,
    :func:`~sage.combinat.words.alphabet.Alphabet`,
    :class:`~sage.combinat.words.morphism.WordMorphism`,
    :class:`~sage.combinat.words.paths.WordPaths`.

A list of common words can be accessed through ``words.<tab>`` and are listed in
the :ref:`words catalog <sage.combinat.words.word_generators>`.

**Internal representation of words:**

- :ref:`sage.combinat.words.word`
- :ref:`sage.combinat.words.word_char`
- :ref:`sage.combinat.words.word_datatypes`
- :ref:`sage.combinat.words.word_infinite_datatypes`

**Options:**

- :ref:`sage.combinat.words.word_options`

See :func:`~sage.combinat.words.word_options.WordOptions`.
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import

from sage.combinat.words.alphabet import Alphabet, build_alphabet
from sage.combinat.words.morphism import WordMorphism
lazy_import('sage.combinat.words.paths', 'WordPaths')
from sage.combinat.words.word import Word
from sage.combinat.words.word_options import WordOptions
from sage.combinat.words.word_generators import words
from sage.combinat.words.words import Words, FiniteWords, InfiniteWords
from sage.combinat.words.lyndon_word import LyndonWord, LyndonWords, StandardBracketedLyndonWords

del install_doc
del lazy_import
