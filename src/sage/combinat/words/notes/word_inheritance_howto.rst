How to define a new class of words
==================================

To define a new class of words, one should inherit from ``Words_all`` (or any
class that inherits from it) and implement the following methods:

- ``self.alphabet()`` -- return the alphabet

- ``sortkey_letters(x)`` -- a function to compare letters in the alphabet;
  should return something that can be readily compared by
  Python. By default, the letters themselves are used.

By implementing the above methods, almost anything can be used as an
alphabet. The alphabet should have a ``cardinality`` method.

EXAMPLE 1: Using a list for the alphabet::

    self._alphabet = [0,1,2,3]

    def sortkey_letters(self, x):
        return self._alphabet.index(x)

EXAMPLE 2: Using a CombinatorialClass as an alphabet::

    self._alphabet = Partitions(3)

    def sortkey_letters(self, letter1):
        return self._alphabet.rank(letter1)

EXAMPLE 3: Integers::

    self._alphabet = ZZ

    def sortkey_letters(self, letter1):
        return letter1
