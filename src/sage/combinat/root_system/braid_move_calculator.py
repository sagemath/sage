# sage.doctest: needs sage.libs.gap
"""
Braid Move Calculator

AUTHORS:

- Dinakar Muthiah (2014-06-03): initial version
"""
# ****************************************************************************
#       Copyright (C) 2014 Dinakar Muthiah <muthiah at ualberta.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method


class BraidMoveCalculator:
    """
    Helper class to compute braid moves.
    """

    def __init__(self, coxeter_group):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.braid_move_calculator import BraidMoveCalculator
            sage: W = CoxeterGroup(['C',3])
            sage: B = BraidMoveCalculator(W)
            sage: TestSuite(B).run(skip='_test_pickling')
        """
        self.coxeter_matrix = coxeter_group.coxeter_matrix()

    def _apply_put_in_front_recur_step(self, k, input_word, coxeter_matrix_entry):
        """
        Recurrence step for :meth:`put_in_front`.

        EXAMPLES::

            sage: from sage.combinat.root_system.braid_move_calculator import BraidMoveCalculator
            sage: W = CoxeterGroup(['C',3])
            sage: B = BraidMoveCalculator(W)
            sage: B.put_in_front(2, (3, 2, 3, 1, 2, 3, 1, 2, 1)) # indirect doctest
            ((3, 2, 3, 1, 2, 3, 1, 2, 1),
             (3, 2, 3, 1, 2, 1, 3, 2, 1),
             (3, 2, 3, 2, 1, 2, 3, 2, 1),
             (2, 3, 2, 3, 1, 2, 3, 2, 1))
        """
        i = input_word[0]

        def partial_braid_word(length, swap=False, i=i, k=k):
            if swap:
                i, k = k, i
            running_braid_word = [i, k] * (length // 2)
            if length % 2:
                running_braid_word.append(i)
            return tuple(running_braid_word)

        current_last_word = input_word
        current_first_letter = k
        output_word_list = [current_last_word]
        for counter in range(1, coxeter_matrix_entry):
            current_word_list = self.put_in_front(current_first_letter, current_last_word[1:])
            output_word_list += [partial_braid_word(counter) + word
                                 for word in current_word_list[1:]]
            if current_first_letter == k:
                current_first_letter = i
            else:
                current_first_letter = k
            current_last_word = current_word_list[-1]
        if i != k:
            output_word_list += [partial_braid_word(coxeter_matrix_entry, swap=True) +
                                 current_last_word[1:]]
        return tuple(output_word_list)

    def put_in_front(self, k, input_word):
        """
        Return a list of reduced words starting with ``input_word``
        and ending with a reduced word whose first letter  is ``k``.

        There still remains an issue with 0 indices.

        EXAMPLES::

            sage: from sage.combinat.root_system.braid_move_calculator import BraidMoveCalculator
            sage: W = CoxeterGroup(['C',3])
            sage: B = BraidMoveCalculator(W)
            sage: B.put_in_front(2, (3, 2, 3, 1, 2, 3, 1, 2, 1))
            ((3, 2, 3, 1, 2, 3, 1, 2, 1),
             (3, 2, 3, 1, 2, 1, 3, 2, 1),
             (3, 2, 3, 2, 1, 2, 3, 2, 1),
             (2, 3, 2, 3, 1, 2, 3, 2, 1))
            sage: B.put_in_front(1, (3, 2, 3, 1, 2, 3, 1, 2, 1))
            ((3, 2, 3, 1, 2, 3, 1, 2, 1),
             (3, 2, 1, 3, 2, 3, 1, 2, 1),
             (3, 2, 1, 3, 2, 3, 2, 1, 2),
             (3, 2, 1, 2, 3, 2, 3, 1, 2),
             (3, 1, 2, 1, 3, 2, 3, 1, 2),
             (1, 3, 2, 1, 3, 2, 3, 1, 2))
            sage: B.put_in_front(1, (1, 3, 2, 3, 2, 1, 2, 3, 2))
            ((1, 3, 2, 3, 2, 1, 2, 3, 2),)
        """
        i = input_word[0]
        if i == 0 or k == 0:  # Is this for affine types? - Travis
            raise NotImplementedError
        entry = self.coxeter_matrix[i, k]
        return self._apply_put_in_front_recur_step(k, input_word, entry)

    @cached_method
    def chain_of_reduced_words(self, start_word, end_word):
        """
        Compute the chain of reduced words from ``stard_word``
        to ``end_word``.

        INPUT:

        - ``start_word``, ``end_word`` -- two reduced expressions
          for the long word

        EXAMPLES::

            sage: from sage.combinat.root_system.braid_move_calculator import BraidMoveCalculator
            sage: W = CoxeterGroup(['A',5])
            sage: B = BraidMoveCalculator(W)
            sage: B.chain_of_reduced_words((1,2,1,3,2,1,4,3,2,1,5,4,3,2,1), # not tested
            ....:                          (5,4,5,3,4,5,2,3,4,5,1,2,3,4,5))
        """
        if start_word == end_word:
            return (start_word,)
        k = end_word[0]
        first_word_list = self.put_in_front(k, start_word)
        first_last_word = first_word_list[-1]
        return (first_word_list[:-1] +
                tuple([(k,) + word for word in
                       self.chain_of_reduced_words(first_last_word[1:],
                                                   end_word[1:])]))
