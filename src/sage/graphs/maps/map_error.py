# This file contain some custom error used in our project

class InvalidMapPermutationArgument(Exception):
    """
    This class represent an invalid argument error in MapPermutation.
    """

    def __init__(self):
        """
        Initialise the InvalidMapPermutationArgument

        EXAMPLES::

            sage: InvalidMapPermutationArgument()
            InvalidMapPermutationArgument('Invalid argument: The argument given must be Permutation or a non empty list of integers representing the permutation or a non empty list of tuples representing the cycles of the permutations or a positive integer.')

        """
        super().__init__("Invalid argument: The argument given must be Permutation or a non empty list of integers representing the permutation or a non empty list of tuples representing the cycles of the permutations or a positive integer.")


class InvalidSwapPermutationArgument(Exception):
    """
    This class represent an invalid argument error in SwapPermutation.
    """

    def __init__(self):
        """
        Initialise the InvalidSwapPermutationArgument

        EXAMPLES::

            sage: InvalidSwapPermutationArgument()
            InvalidSwapPermutationArgument('Invalid argument for swap permutation')

        """
        super().__init__("Invalid argument for swap permutation")


class NotImplemented(Exception):
    """

    This class represent a not implemented error.

    """

    def __init__(self, x):
        """

        Initialise the NotImplemented 

        INPUT:
        x the object on which the method isn't defined

        EXAMPLES::

            sage: NotImplemented(2)
            NotImplemented("This  method isn't implemented for the class <class 'sage.rings.integer.Integer'> ")


        """
        super().__init__(
            f"This  method isn't implemented for the class {x.__class__} ")


class TODO(Exception):
    """

    This class represent a TODO error

    """

    def __init__(self):
        """

        Initialise the TODO error.

        EXAMPLES::

            sage: TODO()
            TODO('Todo')

        """
        super().__init__(
            f"Todo")
