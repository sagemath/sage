"""Define some custom errors used in our library."""


class InvalidMapPermutationArgumentError(Exception):
    """
    This class represent an invalid argument error in MapPermutation.
    """

    def __init__(self):
        """
        Initialise the InvalidMapPermutationArgument

        EXAMPLES::

            sage: from sage.graphs.maps.map_error import InvalidMapPermutationArgumentError
            sage: InvalidMapPermutationArgumentError()
            InvalidMapPermutationArgumentError('Invalid argument: The argument given must be Permutation or a non empty list of integers representing the permutation or a non empty list of tuples representing the cycles of the permutations or a positive integer.')
        """
        super().__init__("Invalid argument: The argument given must be Permutation or a non empty list of integers representing the permutation or a non empty list of tuples representing the cycles of the permutations or a positive integer.")


class InvalidSwapPermutationArgumentError(Exception):
    """
    This class represent an invalid argument error in SwapPermutation.
    """

    def __init__(self):
        """
        Initialise the InvalidSwapPermutationArgumentError

        EXAMPLES::

            sage: from sage.graphs.maps.map_error import InvalidSwapPermutationArgumentError
            sage: InvalidSwapPermutationArgumentError()
            InvalidSwapPermutationArgumentError('Invalid argument for swap permutation')
        """
        super().__init__("Invalid argument for swap permutation")


class NotImplementedError(Exception):
    """

    This class represent a not implemented error.

    """

    def __init__(self, x):
        """

        Initialise the NotImplementedError.

        INPUT:

        -``x`` ; the object on which the method isn't defined

        EXAMPLES::

            sage: from sage.graphs.maps.map_error import NotImplementedError
            sage: NotImplementedError(2)
            NotImplementedError("This  method isn't implemented for the class <class 'sage.rings.integer.Integer'> ")
        """
        super().__init__(
            f"This  method isn't implemented for the class {x.__class__} ")


class TODOError(Exception):
    """

    This class represent a TODO error

    """

    def __init__(self):
        """

        Initialise the TODO error.

        EXAMPLES::

            sage: from sage.graphs.maps.map_error import TODOError
            sage: TODOError()
            TODOError('Todo')
        """
        super().__init__(
            "Todo")
