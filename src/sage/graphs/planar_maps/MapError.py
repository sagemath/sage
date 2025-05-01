# This file contain some custom error used in our project

class InvalidMapPermutationArgument(Exception):
    def __init__(self):
        super().__init__("Invalid argument: The argument given must be Permutation or a non empty list of integers representing the permutation or a non empty list of tuples representing the cycles of the permutations or a positive integer.")


class InvalidSwapPermutationArgument(Exception):
    def __init__(self):
        super().__init__("Invalid argument for swap permutation")


class NotImplementedErrorError(Exception):
    def __init__(self, x):
        super().__init__(
            f"This  method isn't implemented for the class {x.__class__} ")

