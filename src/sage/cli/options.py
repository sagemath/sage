from dataclasses import dataclass


@dataclass
class CliOptions:
    """
    A TypedDict for command-line interface options.

    Attributes:
        verbose (bool): Indicates whether verbose output is enabled.
    """

    verbose: bool = False
