# sage_setup: distribution = sagemath-repl
import doctest
import sys
import textwrap

import pytest

from sage.doctest.marked_output import MarkedOutput
from sage.doctest.parsing import (
    SageDocTestParser,
    SageOutputChecker,
    bitness_value,
    parse_optional_tags,
)

onlyLinux = pytest.mark.skipif(
    sys.platform != "linux",
    reason="No Linux system",
)
"""A decorator to specify that this function should only execute on Linux systems.
"""


def test_parse_optional_tags_known_bug_returns_bug():
    tags = parse_optional_tags("sage: # known bug")
    assert tags == {"bug": None}


def test_parse_optional_tags_known_bug_with_value_returns_bug_and_value():
    tags = parse_optional_tags("sage: # known bug: linux")
    assert tags == {"bug": "linux"}


def test_parse_optional_tags_known_bug_with_description_returns_bug():
    tags = parse_optional_tags("sage: # known bug, #34506")
    assert tags == {"bug": None}


def test_parse_optional_tags_known_bug_with_description_in_parentheses_returns_bug():
    tags = parse_optional_tags("sage: # known bug (#34506)")
    assert tags == {"bug": None}


def test_parse_optional_tags_known_bug_with_value_and_description_returns_bug_and_value():
    tags = parse_optional_tags("sage: # known bug: linux (#34506)")
    assert tags == {"bug": "linux"}


def test_parse_known_bug_returns_empty():
    parser = SageDocTestParser(("sage",))
    parsed = parser.parse("sage: x = int('1'*4301) # known bug")
    assert parsed == ["", ""]


def test_parse_known_bug_returns_code_if_requested():
    parser = SageDocTestParser(("sage", "bug"))
    parsed = parser.parse("sage: x = int('1'*4301) # known bug")
    assert len(parsed) == 3
    assert parsed[1].sage_source == "x = int('1'*4301) # known bug\n"


@onlyLinux
def test_parse_known_bug_returns_code_if_requested_even_on_affected_os():
    parser = SageDocTestParser(("sage", "bug"))
    parsed = parser.parse("sage: x = int('1'*4301) # known bug: macos")
    assert len(parsed) == 3
    assert parsed[1].sage_source == "x = int('1'*4301) # known bug: macos\n"


@onlyLinux
def test_parse_known_bug_returns_code_on_not_affected_os():
    parser = SageDocTestParser(("sage",))
    parsed = parser.parse("sage: x = int('1'*4301) # known bug: macos")
    assert len(parsed) == 3
    assert parsed[1].sage_source == "x = int('1'*4301) # known bug: macos\n"


@onlyLinux
def test_parse_known_bug_returns_empty_on_affected_os():
    parser = SageDocTestParser(("sage",))
    parsed = parser.parse("sage: x = int('1'*4301) # known bug: linux")
    assert parsed == ["", ""]


def test_parse_known_bug_with_description_returns_empty():
    parser = SageDocTestParser(("sage",))
    parsed = parser.parse("sage: x = int('1'*4301) # known bug, #34506")
    assert parsed == ["", ""]


def test_parse_bitness():
    parser = SageDocTestParser(("sage",))
    input = textwrap.dedent(
        """
        sage: sys.maxsize > (1 << 32)
        True # 64-bit
        False # 32-bit
        """
    )
    parsed = parser.parse(input)[1]
    assert isinstance(parsed, doctest.Example)
    assert parsed.want == MarkedOutput("True # 64-bit\nFalse # 32-bit\n").update(
        bitness_32="False \n", bitness_64="True \n"
    )


def test_check_output_bitness():
    checker = SageOutputChecker()
    expected = MarkedOutput("True # 64-bit\nFalse # 32-bit\n").update(
        bitness_32="False \n", bitness_64="True \n"
    )
    assert checker.check_output(expected, str(bitness_value == 64) + " \n", 0)
