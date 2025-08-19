import pytest

from .sage_transformer import DoctestTransformer


def test_consume_field_simple():
    lines = ["- param -- description"]
    dt = DoctestTransformer(lines)
    name, descs = dt._consume_field()
    assert name == "param"
    assert descs == ["description"]


def test_consume_field_multiple_params():
    lines = ["- param1, param2 -- description"]
    dt = DoctestTransformer(lines)
    name, descs = dt._consume_field()
    assert name == "param1, param2"
    assert descs == ["description"]


def test_consume_field_with_blank_lines_in_description():
    lines = [
        "- param -- first line",
        "    ",
        "    second line",
        "- another -- x",
    ]
    dt = DoctestTransformer(lines)
    name, descs = dt._consume_field()
    assert name == "param"
    assert descs == ["first line", "", "second line"]


def test_consume_field_escaping_stars():
    lines = ["- *args -- star args"]
    dt = DoctestTransformer(lines)
    name, descs = dt._consume_field()
    assert name == r"\*args"
    assert descs == ["star args"]

    lines = ["- **kwargs -- kw args"]
    dt = DoctestTransformer(lines)
    name, descs = dt._consume_field()
    assert name == r"\*\*kwargs"
    assert descs == ["kw args"]


def test_consume_field_invalid_line_raises():
    lines = ["param -- missing leading dash"]
    dt = DoctestTransformer(lines)
    with pytest.raises(ValueError):
        dt._consume_field()


def test_consume_fields_simple():
    lines = [
        "- param1 -- first param",
        "- param2 -- second param",
    ]
    dt = DoctestTransformer(lines)
    fields = dt._consume_fields()
    assert fields == [
        ("param1", ["first param"]),
        ("param2", ["second param"]),
    ]


def test_consume_fields_multiple_flag_splits_names():
    lines = ["- x, y -- coordinate values"]
    dt = DoctestTransformer(lines)
    fields = dt._consume_fields(multiple=True)
    assert fields == [
        ("x", ["coordinate values"]),
        ("y", ["coordinate values"]),
    ]
    # Both entries share the same description list object
    assert fields[0][1] is fields[1][1]


def test_consume_fields_multiple_flag_false_no_split():
    lines = ["- x, y -- coordinate values"]
    dt = DoctestTransformer(lines)
    fields = dt._consume_fields(multiple=False)
    assert fields == [
        ("x, y", ["coordinate values"]),
    ]


def test_consume_fields_trims_backticks_and_splits():
    lines = ["- ``x``, ``y`` -- desc"]
    dt = DoctestTransformer(lines)
    fields = dt._consume_fields(multiple=True)
    assert fields == [
        ("x", ["desc"]),
        ("y", ["desc"]),
    ]


def test_consume_fields_skips_leading_blank_lines():
    lines = [
        "",
        "   ",
        "- param -- value",
    ]
    dt = DoctestTransformer(lines)
    fields = dt._consume_fields()
    assert fields == [("param", ["value"])]


def test_consume_fields_stops_at_section_header():
    lines = [
        "- p1 -- first",
        "INPUT:",
        "- p2 -- second",  # should remain untouched
    ]
    dt = DoctestTransformer(lines)
    fields = dt._consume_fields()
    assert fields == [("p1", ["first"])]
    # Header not consumed
    assert dt._lines.get(0) == "INPUT:"


def test_consume_returns_section_basic():
    lines = [
        "result line 1",
        "",
        "OUTPUT:",
    ]
    dt = DoctestTransformer(lines)
    out = dt._consume_returns_section()
    assert out == [("", ["result line 1", ""])]
    # Header not consumed
    assert dt._lines.get(0) == "OUTPUT:"


def test_consume_returns_section_indented_dedents():
    lines = [
        "    first line",
        "        second deeper",
        "OUTPUT:",
    ]
    dt = DoctestTransformer(lines)
    out = dt._consume_returns_section()
    assert out == [("", ["first line", "    second deeper"])]
    assert dt._lines.get(0) == "OUTPUT:"


def test_consume_returns_section_empty():
    lines = [
        "",
        "",
        "OUTPUT:",
    ]
    dt = DoctestTransformer(lines)
    out = dt._consume_returns_section()
    assert out == []
    assert dt._lines.get(0) == "OUTPUT:"


def test_consume_returns_section_trailing_blank_lines_preserved():
    lines = [
        "desc",
        "",
        "",
        "OUTPUT:",
    ]
    dt = DoctestTransformer(lines)
    out = dt._consume_returns_section()
    assert out == [("", ["desc", "", ""])]
    assert dt._lines.get(0) == "OUTPUT:"
