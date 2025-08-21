from .sage_transformer import DoctestTransformer


def test_consume_field_simple():
    lines = ["- param -- description"]
    dt = DoctestTransformer(lines)
    name, descs = dt._consume_field()
    assert name == "param"
    assert descs == ["description"]


def test_consume_field_encoded_param():
    lines = ["- ``param`` -- description"]
    dt = DoctestTransformer(lines)
    name, descs = dt._consume_field()
    assert name == "``param``"
    assert descs == ["description"]


def test_consume_field_multiple_dashes():
    lines = ["- param -- description with -- multiple dashes"]
    dt = DoctestTransformer(lines)
    name, descs = dt._consume_field()
    assert name == "param"
    assert descs == ["description with -- multiple dashes"]


def test_consume_field_multiple_quotes():
    lines = ["- ``param`` -- description with ``inline code``"]
    dt = DoctestTransformer(lines)
    name, descs = dt._consume_field()
    assert name == "``param``"
    assert descs == ["description with ``inline code``"]


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


def test_consume_fields_escaping_stars():
    lines = ["- ``*args`` -- star args"]
    dt = DoctestTransformer(lines)
    fields = dt._consume_fields()
    assert fields == [
        ("\*args", ["star args"]),
    ]

    lines = ["- ``**kwargs`` -- kw args"]
    dt = DoctestTransformer(lines)
    fields = dt._consume_fields()
    assert fields == [
        ("\*\*kwargs", ["kw args"]),
    ]


# def test_consume_field_invalid_line_raises():
#     lines = ["param -- missing leading dash"]
#     dt = DoctestTransformer(lines)
#     with pytest.raises(ValueError):
#         dt._consume_field()


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


def test_consume_fields_preserves_embedded_lists():
    lines = [
        "- param1 -- value1",
        "",
        "    - item1",
        "    - item2",
        "",
        "- param2 -- value2",
    ]
    dt = DoctestTransformer(lines)
    fields = dt._consume_fields()
    assert fields == [
        ("param1", ["value1", "", "- item1", "- item2", ""]),
        ("param2", ["value2"]),
    ]


def test_transform_handles_params_with_embedded_lists():
    lines = [
        "INPUT:",
        "",
        "- param1 -- value1",
        "",
        "    - item1",
        "    - item2",
        "",
        "- param2 -- value2",
    ]
    dt = DoctestTransformer(lines)
    result = dt.transform()
    assert result == [
        ":param param1: value1",
        "",
        "               - item1",
        "               - item2",
        "",
        ":param param2: value2",
        "",
    ]


def test_transform_input_section_followed_by_text_doesnot_convert_params():
    # TODO: We probably want to raise an exception instead in the future
    lines = ["INPUT:", "some text", "- param1 -- value1", ""]
    dt = DoctestTransformer(lines)
    result = dt.transform()
    assert result == ["", "some text", "- param1 -- value1", ""]


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
    lines = ["result line 1", "", "result line 2"]
    dt = DoctestTransformer(lines)
    out = dt._consume_returns_section()
    assert out == [("", ["result line 1", "", "result line 2"])]


def test_consume_returns_section_indented_dedents():
    lines = [
        "    first line",
        "        second deeper",
    ]
    dt = DoctestTransformer(lines)
    out = dt._consume_returns_section()
    assert out == [("", ["first line", "    second deeper"])]


def test_consume_returns_section_empty():
    lines = [
        "",
        "",
    ]
    dt = DoctestTransformer(lines)
    out = dt._consume_returns_section()
    assert out == []


def test_consume_returns_section_trailing_blank_lines_preserved():
    lines = [
        "desc",
        "",
        "",
    ]
    dt = DoctestTransformer(lines)
    out = dt._consume_returns_section()
    assert out == [("", ["desc", "", ""])]


def test_consume_returns_section_closed_by_other_section():
    lines = [
        "result line 1",
        "",
        "INPUT:",  # This should close the OUTPUT section
        "input line",
    ]
    dt = DoctestTransformer(lines)
    out = dt._consume_returns_section()
    assert out == [("", ["result line 1", ""])]
    # Header not consumed
    assert dt._lines.get(0) == "INPUT:"


def test_transform_returns_section():
    lines = [
        "OUTPUT:",
        "result",
    ]
    dt = DoctestTransformer(lines)
    result = dt.transform()
    assert result == [":returns: result", ""]


def test_transform_returns_section_closed_by_other_section():
    lines = [
        "OUTPUT:",
        "result",
        "ALGORITHM:",
    ]
    dt = DoctestTransformer(lines)
    result = dt.transform()
    assert result == [":returns: result", "", ".. rubric:: algorithm", ""]


def test_transform_examples_section():
    lines = [
        "EXAMPLES:",
        "",
        "    sage: foo1",
        "    sage: foo2",
        "",
    ]
    dt = DoctestTransformer(lines)
    result = dt.transform()
    assert result == [
        ".. rubric:: Examples",
        ".. code-block::",
        "",
        "    sage: foo1",
        "    sage: foo2",
        "",
    ]
