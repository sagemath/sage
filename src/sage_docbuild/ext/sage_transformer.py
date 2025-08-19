"""Classes for docstring parsing and formatting."""
# Heavily inspired by https://github.com/sphinx-doc/sphinx/tree/master/sphinx/ext/napoleon

from __future__ import annotations

import collections
import contextlib
import inspect
import re
from typing import TYPE_CHECKING, Any

from sphinx.util import logging
from sphinx.util.typing import get_type_hints, stringify_annotation

if TYPE_CHECKING:
    from collections.abc import Callable
from typing import TypeVar

logger = logging.getLogger(__name__)

# " <section>: "
_section_regex = re.compile(r"^\s*([A-Z]+):+\s*$")
# - <name> -- <rest>
_field_regex = r"^\s*-\s*(?P<name>.+)\s*--\s*(?P<rest>.*)$"

_single_colon_regex = re.compile(r"(?<!:):(?!:)")
_xref_or_code_regex = re.compile(
    r"((?::(?:[a-zA-Z0-9]+[\-_+:.])*[a-zA-Z0-9]+:`.+?`)|"
    r"(?:``.+?``)|"
    r"(?::meta .+:.*)|"
    r"(?:`.+?\s*(?<!\x00)<.*?>`))"
)
_bullet_list_regex = re.compile(r"^(\*|\+|\-)(\s+\S|\s*$)")
_enumerated_list_regex = re.compile(
    r"^(?P<paren>\()?"
    r"(\d+|#|[ivxlcdm]+|[IVXLCDM]+|[a-zA-Z])"
    r"(?(paren)\)|\.)(\s+\S|\s*$)"
)


T = TypeVar("T")


class Deque(collections.deque[T]):
    """A subclass of deque that mimics ``pockets.iterators.modify_iter``.

    The `.Deque.get` and `.Deque.next` methods are added.
    """

    sentinel = object()

    def get(self, n: int) -> T:
        """Return the nth element of the stack, or ``self.sentinel`` if n is
        greater than the stack size.
        """
        return self[n] if n < len(self) else self.sentinel  # type: ignore[return-value]

    def next(self) -> T:
        if self:
            return super().popleft()
        else:
            raise StopIteration


class DoctestTransformer:
    """Convert Sage docstrings to reStructuredText."""

    _name_rgx = re.compile(
        r"^\s*((?::(?P<role>\S+):)?`(?P<name>~?[a-zA-Z0-9_.-]+)`|"
        r" (?P<name2>~?[a-zA-Z0-9_.-]+))\s*",
        re.VERBOSE,
    )

    def __init__(
        self,
        lines: list[str],
        what: str = "",
        name: str = "",
        obj: Any = None,
    ) -> None:
        if not what:
            if inspect.isclass(obj):
                what = "class"
            elif inspect.ismodule(obj):
                what = "module"
            elif callable(obj):
                what = "function"
            else:
                what = "object"

        self._what = what
        self._name = name
        self._obj = obj
        self._lines: Deque[str] = Deque(map(str.rstrip, lines))
        self._is_in_section = False
        self._section_indent = 0
        self._sections: dict[str, Callable[..., list[str]]] = {
            "input": self._parse_parameters_section,
            "output": self._parse_returns_section,
            #'examples': self._parse_examples_section,
            #'algorithm': partial(self._parse_admonition, 'algorithm'),
            #'references': partial(self._parse_admonition, 'references'),
            #'authors': partial(self._parse_admonition, 'authors'),
        }

    def _get_location(self) -> str | None:
        try:
            filepath = inspect.getfile(self._obj) if self._obj is not None else None
        except TypeError:
            filepath = None
        name = self._name

        if filepath is None and name is None:
            return None
        elif filepath is None:
            filepath = ""

        return f"{filepath}:docstring of {name}"

    def _consume_indented_block(self, indent: int = 1) -> list[str]:
        lines = []
        line = self._lines.get(0)
        while not self._is_section_break() and (
            not line or self._is_indented(line, indent)
        ):
            lines.append(self._lines.next())
            line = self._lines.get(0)
        return lines

    def _consume_contiguous(self) -> list[str]:
        lines = []
        while self._lines and self._lines.get(0) and not self._is_section_header():
            lines.append(self._lines.next())
        return lines

    def _consume_empty(self) -> list[str]:
        lines = []
        line = self._lines.get(0)
        while self._lines and not line:
            lines.append(self._lines.next())
            line = self._lines.get(0)
        return lines

    def _consume_field(self) -> tuple[str, list[str]]:
        """
        Consume a single field/parameter from the docstring.
        """
        line = self._lines.next()
        match = re.match(_field_regex, line)
        if match:
            _name = self._escape_args_and_kwargs(match.group("name").strip())
            _desc = match.group("rest").strip()
            indent = self._get_indent(line) + 1
            _descs = [_desc, *self._dedent(self._consume_indented_block(indent))]
            _descs = self.__class__(_descs).transform()
            return _name, _descs
        else:
            raise ValueError(f"Invalid field line: {line}")

    def _consume_fields(self, multiple: bool = False) -> list[tuple[str, list[str]]]:
        """
        Consume all fields/parameters from the docstring,
        until a section break is encountered.
        """
        self._consume_empty()
        fields: list[tuple[str, list[str]]] = []
        while not self._is_section_break():
            _name, _desc = self._consume_field()
            if multiple and _name:
                fields.extend(
                    (name.strip().strip(r"`").strip(), _desc)
                    for name in _name.split(",")
                )
            elif _name or _desc:
                fields.append((_name.strip().strip(r"`").strip(), _desc))
        return fields

    def _consume_inline_attribute(self) -> tuple[str, list[str]]:
        line = self._lines.next()
        _type, colon, _desc = self._partition_field_on_colon(line)
        if not colon or not _desc:
            _type, _desc = _desc, _type
            _desc += colon
        _descs = [_desc, *self._dedent(self._consume_to_end())]
        _descs = self.__class__(_descs).transform()
        return _type, _descs

    def _consume_returns_section(self) -> list[tuple[str, list[str]]]:
        lines = self._dedent(self._consume_to_next_section())
        if not lines:
            return []
        transformed = self.__class__(lines).transform()
        return [("", transformed)]

    def _consume_section_header(self) -> str:
        section = self._lines.next()
        return section.strip().strip(":").lower()

    def _consume_to_end(self) -> list[str]:
        lines = []
        while self._lines:
            lines.append(self._lines.next())
        return lines

    def _consume_to_next_section(self) -> list[str]:
        self._consume_empty()
        lines = []
        while not self._is_section_break():
            lines.append(self._lines.next())
        return lines + self._consume_empty()

    def _dedent(self, lines: list[str], full: bool = False) -> list[str]:
        if full:
            return [line.lstrip() for line in lines]
        else:
            min_indent = self._get_min_indent(lines)
            return [line[min_indent:] for line in lines]

    def _escape_args_and_kwargs(self, name: str) -> str:
        if name[:2] == "**":
            return r"\*\*" + name[2:]
        elif name[:1] == "*":
            return r"\*" + name[1:]
        else:
            return name

    def _fix_field_desc(self, desc: list[str]) -> list[str]:
        if self._is_list(desc):
            desc = ["", *desc]
        elif desc[0].endswith("::"):
            desc_block = desc[1:]
            indent = self._get_indent(desc[0])
            block_indent = self._get_initial_indent(desc_block)
            if block_indent > indent:
                desc = ["", *desc]
            else:
                desc = ["", desc[0], *self._indent(desc_block, 4)]
        return desc

    def _format_admonition(self, admonition: str, lines: list[str]) -> list[str]:
        lines = self._strip_empty(lines)
        if len(lines) == 1:
            return [f".. {admonition}:: {lines[0].strip()}", ""]
        elif lines:
            lines = self._indent(self._dedent(lines), 3)
            return [f".. {admonition}::", "", *lines, ""]
        else:
            return [f".. {admonition}::", ""]

    def _format_block(
        self,
        prefix: str,
        lines: list[str],
        padding: str | None = None,
    ) -> list[str]:
        if lines:
            if padding is None:
                padding = " " * len(prefix)
            result_lines = []
            for i, line in enumerate(lines):
                if i == 0:
                    result_lines.append((prefix + line).rstrip())
                elif line:
                    result_lines.append(padding + line)
                else:
                    result_lines.append("")
            return result_lines
        else:
            return [prefix]

    def _format_docutils_params(
        self,
        fields: list[tuple[str, list[str]]],
        field_role: str = "param",
    ) -> list[str]:
        lines = []
        for _name, _desc in fields:
            _desc = self._strip_empty(_desc)
            if any(_desc):
                _desc = self._fix_field_desc(_desc)
                field = f":{field_role} {_name}: "
                lines.extend(self._format_block(field, _desc))
            else:
                lines.append(f":{field_role} {_name}:")
        return [*lines, ""]

    def _format_field(self, _name: str, _desc: list[str]) -> list[str]:
        _desc = self._strip_empty(_desc)
        has_desc = any(_desc)
        separator = " -- " if has_desc else ""
        if _name:
            field = f"**{_name}**{separator}"
        else:
            field = ""

        if has_desc:
            _desc = self._fix_field_desc(_desc)
            if _desc[0]:
                return [field + _desc[0], *_desc[1:]]
            else:
                return [field, *_desc]
        else:
            return [field]

    def _format_fields(
        self,
        field_type: str,
        fields: list[tuple[str, list[str]]],
    ) -> list[str]:
        field_type = f":{field_type.strip()}:"
        padding = " " * len(field_type)
        multi = len(fields) > 1
        lines: list[str] = []
        for _name, _desc in fields:
            field = self._format_field(_name, _desc)
            if multi:
                if lines:
                    lines.extend(self._format_block(padding + " * ", field))
                else:
                    lines.extend(self._format_block(field_type + " * ", field))
            else:
                lines.extend(self._format_block(field_type + " ", field))
        if lines and lines[-1]:
            lines.append("")
        return lines

    def _get_current_indent(self, peek_ahead: int = 0) -> int:
        line = self._lines.get(peek_ahead)
        while line is not self._lines.sentinel:
            if line:
                return self._get_indent(line)
            peek_ahead += 1
            line = self._lines.get(peek_ahead)
        return 0

    def _get_indent(self, line: str) -> int:
        for i, s in enumerate(line):
            if not s.isspace():
                return i
        return len(line)

    def _get_initial_indent(self, lines: list[str]) -> int:
        for line in lines:
            if line:
                return self._get_indent(line)
        return 0

    def _get_min_indent(self, lines: list[str]) -> int:
        min_indent = None
        for line in lines:
            if line:
                indent = self._get_indent(line)
                if min_indent is None or indent < min_indent:
                    min_indent = indent
        return min_indent or 0

    def _indent(self, lines: list[str], n: int = 4) -> list[str]:
        return [(" " * n) + line for line in lines]

    def _is_indented(self, line: str, indent: int = 1) -> bool:
        for i, s in enumerate(line):
            if i >= indent:
                return True
            elif not s.isspace():
                return False
        return False

    def _is_list(self, lines: list[str]) -> bool:
        if not lines:
            return False
        if _bullet_list_regex.match(lines[0]):
            return True
        if _enumerated_list_regex.match(lines[0]):
            return True
        if len(lines) < 2 or lines[0].endswith("::"):
            return False
        indent = self._get_indent(lines[0])
        next_indent = indent
        for line in lines[1:]:
            if line:
                next_indent = self._get_indent(line)
                break
        return next_indent > indent

    def _is_section_header(self) -> bool:
        line = self._lines.get(0)
        match = _section_regex.match(line)
        if match:
            section = match.group(1).lower()
            return section in self._sections
        return False

    def _is_section_break(self) -> bool:
        line = self._lines.get(0)
        return (
            not self._lines
            or self._is_section_header()
            or (
                self._is_in_section
                and line is not self._lines.sentinel
                and not self._is_indented(line, self._section_indent)
            )
        )

    def transform(self) -> list[str]:
        """
        Return the parsed lines of the docstring in reStructuredText format.
        """
        _parsed_lines = self._consume_empty()

        if self._name and self._what in {"attribute", "data", "property"}:
            res: list[str] = []
            with contextlib.suppress(StopIteration):
                res = self._parse_attribute_docstring()

            _parsed_lines.extend(res)
            return _parsed_lines

        while self._lines:
            if self._is_section_header():
                try:
                    section = self._consume_section_header()
                    self._is_in_section = True
                    self._section_indent = self._get_current_indent()
                    lines = self._sections[section](section)
                finally:
                    self._is_in_section = False
                    self._section_indent = 0
            elif not _parsed_lines:
                lines = self._consume_contiguous() + self._consume_empty()
            else:
                lines = self._consume_to_next_section()
            _parsed_lines.extend(lines)
        return _parsed_lines

    def _parse_admonition(self, admonition: str, section: str) -> list[str]:
        lines = self._consume_to_next_section()
        return self._format_admonition(admonition, lines)

    def _parse_attribute_docstring(self) -> list[str]:
        _type, _desc = self._consume_inline_attribute()
        lines = self._format_field("", _desc)
        if _type:
            lines.extend(["", f":type: {_type}"])
        return lines

    def _parse_examples_section(self, section: str) -> list[str]:
        return self._parse_generic_section("Examples", False)

    def _parse_custom_generic_section(self, section: str) -> list[str]:
        # for now, no admonition for simple custom sections
        return self._parse_generic_section(section, False)

    def _parse_custom_params_style_section(self, section: str) -> list[str]:
        return self._format_fields(section, self._consume_fields())

    def _parse_custom_returns_style_section(self, section: str) -> list[str]:
        fields = self._consume_returns_section()
        return self._format_fields(section, fields)

    def _parse_generic_section(self, section: str, use_admonition: bool) -> list[str]:
        lines = self._strip_empty(self._consume_to_next_section())
        lines = self._dedent(lines)
        if use_admonition:
            header = f".. admonition:: {section}"
            lines = self._indent(lines, 3)
        else:
            header = f".. rubric:: {section}"
        if lines:
            return [header, "", *lines, ""]
        else:
            return [header, ""]

    def _parse_parameters_section(self, section: str) -> list[str]:
        fields = self._consume_fields(multiple=True)
        return self._format_docutils_params(fields)

    def _parse_returns_section(self, section: str) -> list[str]:
        fields = self._consume_returns_section()
        multi = len(fields) > 1
        lines: list[str] = []

        for _name, _desc in fields:
            field = self._format_field(_name, _desc)
            if multi:
                if lines:
                    lines.extend(self._format_block("          * ", field))
                else:
                    lines.extend(self._format_block(":returns: * ", field))
            elif any(field):  # only add :returns: if there's something to say
                lines.extend(self._format_block(":returns: ", field))
        if lines and lines[-1]:
            lines.append("")
        return lines

    def _partition_field_on_colon(self, line: str) -> tuple[str, str, str]:
        before_colon = []
        after_colon = []
        colon = ""
        found_colon = False
        for i, source in enumerate(_xref_or_code_regex.split(line)):
            if found_colon:
                after_colon.append(source)
            else:
                m = _single_colon_regex.search(source)
                if (i % 2) == 0 and m:
                    found_colon = True
                    colon = source[m.start() : m.end()]
                    before_colon.append(source[: m.start()])
                    after_colon.append(source[m.end() :])
                else:
                    before_colon.append(source)

        return "".join(before_colon).strip(), colon, "".join(after_colon).strip()

    def _strip_empty(self, lines: list[str]) -> list[str]:
        if lines:
            start = -1
            for i, line in enumerate(lines):
                if line:
                    start = i
                    break
            if start == -1:
                lines = []
            end = -1
            for i in reversed(range(len(lines))):
                line = lines[i]
                if line:
                    end = i
                    break
            if start > 0 or end + 1 < len(lines):
                lines = lines[start : end + 1]
        return lines

    def _lookup_annotation(self, _name: str) -> str:
        if False:  # True is default
            if self._what in {"module", "class", "exception"} and self._obj:
                # cache the class annotations
                if not hasattr(self, "_annotations"):
                    localns = getattr(self._config, "autodoc_type_aliases", {})
                    localns.update(
                        getattr(
                            self._config,
                            "napoleon_type_aliases",
                            {},
                        )
                        or {}
                    )
                    self._annotations = get_type_hints(self._obj, None, localns)
                if _name in self._annotations:
                    short_literals = getattr(
                        self._config, "python_display_short_literal_types", False
                    )
                    return stringify_annotation(
                        self._annotations[_name],
                        mode="fully-qualified-except-typing",
                        short_literals=short_literals,
                    )
        # No annotation found
        return ""
