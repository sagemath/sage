"""
Inventory builder

A customized builder which only generates intersphinx "object.inv"
inventory files. The documentation files are not written.
"""
from __future__ import annotations

from os import path
from typing import Any, Iterable
from urllib.parse import quote

from sphinx.application import Sphinx
from sphinx.builders.dummy import DummyBuilder
from sphinx.util.inventory import InventoryFile

INVENTORY_FILENAME = "objects.inv"


class InventoryBuilder(DummyBuilder):
    """
    A customized builder which only generates intersphinx "object.inv"
    inventory files. The documentation files are not written.
    """
    name = "inventory"
    format = "inventory"
    epilog = "The inventory file is in %(outdir)s."

    def get_outdated_docs(self) -> Iterable[str]:
        """
        Return an iterable of output files that are outdated.
        """
        assert self.env is not None

        for docname in self.env.found_docs:
            if docname not in self.env.all_docs:
                yield docname
                continue
            targetname = path.join(self.outdir, INVENTORY_FILENAME)
            try:
                targetmtime = path.getmtime(targetname)
            except Exception:
                targetmtime = 0
            try:
                srcmtime = path.getmtime(self.env.doc2path(docname))
                if srcmtime > targetmtime:
                    yield docname
            except OSError:
                # source doesn't exist anymore
                pass

    def get_target_uri(self, docname: str, typ: str | None = None) -> str:
        """
        Return the target URI for a document name.
        """
        return quote(docname) + ".html"

    def finish(self) -> None:
        """
        Only write the inventory files.
        """
        assert self.env is not None

        InventoryFile.dump(
            path.join(self.outdir, INVENTORY_FILENAME), self.env, self
        )


def setup(app: Sphinx) -> dict[str, Any]:
    app.add_builder(InventoryBuilder)
    return {"parallel_read_safe": True}
