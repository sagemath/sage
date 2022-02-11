# -*- coding: utf-8 -*-
"""
    inventory builder
    ~~~~~~~~~~~~~~~~~

    A customized HTML builder which only generates intersphinx "object.inv"
    inventory files and pickle files. The documentation files are not written.
"""
from sphinx.builders.html import StandaloneHTMLBuilder
from os import path
import shutil


class InventoryBuilder(StandaloneHTMLBuilder):
    """
    A customized HTML builder which only generates intersphinx "object.inv"
    inventory files and pickle files. The documentation files are not written.
    """
    name = "inventory"
    format = "inventory"
    epilog = "The inventory files are in %(outdir)s."

    def get_outdated_docs(self):
        from sphinx.builders.html import BuildInfo
        old = BuildInfo()
        try:
            with open(path.join(self.outdir, '.buildinfo')) as fp:
                old = BuildInfo.load(fp)
        except ValueError:
            self.warn('unsupported build info format in %r, building all' %
                      path.join(self.outdir, '.buildinfo'))
        except Exception:
            pass
        if self.build_info != old:
            for docname in self.env.found_docs:
                yield docname
            return

        if self.templates:
            template_mtime = self.templates.newest_template_mtime()
        else:
            template_mtime = 0
        for docname in self.env.found_docs:
            if docname not in self.env.all_docs:
                yield docname
                continue
            targetname = path.join(self.outdir, 'objects.inv')
            try:
                targetmtime = path.getmtime(targetname)
            except Exception:
                targetmtime = 0
            try:
                srcmtime = max(path.getmtime(self.env.doc2path(docname)),
                               template_mtime)
                if srcmtime > targetmtime:
                    yield docname
            except EnvironmentError:
                # source doesn't exist anymore
                pass

    def write_doc(self, docname, doctree):
        """
        Don't write any doc
        """

    def finish(self):
        """
        Only write the inventory files.
        """
        self.write_buildinfo()
        self.dump_inventory()

    def removed_method_error(self):
        """
        Raise an error if this method is called.

        This is just for making sure that some writer methods are indeed
        deactivated.
        """
        raise RuntimeError("This function shouldn't be called in \"%s\" builder"%(self.name))

    def cleanup(self):
        """
        Remove the '_static' directory.

        This directory is unnecessary for the inventory build, but it
        may be created by the graphviz extension. Its presence will
        break the docbuild later on, so remove it.
        """
        if path.isdir(path.join(self.outdir, '_static')):
            shutil.rmtree(path.join(self.outdir, '_static'))


    copy_image_files = removed_method_error
    copy_download_files = removed_method_error
    copy_static_files = removed_method_error
    handle_finish = removed_method_error


def setup(app):
    app.add_builder(InventoryBuilder)
    return {'parallel_read_safe': True} 
