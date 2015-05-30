# -*- coding: utf-8 -*-
"""
    inventory builder
    ~~~~~~~~~~~~~~~~~

    A customized HTML builder which only generates intersphinx "object.inv"
    inventory files and pickle files. The documentation files are not written.
"""
from sphinx.builders.html import StandaloneHTMLBuilder
from sphinx.util.console import bold
from os import path
try:
    from hashlib import md5
except ImportError:
    # 2.4 compatibility
    from md5 import md5

class InventoryBuilder(StandaloneHTMLBuilder):
    """
    A customized HTML builder which only generates intersphinx "object.inv"
    inventory files and pickle files. The documentation files are not written.
    """
    name = "inventory"

    def get_outdated_docs(self):
        cfgdict = dict((name, self.config[name])
                       for (name, desc) in self.config.values.iteritems()
                       if desc[1] == 'html')
        self.config_hash = md5(unicode(cfgdict).encode('utf-8')).hexdigest()
        self.tags_hash = md5(unicode(sorted(self.tags)).encode('utf-8')) \
                .hexdigest()
        old_config_hash = old_tags_hash = ''
        try:
            fp = open(path.join(self.outdir, '.buildinfo'))
            try:
                version = fp.readline()
                if version.rstrip() != '# Sphinx build info version 1':
                    raise ValueError
                fp.readline()  # skip commentary
                cfg, old_config_hash = fp.readline().strip().split(': ')
                if cfg != 'config':
                    raise ValueError
                tag, old_tags_hash = fp.readline().strip().split(': ')
                if tag != 'tags':
                    raise ValueError
            finally:
                fp.close()
        except ValueError:
            self.warn('unsupported build info format in %r, building all' %
                      path.join(self.outdir, '.buildinfo'))
        except Exception:
            pass
        if old_config_hash != self.config_hash or \
               old_tags_hash != self.tags_hash:
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

    copy_image_files = removed_method_error
    copy_download_files = removed_method_error
    copy_static_files = removed_method_error
    handle_finish = removed_method_error

def setup(app):
    app.add_builder(InventoryBuilder)

