# -*- coding: utf-8 -*-
"""
Utility to manage lists of packages
"""

# ****************************************************************************
#       Copyright (C) 2016      Volker Braun <vbraun.name@gmail.com>
#                     2020-2024 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import logging
from sage_bootstrap.package import Package

log = logging.getLogger()


class PackageClass(object):

    def __init__(self, *package_names_or_classes, **filters):
        self.__names = set()
        exclude = filters.pop('exclude', ())
        include_dependencies = filters.pop('include_dependencies', False)
        exclude_dependencies = filters.pop('exclude_dependencies', False)
        filenames = filters.pop('has_files', [])
        no_filenames = filters.pop('no_files', [])
        excluded = []
        for package_names in exclude:
            excluded.extend(package_names)
        if filters:
            raise ValueError('filter not supported')

        def included_in_filter(pkg):
            if not all(any(pkg.has_file(filename)
                           for filename in filename_disjunction.split('|'))
                       for filename_disjunction in filenames):
                return False
            return not any(any(pkg.has_file(filename)
                               for filename in no_filename_disjunction.split('|'))
                           for no_filename_disjunction in no_filenames)

        for package_name_or_class in package_names_or_classes:
            if package_name_or_class == ':all:':
                self._init_all(predicate=included_in_filter)
            elif package_name_or_class == ':standard:':
                self._init_standard(predicate=included_in_filter)
            elif package_name_or_class == ':optional:':
                self._init_optional(predicate=included_in_filter)
            elif package_name_or_class == ':experimental:':
                self._init_experimental(predicate=included_in_filter)
            elif any(package_name_or_class.startswith(prefix)
                     for prefix in ["pkg:", "pypi/", "generic"]):
                self.__names.add(Package(package_name_or_class).name)
            else:
                if ':' in package_name_or_class:
                    raise ValueError('a colon may only appear in a PURL such as '
                                     'pkg:pypi/DISTRIBUTION-NAME '
                                     'and in designators of package types, '
                                     'which must be one of '
                                     ':all:, :standard:, :optional:, or :experimental:'
                                     'got {}'.format(package_name_or_class))
                if '-' in package_name_or_class:
                    raise ValueError('dashes may only appear in a PURL such as '
                                     'pkg:pypi/DISTRIBUTION-NAME; '
                                     'SPKG names use underscores')
                self.__names.add(package_name_or_class)

        def include_recursive_dependencies(names, package_name):
            if package_name in names:
                return
            try:
                pkg = Package(package_name)
            except FileNotFoundError:
                # Silently ignore unknown packages,
                # substitutions such as $(BLAS) $(PYTHON),
                # and optional dependencies of the form $(find-string ...).
                return
            names.add(package_name)
            for dependency in pkg.dependencies:
                include_recursive_dependencies(names, dependency)

        if include_dependencies:
            package_names = set()
            for name in self.__names:
                include_recursive_dependencies(package_names, name)
            self.__names = package_names

        def exclude_recursive_dependencies(names, package_name):
            try:
                pkg = Package(package_name)
            except FileNotFoundError:
                return
            for dependency in pkg.dependencies:
                names.discard(dependency)
                exclude_recursive_dependencies(names, dependency)

        if exclude_dependencies:
            for name in list(self.__names):
                exclude_recursive_dependencies(self.__names, name)

        self.__names.difference_update(excluded)

    @property
    def names(self):
        return sorted(self.__names)

    def _init_all(self, predicate):
        self.__names.update(pkg.name for pkg in Package.all() if predicate(pkg))

    def _init_standard(self, predicate):
        self.__names.update(pkg.name for pkg in Package.all() if pkg.type == 'standard' and predicate(pkg))

    def _init_optional(self, predicate):
        self.__names.update(pkg.name for pkg in Package.all() if pkg.type == 'optional' and predicate(pkg))

    def _init_experimental(self, predicate):
        self.__names.update(pkg.name for pkg in Package.all() if pkg.type == 'experimental' and predicate(pkg))

    def apply(self, function, *args, **kwds):
        for package_name in self.names:
            function(package_name, *args, **kwds)
