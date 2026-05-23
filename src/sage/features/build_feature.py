r"""
Features that can be explicitly enabled or disabled at build time

These features are unique in that, if they are enabled or disabled at
build-time, then we usually do not want to detect them on-the-fly.
Instead we trust the value supplied (or detected) at build-time. There
is however an overrride to defer these checks to run-time (the classic
behavior) for use on binary distros or anywhere it is desirable to
enable/disable features without rebuilding sage.

This is an implementation of Option 3 in `Github discussion 41067
<https://github.com/sagemath/sage/discussions/41067>`__.
"""

from sage.features import Feature, FeatureTestResult

class BuildFeature(Feature):
    r"""
    A class for features that can be enabled or disabled at
    build-time.

    The current implementation refers to build features that are
    configurable in meson. For example::

        option(
          'foo',
          type: 'feature',
          value: 'auto',
          description: 'support for foo'
        )

    At build time, support for this "foo" will be automatically
    detected, and either enabled or disabled depending on whether or
    not its requirements are met. Alternatively, users may pass either
    ``-Dfoo=enabled`` or ``-Dfoo=disabled`` to explicitly enable or
    disable the feature. Features may be disabled regardless of
    whether or not they are installed, but usually features may only
    be enabled if their dependencies are present and usable.

    In any event, after ``meson setup``, support for "foo" is either
    enabled or disabled, and a boolean variable called something like
    ``foo_enabled`` is written to ``sage.config``. In your subclass,
    you should set the member variable ``_enabled_in_build`` to the
    value of that config variable.

    The :meth:`_is_present` method for this class will return the
    value of that config variable unless ``defer_feature_checks`` is
    set to ``True`` in ``sage.config``. If checks are deferred, the
    :meth:`_is_present` method will try to return the value of
    :meth:`is_present_at_runtime` instead. If your feature can be
    detected at run-time, you should implement that check in
    :meth:`is_present_at_runtime`. Otherwise, leave it unimplemented;
    and :meth:`_is_present` will return ``False``.

    EXAMPLES::

        sage: from sage.features.build_feature import BuildFeature
        sage: BuildFeature("foo")
        Feature('foo')

    """

    # Set this in subclasses.
    _enabled_in_build = None

    # Implement this method if your feature is detectable at run-time.
    # Your test should only return True if the feature meets Sage's
    # requirements; for example, if there are doctests for gzipped foo
    # data files hidden behind "needs foo", then you should ensure
    # that foo was compiled with (say) --enable-zlib in your check.
    #
    # def is_present_at_runtime(self):
    #     pass

    def is_runtime_detectable(self):
        r"""
        Return whether or not this feature can (and should) be
        detected at runtime.

        A feature is runtime detectable if both of the following hold:

        - Deferred feature checks have been enabled globally by
          passing ``-Ddefer_feature_checks=true`` to ``meson setup``.

        - An ``is_present_at_runtime`` method has been implemented for
          the feature.

        EXAMPLES:

        The method returns ``False`` if you have not implemented
        ``is_present_at_runtime``::

            sage: from sage.features.build_feature import BuildFeature
            sage: bf = BuildFeature("example")
            sage: bf.is_runtime_detectable()
            False

        """
        from sage.config import defer_feature_checks
        if not defer_feature_checks:
            return False
        elif hasattr(self, "is_present_at_runtime"):
            return True
        else:
            return False

    def _is_present(self):
        r"""
        Default presence check for build features.

        If this feature :meth:`is_runtime_detectable`, we return the
        result of that method. Otherwise, we use the value of
        ``self._enabled_in_build``.

        EXAMPLES:

        When feature checks are deferred, runtime-detectable features
        can be detected without ``self._enabled_in_build`` being set,
        but this will fail by surprise when they are un-deferred::

            sage: from sage.config import defer_feature_checks
            sage: from sage.features.build_feature import BuildFeature
            sage: bf = BuildFeature("example")
            sage: const_True = lambda s: True
            sage: bf.is_present_at_runtime = const_True.__get__(bf)
            sage: (not defer_feature_checks) or bf.is_present().is_present
            True

        """
        if self.is_runtime_detectable():
            return self.is_present_at_runtime()
        else:
            import sage.config
            # Wrap with bool() so that we can be lazy and use meson's
            # set10() rather than painstakingly writing "True" and
            # "False" to the config file.
            result =  bool(self._enabled_in_build)
            return FeatureTestResult(self, result)
