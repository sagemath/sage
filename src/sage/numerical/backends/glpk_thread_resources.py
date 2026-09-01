"""
Shared thread-local GLPK resource tracking

GLPK problem and graph objects use the same thread-local GLPK allocator state,
so pending releases from the LP and graph backends must be swept together.

TESTS:

Test the local resource collection helper with small fake resources::

    sage: from sage.numerical.backends.glpk_thread_resources import GLPKThreadResources
    sage: class RecordingResource:
    ....:     def __init__(self, name, events, released_error=None):
    ....:         self.name = name
    ....:         self.events = events
    ....:         self.released_error = released_error
    ....:     def _free_from_owner_thread(self, only_if_released=False):
    ....:         kind = "released" if only_if_released else "final"
    ....:         self.events.append((kind, self.name))
    ....:         if only_if_released and self.released_error is not None:
    ....:             raise self.released_error

Missing resources can be discarded, and a collection with no pending releases
does not touch registered resources::

    sage: resources = GLPKThreadResources()
    sage: events = []
    sage: a = RecordingResource("a", events)
    sage: b = RecordingResource("b", events)
    sage: resources.add(a)
    sage: resources.add(b)
    sage: resources.discard(a)
    sage: resources.discard(a)
    sage: resources.resources == [b]
    True
    sage: resources.collect_released()
    sage: events
    []

Successful sweeps free pending resources and clear the retry flag::

    sage: resources.has_released = True
    sage: resources.collect_released()
    sage: events
    [('released', 'b')]
    sage: resources.has_released
    False

The destructor frees all resources that are still registered and clears the
registry::

    sage: resources = GLPKThreadResources()
    sage: events = []
    sage: resources.add(RecordingResource("a", events))
    sage: resources.__del__()
    sage: events
    [('final', 'a')]
    sage: resources.resources
    []

Released resources unregister themselves while they are being swept in the
real backends, so collection iterates over a snapshot::

    sage: class SelfDiscardingResource:
    ....:     def __init__(self, resources, name, events):
    ....:         self.resources = resources
    ....:         self.name = name
    ....:         self.events = events
    ....:     def _free_from_owner_thread(self, only_if_released=False):
    ....:         if only_if_released:
    ....:             self.events.append(self.name)
    ....:             self.resources.discard(self)
    sage: resources = GLPKThreadResources()
    sage: events = []
    sage: resources.add(SelfDiscardingResource(resources, "a", events))
    sage: resources.add(SelfDiscardingResource(resources, "b", events))
    sage: resources.has_released = True
    sage: resources.collect_released()
    sage: events
    ['a', 'b']
    sage: resources.resources
    []
    sage: resources.has_released
    False

If a release sweep is interrupted, it leaves the retry flag armed so the
owner thread can try again later::

    sage: resources = GLPKThreadResources()
    sage: resources.add(RecordingResource("a", [], KeyboardInterrupt()))
    sage: resources.has_released = True
    sage: try:
    ....:     resources.collect_released()
    ....: except KeyboardInterrupt:
    ....:     pass
    sage: resources.has_released
    True

Resource tracking is thread-local::

    sage: from sage.numerical.backends import glpk_thread_resources as resources_module
    sage: import queue
    sage: import threading
    sage: main_resources = resources_module.glpk_thread_resources()
    sage: q = queue.Queue()
    sage: def get_thread_resource_identity():
    ....:     q.put(resources_module.glpk_thread_resources() is main_resources)
    sage: t = threading.Thread(target=get_thread_resource_identity)
    sage: t.start(); t.join()
    sage: q.get()
    False
    sage: resources_module.glpk_thread_resources() is main_resources
    True
"""

import threading

_glpk_thread_data = threading.local()


class GLPKThreadResources:
    def __init__(self):
        """
        Create an empty registry for one thread's GLPK resources.

        TESTS::

            sage: from sage.numerical.backends.glpk_thread_resources import GLPKThreadResources
            sage: resources = GLPKThreadResources()
            sage: (resources.resources, resources.has_released,
            ....:  resources.environment_cleaner)
            ([], False, None)
        """
        self.resources = []
        self.has_released = False
        self.environment_cleaner = None

    def add(self, resource):
        """
        Register ``resource`` for owner-thread cleanup.

        TESTS::

            sage: from sage.numerical.backends.glpk_thread_resources import GLPKThreadResources
            sage: resources = GLPKThreadResources()
            sage: resource = object(); resources.add(resource)
            sage: resources.resources == [resource]
            True
        """
        self.resources.append(resource)

    def discard(self, resource):
        """
        Discard ``resource`` if it is still registered.

        TESTS::

            sage: from sage.numerical.backends.glpk_thread_resources import GLPKThreadResources
            sage: resources = GLPKThreadResources()
            sage: resource = object(); resources.add(resource)
            sage: resources.discard(resource); resources.discard(resource)
            sage: resources.resources
            []
        """
        try:
            self.resources.remove(resource)
        except ValueError:
            pass

    def set_environment_cleaner(self, resource):
        """
        Retain one resource wrapper that can call ``glp_free_env`` at exit.

        The cleaner is kept even after its backend is released because GLPK's
        environment retains internal buffers until ``glp_free_env`` is called.

        TESTS::

            sage: from sage.numerical.backends.glpk_thread_resources import GLPKThreadResources
            sage: resources = GLPKThreadResources()
            sage: first = object()
            sage: resources.set_environment_cleaner(first)
            sage: resources.set_environment_cleaner(object())
            sage: resources.environment_cleaner is first
            True
        """
        if self.environment_cleaner is None:
            self.environment_cleaner = resource

    def collect_released(self):
        """
        Free resources released from another thread, if there are any.

        TESTS::

            sage: from unittest.mock import Mock
            sage: from sage.numerical.backends.glpk_thread_resources import GLPKThreadResources
            sage: resources = GLPKThreadResources()
            sage: resource = Mock(); resources.add(resource)
            sage: resources.has_released = True; resources.collect_released()
            sage: resource._free_from_owner_thread.assert_called_once_with(only_if_released=True)
            sage: resources.has_released
            False
        """
        if not self.has_released:
            return
        # Clear the flag *before* sweeping: if a cross-thread release lands
        # mid-sweep it sets ``has_released`` back to True, and we must not
        # overwrite that (otherwise that pending free is deferred all the way
        # to owner-thread exit).  If the sweep is interrupted, re-arm the flag
        # so the owner thread retries the remaining resources later.
        self.has_released = False
        try:
            for resource in list(self.resources):
                resource._free_from_owner_thread(only_if_released=True)
        except BaseException:
            self.has_released = True
            raise

    def _cleanup_from_owner_thread(self):
        """
        Free every registered resource while its GLPK environment is alive.

        TESTS::

            sage: from unittest.mock import Mock
            sage: from sage.numerical.backends.glpk_thread_resources import GLPKThreadResources
            sage: resources = GLPKThreadResources()
            sage: resource = Mock(); resources.add(resource)
            sage: resources._cleanup_from_owner_thread()
            sage: resource._free_from_owner_thread.assert_called_once_with()
            sage: resources.resources
            []
        """
        # Do not let one unexpected Python exception prevent the remaining
        # resources from being offered for cleanup.  Each Cython resource
        # independently verifies the CPython interpreter/thread-state IDs
        # before touching GLPK, so delayed finalization on another thread can
        # only leak, never free from the wrong allocator environment.
        resources = self.resources
        self.resources = []
        self.has_released = False
        for resource in resources:
            try:
                resource._free_from_owner_thread()
            except BaseException:
                pass

        cleaner = self.environment_cleaner
        self.environment_cleaner = None
        if cleaner is not None:
            try:
                cleaner._free_environment_from_owner_thread()
            except BaseException:
                pass

    def __del__(self):
        """
        Run owner-thread cleanup when this registry is finalized.

        TESTS::

            sage: from unittest.mock import Mock
            sage: from sage.numerical.backends.glpk_thread_resources import GLPKThreadResources
            sage: resources = GLPKThreadResources()
            sage: resource = Mock(); resources.add(resource)
            sage: resources.__del__()
            sage: resource._free_from_owner_thread.assert_called_once_with()
        """
        self._cleanup_from_owner_thread()


class _GLPKThreadCleanup:
    """
    Private thread-local cleanup guard.

    The public registry can be retained after its owner exits.  Keeping the
    cleanup hook in a separate, unexposed object ensures that such a reference
    does not postpone normal owner-thread cleanup.
    """

    def __init__(self, resources):
        """
        Store the public registry without exposing this cleanup guard.

        TESTS::

            sage: from sage.numerical.backends.glpk_thread_resources import GLPKThreadResources, _GLPKThreadCleanup
            sage: resources = GLPKThreadResources()
            sage: cleanup = _GLPKThreadCleanup(resources)
            sage: cleanup.resources is resources
            True
        """
        self.resources = resources

    def __del__(self):
        """
        Clean the public registry when the owning thread exits.

        TESTS::

            sage: from unittest.mock import Mock
            sage: from sage.numerical.backends.glpk_thread_resources import _GLPKThreadCleanup
            sage: resources = Mock(); cleanup = _GLPKThreadCleanup(resources)
            sage: cleanup.__del__()
            sage: resources._cleanup_from_owner_thread.assert_called_once_with()
        """
        self.resources._cleanup_from_owner_thread()


def glpk_thread_resources():
    """
    Return the GLPK resource registry for the current thread.

    TESTS::

        sage: from sage.numerical.backends.glpk_thread_resources import glpk_thread_resources
        sage: glpk_thread_resources() is glpk_thread_resources()
        True
    """
    try:
        return _glpk_thread_data.cleanup.resources
    except AttributeError:
        resources = GLPKThreadResources()
        _glpk_thread_data.cleanup = _GLPKThreadCleanup(resources)
        return resources
