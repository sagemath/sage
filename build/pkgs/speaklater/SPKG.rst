speaklater: Lazy strings for Python
===================================

Description
-----------

Implements a lazy string for python useful for use with gettext

A module that provides lazy strings for translations. Basically you get
an object that appears to be a string but changes the value every time
the value is evaluated based on a callable you provide.

For example you can have a global lazy_gettext function that returns a
lazy string with the value of the current set language.
