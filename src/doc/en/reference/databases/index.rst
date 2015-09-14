Databases
=========

There are numerous specific mathematical databases either included
in Sage or available as optional packages. Also, Sage includes two
powerful general database packages.

Sage includes the ZOPE object oriented database ZODB, which
"is a Python object persistence system. It provides transparent object-oriented persistency."

Sage also includes the powerful relational database SQLite, along
with a Python interface to SQLite. SQlite is a small C library that
implements a self-contained, embeddable, zero-configuration SQL
database engine.


-  Transactions are atomic, consistent, isolated, and durable
   (ACID) even after system crashes and power failures.

-  Zero-configuration - no setup or administration needed.

-  Implements most of SQL92. (Features not supported)

-  A complete database is stored in a single disk file.

-  Database files can be freely shared between machines with
   different byte orders.

-  Supports databases up to 2 tebibytes (2^41 bytes) in size.

-  Strings and BLOBs up to 2 gibibytes (2^31 bytes) in size.

-  Small code footprint: less than 250KiB fully configured or less
   than 150KiB with optional features omitted.

-  Faster than popular client/server database engines for most
   common operations.

-  Simple, easy to use API.

-  TCL bindings included. Bindings for many other languages
   available separately.

-  Well-commented source code with over 95% test coverage.

-  Self-contained: no external dependencies.

-  Sources are in the public domain. Use for any purpose.

.. toctree::
   :maxdepth: 1

   sage/databases/cremona
   sage/databases/stein_watkins
   sage/databases/jones
   sage/databases/oeis
   sage/databases/sloane
   sage/databases/findstat
   sage/databases/conway
   sage/databases/odlyzko
   sage/databases/symbolic_data
   sage/databases/cunningham_tables
   sage/databases/db_class_polynomials
   sage/databases/db_modular_polynomials

.. include:: ../footer.txt
