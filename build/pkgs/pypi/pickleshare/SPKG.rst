pickleshare: A 'shelve' like datastore with concurrency support
===============================================================

Description
-----------

PickleShare - a small 'shelve' like datastore with concurrency support

Like shelve, a PickleShareDB object acts like a normal dictionary.
Unlike shelve, many processes can access the database simultaneously.
Changing a value in database is immediately visible to other processes
accessing the same database.

Concurrency is possible because the values are stored in separate files.
Hence the "database" is a directory where all files are governed by
PickleShare.
