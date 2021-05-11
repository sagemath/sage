scandir: Fast file system iteration for Python
==============================================

Description
-----------

scandir, a better directory iterator and faster os.walk()

scandir() is a directory iteration function like os.listdir(), except
that instead of returning a list of bare filenames, it yields DirEntry
objects that include file type and stat information along with the name.
Using scandir() increases the speed of os.walk() by 2-20 times
(depending on the platform and file system) by avoiding unnecessary
calls to os.stat() in most cases.
