send2trash: Send file to trash natively under Mac OS X, Windows and Linux
=========================================================================

Description
-----------

Send file to trash natively under Mac OS X, Windows and Linux.

Send2Trash is a small package that sends files to the Trash (or Recycle
Bin) natively and on all platforms. On OS X, it uses native
FSMoveObjectToTrashSync Cocoa calls, on Windows, it uses native (and
ugly) SHFileOperation win32 calls. On other platforms, if PyGObject and
GIO are available, it will use this. Otherwise, it will fallback to its
own implementation of the trash specifications from freedesktop.org.

ctypes is used to access native libraries, so no compilation is
necessary.

Send2Trash supports Python 2.7 and up (Python 3 is supported).
