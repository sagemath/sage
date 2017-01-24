import time


def is_url(url):
    """
    Test whether argument is url
    """
    url = url.rstrip()
    if len(url.splitlines()) > 1:
        return False
    if url.find(' ') >= 0:
        return False
    return (
        url.startswith('http://') or
        url.startswith('https://') or
        url.startswith('ftp://')
    )


def retry(func, exc=Exception, tries=3, delay=1):
    """
    Call ``func()`` up to ``tries`` times, exiting only if the function
    returns without an exception.  If the function raises an exception on
    the final try that exception is raised.

    If given, ``exc`` can be either an `Exception` or a tuple of `Exception`s
    in which only those exceptions result in a retry, and all other exceptions
    are raised.  ``delay`` is the time in seconds between each retry, and
    doubles after each retry.
    """

    while True:
        try:
            return func()
        except exc:
            tries -= 1
            if tries == 0:
                raise

            time.sleep(delay)
            delay *= 2
