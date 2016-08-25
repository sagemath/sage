


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

