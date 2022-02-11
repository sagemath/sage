"""
Message delivery

Various interfaces to messaging services. Currently:

- ``pushover`` - a platform for sending and receiving push notifications

is supported.

AUTHORS:

- Martin Albrecht (2012) - initial implementation
"""

import http.client as httplib
from urllib.parse import urlencode
from ssl import SSLContext

pushover_defaults = {"token": "Eql67F14ohOZJ0AtEBJJU7FiLAk8wK"}


def pushover(message, **kwds):
    """
    Send a push notification with ``message`` to ``user`` using https://pushover.net/.

    Pushover is a platform for sending and receiving push notifications. On the server side, it
    provides an HTTP API for queueing messages to deliver to devices. On the device side, iOS and
    Android clients receive those push notifications, show them to the user, and store them for
    offline viewing.

    An account on https://pushover.net is required and the Pushover app must be installed on your
    phone for this function to be able to deliver messages to you.

    INPUT:

      - ``message`` - your message

      - ``user`` - the user key (not e-mail address) of your user (or you), viewable when logged
        into the Pushover dashboard. (default: ``None``)

      - ``device`` - your user's device identifier to send the message directly to that device,
        rather than all of the user's devices (default: ``None``)

      - ``title`` - your message's title, otherwise uses your app's name (default: ``None``)

      - ``url`` - a supplementary URL to show with your message (default: ``None``)

      - ``url_title`` - a title for your supplementary URL (default: ``None``)

      - ``priority`` - set to 1 to display as high-priority and bypass quiet hours, or -1 to always
        send as a quiet notification (default: ``0``)

      - ``timestamp`` - set to a unix timestamp to have your message show with a particular time,
        rather than now (default: ``None``)

      - ``sound`` - set to the name of one of the sounds supported by device clients to override the
        user's default sound choice (default: ``None``)

      - ``token`` - your application's API token (default: Sage's default App token)

    EXAMPLES::

        sage: import sage.misc.messaging
        sage: sage.misc.messaging.pushover("Hi, how are you?", user="XXX") # not tested

    To set default values populate ``pushover_defaults``::

        sage: sage.misc.messaging.pushover_defaults["user"] = "USER_TOKEN"
        sage: sage.misc.messaging.pushover("Hi, how are you?") # not tested

    .. note::

        You may want to populate ``sage.misc.messaging.pushover_defaults`` with default values such
        as the default user in ``$HOME/.sage/init.sage``.
    """
    request = {"message": message}
    request.update(pushover_defaults)
    request.update(kwds)

    conn = httplib.HTTPSConnection("api.pushover.net:443", context=SSLContext())
    conn.request("POST", "/1/messages.json",
                 urlencode(request),
                 {"Content-type": "application/x-www-form-urlencoded"})
    return conn.getresponse().status == 200
