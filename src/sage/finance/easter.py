"""
Copyright (c) 2003-2005  Gustavo Niemeyer <gustavo@niemeyer.net>

This module offers extensions to the standard python 2.3+
datetime module.
"""
__author__ = "Gustavo Niemeyer <gustavo@niemeyer.net>"
__license__ = "PSF License"

import datetime

__all__ = ["easter", "EASTER_JULIAN", "EASTER_ORTHODOX", "EASTER_WESTERN"]

EASTER_JULIAN   = 1
EASTER_ORTHODOX = 2
EASTER_WESTERN  = 3


def easter(year, algorithm=EASTER_WESTERN):
    """
    This function was ported from the work done by GM Arts,
    on top of the algorithm by Claus Tondering, which was
    based in part on the algorithm of Ouding (1940), as
    quoted in "Explanatory Supplement to the Astronomical
    Almanac", P.  Kenneth Seidelmann, editor.

    This function implements three different easter
    calculation algorithms:

    1 - Original calculation in Julian calendar, valid in
        dates after 326 AD
    2 - Original algorithm, with date converted to Gregorian
        calendar, valid in years 1583 to 4099
    3 - Revised algorithm, in Gregorian calendar, valid in
        years 1583 to 4099 as well

    These algorithms are represented by the constants:

    EASTER_JULIAN   = 1
    EASTER_ORTHODOX = 2
    EASTER_WESTERN  = 3

    The default algorithm is algorithm 3.

    More about the algorithm may be found at:

    http://users.chariot.net.au/~gmarts/eastalg.htm

    and

    https://www.tondering.dk/claus/calendar.html

    EXAMPLES::

        sage: import sage.finance.easter
        doctest:warning...
        DeprecationWarning: the package sage.finance is deprecated...

        sage: sage.finance.easter.easter(2009, 1)
        datetime.date(2009, 4, 6)
        sage: sage.finance.easter.easter(2009, 2)
        datetime.date(2009, 4, 19)
        sage: sage.finance.easter.easter(2009, 3)
        datetime.date(2009, 4, 12)
        sage: sage.finance.easter.easter(2010, 1)
        datetime.date(2010, 3, 22)
        sage: sage.finance.easter.easter(2010, 2)
        datetime.date(2010, 4, 4)
        sage: sage.finance.easter.easter(2010, 3)
        datetime.date(2010, 4, 4)

    AUTHORS:

     - Gustavo Niemeyer (author of function)
     - Phaedon Sinis (adapted code for sage)
    """
    year = int(year)
    algorithm = int(algorithm)

    if not (1 <= algorithm <= 3):
        raise ValueError("invalid algorithm")

    # g - Golden year - 1
    # c - Century
    # h - (23 - Epact) mod 30
    # i - Number of days from March 21 to Paschal Full Moon
    # j - Weekday for PFM (0=Sunday, etc)
    # p - Number of days from March 21 to Sunday on or before PFM
    #     (-6 to 28 algorithms 1 & 3, to 56 for algorithm 2)
    # e - Extra days to add for algorithm 2 (converting Julian
    #     date to Gregorian date)

    y = year
    g = y % 19
    e = 0
    if algorithm < 3:
        # Old algorithm
        i = (19*g+15)%30
        j = (y+y//4+i)%7
        if algorithm == 2:
            # Extra dates to convert Julian to Gregorian date
            e = 10
            if y > 1600:
                e = e+y//100-16-(y//100-16)//4
    else:
        # New algorithm
        c = y//100
        h = (c-c//4-(8*c+13)//25+19*g+15)%30
        i = h-(h//28)*(1-(h//28)*(29//(h+1))*((21-g)//11))
        j = (y+y//4+i+2-c+c//4)%7

    # p can be from -6 to 56 corresponding to dates 22 March to 23 May
    # (later dates apply to algorithm 2, although 23 May never actually occurs)
    p = i-j+e
    d = 1+(p+27+(p+6)//40)%31
    m = 3+(p+26)//30
    return datetime.date(y, m, d)
