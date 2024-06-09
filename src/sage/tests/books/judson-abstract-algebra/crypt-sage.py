# sage_setup: distribution = sagemath-repl
##          Sage Doctest File         ##
#**************************************#
#*    Generated from PreTeXt source   *#
#*    on 2017-08-24T11:43:34-07:00    *#
#*                                    *#
#*   http://mathbook.pugetsound.edu   *#
#*                                    *#
#**************************************#
##
"""
Please contact Rob Beezer (beezer@ups.edu) with
any test failures here that need to be changed
as a result of changes accepted into Sage.  You
may edit/change this file in any sensible way, so
that development work may procede.  Your changes
may later be replaced by the authors of "Abstract
Algebra: Theory and Applications" when the text is
updated, and a replacement of this file is proposed
for review.
"""
##
## To execute doctests in these files, run
##   $ $SAGE_ROOT/sage -t <directory-of-these-files>
## or
##   $ $SAGE_ROOT/sage -t <a-single-file>
##
## Replace -t by "-tp n" for parallel testing,
##   "-tp 0" will use a sensible number of threads
##
## See: http://www.sagemath.org/doc/developer/doctesting.html
##   or run  $ $SAGE_ROOT/sage --advanced  for brief help
##
## Generated at 2017-08-24T11:43:34-07:00
## From "Abstract Algebra"
## At commit 26d3cac0b4047f4b8d6f737542be455606e2c4b4
##
## Section 7.6 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p_a = next_prime(10^10)
    sage: q_a = next_prime(p_a)
    sage: p_b = next_prime((3/2)*10^10)
    sage: q_b = next_prime(p_b)
    sage: n_a = p_a * q_a
    sage: n_b = p_b * q_b
    sage: n_a, n_b
    (100000000520000000627, 225000000300000000091)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: m_a = euler_phi(n_a)
    sage: m_b = euler_phi(n_b)
    sage: m_a, m_b
    (100000000500000000576, 225000000270000000072)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: factor(m_a)
    2^6 * 3 * 11 * 17 * 131 * 521 * 73259 * 557041

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: E_a = 5*23
    sage: D_a = inverse_mod(E_a, m_a)
    sage: D_a
    20869565321739130555

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: factor(m_b)
    2^3 * 3^4 * 107 * 1298027 * 2500000001

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: E_b = 7*29
    sage: D_b = inverse_mod(E_b, m_b)
    sage: D_b
    24384236482463054195

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: print("Alice's public key, n:", n_a, "E:", E_a)
    Alice's public key, n: 100000000520000000627 E: 115

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: print("Alice's private key, D:", D_a)
    Alice's private key, D: 20869565321739130555

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: print("Bob's public key, n:", n_b, "E:", E_b)
    Bob's public key, n: 225000000300000000091 E: 203

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: print("Bob's private key, D:", D_b)
    Bob's private key, D: 24384236482463054195

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: word = 'Sage'
    sage: digits = [ord(letter) for letter in word]
    sage: digits
    [83, 97, 103, 101]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: message = ZZ(digits, 128)
    sage: message
    213512403

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: signed = power_mod(message, D_a, n_a)
    sage: signed
    47838774644892618423

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: encrypted = power_mod(signed, E_b, n_b)
    sage: encrypted
    111866209291209840488

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: decrypted = power_mod(encrypted, D_b, n_b)
    sage: decrypted
    47838774644892618423

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: received = power_mod(decrypted, E_a, n_a)
    sage: received
    213512403

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: digits = received.digits(base=128)
    sage: letters = [chr(ascii) for ascii in digits]
    sage: letters
    ['S', 'a', 'g', 'e']

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: ''.join(letters)
    'Sage'
"""
