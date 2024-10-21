# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
DES

The Data Encryption Standard.

This file implements the Data Encryption Standard and the corresponding key
schedule as described in [U.S1999]_.

This implementation is meant for experimental and educational usage only,
do not use it in production code!

EXAMPLES:

Encrypt a message::

    sage: from sage.crypto.block_cipher.des import DES
    sage: des = DES()
    sage: P = 0x01A1D6D039776742
    sage: K = 0x7CA110454A1A6E57
    sage: C = des.encrypt(plaintext=P, key=K); C.hex()
    '690f5b0d9a26939b'

And decrypt it again::

    sage: des.decrypt(ciphertext=C, key=K).hex()
    '1a1d6d039776742'

Have a look at the used round keys::

    sage: from sage.crypto.block_cipher.des import DES_KS
    sage: ks = DES_KS()
    sage: [k.hex() for k in ks(0x1F08260D1AC2465E)]
    ['103049bfb90e',
     '808d40f07bf',
      ...
     '231000f2dd97']

Validate the Sample Round Outputs for DES (cf. [KeSm1998]_ p. 124)::

    sage: from sage.crypto.block_cipher.des import DES
    sage: P = 0
    sage: K = 0x10316E028C8F3B4A
    sage: for r in range(1, 17):
    ....:     DES(rounds=r, doFinalRound=False).encrypt(P, K).hex()
    '47092b5b'
    '47092b5b53f372af'
    '53f372af9f1d158b'
    ...
    '3f6c3efd5a1e5228'
    sage: DES().encrypt(P, K).hex()
    '82dcbafbdeab6602'

Change cipher internals::

    sage: from sage.crypto.sbox import SBox
    sage: cipher = DES(rounds=1, doFinalRound=False)
    sage: cipher.sboxes = [[SBox(range(16))]*4]*8
    sage: cipher.keySchedule = lambda x: [0]  # return the 0 key as round key
    sage: cipher.encrypt(plaintext=0x1234, key=0x0).hex()
    '80004000d08100'

AUTHORS:

- Lukas Stennes (2019-03-29): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013 Lukas Stennes <lukas.stennes@rub.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.structure.element import Vector
from sage.rings.integer import Integer
from sage.crypto.sboxes import DES_S1_1, DES_S1_2, DES_S1_3, DES_S1_4
from sage.crypto.sboxes import DES_S2_1, DES_S2_2, DES_S2_3, DES_S2_4
from sage.crypto.sboxes import DES_S3_1, DES_S3_2, DES_S3_3, DES_S3_4
from sage.crypto.sboxes import DES_S4_1, DES_S4_2, DES_S4_3, DES_S4_4
from sage.crypto.sboxes import DES_S5_1, DES_S5_2, DES_S5_3, DES_S5_4
from sage.crypto.sboxes import DES_S6_1, DES_S6_2, DES_S6_3, DES_S6_4
from sage.crypto.sboxes import DES_S7_1, DES_S7_2, DES_S7_3, DES_S7_4
from sage.crypto.sboxes import DES_S8_1, DES_S8_2, DES_S8_3, DES_S8_4
from itertools import chain


sboxes = [[DES_S1_1, DES_S1_2, DES_S1_3, DES_S1_4],
          [DES_S2_1, DES_S2_2, DES_S2_3, DES_S2_4],
          [DES_S3_1, DES_S3_2, DES_S3_3, DES_S3_4],
          [DES_S4_1, DES_S4_2, DES_S4_3, DES_S4_4],
          [DES_S5_1, DES_S5_2, DES_S5_3, DES_S5_4],
          [DES_S6_1, DES_S6_2, DES_S6_3, DES_S6_4],
          [DES_S7_1, DES_S7_2, DES_S7_3, DES_S7_4],
          [DES_S8_1, DES_S8_2, DES_S8_3, DES_S8_4]]


class DES(SageObject):
    r"""
    This class implements DES described in [U.S1999]_.

    EXAMPLES:

    You can invoke DES encryption/decryption either by calling DES with an
    appropriate flag::

        sage: from sage.crypto.block_cipher.des import DES
        sage: des = DES()
        sage: P = 0x8000000000000000
        sage: K = 0x0
        sage: C = des(P, K, 'encrypt'); C.hex()
        '95f8a5e5dd31d900'
        sage: des(C, K, 'decrypt').hex()
        '8000000000000000'

    Or by calling encryption/decryption methods directly::

        sage: C = des.encrypt(P, K)
        sage: P == des.decrypt(C, K)
        True

    The number of rounds can be reduced easily::

        sage: des = DES(rounds=15)
        sage: des(des(P, K, 'encrypt'), K, 'decrypt') == P
        True

    You can use hex (i.e. integers) or a list-like bit representation for the
    inputs. If the input is an integer the output will be too. If it is
    list-like the output will be a bit vector::

        sage: des = DES()
        sage: P = vector(GF(2), 64, [1] + [0]*63)
        sage: K = vector(GF(2), 64, [0,0,0,0,0,0,0,1]*8)
        sage: des.encrypt(P, K)
        (1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0,
         1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0,
         0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
        sage: P = 0x8000000000000000
        sage: K = 0x0101010101010101
        sage: C = des.encrypt(P, K); C; C.hex()
        10806569712552630528
        '95f8a5e5dd31d900'

    .. SEEALSO::

        :class:`DES_KS`
        :mod:`sage.crypto.sboxes`

    TESTS:

    Test test vectors from [KeSm1998]_ pp. 125-136::

        sage: from sage.crypto.block_cipher.des import DES
        sage: test = \
        ....: [[0x0101010101010101, 0x8000000000000000, 0x95F8A5E5DD31D900],
        ....:  [0x0101010101010101, 0x4000000000000000, 0xDD7F121CA5015619],
        ....:  [0x0101010101010101, 0x2000000000000000, 0x2E8653104F3834EA],
        ....:  [0x0101010101010101, 0x1000000000000000, 0x4BD388FF6CD81D4F],
        ....:  [0x0101010101010101, 0x0800000000000000, 0x20B9E767B2FB1456],
        ....:  [0x0101010101010101, 0x0400000000000000, 0x55579380D77138EF],
        ....:  [0x0101010101010101, 0x0200000000000000, 0x6CC5DEFAAF04512F],
        ....:  [0x0101010101010101, 0x0100000000000000, 0x0D9F279BA5D87260],
        ....:  [0x0101010101010101, 0x0080000000000000, 0xD9031B0271BD5A0A],
        ....:  [0x0101010101010101, 0x0040000000000000, 0x424250B37C3DD951],
        ....:  [0x0101010101010101, 0x0020000000000000, 0xB8061B7ECD9A21E5],
        ....:  [0x0101010101010101, 0x0010000000000000, 0xF15D0F286B65BD28],
        ....:  [0x0101010101010101, 0x0008000000000000, 0xADD0CC8D6E5DEBA1],
        ....:  [0x0101010101010101, 0x0004000000000000, 0xE6D5F82752AD63D1],
        ....:  [0x0101010101010101, 0x0002000000000000, 0xECBFE3BD3F591A5E],
        ....:  [0x0101010101010101, 0x0001000000000000, 0xF356834379D165CD],
        ....:  [0x0101010101010101, 0x0000800000000000, 0x2B9F982F20037FA9],
        ....:  [0x0101010101010101, 0x0000400000000000, 0x889DE068A16F0BE6],
        ....:  [0x0101010101010101, 0x0000200000000000, 0xE19E275D846A1298],
        ....:  [0x0101010101010101, 0x0000100000000000, 0x329A8ED523D71AEC],
        ....:  [0x0101010101010101, 0x0000080000000000, 0xE7FCE22557D23C97],
        ....:  [0x0101010101010101, 0x0000040000000000, 0x12A9F5817FF2D65D],
        ....:  [0x0101010101010101, 0x0000020000000000, 0xA484C3AD38DC9C19],
        ....:  [0x0101010101010101, 0x0000010000000000, 0xFBE00A8A1EF8AD72],
        ....:  [0x0101010101010101, 0x0000008000000000, 0x750D079407521363],
        ....:  [0x0101010101010101, 0x0000004000000000, 0x64FEED9C724C2FAF],
        ....:  [0x0101010101010101, 0x0000002000000000, 0xF02B263B328E2B60],
        ....:  [0x0101010101010101, 0x0000001000000000, 0x9D64555A9A10B852],
        ....:  [0x0101010101010101, 0x0000000800000000, 0xD106FF0BED5255D7],
        ....:  [0x0101010101010101, 0x0000000400000000, 0xE1652C6B138C64A5],
        ....:  [0x0101010101010101, 0x0000000200000000, 0xE428581186EC8F46],
        ....:  [0x0101010101010101, 0x0000000100000000, 0xAEB5F5EDE22D1A36],
        ....:  [0x0101010101010101, 0x0000000080000000, 0xE943D7568AEC0C5C],
        ....:  [0x0101010101010101, 0x0000000040000000, 0xDF98C8276F54B04B],
        ....:  [0x0101010101010101, 0x0000000020000000, 0xB160E4680F6C696F],
        ....:  [0x0101010101010101, 0x0000000010000000, 0xFA0752B07D9C4AB8],
        ....:  [0x0101010101010101, 0x0000000008000000, 0xCA3A2B036DBC8502],
        ....:  [0x0101010101010101, 0x0000000004000000, 0x5E0905517BB59BCF],
        ....:  [0x0101010101010101, 0x0000000002000000, 0x814EEB3B91D90726],
        ....:  [0x0101010101010101, 0x0000000001000000, 0x4D49DB1532919C9F],
        ....:  [0x0101010101010101, 0x0000000000800000, 0x25EB5FC3F8CF0621],
        ....:  [0x0101010101010101, 0x0000000000400000, 0xAB6A20C0620D1C6F],
        ....:  [0x0101010101010101, 0x0000000000200000, 0x79E90DBC98F92CCA],
        ....:  [0x0101010101010101, 0x0000000000100000, 0x866ECEDD8072BB0E],
        ....:  [0x0101010101010101, 0x0000000000080000, 0x8B54536F2F3E64A8],
        ....:  [0x0101010101010101, 0x0000000000040000, 0xEA51D3975595B86B],
        ....:  [0x0101010101010101, 0x0000000000020000, 0xCAFFC6AC4542DE31],
        ....:  [0x0101010101010101, 0x0000000000010000, 0x8DD45A2DDF90796C],
        ....:  [0x0101010101010101, 0x0000000000008000, 0x1029D55E880EC2D0],
        ....:  [0x0101010101010101, 0x0000000000004000, 0x5D86CB23639DBEA9],
        ....:  [0x0101010101010101, 0x0000000000002000, 0x1D1CA853AE7C0C5F],
        ....:  [0x0101010101010101, 0x0000000000001000, 0xCE332329248F3228],
        ....:  [0x0101010101010101, 0x0000000000000800, 0x8405D1ABE24FB942],
        ....:  [0x0101010101010101, 0x0000000000000400, 0xE643D78090CA4207],
        ....:  [0x0101010101010101, 0x0000000000000200, 0x48221B9937748A23],
        ....:  [0x0101010101010101, 0x0000000000000100, 0xDD7C0BBD61FAFD54],
        ....:  [0x0101010101010101, 0x0000000000000080, 0x2FBC291A570DB5C4],
        ....:  [0x0101010101010101, 0x0000000000000040, 0xE07C30D7E4E26E12],
        ....:  [0x0101010101010101, 0x0000000000000020, 0x0953E2258E8E90A1],
        ....:  [0x0101010101010101, 0x0000000000000010, 0x5B711BC4CEEBF2EE],
        ....:  [0x0101010101010101, 0x0000000000000008, 0xCC083F1E6D9E85F6],
        ....:  [0x0101010101010101, 0x0000000000000004, 0xD2FD8867D50D2DFE],
        ....:  [0x0101010101010101, 0x0000000000000002, 0x06E7EA22CE92708F],
        ....:  [0x0101010101010101, 0x0000000000000001, 0x166B40B44ABA4BD6],
        ....:  [0x8001010101010101, 0x0000000000000000, 0x95A8D72813DAA94D],
        ....:  [0x4001010101010101, 0x0000000000000000, 0x0EEC1487DD8C26D5],
        ....:  [0x2001010101010101, 0x0000000000000000, 0x7AD16FFB79C45926],
        ....:  [0x1001010101010101, 0x0000000000000000, 0xD3746294CA6A6CF3],
        ....:  [0x0801010101010101, 0x0000000000000000, 0x809F5F873C1FD761],
        ....:  [0x0401010101010101, 0x0000000000000000, 0xC02FAFFEC989D1FC],
        ....:  [0x0201010101010101, 0x0000000000000000, 0x4615AA1D33E72F10],
        ....:  [0x0180010101010101, 0x0000000000000000, 0x2055123350C00858],
        ....:  [0x0140010101010101, 0x0000000000000000, 0xDF3B99D6577397C8],
        ....:  [0x0120010101010101, 0x0000000000000000, 0x31FE17369B5288C9],
        ....:  [0x0110010101010101, 0x0000000000000000, 0xDFDD3CC64DAE1642],
        ....:  [0x0108010101010101, 0x0000000000000000, 0x178C83CE2B399D94],
        ....:  [0x0104010101010101, 0x0000000000000000, 0x50F636324A9B7F80],
        ....:  [0x0102010101010101, 0x0000000000000000, 0xA8468EE3BC18F06D],
        ....:  [0x0101800101010101, 0x0000000000000000, 0xA2DC9E92FD3CDE92],
        ....:  [0x0101400101010101, 0x0000000000000000, 0xCAC09F797D031287],
        ....:  [0x0101200101010101, 0x0000000000000000, 0x90BA680B22AEB525],
        ....:  [0x0101100101010101, 0x0000000000000000, 0xCE7A24F350E280B6],
        ....:  [0x0101080101010101, 0x0000000000000000, 0x882BFF0AA01A0B87],
        ....:  [0x0101040101010101, 0x0000000000000000, 0x25610288924511C2],
        ....:  [0x0101020101010101, 0x0000000000000000, 0xC71516C29C75D170],
        ....:  [0x0101018001010101, 0x0000000000000000, 0x5199C29A52C9F059],
        ....:  [0x0101014001010101, 0x0000000000000000, 0xC22F0A294A71F29F],
        ....:  [0x0101012001010101, 0x0000000000000000, 0xEE371483714C02EA],
        ....:  [0x0101011001010101, 0x0000000000000000, 0xA81FBD448F9E522F],
        ....:  [0x0101010801010101, 0x0000000000000000, 0x4F644C92E192DFED],
        ....:  [0x0101010401010101, 0x0000000000000000, 0x1AFA9A66A6DF92AE],
        ....:  [0x0101010201010101, 0x0000000000000000, 0xB3C1CC715CB879D8],
        ....:  [0x0101010180010101, 0x0000000000000000, 0x19D032E64AB0BD8B],
        ....:  [0x0101010140010101, 0x0000000000000000, 0x3CFAA7A7DC8720DC],
        ....:  [0x0101010120010101, 0x0000000000000000, 0xB7265F7F447AC6F3],
        ....:  [0x0101010110010101, 0x0000000000000000, 0x9DB73B3C0D163F54],
        ....:  [0x0101010108010101, 0x0000000000000000, 0x8181B65BABF4A975],
        ....:  [0x0101010104010101, 0x0000000000000000, 0x93C9B64042EAA240],
        ....:  [0x0101010102010101, 0x0000000000000000, 0x5570530829705592],
        ....:  [0x0101010101800101, 0x0000000000000000, 0x8638809E878787A0],
        ....:  [0x0101010101400101, 0x0000000000000000, 0x41B9A79AF79AC208],
        ....:  [0x0101010101200101, 0x0000000000000000, 0x7A9BE42F2009A892],
        ....:  [0x0101010101100101, 0x0000000000000000, 0x29038D56BA6D2745],
        ....:  [0x0101010101080101, 0x0000000000000000, 0x5495C6ABF1E5DF51],
        ....:  [0x0101010101040101, 0x0000000000000000, 0xAE13DBD561488933],
        ....:  [0x0101010101020101, 0x0000000000000000, 0x024D1FFA8904E389],
        ....:  [0x0101010101018001, 0x0000000000000000, 0xD1399712F99BF02E],
        ....:  [0x0101010101014001, 0x0000000000000000, 0x14C1D7C1CFFEC79E],
        ....:  [0x0101010101012001, 0x0000000000000000, 0x1DE5279DAE3BED6F],
        ....:  [0x0101010101011001, 0x0000000000000000, 0xE941A33F85501303],
        ....:  [0x0101010101010801, 0x0000000000000000, 0xDA99DBBC9A03F379],
        ....:  [0x0101010101010401, 0x0000000000000000, 0xB7FC92F91D8E92E9],
        ....:  [0x0101010101010201, 0x0000000000000000, 0xAE8E5CAA3CA04E85],
        ....:  [0x0101010101010180, 0x0000000000000000, 0x9CC62DF43B6EED74],
        ....:  [0x0101010101010140, 0x0000000000000000, 0xD863DBB5C59A91A0],
        ....:  [0x0101010101010120, 0x0000000000000000, 0xA1AB2190545B91D7],
        ....:  [0x0101010101010110, 0x0000000000000000, 0x0875041E64C570F7],
        ....:  [0x0101010101010108, 0x0000000000000000, 0x5A594528BEBEF1CC],
        ....:  [0x0101010101010104, 0x0000000000000000, 0xFCDB3291DE21F0C0],
        ....:  [0x0101010101010102, 0x0000000000000000, 0x869EFD7F9F265A09],
        ....:  [0x1046913489980131, 0x0000000000000000, 0x88D55E54F54C97B4],
        ....:  [0x1007103489988020, 0x0000000000000000, 0x0C0CC00C83EA48FD],
        ....:  [0x10071034C8980120, 0x0000000000000000, 0x83BC8EF3A6570183],
        ....:  [0x1046103489988020, 0x0000000000000000, 0xDF725DCAD94EA2E9],
        ....:  [0x1086911519190101, 0x0000000000000000, 0xE652B53B550BE8B0],
        ....:  [0x1086911519580101, 0x0000000000000000, 0xAF527120C485CBB0],
        ....:  [0x5107B01519580101, 0x0000000000000000, 0x0F04CE393DB926D5],
        ....:  [0x1007B01519190101, 0x0000000000000000, 0xC9F00FFC74079067],
        ....:  [0x3107915498080101, 0x0000000000000000, 0x7CFD82A593252B4E],
        ....:  [0x3107919498080101, 0x0000000000000000, 0xCB49A2F9E91363E3],
        ....:  [0x10079115B9080140, 0x0000000000000000, 0x00B588BE70D23F56],
        ....:  [0x3107911598090140, 0x0000000000000000, 0x406A9A6AB43399AE],
        ....:  [0x1007D01589980101, 0x0000000000000000, 0x6CB773611DCA9ADA],
        ....:  [0x9107911589980101, 0x0000000000000000, 0x67FD21C17DBB5D70],
        ....:  [0x9107D01589190101, 0x0000000000000000, 0x9592CB4110430787],
        ....:  [0x1007D01598980120, 0x0000000000000000, 0xA6B7FF68A318DDD3],
        ....:  [0x1007940498190101, 0x0000000000000000, 0x4D102196C914CA16],
        ....:  [0x0107910491190401, 0x0000000000000000, 0x2DFA9F4573594965],
        ....:  [0x0107910491190101, 0x0000000000000000, 0xB46604816C0E0774],
        ....:  [0x0107940491190401, 0x0000000000000000, 0x6E7E6221A4F34E87],
        ....:  [0x19079210981A0101, 0x0000000000000000, 0xAA85E74643233199],
        ....:  [0x1007911998190801, 0x0000000000000000, 0x2E5A19DB4D1962D6],
        ....:  [0x10079119981A0801, 0x0000000000000000, 0x23A866A809D30894],
        ....:  [0x1007921098190101, 0x0000000000000000, 0xD812D961F017D320],
        ....:  [0x100791159819010B, 0x0000000000000000, 0x055605816E58608F],
        ....:  [0x1004801598190101, 0x0000000000000000, 0xABD88E8B1B7716F1],
        ....:  [0x1004801598190102, 0x0000000000000000, 0x537AC95BE69DA1E1],
        ....:  [0x1004801598190108, 0x0000000000000000, 0xAED0F6AE3C25CDD8],
        ....:  [0x1002911598100104, 0x0000000000000000, 0xB3E35A5EE53E7B8D],
        ....:  [0x1002911598190104, 0x0000000000000000, 0x61C79C71921A2EF8],
        ....:  [0x1002911598100201, 0x0000000000000000, 0xE2F5728F0995013C],
        ....:  [0x1002911698100101, 0x0000000000000000, 0x1AEAC39A61F0A464],
        ....:  [0x7CA110454A1A6E57, 0x01A1D6D039776742, 0x690F5B0D9A26939B],
        ....:  [0x0131D9619DC1376E, 0x5CD54CA83DEF57DA, 0x7A389D10354BD271],
        ....:  [0x07A1133E4A0B2686, 0x0248D43806F67172, 0x868EBB51CAB4599A],
        ....:  [0x3849674C2602319E, 0x51454B582DDF440A, 0x7178876E01F19B2A],
        ....:  [0x04B915BA43FEB5B6, 0x42FD443059577FA2, 0xAF37FB421F8C4095],
        ....:  [0x0113B970FD34F2CE, 0x059B5E0851CF143A, 0x86A560F10EC6D85B],
        ....:  [0x0170F175468FB5E6, 0x0756D8E0774761D2, 0x0CD3DA020021DC09],
        ....:  [0x43297FAD38E373FE, 0x762514B829BF486A, 0xEA676B2CB7DB2B7A],
        ....:  [0x07A7137045DA2A16, 0x3BDD119049372802, 0xDFD64A815CAF1A0F],
        ....:  [0x04689104C2FD3B2F, 0x26955F6835AF609A, 0x5C513C9C4886C088],
        ....:  [0x37D06BB516CB7546, 0x164D5E404F275232, 0x0A2AEEAE3FF4AB77],
        ....:  [0x1F08260D1AC2465E, 0x6B056E18759F5CCA, 0xEF1BF03E5DFA575A],
        ....:  [0x584023641ABA6176, 0x004BD6EF09176062, 0x88BF0DB6D70DEE56],
        ....:  [0x025816164629B007, 0x480D39006EE762F2, 0xA1F9915541020B56],
        ....:  [0x49793EBC79B3258F, 0x437540C8698F3CFA, 0x6FBF1CAFCFFD0556],
        ....:  [0x4FB05E1515AB73A7, 0x072D43A077075292, 0x2F22E49BAB7CA1AC],
        ....:  [0x49E95D6D4CA229BF, 0x02FE55778117F12A, 0x5A6B612CC26CCE4A],
        ....:  [0x018310DC409B26D6, 0x1D9D5C5018F728C2, 0x5F4C038ED12B2E41],
        ....:  [0x1C587F1C13924FEF, 0x305532286D6F295A, 0x63FAC0D034D9F793]]
        sage: des = DES()
        sage: for K, P, C in test: # long time
        ....:    if des.encrypt(P, K) != C or des.decrypt(C, K) != P:
        ....:        print("DES tests failed for K=0x%s, P=0x%s, C=0x%s" %
        ....:              (K.hex(), P.hex(), C.hex()))

    .. automethod:: __init__
    .. automethod:: __call__
    """

    def __init__(self, rounds=None, keySchedule='DES_KS', keySize=64, doFinalRound=True):
        r"""
        Construct an instance of DES.

        INPUT:

        - ``rounds`` -- integer (default: ``None``); the number of rounds. If
          ``None`` the number of rounds of the key schedule is used.

        - ``keySchedule`` -- (default: ``'DES_KS'``) the key schedule that
          will be used for encryption and decryption. If ``'DES_KS'`` the
          default DES key schedule is used.

        - ``keySize`` -- (default: ``64``) the key length in bits. Must be
          ``56`` of ``64``. In the latter case the key contains 8 parity bits.

        - ``doFinalRound`` -- boolean (default: ``True``); if ``False`` a swap
          takes places but the inverse initial permutation is omitted (i.e. you
          can get the state after ``rounds``). This only effects encryption.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: DES() # indirect doctest
            DES block cipher with 16 rounds and the following key schedule:
            Original DES key schedule with 16 rounds

        Reducing the number of rounds is simple. But increasing it is only
        possible if the key schedule can produce enough round keys::

            sage: DES(rounds=11) # indirect doctest
            DES block cipher with 11 rounds and the following key schedule:
            Original DES key schedule with 16 rounds
            sage: DES(rounds=42) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: number of rounds must be less or equal to the number
            of rounds of the key schedule

        You can use arbitrary key schedules. Since it is the only one
        implemented here the original key schedule is used for demonstration::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: DES(keySchedule=DES_KS(11)) # indirect doctest
            DES block cipher with 11 rounds and the following key schedule:
            Original DES key schedule with 11 rounds
        """
        self.keySchedule = DES_KS() if keySchedule == 'DES_KS' else keySchedule
        self._rounds = self.keySchedule._rounds if rounds is None else rounds
        if self._rounds > self.keySchedule._rounds:
            raise ValueError('number of rounds must be less or equal to the '
                             'number of rounds of the key schedule')
        self._keySize = keySize
        if keySize not in (56, 64):
            raise ValueError('key size must be 56 or 64')
        self.sboxes = sboxes
        self._doFinalRound = doFinalRound
        self._blocksize = 64

    def __call__(self, block, key, algorithm='encrypt'):
        r"""
        Apply DES encryption or decryption on ``block`` using ``key``. The flag
        ``algorithm`` controls what action is to be performed on ``block``.

        INPUT:

        - ``block`` -- integer or bit list-like; the plaintext or ciphertext

        - ``key`` -- integer or bit list-like; the key

        - ``algorithm`` -- string (default: ``'encrypt'``); a flag to signify
          whether encryption or decryption is to be applied to ``block``. The
          encryption flag is ``'encrypt'`` and the decryption flag is
          ``'decrypt'``

        OUTPUT:

        - The plaintext or ciphertext corresponding to ``block``, obtained
          using ``key``. If ``block`` is an integer the output will be too. If
          ``block`` is list-like the output will be a bit vector.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: P = 0x480D39006EE762F2
            sage: K = 0x025816164629B007
            sage: des(P, K, 'encrypt').hex()
            'a1f9915541020b56'
        """
        if algorithm == 'encrypt':
            return self.encrypt(block, key)
        elif algorithm == 'decrypt':
            return self.decrypt(block, key)
        else:
            raise ValueError('Algorithm must be \'encrypt\' or \'decrypt\' and'
                             ' not \'%s\'' % algorithm)

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        DES objects are the same if all attributes are the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des1 = DES()
            sage: des2 = DES()
            sage: des3 = DES(rounds=11)
            sage: des1 == des2
            True
            sage: des1 == des3
            False
            sage: des2 == 42
            False
        """
        if not isinstance(other, DES):
            return False
        else:
            return self.__dict__ == other.__dict__

    def __repr__(self):
        r"""
        A string representation of this DES.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: DES() # indirect doctest
            DES block cipher with 16 rounds and the following key schedule:
            Original DES key schedule with 16 rounds
        """
        return ('DES block cipher with %s rounds and the following key '
                'schedule:\n%s' % (self._rounds, self.keySchedule.__repr__()))

    def encrypt(self, plaintext, key):
        r"""
        Return the ciphertext corresponding to ``plaintext``, using DES
        encryption with ``key``.

        INPUT:

        - ``plaintext`` -- integer or bit list-like; the plaintext that will be
          encrypted

        - ``key`` -- integer or bit list-like; the key

        OUTPUT:

        - The ciphertext corresponding to ``plaintext``, obtained using
          ``key``. If ``plaintext`` is an integer the output will be too. If
          ``plaintext`` is list-like the output will be a bit vector.

        EXAMPLES:

        Encrypt a message::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: K64 = 0x133457799BBCDFF1
            sage: P = 0x0123456789ABCDEF
            sage: C = des.encrypt(P, K64); C.hex()
            '85e813540f0ab405'

        You can also use 56 bit keys i.e. you can leave out the parity bits::

            sage: K56 = 0x12695BC9B7B7F8
            sage: des = DES(keySize=56)
            sage: des.encrypt(P, K56) == C
            True
        """
        if isinstance(plaintext, (list, tuple, Vector)):
            inputType = 'vector'
        elif isinstance(plaintext, (Integer, int)):
            inputType = 'integer'
        state = convert_to_vector(plaintext, self._blocksize)
        key = convert_to_vector(key, self._keySize)
        if self._keySize == 56:
            # insert 'parity' bits
            key = list(key)
            for i in range(7, 64, 8):
                key.insert(i, 0)
            key = vector(GF(2), 64, key)
        roundKeys = self.keySchedule(key)
        state = self._ip(state)
        for k in roundKeys[:self._rounds]:
            state = self.round(state, k)
        if self._doFinalRound:
            state = vector(GF(2), 64, list(state[32:64])+list(state[0:32]))
            state = self._inv_ip(state)
        return state if inputType == 'vector' else ZZ(list(state)[::-1], 2)

    def decrypt(self, ciphertext, key):
        r"""
        Return the plaintext corresponding to ``ciphertext``, using DES
        decryption with ``key``.

        INPUT:

        - ``ciphertext`` -- integer or bit list-like; the ciphertext that will
          be decrypted

        - ``key`` -- integer or bit list-like; the key

        OUTPUT:

        - The plaintext corresponding to ``ciphertext``, obtained using
          ``key``. If ``ciphertext`` is an integer the output will be too. If
          ``ciphertext`` is list-like the output will be a bit vector.

        EXAMPLES:

        Decrypt a message::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: K64 = 0x7CA110454A1A6E57
            sage: C = 0x690F5B0D9A26939B
            sage: P = des.decrypt(C, K64).hex(); P
            '1a1d6d039776742'

        You can also use 56 bit keys i.e. you can leave out the parity bits::

            sage: K56 = 0x7D404224A35BAB
            sage: des = DES(keySize=56)
            sage: des.decrypt(C, K56).hex() == P
            True
        """
        if isinstance(ciphertext, (list, tuple, Vector)):
            inputType = 'vector'
        elif isinstance(ciphertext, (Integer, int)):
            inputType = 'integer'
        state = convert_to_vector(ciphertext, 64)
        key = convert_to_vector(key, self._keySize)
        if self._keySize == 56:
            # insert 'parity' bits
            key = list(key)
            for i in range(7, 64, 8):
                key.insert(i, 0)
            key = vector(GF(2), 64, key)
        roundKeys = self.keySchedule(key)
        state = self._ip(state)
        for k in roundKeys[:self._rounds][::-1]:
            state = self.round(state, k)
        state = vector(GF(2), 64, list(state[32:64])+list(state[0:32]))
        state = self._inv_ip(state)
        return state if inputType == 'vector' else ZZ(list(state)[::-1], 2)

    def _ip(self, block):
        r"""
        Return the initial permutation of ``block``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: B = vector(GF(2), 64, [0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,
            ....:                        0,0,1,0,1,0,1,1,0,0,1,1,1,1,0,0,0,1,0,
            ....:                        0,1,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,1,1,
            ....:                        1,1,0,1,1,1,1])
            sage: des._ip(B)
            (1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1,
             0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0,
             1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0)
        """
        IP = [58, 50, 42, 34, 26, 18, 10, 2,
              60, 52, 44, 36, 28, 20, 12, 4,
              62, 54, 46, 38, 30, 22, 14, 6,
              64, 56, 48, 40, 32, 24, 16, 8,
              57, 49, 41, 33, 25, 17,  9, 1,
              59, 51, 43, 35, 27, 19, 11, 3,
              61, 53, 45, 37, 29, 21, 13, 5,
              63, 55, 47, 39, 31, 23, 15, 7]
        return vector(GF(2), 64, [block[i-1] for i in IP])

    def round(self, state, round_key):
        r"""
        Apply one round of DES to ``state`` and return the result.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: from sage.crypto.block_cipher.des import convert_to_vector
            sage: des = DES()
            sage: k1 = convert_to_vector(0xFFFFFFFFFFFF, 48)
            sage: state = convert_to_vector(0xFFFFFFFF11111111, 64)
            sage: ZZ(list(des.round(state, k1))[::-1], 2).hex()
            '11111111b29684b8'
        """
        L, R = state[0:32], state[32:64]
        L, R = R, L + self._f(R, round_key)
        state = vector(GF(2), 64, list(L)+list(R))
        return state

    def _f(self, right, subkey):
        r"""
        Apply the cipher function to ``right`` and ``subkey``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: R = vector(GF(2), 32, [1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0,1,1,1,
            ....:                        1,0,0,0,0,1,0,1,0,1,0,1,0])
            sage: K = vector(GF(2), 48, [0,0,0,1,1,0,1,1,0,0,0,0,0,0,1,0,1,1,1,
            ....:                        0,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,0,0,
            ....:                        0,0,0,1,1,1,0,0,1,0])
            sage: des._f(R, K)
            (0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
             0, 1, 1, 0, 1, 1, 1, 0, 1, 1)
        """
        return self._permutaion(self.sbox_layer(self._expand(right)+subkey))

    def _expand(self, right):
        r"""
        Apply the expansion function to ``right``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: R = vector(GF(2), 32, [1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0,1,1,1,
            ....:                        1,0,0,0,0,1,0,1,0,1,0,1,0])
            sage: des._expand(R)
            (0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
             0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
             0, 1, 0, 1)
        """
        E = [32,  1,  2,  3,  4,  5,
              4,  5,  6,  7,  8,  9,
              8,  9, 10, 11, 12, 13,
             12, 13, 14, 15, 16, 17,
             16, 17, 18, 19, 20, 21,
             20, 21, 22, 23, 24, 25,
             24, 25, 26, 27, 28, 29,
             28, 29, 30, 31, 32,  1]
        return vector(GF(2), 48, [right[i-1] for i in E])

    def sbox_layer(self, block):
        r"""
        Apply the Sboxes to ``block``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: B = vector(GF(2), 48, [0,1,1,0,0,0,0,1,0,0,0,1,0,1,1,1,1,0,1,
            ....:                        1,1,0,1,0,1,0,0,0,0,1,1,0,0,1,1,0,0,1,
            ....:                        0,1,0,0,1,0,0,1,1,1])
            sage: des.sbox_layer(B)
            (0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1,
             0, 1, 1, 0, 0, 1, 0, 1, 1, 1)

        .. SEEALSO::

            :mod:`sage.crypto.sboxes`
        """
        s = self.sboxes
        block = [block[i:i+6] for i in range(0, 48, 6)]
        block = list(chain.from_iterable([s[i][ZZ([b[5], b[0]], 2)](b[1:5])
                                          for i, b in enumerate(block)]))
        return vector(GF(2), 32, block)

    def _permutaion(self, block):
        r"""
        Apply the permutation function to ``block``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: B = vector(GF(2), 32, [0,1,0,1,1,1,0,0,1,0,0,0,0,0,1,0,1,0,1,
            ....:                        1,0,1,0,1,1,0,0,1,0,1,1,1])
            sage: des._permutaion(B)
            (0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
             0, 1, 1, 0, 1, 1, 1, 0, 1, 1)
        """
        P = [16,  7, 20, 21,
             29, 12, 28, 17,
              1, 15, 23, 26,
              5, 18, 31, 10,
              2,  8, 24, 14,
             32, 27,  3,  9,
             19, 13, 30,  6,
             22, 11,  4, 25]
        return vector(GF(2), 32, [block[i-1] for i in P])

    def _inv_ip(self, block):
        r"""
        Apply the inverse permutation function to ``block``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: B = vector(GF(2), 64, [0,0,0,0,1,0,1,0,0,1,0,0,1,1,0,0,1,1,0,
            ....:                        1,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,
            ....:                        1,1,0,1,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,
            ....:                        0,1,1,0,1,0,0])
            sage: des._inv_ip(B)
            (1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
             1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
             1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1)
        """
        invIP = [40, 8, 48, 16, 56, 24, 64, 32,
                 39, 7, 47, 15, 55, 23, 63, 31,
                 38, 6, 46, 14, 54, 22, 62, 30,
                 37, 5, 45, 13, 53, 21, 61, 29,
                 36, 4, 44, 12, 52, 20, 60, 28,
                 35, 3, 43, 11, 51, 19, 59, 27,
                 34, 2, 42, 10, 50, 18, 58, 26,
                 33, 1, 41,  9, 49, 17, 57, 25]
        return vector(GF(2), 64, [block[i-1] for i in invIP])


class DES_KS(SageObject):
    r"""
    This class implements the DES key schedules described in [U.S1999]_.

    EXAMPLES:

    Initialise the key schedule with a `masterKey` to use it as an iterable::

        sage: from sage.crypto.block_cipher.des import DES_KS
        sage: ks = DES_KS(masterKey=0)
        sage: ks[0]
        0
        sage: ks[15]
        0

    Or omit the `masterKey` and pass a key when calling the key schedule::

        sage: ks = DES_KS()
        sage: K = ks(0x584023641ABA6176)
        sage: K[0].hex()
        'd0a2ed2fa124'
        sage: K[15].hex()
        '43b42af81183'

    .. SEEALSO::

        :class:`DES`

    .. automethod:: __init__
    .. automethod:: __call__
    """

    def __init__(self, rounds=16, masterKey=None):
        r"""
        Construct an instance of DES_KS.

        INPUT:

        - ``rounds`` -- integer (default: `16`); the number of rounds
          ``self`` can create keys for

        - ``masterKey`` -- integer or bit list-like (default: ``None``); the
          64-bit key that will be used

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: DES_KS()
            Original DES key schedule with 16 rounds

        .. NOTE::

            If you want to use a DES_KS object as an iterable you have to
            pass a ``masterKey`` value on initialisation. Otherwise you can
            omit ``masterKey`` and pass a key when you call the object.
        """
        self._rounds = rounds
        self._masterKey = masterKey
        self._keySize = 64

    def __call__(self, key):
        r"""
        Return all round keys in a list.

        INPUT:

        - ``key`` -- integer or bit list-like; the 64-bit key

        OUTPUT:

        - A list containing the round keys. If ``key`` is an integer the
          elements of the output list will be too. If ``key`` is list-like the
          element of the output list will be  bit vectors.

        EXAMPLES:

        This implementation is using bit vectors for all internal
        representations. So you can invoke the key schedule with a bit
        vector::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: K = vector(GF(2),[0,0,0,1,0,0,1,1,0,0,1,1,0,1,0,0,0,1,0,1,0,
            ....:                   1,1,1,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,1,1,0,
            ....:                   1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,1])
            sage: ks = DES_KS(16, K)
            sage: [k for k in ks]
            [(0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1,
              1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1,
              0, 0, 1, 0),
             (0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0,
              0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0,
              0, 1, 0, 1),
             ...
             (1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0,
              1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1,
              0, 1, 0, 1)]

        But of course you can invoke it with hex representation as well::

            sage: K = 0x133457799bbcdff1
            sage: ks = DES_KS(16, K)
            sage: [k.hex() for k in ks]
            ['1b02effc7072',
             '79aed9dbc9e5',
             ...
             'cb3d8b0e17f5']

        .. NOTE::

            If you want to use a DES_KS object as an iterable you have to
            pass a ``masterKey`` value on initialisation. Otherwise you can
            omit ``masterKey`` and pass a key when you call the object.
        """
        if isinstance(key, (list, tuple, Vector)):
            inputType = 'vector'
        elif isinstance(key, (Integer, int)):
            inputType = 'integer'
        key = convert_to_vector(key, self._keySize)
        roundKeys = []
        C, D = self._pc1(key)
        for i in range(16):
            C, D = self._left_shift(C, i), self._left_shift(D, i)
            roundKeys.append(self._pc2(list(C)+list(D)))
        return roundKeys if inputType == 'vector' else [ZZ(list(k)[::-1], 2)
                                                        for k in roundKeys]

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        DES_KS objects are the same if all attributes are the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: DES_KS() == DES_KS() # indirect doctest
            True
            sage: DES_KS() == DES_KS(11) # indirect doctest
            False
        """
        if not isinstance(other, DES_KS):
            return False
        else:
            return self.__dict__ == other.__dict__

    def __repr__(self):
        r"""
        A string representation of this DES_KS.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: DES_KS() # indirect doctest
            Original DES key schedule with 16 rounds
        """
        return ('Original DES key schedule with %s rounds' % (self._rounds))

    def __getitem__(self, r):
        r"""
        Compute the sub key for round ``r`` derived from initial master key.

        The key schedule object has to have been initialised with the
        `masterKey` argument.

        INPUT:

        - ``r`` -- integer; the round for which the sub key is computed

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS(masterKey=0x1F08260D1AC2465E)
            sage: ks[0].hex() # indirect doctest
            '103049bfb90e'
            sage: ks[15].hex() # indirect doctest
            '231000f2dd97'
        """
        if self._masterKey is None:
            raise ValueError('Key not set during initialisation')
        return self(self._masterKey)[r]

    def __iter__(self):
        r"""
        Iterate over the DES round keys, derived from `masterKey`.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: K = [k for k in DES_KS(masterKey=0x0113B970FD34F2CE)]
            sage: K[0].hex() # indirect doctest
            '6f26cc480fc6'
            sage: K[15].hex() # indirect doctest
            '9778f17524a'
        """
        if self._masterKey is None:
            raise ValueError('Key not set during initialisation')
        return iter(self(self._masterKey))

    def _pc1(self, key):
        r"""
        Return Permuted Choice 1 of ``key``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS()
            sage: K = vector(GF(2),[0,0,0,1,0,0,1,1,0,0,1,1,0,1,0,0,0,1,0,1,0,
            ....:                   1,1,1,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,1,1,0,
            ....:                   1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,1])
            sage: C, D = ks._pc1(K)
            sage: C
            (1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0,
             1, 0, 1, 1, 1, 1)
            sage: D
            (0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0,
             0, 0, 1, 1, 1, 1)
        """
        PC1_C = [57, 49, 41, 33, 25, 17,  9,
                  1, 58, 50, 42, 34, 26, 18,
                 10,  2, 59, 51, 43, 35, 27,
                 19, 11,  3, 60, 52, 44, 36]
        PC1_D = [63, 55, 47, 39, 31, 23, 15,
                  7, 62, 54, 46, 38, 30, 22,
                 14,  6, 61, 53, 45, 37, 29,
                 21, 13,  5, 28, 20, 12,  4]
        C = vector(GF(2), 28, [key[i-1] for i in PC1_C])
        D = vector(GF(2), 28, [key[i-1] for i in PC1_D])
        return C, D

    def _pc2(self, key):
        r"""
        Return Permuted Choice 2 of ``key``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS()
            sage: K = vector(GF(2),[1,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,1,0,1,0,
            ....:                   1,0,1,1,1,1,1,1,0,1,0,1,0,1,0,1,1,0,0,1,1,
            ....:                   0,0,1,1,1,1,0,0,0,1,1,1,1,0])
            sage: ks._pc2(K)
            (0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1,
             1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1,
             0, 0, 1, 0)
        """
        PC2 = [14, 17, 11, 24,  1,  5,
                3, 28, 15,  6, 21, 10,
               23, 19, 12,  4, 26,  8,
               16,  7, 27, 20, 13,  2,
               41, 52, 31, 37, 47, 55,
               30, 40, 51, 45, 33, 48,
               44, 49, 39, 56, 34, 53,
               46, 42, 50, 36, 29, 32]
        return vector(GF(2), 48, [key[i-1] for i in PC2])

    def _left_shift(self, half, i):
        r"""
        Shift ``half`` one or two positions to the left depending on the
        iteration number ``i``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS()
            sage: bits = vector(GF(2), 6, [1,0,1,0,1,0])
            sage: ks._left_shift(bits, 0)
            (0, 1, 0, 1, 0, 1)
            sage: bits
            (1, 0, 1, 0, 1, 0)
            sage: ks._left_shift(bits, 2)
            (1, 0, 1, 0, 1, 0)
        """
        amount = [1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1]
        return vector(GF(2),
                      list(half[amount[i]:]) + list(half[0:amount[i]]))


def convert_to_vector(I, L):
    r"""
    Convert ``I`` to a bit vector of length ``L``.

    INPUT:

    - ``I`` -- integer or bit list-like

    - ``L`` -- integer; the desired bit length of the ouput

    OUTPUT: the ``L``-bit vector representation of ``I``

    EXAMPLES::

        sage: from sage.crypto.block_cipher.des import convert_to_vector
        sage: convert_to_vector(0x1F, 8)
        (0, 0, 0, 1, 1, 1, 1, 1)
        sage: v = vector(GF(2), 4, [1,0,1,0])
        sage: convert_to_vector(v, 4)
        (1, 0, 1, 0)
    """
    try:
        state = vector(GF(2), L, ZZ(I).digits(2, padto=L)[::-1])
        return state
    except TypeError:
        # ignore the error and try list-like types
        pass
    state = vector(GF(2), L, [0]*(L-len(I))+list(I))
    return state
