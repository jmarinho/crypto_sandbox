#
# Implement Secure Hash functions defined in FIPS 180-4
#

import sys
import numpy as np
import logging as log

sha256_H = [0x6a09e667, 0xbb67ae85,
            0x3c6ef372, 0xa54ff53a,
            0x510e527f, 0x9b05688c,
            0x1f83d9ab, 0x5be0cd19]

K256 = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2]

def W_sha256 (t, Mi):

    if t<=15:
        return int.from_bytes(Mi[t*4:(t*4)+4], byteorder="big")

    elif t>= 16 and t<= 63:
        return np.uint32(sigma_1_256(W_sha256(t-2, Mi)) + W_sha256(t-7, Mi) + sigma_0_256(W_sha256(t-15, Mi)) + W_sha256(t-16, Mi))
    else:
        log.error("w_sha256: t our of range ({})".format(t))

def sha256(M):
    block_size = 512 #512 bits

    M = padding(M, block_size)

    # M is an array of 32 bit integers
    M_bit_len = len(M) * 8

    N = int(M_bit_len/block_size)

    a = sha256_H[0]
    b = sha256_H[1]
    c = sha256_H[2]
    d = sha256_H[3]
    e = sha256_H[4]
    f = sha256_H[5]
    g = sha256_H[6]
    h = sha256_H[7]

    for i in range(0, N):

        block_bytes = int(512/8)
        Mi = M[i*block_bytes : ((i+1)*block_bytes)]

        for t in range(0, 64):
            T1 = np.uint32(h + SIGMA_1_256(e) + ch(e, f, g) + K256[t] + W_sha256(t, Mi))
            T2 = np.uint32(SIGMA_0_256(a) + maj(a, b, c))
            h = np.uint32(g)
            g = np.uint32(f)
            f = np.uint32(e)
            e = np.uint32(d + T1)
            d = np.uint32(c)
            c = np.uint32(b)
            b = np.uint32(a)
            a = np.uint32(T1 + T2)

    return [np.uint32(sha256_H[0]+a),
            np.uint32(sha256_H[1]+b),
            np.uint32(sha256_H[2]+c),
            np.uint32(sha256_H[3]+d),
            np.uint32(sha256_H[4]+e),
            np.uint32(sha256_H[5]+f),
            np.uint32(sha256_H[6]+g),
            np.uint32(sha256_H[7]+h)]

def ch (x :np.uint32, y: np.uint32, z: np.uint32) -> np.uint32:
    neg_x = 0xffffffff ^ x #xor with a 32-bit string of 1's is effectively bitwise not
                           # hack for python, since ~ does no bit-wise negate 32 bit integers.
    return (x&y)^(neg_x & z)

def maj (x :np.uint32, y: np.uint32, z: np.uint32) -> np.uint32:
    return ((x&y) ^ (x&z) ^ (y&z))

def rotr(x, i):
    x32 = np.uint32(x)
    i = np.uint32(i)
    bottom = np.right_shift(x32, i)
    top = np.left_shift(x32, np.uint32(32-i))

    return top | bottom

def padding(x: list, pad_to):

    # get number of words to pad to "pad_to" bits
    pad_to_bytes = pad_to/8

    # The first bit after the message end is set to 1
    pad_size = pad_to_bytes - (len(x) + 1)

    ret_list = x + [128] + [0] * int(pad_size)

    ret_list[-1] = len(x) * 8
    return ret_list

def SIGMA_0_256 (x):
    x = np.uint32(x)
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def SIGMA_1_256 (x):
    x = np.uint32(x)
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def sigma_0_256 (x):
    x = np.uint32(x)
    return rotr(x, 7) ^ rotr(x, 18) ^ np.right_shift(x, np.uint32(3))

def sigma_1_256 (x):
    x = np.uint32(x)
    rotr17 = rotr(x, 17)
    rotr19 = rotr(x, 19)
    rshift10 = np.right_shift(x, np.uint32(10))

    return rotr17 ^ rotr19 ^ rshift10

def string_to_int_list(in_str):
    in_len = len(in_str)
    out_list = []

    for i in range(0, in_len, 2):
        out_list = out_list + [int(in_str[i : i+2], base=16)]

    return out_list

if __name__ == "__main__":
    data_in = sys.argv[1]
    data_in = "68656c6c6f20776f726c64"
    array_32bit_int = string_to_int_list(data_in)

    sha256_out = sha256(array_32bit_int);
    out = ''.join('{:04x} '.format(x) for x in sha256_out)
    print("{}".format(out));