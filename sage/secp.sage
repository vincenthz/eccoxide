# create a SEC2 prime curve "object"
def sec(name, p, ps, a, b, gx, gy, size):
    F = FiniteField (p)
    a = F(a)
    b = F(b)
    E = EllipticCurve([a, b])
    G = E(gx, gy)
    #order = E.order()
    # scalar finite field
    #S = FiniteField(order)
    return { 'name': name, 'p' : p, 'ps': ps, 'E': E, 'G': G, 'size': size }

def secp192r1():
    p = 2^192 - 2^64 - 1
    ps = "2^192 - 2^64 - 1"
    a = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFC
    b = 0x64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1
    gx = 0x188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012
    gy = 0x07192B95FFC8DA78631011ED6B24CDD573F977A11E794811
    return sec("p192r1", p, ps, a, b, gx, gy, 24)

def secp224k1():
    p = 2^224 - 2^32 - 2^12 - 2^11 - 2^9 - 2^7 - 2^4 - 2 - 1
    ps = "2^224 - 2^32 - 2^12 - 2^11 - 2^9 - 2^7 - 2^4 - 2 - 1"
    a = 0
    b = 5
    gx = 0xA1455B334DF099DF30FC28A169A467E9E47075A90F7E650EB6B7A45C
    gy = 0x7E089FED7FBA344282CAFBD6F7E319F7C0B0BD59E2CA4BDB556D61A5
    return sec("p224k1", p, ps, a, b, gx, gy, 28)

def secp224r1():
    p = 2^224 - 2^96 + 1
    ps = "2^224 - 2^96 + 1"
    a = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFE
    b = 0xB4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4
    gx = 0xB70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21
    gy = 0xBD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34
    return sec("p224r1", p, ps, a, b, gx, gy, 28)

def secp256k1():
    p = 2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1
    ps = "2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1"
    a = 0
    b = 7
    gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
    gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
    return sec("p256k1", p, ps, a, b, gx, gy, 32)

def secp256r1():
    # p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
    p = 2^256 - 2^224 + 2^192 + 2^96 - 1
    ps = "2^256 - 2^224 + 2^192 + 2^96 - 1"
    a = 0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc
    b = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
    gx = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
    gy = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
    return sec("p256r1", p, ps, a, b, gx, gy, 32)

def secp384r1():
    p = 2^384 - 2^128 - 2^96 + 2^32 - 1
    ps = "2^384 - 2^128 - 2^96 + 2^32 - 1"
    a = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000fffffffc
    b = 0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef
    gx = 0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7
    gy = 0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f
    return sec("p384r1", p, ps, a, b, gx, gy, 48)

def secp521r1():
    p = 2^521 - 1
    ps = "2^521 - 1"
    a = 0x01fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffc
    b = 0x0051953eb9618e1c9a1f929a21a0b68540eea2da725b99b315f3b8b489918ef109e156193951ec7e937b1652c0bd3bb1bf073573df883d2c34f1ef451fd46b503f00
    gx = 0x00c6858e06b70404e9cd9e3ecb662395b4429c648139053fb521f828af606b4d3dbaa14b5e77efe75928fe1dc127a2ffa8de3348b3c1856a429bf97e7e31c2e5bd66
    gy = 0x011839296a789a3bc0045c8a5fb42c7d1bd998f54449579b446817afbd17273e662c97ee72995ef42640c550b9013fad0761353c7086a272c24088be94769fd16650
    return sec("p521r1", p, ps, a, b, gx, gy, 66)


# get a generator of E
#x = E.gens()

def print_n(n, z):
    s = format('%x' % n);
    while len(s) < z * 2:
        s = "0" + s
    x = "0x" + ", 0x".join([ s[i*2:i*2+2] for i in range(0, len(s) / 2) ])
    return "[" + x + "]"

def print_compress(p, sz):
    x, y = p.xy()
    print("x: %s," % print_n(x, sz))
    print("y: %s," % print_n(y, sz))

def print_kats(d):
    G = d['G'];
    size = d['size'];
    print("struct KAT { n: u64, x: [u8;%d], y: [u8;%d] }" % (size, size))
    print("const KATS : [KAT; 100] = [")
    for x in range(1, 10):
        p = x * G
        print("KAT {")
        print("n: %d," % x)
        print_compress(p, size)
        print("}, ")
    print("];")

#print_kats(secp256k1())
#print_kats(secp256r1())
#print_kats(secp384r1())
#print_kats(secp521r1())

def print_docs(d):
    G = d['G'];
    size = d['size'];
    print("//! Curve %s as defined over the prime field of order %s" % (d['name'], d['ps']))

print_docs(secp256k1())
print_docs(secp256r1())
print_docs(secp384r1())
print_docs(secp521r1())
