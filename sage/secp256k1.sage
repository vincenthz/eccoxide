p = 2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1
a = 0
b = 7

F = FiniteField (p)
E = EllipticCurve([F(a), F(b)])
order = E.order()

# print('%x' % order)

gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
G = E(gx, gy)

M = FiniteField(order);

x = E.gens()

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
    #print(print_n(x, 32))
    #print('%x%x' % (x,y))


print("struct KAT { n: u64, x: [u8;32], y: [u8;32] }")
print("const KATS : [KAT; 100] = [")
for x in range(1, 101):
    p = x * G
    print("KAT {")
    print("n: %d," % x)
    print_compress(p, 32)
    print("}, ")
print("];")
