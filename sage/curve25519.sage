# prime field = 2^255 - 19
# order = (2^252 + 27742317777372353535851937790883648493) * 8
#       = 2^255 + 221938542218978828286815502327069187944

f = GF(2**255-19)
ec = EllipticCurve(GF(2**255-19), [0,486662,0,1,0])
base_point = ec.lift_x(9)
point_at_infinity = ec(0)


print(f.order() % 8)

p = 2**255-19
r = (p + 3) / 8

def curve25519_sqrt(a):
    p = 2**255 - 19
    # Since p % 8 == 5, we use this formula
    r = pow(a, (p + 3) // 8, p)
    if pow(r, 2, p) != a % p:
        # Multiply by sqrt(-1)
        r = (r * pow(2, (p - 1) // 4, p)) % p

    # Final check if root is valid
    if pow(r, 2, p) != a % p:
        return None  # No root
    return r

#print(p)
#print(r)
