# Derive the BLS12-381 constants behind the G1 / G2 subgroup membership checks,
# cofactor clearing and the hash-to-curve square roots.
#
# Run with:   sage sage/bls12_381.sage        (plain `python3` works too)
#
# Prints, as ready-to-paste Rust, the byte constants of src/params/bls12_381.rs
# that `is_in_subgroup` needs:
#
#     BETA_BYTES              beta, a primitive cube root of unity in Fp
#     g2::PSI_X_COEFF_BYTES   xi^-((p-1)/3)
#     g2::PSI_Y_COEFF_BYTES   xi^-((p-1)/2)
#
# together with the xi powers those two invert, already in the tree as the
# Frobenius coefficients FROBENIUS_FP12_C1_BYTES (xi^((p-1)/6)),
# FROBENIUS_FP6_C1_BYTES (xi^((p-1)/3)) and FROBENIUS_FP6_C2_BYTES
# (xi^(2(p-1)/3)), so the whole family can be eyeballed at once; the ones of
# src/params/bls12_381_h2c.rs that the optimized `sqrt_ratio` of RFC 9380
# appendix F.2.1 needs, and which the RFC does not publish:
#
#     g1::SQRT_RATIO_C2_BYTES  sqrt(-Z)               (appendix F.2.1.2)
#     g2::SQRT_RATIO_C3_BYTES  ((q-1)/2^c1 - 1)/2     (appendix F.2.1.1)
#     g2::SQRT_RATIO_C6_BYTES  Z^((q-1)/2^c1)
#     g2::SQRT_RATIO_C7_BYTES  Z^(((q-1)/2^c1 + 1)/2)
#
# and the test vectors of the Rust `mod tests::subgroup`: points off the
# subgroups, and the G2 h_eff. `cargo fmt` re-wraps the emitted arrays.
#
# Nothing is taken on faith: this asserts the seed relations
# (r = x^4 - x^2 + 1, p == x mod r, the two cofactors), that beta is the cube
# root of unity for which sigma(x, y) = (beta*x, y) acts on G1 as [-x^2] -- the
# other one gives [x^2 - 1], which is why the choice matters -- that
# psi(x, y) = (x^p * cx, y^p * cy) maps the twist to itself and acts on G2 as
# [x], and that both endomorphism checks agree with the definition [r]P == O,
# on points of the subgroups and on points outside them. The cofactor-clearing
# maps are checked the same way: that each one lands in its subgroup for any
# input point, and that it is multiplication by the h_eff of RFC 9380 section
# 8.8. Each sqrt_ratio is run against its own definition, on inputs of either
# quadratic character. The Rust test suite re-checks the same relations against
# the bytes as committed (`g1::tests::subgroup`, `g2::tests::subgroup`,
# `hash_to_curve::tests::sqrt_ratio_matches_the_definition`).
#
# Only +, -, *, // and pow are used, so this runs under plain python3 as well as
# sage, like the pure-python table generators in comb.sage.

# --- parameters -------------------------------------------------------------

# base field prime and group order, both derived from the seed x below
p = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
r = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001

# the seed is negative; the Rust side carries it as BLS_X = -x with
# BLS_X_IS_NEGATIVE, since every use iterates over the bits of the magnitude
x = -0xD201000000010000

# G1: y^2 = x^3 + 4 over Fp
g1_b = 4
g1_gx = 0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB
g1_gy = 0x08B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1

# G2: y^2 = x^3 + 4(1 + u) over Fp2 = Fp[u]/(u^2 + 1), elements as (c0, c1)
g2_b = (4, 4)
g2_gx = (0x024AA2B2F08F0A91260805272DC51051C6E47AD4FA403B02B4510B647AE3D1770BAC0326A805BBEFD48056C8C121BDB8,
         0x13E02B6052719F607DACD3A088274F65596BD0D09920B61AB5DA61BBDC7F5049334CF11213945D57E5AC7D055D042B7E)
g2_gy = (0x0CE5D527727D6E118CC9CDC6DA2E351AADFD9BAA8CBDD3A76D429A695160D12C923AC9CC3BACA289E193548608B82801,
         0x0606C4A02EA734CC32ACD2B02BC28B99CB3E287E85A763AF267492AB572E99AB3F370D275CEC1DA1AAA9075FF05F79BE)

# the tower non-residue xi = 1 + u, of Fp6 = Fp2[v]/(v^3 - xi), Fp12 = Fp6[w]/(w^2 - v)
xi = (1, 1)

FIELD_BYTES = 48

# --- Fp and Fp2 -------------------------------------------------------------

def inv(a):
    return pow(a, p - 2, p)


def f2_add(a, b):
    return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)


def f2_sub(a, b):
    return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)


def f2_mul(a, b):
    return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)


def f2_inv(a):
    t = inv((a[0] * a[0] + a[1] * a[1]) % p)
    return ((a[0] * t) % p, ((-a[1]) * t) % p)


def f2_pow(a, e):
    acc, base = (1, 0), a
    while e > 0:
        if e % 2 == 1:
            acc = f2_mul(acc, base)
        base = f2_mul(base, base)
        e = e // 2
    return acc


def f2_frobenius(a):
    # a^p is the conjugation, the non-trivial automorphism of Fp2/Fp
    return (a[0], (-a[1]) % p)


# The two fields as a uniform interface, so the curve code below serves both.
class Fp:
    zero, one = 0, 1
    add = staticmethod(lambda a, b: (a + b) % p)
    sub = staticmethod(lambda a, b: (a - b) % p)
    mul = staticmethod(lambda a, b: (a * b) % p)
    inv = staticmethod(inv)


class Fp2:
    zero, one = (0, 0), (1, 0)
    add = staticmethod(f2_add)
    sub = staticmethod(f2_sub)
    mul = staticmethod(f2_mul)
    inv = staticmethod(f2_inv)


# --- affine short-Weierstrass arithmetic (a = 0), None being the identity ----

def curve(F, b):
    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if F.add(y1, y2) == F.zero:
                return None
            xx = F.mul(x1, x1)
            num = F.add(F.add(xx, xx), xx)  # 3x^2, the curve having a = 0
            den = F.add(y1, y1)
        else:
            num = F.sub(y2, y1)
            den = F.sub(x2, x1)
        l = F.mul(num, F.inv(den))
        x3 = F.sub(F.sub(F.mul(l, l), x1), x2)
        y3 = F.sub(F.mul(l, F.sub(x1, x3)), y1)
        return (x3, y3)

    def neg(P):
        return None if P is None else (P[0], F.sub(F.zero, P[1]))

    def mul(k, P):
        negate, k = k < 0, abs(k)
        acc, base = None, P
        while k > 0:
            if k % 2 == 1:
                acc = add(acc, base)
            base = add(base, base)
            k = k // 2
        return neg(acc) if negate else acc

    def on_curve(P):
        if P is None:
            return True
        cx, cy = P
        x3 = F.mul(F.mul(cx, cx), cx)
        return F.mul(cy, cy) == F.add(x3, b)

    return add, mul, on_curve


g1_add, g1_mul, g1_on_curve = curve(Fp, g1_b)
g2_add, g2_mul, g2_on_curve = curve(Fp2, g2_b)

G1 = (g1_gx, g1_gy)
G2 = (g2_gx, g2_gy)

# --- points off the subgroups, to check the rejecting side of each test -----

def g1_points_off_subgroup(count):
    # decompress small x-coordinates: the 126-bit cofactor makes it a
    # certainty that none of them lands in G1
    out, cx = [], 1
    while len(out) < count:
        yy = (cx * cx * cx + g1_b) % p
        if pow(yy, (p - 1) // 2, p) == 1:
            P = (cx, pow(yy, (p + 1) // 4, p))  # p = 3 mod 4
            assert g1_on_curve(P)
            if g1_mul(r, P) is not None:
                out.append(P)
        cx += 1
    return out


def f2_sqrt(a):
    # the complex method (Alg. 9 of eprint 2012/685), valid for p = 3 mod 4
    a1 = f2_pow(a, (p - 3) // 4)
    alpha = f2_mul(f2_mul(a1, a1), a)
    x0 = f2_mul(a1, a)
    if alpha == ((p - 1) % p, 0):
        root = f2_mul(x0, (0, 1))
    else:
        root = f2_mul(f2_pow(f2_add((1, 0), alpha), (p - 1) // 2), x0)
    return root if f2_mul(root, root) == a else None


def g2_points_off_subgroup(count):
    # likewise on the twist, whose cofactor is 507 bits wide
    out, i = [], 1
    while len(out) < count:
        cx = (i, i + 1)
        yy = f2_add(f2_mul(f2_mul(cx, cx), cx), g2_b)
        cy = f2_sqrt(yy)
        if cy is not None:
            Q = (cx, cy)
            assert g2_on_curve(Q)
            if g2_mul(r, Q) is not None:
                out.append(Q)
        i += 1
    return out


# --- the seed relations -----------------------------------------------------

assert r == x**4 - x**2 + 1, "r is not the seed polynomial"
assert (p - x) % r == 0, "p != x (mod r), which both checks rely on"

t = x + 1  # trace of Frobenius
h1 = (p + 1 - t) // r  # #E(Fp) = p + 1 - t
assert h1 * r == p + 1 - t
assert h1 == (x - 1)**2 // 3, "G1 cofactor is not (x-1)^2/3"
# The G2 cofactor is that of the sextic *twist*, not of E(Fp2): the two curves
# have different orders, and (p^2 + 1 - (t^2 - 2p))/r -- the cofactor of E(Fp2),
# whose leading digits are deceptively close -- does not annihilate the twist.
# This is the standard closed form, checked against the twist below.
h2 = (x**8 - 4 * x**7 + 5 * x**6 - 4 * x**4 + 6 * x**3 - 4 * x**2 - 4 * x + 13) // 9
assert 9 * h2 == x**8 - 4 * x**7 + 5 * x**6 - 4 * x**4 + 6 * x**3 - 4 * x**2 - 4 * x + 13

assert g1_on_curve(G1) and g1_mul(r, G1) is None, "G1 generator"
assert g2_on_curve(G2) and g2_mul(r, G2) is None, "G2 generator"
# the group order kills every point of the twist, cofactor included
for Q in g2_points_off_subgroup(3):
    assert g2_mul(h2 * r, Q) is None, "h2 is not the cofactor of the twist"
assert g2_mul(h2, G2) is not None, "h2 must not kill G2 itself"

print("p           %d bits" % p.bit_length())
print("r           %d bits, = x^4 - x^2 + 1" % r.bit_length())
print("x           -0x%x (negative), p == x (mod r)" % (-x))
print("h1          %d bits, = (x-1)^2/3 = 0x%x" % (h1.bit_length(), h1))
print("h2          %d bits, twist cofactor = 0x%x" % (h2.bit_length(), h2))

# --- beta: the G1 endomorphism sigma(x, y) = (beta*x, y) --------------------

# the two non-trivial cube roots of unity in Fp, beta and beta^2
g = 2
while pow(g, (p - 1) // 3, p) == 1:
    g += 1
roots = [pow(g, (p - 1) // 3, p)]
roots.append((roots[0] * roots[0]) % p)
for b in roots:
    assert b != 1 and pow(b, 3, p) == 1
assert (roots[0] + roots[1] + 1) % p == 0, "1 + beta + beta^2 != 0"


def sigma(P, beta):
    return None if P is None else ((beta * P[0]) % p, P[1])


# sigma acts on G1 as one of the two primitive cube roots of unity mod r, and
# which one depends on the root picked: -x^2 and x^2 - 1 are the two, being
# each other's square mod r. Only [-x^2] is what the Rust check compares to.
lambda_minus_x2 = (-x * x) % r
lambda_x2_minus_1 = (x * x - 1) % r
assert (lambda_minus_x2 * lambda_minus_x2) % r == lambda_x2_minus_1
assert (lambda_minus_x2**2 + lambda_minus_x2 + 1) % r == 0

eigenvalues = {}
for b in roots:
    S = sigma(G1, b)
    assert g1_on_curve(S)
    if S == g1_mul(lambda_minus_x2, G1):
        eigenvalues[b] = "[-x^2]"
    elif S == g1_mul(lambda_x2_minus_1, G1):
        eigenvalues[b] = "[x^2-1]"
    else:
        raise AssertionError("sigma is not a cube root of unity on G1")
assert sorted(eigenvalues.values()) == ["[-x^2]", "[x^2-1]"]
BETA = [b for b, e in eigenvalues.items() if e == "[-x^2]"][0]

# on the whole subgroup, not just on the generator
for k in (1, 2, 3, 12345, 0xDEADBEEF, r - 1):
    P = g1_mul(k, G1)
    assert sigma(P, BETA) == g1_mul(lambda_minus_x2, P), "sigma != [-x^2] on [%d]G1" % k

# and the check rejects exactly the points the definition rejects
for P in g1_points_off_subgroup(8):
    in_g1 = g1_mul(r, P) is None
    assert not in_g1
    assert (sigma(P, BETA) == g1_mul(lambda_minus_x2, P)) == in_g1, "G1 check accepted a stray point"
    # a G1 point with a stray component added -- what an attacker would send
    Q = g1_add(G1, P)
    assert g1_mul(r, Q) is not None
    assert sigma(Q, BETA) != g1_mul(lambda_minus_x2, Q), "G1 check accepted G + stray"

print("beta        0x%x  (the other root gives [x^2-1])" % BETA)

# --- psi: the G2 untwist-Frobenius-twist endomorphism ----------------------

# psi = phi . frobenius . phi^-1 for the twist map phi(x, y) = (x*w^2, y*w^3),
# and w^(2(p-1)) = xi^((p-1)/3), w^(3(p-1)) = xi^((p-1)/2) since w^6 = xi
PSI_X_COEFF = f2_inv(f2_pow(xi, (p - 1) // 3))
PSI_Y_COEFF = f2_inv(f2_pow(xi, (p - 1) // 2))


def psi(Q):
    if Q is None:
        return None
    cx, cy = Q
    return (f2_mul(f2_frobenius(cx), PSI_X_COEFF), f2_mul(f2_frobenius(cy), PSI_Y_COEFF))


# psi acts as [p] on G2, and p == x (mod r), which is the whole point
for k in (1, 2, 3, 12345, 0xDEADBEEF, r - 1):
    Q = g2_mul(k, G2)
    assert g2_on_curve(psi(Q)), "psi([%d]G2) fell off the twist" % k
    assert psi(Q) == g2_mul(p % r, Q), "psi != [p] on [%d]G2" % k
    assert psi(Q) == g2_mul(x % r, Q), "psi != [x] on [%d]G2" % k

for Q in g2_points_off_subgroup(8):
    in_g2 = g2_mul(r, Q) is None
    assert not in_g2
    assert g2_on_curve(psi(Q)), "psi is not an endomorphism of the twist"
    assert (psi(Q) == g2_mul(x % r, Q)) == in_g2, "G2 check accepted a stray point"
    S = g2_add(G2, Q)
    assert g2_mul(r, S) is not None
    assert psi(S) != g2_mul(x % r, S), "G2 check accepted G + stray"

# the coefficients are the inverses of the Frobenius ones already in the tree
FROBENIUS_FP12_C1 = f2_pow(xi, (p - 1) // 6)
FROBENIUS_FP6_C1 = f2_pow(xi, (p - 1) // 3)
FROBENIUS_FP6_C2 = f2_pow(xi, 2 * (p - 1) // 3)
assert f2_mul(PSI_X_COEFF, FROBENIUS_FP6_C1) == (1, 0)
assert f2_mul(PSI_Y_COEFF, f2_pow(FROBENIUS_FP12_C1, 3)) == (1, 0)

print("psi_x       0x%x + 0x%x*u" % PSI_X_COEFF)
print("psi_y       0x%x + 0x%x*u" % PSI_Y_COEFF)

# --- cofactor clearing ------------------------------------------------------
#
# Both maps send an arbitrary point of the curve (resp. the twist) into the
# prime-order subgroup, which is what hash-to-curve needs after the map-to-curve
# step. Neither multiplies by the cofactor itself: G1 uses the much shorter
# 1 - x (Wahby-Boneh, eprint 2019/403), which works because E(Fp) is not cyclic,
# and G2 uses the Budroni-Pintore endomorphism chain (eprint 2017/419), whose
# result is [3*h2]Q. Both are the h_eff of RFC 9380 section 8.8.

G1_H_EFF = 1 - x  # = 0xd201000000010001, the RFC 9380 h_eff for G1
assert G1_H_EFF == 0xD201000000010001


def g1_clear_cofactor(P):
    return g1_mul(G1_H_EFF, P)


def g2_clear_cofactor(Q):
    # [x^2 - x - 1]Q + [x - 1]psi(Q) + [2]psi^2(Q), evaluated as
    # psi^2([2]Q) + [x]([x]Q + psi(Q)) - [x]Q - psi(Q) - Q
    t1 = g2_mul(x, Q)
    t2 = psi(Q)
    acc = psi(psi(g2_mul(2, Q)))
    acc = g2_add(acc, g2_mul(x, g2_add(t1, t2)))
    acc = g2_add(acc, g2_mul(-1, t1))
    acc = g2_add(acc, g2_mul(-1, t2))
    return g2_add(acc, g2_mul(-1, Q))


# G1: [1 - x] lands in the subgroup whatever the input point, and on a point
# that is already in G1 it is plain multiplication by 1 - x
for P in g1_points_off_subgroup(8) + [G1, g1_mul(12345, G1), g1_add(G1, g1_points_off_subgroup(1)[0])]:
    cleared = g1_clear_cofactor(P)
    assert g1_on_curve(cleared)
    assert g1_mul(r, cleared) is None, "G1 cofactor clearing left the subgroup"
assert g1_clear_cofactor(G1) == g1_mul(G1_H_EFF % r, G1)
# it is not the identity map on G1, and it does not collapse stray points to O
assert g1_clear_cofactor(G1) != G1
assert g1_clear_cofactor(g1_points_off_subgroup(1)[0]) is not None

# G2: the endomorphism chain is multiplication by 3(x^2 - 1)*h2 -- a multiple of
# the cofactor, not the cofactor itself, which is just as good and is the h_eff
# RFC 9380 section 8.8.2 specifies
G2_H_EFF = 3 * (x * x - 1) * h2
assert G2_H_EFF == 0xBC69F08F2EE75B3584C6A0EA91B352888E2A8E9145AD7689986FF031508FFE1329C2F178731DB956D82BF015D1212B02EC0EC69D7477C1AE954CBC06689F6A359894C0ADEBBF6B4E8020005AAA95551, "not the RFC 9380 h_eff"
for Q in g2_points_off_subgroup(4) + [G2, g2_mul(12345, G2), g2_add(G2, g2_points_off_subgroup(1)[0])]:
    cleared = g2_clear_cofactor(Q)
    assert g2_on_curve(cleared)
    assert cleared == g2_mul(G2_H_EFF, Q), "the psi chain is not [h_eff]"
    assert g2_mul(r, cleared) is None, "G2 cofactor clearing left the subgroup"
assert g2_clear_cofactor(G2) is not None and g2_clear_cofactor(G2) != G2

# on the subgroup psi is [x], so the chain collapses to a scalar the Rust tests
# can rebuild from the seed alone
assert G2_H_EFF % r == (4 * x * x - 2 * x - 1) % r

print("h_eff(G1)   0x%x = 1 - x" % G1_H_EFF)
print("h_eff(G2)   %d bits = 3(x^2-1)*h2, == 4x^2 - 2x - 1 (mod r) on the subgroup"
      % G2_H_EFF.bit_length())

# --- hash-to-curve: the sqrt_ratio constants of RFC 9380 appendix F.2.1 ------
#
# Simplified SWU needs sqrt_ratio(u, v) = (is_square(u/v), sqrt(u/v)), answering
# with sqrt(Z*(u/v)) when u/v is not a square. Spelled out, that is an inversion
# and two square roots; the RFC's shortcuts each cost a single exponentiation
# plus a handful of multiplications, in exchange for the constants derived here.
#
#   Fp,  q = p = 3 mod 4      appendix F.2.1.2, needing c1 = (p-3)/4 -- already
#                             in the tree as P_MINUS3_DIV4_BYTES -- and
#                             c2 = sqrt(-Z)
#   Fp2, q = p^2 = 9 mod 16   appendix F.2.1.1, the any-field variant, whose
#                             constants follow from c1 = v_2(q-1) = 3
#
# Z is the non-square the map is built around, 11 for G1 and -(2 + u) for G2
# (RFC section 8.8); it has to be the same Z on both sides, so the Rust code
# takes all of these from one place per field.

G1_Z = 11
G2_Z = ((-2) % p, (-1) % p)
q2 = p * p

assert p % 4 == 3, "Fp is not the 3 mod 4 case"
assert q2 % 16 == 9, "Fp2 is not the 9 mod 16 case"
assert pow(G1_Z, (p - 1) // 2, p) == p - 1, "Z is a square in Fp"
assert f2_sqrt(G2_Z) is None, "Z is a square in Fp2"

# Fp: c2 is either square root of -Z (-Z is a square because both -1 and Z are
# not); the two differ by a sign that sqrt_ratio's caller normalizes away with
# sgn0, so pin the even one for reproducibility
SQRT_RATIO_C2 = pow((-G1_Z) % p, (p + 1) // 4, p)
assert SQRT_RATIO_C2 * SQRT_RATIO_C2 % p == (-G1_Z) % p
if SQRT_RATIO_C2 % 2 == 1:
    SQRT_RATIO_C2 = p - SQRT_RATIO_C2

# Fp2: c1 is the 2-adic valuation of q - 1 and the rest follows from it. c4, c5
# and the loop trip count are small enough to sit in the Rust source as
# literals; c3, c6 and c7 become byte constants.
F2_C1 = 0
while (q2 - 1) % (1 << (F2_C1 + 1)) == 0:
    F2_C1 += 1
F2_C2 = (q2 - 1) // (1 << F2_C1)
SQRT_RATIO_C3 = (F2_C2 - 1) // 2
F2_C4 = (1 << F2_C1) - 1
F2_C5 = 1 << (F2_C1 - 1)
SQRT_RATIO_C6 = f2_pow(G2_Z, F2_C2)
SQRT_RATIO_C7 = f2_pow(G2_Z, (F2_C2 + 1) // 2)
assert (F2_C1, F2_C4, F2_C5) == (3, 7, 4), "the Rust literals no longer match"
# c6 is a primitive 2^c1-th root of unity, which is what makes the loop of steps
# 17-26 converge, and c7 squares to it
assert f2_pow(SQRT_RATIO_C6, 1 << F2_C1) == (1, 0)
assert f2_pow(SQRT_RATIO_C6, 1 << (F2_C1 - 1)) != (1, 0)
assert f2_mul(SQRT_RATIO_C7, SQRT_RATIO_C7) == f2_mul(SQRT_RATIO_C6, G2_Z)


def fp_sqrt_ratio(u, v):
    # appendix F.2.1.2, step for step
    tv1 = v * v % p
    tv2 = u * v % p
    tv1 = tv1 * tv2 % p
    y1 = pow(tv1, (p - 3) // 4, p)
    y1 = y1 * tv2 % p
    y2 = y1 * SQRT_RATIO_C2 % p
    tv3 = y1 * y1 % p
    tv3 = tv3 * v % p
    is_qr = tv3 == u
    return is_qr, (y1 if is_qr else y2)


def fp2_sqrt_ratio(u, v):
    # appendix F.2.1.1, step for step, with the u = 0 guard noted below
    tv1 = SQRT_RATIO_C6
    tv2 = f2_pow(v, F2_C4)
    tv3 = f2_mul(tv2, tv2)
    tv3 = f2_mul(tv3, v)
    tv5 = f2_mul(u, tv3)
    tv5 = f2_pow(tv5, SQRT_RATIO_C3)
    tv5 = f2_mul(tv5, tv2)
    tv2 = f2_mul(tv5, v)
    tv3 = f2_mul(tv5, u)
    tv4 = f2_mul(tv3, tv2)
    tv5 = f2_pow(tv4, F2_C5)
    # the published steps read `isQR = tv5 == 1`, which answers False for u = 0
    # even though 0 is a square; every intermediate is 0 there, so the test has
    # nothing to go on and the zero case has to be added back
    is_qr = tv5 == (1, 0) or u == (0, 0)
    tv2 = f2_mul(tv3, SQRT_RATIO_C7)
    tv5 = f2_mul(tv4, tv1)
    tv3 = tv3 if is_qr else tv2
    tv4 = tv4 if is_qr else tv5
    for i in range(F2_C1, 1, -1):
        tv5 = f2_pow(tv4, 1 << (i - 2))
        e1 = tv5 == (1, 0)
        tv2 = f2_mul(tv3, tv1)
        tv1 = f2_mul(tv1, tv1)
        tv5 = f2_mul(tv4, tv1)
        tv3 = tv3 if e1 else tv2
        tv4 = tv4 if e1 else tv5
    return is_qr, tv3


# both agree with the definition, on inputs of either quadratic character and on
# the u = 0 corner; the SSWU caller never hands over v = 0
squares = 0
for i in range(1, 65):
    u, v = pow(7, i, p), pow(11, i + 3, p)
    is_qr, y = fp_sqrt_ratio(u, v)
    ratio = u * inv(v) % p
    assert is_qr == (pow(ratio, (p - 1) // 2, p) == 1), "Fp: wrong verdict"
    assert y * y % p == (ratio if is_qr else ratio * G1_Z % p), "Fp: wrong root"

    u2 = (pow(5, i, p), pow(13, i + 1, p))
    v2 = (pow(3, i + 2, p), pow(17, i, p))
    is_qr2, y2 = fp2_sqrt_ratio(u2, v2)
    ratio2 = f2_mul(u2, f2_inv(v2))
    assert is_qr2 == (f2_sqrt(ratio2) is not None), "Fp2: wrong verdict"
    want = ratio2 if is_qr2 else f2_mul(G2_Z, ratio2)
    assert f2_mul(y2, y2) == want, "Fp2: wrong root"
    squares += is_qr + is_qr2
assert 20 < squares < 108, "the sample missed one of the two branches"
assert fp_sqrt_ratio(0, 7) == (True, 0)
assert fp2_sqrt_ratio((0, 0), (7, 3)) == (True, (0, 0))

print("sqrt_ratio  Fp: F.2.1.2 (3 mod 4); Fp2: F.2.1.1 with c1 = %d" % F2_C1)
print("all checks passed\n")

# --- emit the Rust byte arrays ---------------------------------------------

def rust_array(name, raw, indent, doc):
    lines = ["%s/// %s" % (indent, doc)]
    lines.append("%spub const %s: [u8; %d] = [" % (indent, name, len(raw)))
    per_line = 16
    for i in range(0, len(raw), per_line):
        chunk = ", ".join("0x%02x" % b for b in raw[i:i + per_line])
        lines.append("%s    %s," % (indent, chunk))
    lines.append("%s];" % indent)
    return "\n".join(lines)


def fp_bytes(a):
    return int(a).to_bytes(FIELD_BYTES, "big")


def fp2_bytes(a):
    # the imaginary component leads, as `Fp2::to_bytes` does
    return int(a[1]).to_bytes(FIELD_BYTES, "big") + int(a[0]).to_bytes(FIELD_BYTES, "big")


print("// src/params/bls12_381.rs, top level:\n")
print(rust_array("BETA_BYTES", fp_bytes(BETA), "", "beta, a primitive cube root of unity in Fp"))
print("\n// src/params/bls12_381.rs, inside `pub mod g2`:\n")
print(rust_array("PSI_X_COEFF_BYTES", fp2_bytes(PSI_X_COEFF), "    ", "xi^-((p-1)/3)"))
print()
print(rust_array("PSI_Y_COEFF_BYTES", fp2_bytes(PSI_Y_COEFF), "    ", "xi^-((p-1)/2)"))
print("\n// already in the tree, from the same xi powers:\n")
print(rust_array("FROBENIUS_FP12_C1_BYTES", fp2_bytes(FROBENIUS_FP12_C1), "", "xi^((p-1)/6)"))
print()
print(rust_array("FROBENIUS_FP6_C1_BYTES", fp2_bytes(FROBENIUS_FP6_C1), "", "xi^((p-1)/3)"))
print()
print(rust_array("FROBENIUS_FP6_C2_BYTES", fp2_bytes(FROBENIUS_FP6_C2), "", "xi^(2(p-1)/3)"))

print("\n// src/params/bls12_381_h2c.rs, inside `pub mod g1`:\n")
print(rust_array("SQRT_RATIO_C2_BYTES", fp_bytes(SQRT_RATIO_C2), "    ",
                 "`c2 = sqrt(-Z)`, the appendix F.2.1.2 constant"))
print("\n// src/params/bls12_381_h2c.rs, inside `pub mod g2`:\n")
print(rust_array("SQRT_RATIO_C3_BYTES", int(SQRT_RATIO_C3).to_bytes(2 * FIELD_BYTES, "big"), "    ",
                 "`c3 = ((q - 1)/2^c1 - 1)/2` with `q = p^2`, `c1 = 3`"))
print()
print(rust_array("SQRT_RATIO_C6_BYTES", fp2_bytes(SQRT_RATIO_C6), "    ",
                 "`c6 = Z^((q - 1)/2^c1)`, a primitive 8th root of unity"))
print()
print(rust_array("SQRT_RATIO_C7_BYTES", fp2_bytes(SQRT_RATIO_C7), "    ",
                 "`c7 = Z^(((q - 1)/2^c1 + 1)/2)`, a square root of `c6*Z`"))

# --- test vectors: on the curve, outside the prime-order subgroup -----------
#
# The standard encodings of a few such points, for the Rust `off_subgroup_kat`
# tests: `from_compressed` / `from_uncompressed` must reject them and the
# `_oncurve_only` decoders must accept them. The flags are computed here from
# the format definition (bit 7 compression, bit 5 "y is the larger root",
# comparing c1 first for an Fp2 coordinate) rather than by the Rust encoder.

COMPRESSION_FLAG = 0x80
SORT_FLAG = 0x20


def fp_is_largest(a):
    return a > (p - 1) // 2


def fp2_is_largest(a):
    return fp_is_largest(a[1]) or (a[1] == 0 and fp_is_largest(a[0]))


def encodings(P, to_bytes, is_largest):
    cx, cy = P
    compressed = bytearray(to_bytes(cx))
    compressed[0] |= COMPRESSION_FLAG
    if is_largest(cy):
        compressed[0] |= SORT_FLAG
    return bytes(compressed), to_bytes(cx) + to_bytes(cy)


def vectors(points, to_bytes, is_largest, mul, on_curve):
    lines = []
    for P in points:
        assert on_curve(P), "vector is not on the curve"
        assert mul(r, P) is not None, "vector is in the prime-order subgroup"
        compressed, uncompressed = encodings(P, to_bytes, is_largest)
        lines.append("    (")
        for raw in (compressed, uncompressed):
            lines.append("        [")
            for i in range(0, len(raw), 16):
                lines.append("            " + ", ".join("0x%02x" % b for b in raw[i:i + 16]) + ",")
            lines.append("        ],")
        lines.append("    ),")
    size = len(to_bytes(points[0][0]))
    return ("/// Points on the curve but outside the prime-order subgroup, as\n"
            "/// `(compressed, uncompressed)` encodings; see `sage/bls12_381.sage`.\n"
            "const OFF_SUBGROUP: &[([u8; %d], [u8; %d])] = &[\n%s\n];" % (size, 2 * size, "\n".join(lines)))


def pick(points, neg, is_largest):
    # Two points from the smallest x-coordinates that decompress off the
    # subgroup, with the second one's y flipped if needed so that the two
    # disagree on the sort flag -- that way the vectors pin the flag both ways.
    first, second = points[0], points[1]
    if is_largest(second[1]) == is_largest(first[1]):
        second = neg(second)
    return [first, second]


def fp_neg(P):
    return (P[0], (-P[1]) % p)


def fp2_neg(P):
    return (P[0], f2_sub((0, 0), P[1]))


g1_off = pick(g1_points_off_subgroup(2), fp_neg, fp_is_largest)
g2_off = pick(g2_points_off_subgroup(2), fp2_neg, fp2_is_largest)
# and one of the shape [1]G + stray, which is what an attacker would send
g1_off.append(g1_add(G1, g1_off[0]))
g2_off.append(g2_add(G2, g2_off[0]))
assert fp_is_largest(g1_off[0][1]) != fp_is_largest(g1_off[1][1])
assert fp2_is_largest(g2_off[0][1]) != fp2_is_largest(g2_off[1][1])

# the endomorphism checks must reject every one of them, as the definition does
for P in g1_off:
    assert sigma(P, BETA) != g1_mul(lambda_minus_x2, P)
for Q in g2_off:
    assert psi(Q) != g2_mul(x % r, Q)

print("\n// src/curve/bls12_381/g1.rs, in `mod tests::subgroup`:\n")
print(vectors(g1_off, fp_bytes, fp_is_largest, g1_mul, g1_on_curve))
print("\n// src/curve/bls12_381/g2.rs, in `mod tests::subgroup`:\n")
print(rust_array("H_EFF", G2_H_EFF.to_bytes(80, "big"), "        ",
                 "h_eff for G2 = 3(x^2-1)*h2, as RFC 9380 section 8.8.2 gives it")
      .replace("pub const", "const"))
print()
print(vectors(g2_off, fp2_bytes, fp2_is_largest, g2_mul, g2_on_curve))
