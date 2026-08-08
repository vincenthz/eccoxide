//! Hashing to G1 and G2, following [RFC 9380].
//!
//! The two suites implemented here are the BLS12-381 ones of section 8.8:
//!
//! | suite | entry point |
//! |-------|-------------|
//! | `BLS12381G1_XMD:SHA-256_SSWU_RO_` | [`g1::Point::hash_to_curve`](super::g1::Point::hash_to_curve) |
//! | `BLS12381G1_XMD:SHA-256_SSWU_NU_` | [`g1::Point::encode_to_curve`](super::g1::Point::encode_to_curve) |
//! | `BLS12381G2_XMD:SHA-256_SSWU_RO_` | [`g2::Point::hash_to_curve`](super::g2::Point::hash_to_curve) |
//! | `BLS12381G2_XMD:SHA-256_SSWU_NU_` | [`g2::Point::encode_to_curve`](super::g2::Point::encode_to_curve) |
//!
//! `hash_to_curve` (the `RO_` suites) is indifferentiable from a random oracle
//! into the group, and is what protocols should use; `encode_to_curve` (`NU_`)
//! is a cheaper non-uniform encoding whose image is only a fraction of the
//! group, suitable when the input is already uniform.
//!
//! # Pipeline
//!
//! Both suites are built from the same four steps, in this module and its
//! callers:
//!
//! 1. [`expand_message_xmd`] stretches the message and the domain separation
//!    tag into uniform bytes with SHA-256 (section 5.3.1).
//! 2. `hash_to_field` reads those bytes as one or two field elements
//!    (section 5.2), 64 bytes per `Fp` component reduced modulo `p` — wide
//!    enough that the bias of the reduction is negligible.
//! 3. `map_to_curve` sends a field element to a point, by the Simplified SWU
//!    map onto an isogenous curve `E'` followed by the isogeny back to `E`
//!    (section 6.6.3). Neither BLS12-381 curve can be mapped onto directly:
//!    SSWU needs `A·B != 0`, and both have `A = 0`.
//! 4. `clear_cofactor` finishes the job, turning a point of the curve into a
//!    point of the prime-order group (see
//!    [`g1::Point::clear_cofactor`](super::g1::Point::clear_cofactor)).
//!
//! The domain separation tag is the caller's: RFC 9380 section 3.1 requires it
//! to be unique per application and recommends embedding the suite name. Tags
//! longer than 255 bytes are hashed down as section 5.3.3 prescribes, so any
//! length is accepted.
//!
//! [RFC 9380]: https://www.rfc-editor.org/rfc/rfc9380.html

use super::fp::Fp;
use super::fp2::Fp2;
use crate::curve::field::Field;
use crate::curve::projective;
use crate::mp::ct::{Choice, CtEqual, CtSelect, CtZero};
use crate::params::bls12_381::h2c;
use crate::params::bls12_381::P_MINUS3_DIV4_BYTES;
use cryptoxide::hashing::sha2::Sha256;

/// Output size of SHA-256, `b / 8` in the language of section 5.3.1.
const B_IN_BYTES: usize = 32;

/// Input block size of SHA-256, `s / 8`.
const S_IN_BYTES: usize = 64;

/// Bytes of uniform randomness per `Fp` element:
/// `L = ceil((ceil(log2(p)) + k) / 8) = ceil((381 + 128) / 8)`, for the
/// suites' security parameter `k = 128`.
const L: usize = 64;

/// The 17-byte prefix section 5.3.3 hashes an oversize tag under.
const OVERSIZE_DST_PREFIX: &[u8] = b"H2C-OVERSIZE-DST-";

/// Expand a message and a domain separation tag into uniform bytes with
/// SHA-256, as `expand_message_xmd` of section 5.3.1 of [RFC 9380].
///
/// The whole of `out` is filled; its length is the `len_in_bytes` of the
/// specification and must be at most `255 * 32 = 8160` bytes (the panic bound
/// below is the tighter of that and the specification's 65535).
///
/// A `dst` longer than the 255 bytes the format allows is replaced by
/// `SHA-256("H2C-OVERSIZE-DST-" || dst)`, which is what section 5.3.3 requires
/// of an implementation rather than an error.
///
/// [RFC 9380]: https://www.rfc-editor.org/rfc/rfc9380.html
pub fn expand_message_xmd(msg: &[u8], dst: &[u8], out: &mut [u8]) {
    // `ell`, the number of blocks, is what the specification bounds by 255;
    // the loop below counts them, so bound the length directly
    assert!(
        out.len() <= 255 * B_IN_BYTES,
        "expand_message_xmd: output too long"
    );

    // an oversize tag is hashed down to 32 bytes (5.3.3), and the tag is
    // always used with its length appended (`DST_prime`)
    let hashed;
    let dst = if dst.len() > 255 {
        hashed = Sha256::new()
            .update(OVERSIZE_DST_PREFIX)
            .update(dst)
            .finalize();
        &hashed[..]
    } else {
        dst
    };
    let dst_len = [dst.len() as u8];

    // b_0 = H(Z_pad || msg || I2OSP(len_in_bytes, 2) || I2OSP(0, 1) || DST_prime)
    let b_0 = Sha256::new()
        .update(&[0u8; S_IN_BYTES])
        .update(msg)
        .update(&(out.len() as u16).to_be_bytes())
        .update(&[0u8])
        .update(dst)
        .update(&dst_len)
        .finalize();

    // b_1 = H(b_0 || I2OSP(1, 1) || DST_prime), then
    // b_i = H((b_0 xor b_(i-1)) || I2OSP(i, 1) || DST_prime)
    let mut b_i = Sha256::new()
        .update(&b_0)
        .update(&[1u8])
        .update(dst)
        .update(&dst_len)
        .finalize();

    for (i, chunk) in out.chunks_mut(B_IN_BYTES).enumerate() {
        if i > 0 {
            let mut xored = [0u8; B_IN_BYTES];
            for (x, (a, b)) in xored.iter_mut().zip(b_0.iter().zip(b_i.iter())) {
                *x = a ^ b;
            }
            b_i = Sha256::new()
                .update(&xored)
                .update(&[(i + 1) as u8])
                .update(dst)
                .update(&dst_len)
                .finalize();
        }
        chunk.copy_from_slice(&b_i[..chunk.len()]);
    }
}

/// The field operations the maps below need, as plain methods.
///
/// `Fp` and `Fp2` both provide all of this already; going through a local trait
/// rather than the generic [`Field`] bounds keeps the mapping code written once
/// and free of the operator lifetime plumbing.
trait H2cField: Clone {
    const ZERO: Self;
    const ONE: Self;

    /// The Simplified SWU non-square, `Z` of section 8.8.
    ///
    /// Appendix F.2.1 requires [`Self::sqrt_ratio`] and the map to agree on it,
    /// which holding it here makes automatic.
    const Z: Self;

    fn add(&self, other: &Self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn neg(&self) -> Self;
    fn inverse(&self) -> Self;

    /// `(is_square(u/v), sqrt(u/v))`, or `(false, sqrt(Z·u/v))` when `u/v` is
    /// not a square — the `sqrt_ratio` of appendix F.2.1, `v != 0` assumed.
    ///
    /// The naive reading of that — an inversion, a quadratic character test and
    /// a square root — costs about twice what the shortcut of each field does,
    /// so the two implementations follow the specialized appendices instead.
    /// Both are free of branches on the field values.
    fn sqrt_ratio(u: &Self, v: &Self) -> (Choice, Self)
    where
        Self: Sized;

    fn is_zero_ct(&self) -> Choice;
    fn select(cond: Choice, a: &Self, b: &Self) -> Self;

    /// `sgn0` of section 4.1: the parity of the element, with the
    /// component-wise rule for the extension field.
    fn sgn0(&self) -> Choice;
}

impl H2cField for Fp {
    const ZERO: Self = <Fp as Field>::ZERO;
    const ONE: Self = <Fp as Field>::ONE;
    /// `Z = 11`, the non-square of the G1 suites
    const Z: Self = Fp::from_u64(11);

    fn add(&self, other: &Self) -> Self {
        self + other
    }
    fn mul(&self, other: &Self) -> Self {
        self * other
    }
    fn neg(&self) -> Self {
        -self
    }
    fn inverse(&self) -> Self {
        Fp::inverse(self)
    }

    /// `sqrt_ratio_3mod4` of appendix F.2.1.2, step for step, which `p = 3 mod
    /// 4` allows: one exponentiation instead of the inversion and two square
    /// roots of the definition.
    ///
    /// It works because `(u·v^3)^((p-3)/4)·u·v` squares to `u/v` when `u·v^3`
    /// is a square and to `-u/v` when it is not, the two cases being told apart
    /// by squaring the candidate back; the second is corrected by `sqrt(-Z)`.
    fn sqrt_ratio(u: &Self, v: &Self) -> (Choice, Self) {
        /// `c2 = sqrt(-Z)`
        const C2: Fp = Fp::from_bytes_unchecked(&h2c::g1::SQRT_RATIO_C2_BYTES);

        let tv1 = v.square(); // 1
        let tv2 = u * v; // 2
        let tv1 = &tv1 * &tv2; // 3
        let y1 = tv1.power(&P_MINUS3_DIV4_BYTES); // 4, c1 = (p-3)/4
        let y1 = &y1 * &tv2; // 5
        let y2 = &y1 * &C2; // 6
        let tv3 = &y1.square() * v; // 7, 8
        let is_qr = CtEqual::ct_eq(&tv3, u); // 9
        (is_qr, Fp::ct_select(is_qr, &y1, &y2)) // 10, 11
    }

    fn is_zero_ct(&self) -> Choice {
        self.ct_zero()
    }
    fn select(cond: Choice, a: &Self, b: &Self) -> Self {
        Fp::ct_select(cond, a, b)
    }
    fn sgn0(&self) -> Choice {
        // the parity of the canonical representative
        let bytes = self.to_bytes();
        Choice::from(bytes[Fp::SIZE_BYTES - 1] & 1 == 1)
    }
}

impl H2cField for Fp2 {
    const ZERO: Self = <Fp2 as Field>::ZERO;
    const ONE: Self = <Fp2 as Field>::ONE;
    /// `Z = -(2 + u)`, the non-square of the G2 suites
    const Z: Self = Fp2::from_bytes_unchecked(&h2c::g2::ISO_Z_BYTES);

    fn add(&self, other: &Self) -> Self {
        self + other
    }
    fn mul(&self, other: &Self) -> Self {
        self * other
    }
    fn neg(&self) -> Self {
        -self
    }
    fn inverse(&self) -> Self {
        Fp2::inverse(self)
    }

    /// The any-field `sqrt_ratio` of appendix F.2.1.1, step for step.
    ///
    /// `q = p^2` is `9 mod 16`, so neither the `3 mod 4` nor the `5 mod 8`
    /// shortcut applies and this is the constant-time Tonelli-Shanks form: one
    /// exponentiation by `c3` (about `q/16`) settles the square root up to a
    /// power of `c6`, and the loop of steps 17-26 walks that power off in
    /// `c1 = 3` fixed rounds. Still half the work of the definition, which
    /// spends two full square roots and an inversion.
    fn sqrt_ratio(u: &Self, v: &Self) -> (Choice, Self) {
        /// `c1`, the 2-adic valuation of `q - 1`
        const C1: usize = 3;
        /// `c4 = 2^c1 - 1`
        const C4: u64 = 7;
        /// `c5 = 2^(c1 - 1)`
        const C5: u64 = 4;
        /// `c6 = Z^((q - 1)/2^c1)`
        const C6: Fp2 = Fp2::from_bytes_unchecked(&h2c::g2::SQRT_RATIO_C6_BYTES);
        /// `c7 = Z^(((q - 1)/2^c1 + 1)/2)`
        const C7: Fp2 = Fp2::from_bytes_unchecked(&h2c::g2::SQRT_RATIO_C7_BYTES);

        let mut tv1 = C6; // 1
        let tv2 = v.power(&C4.to_be_bytes()); // 2
        let tv3 = &tv2.square() * v; // 3, 4
        let tv5 = u * &tv3; // 5
        let tv5 = tv5.power(&h2c::g2::SQRT_RATIO_C3_BYTES); // 6
        let tv5 = &tv5 * &tv2; // 7
        let tv2 = &tv5 * v; // 8
        let mut tv3 = &tv5 * u; // 9
        let mut tv4 = &tv3 * &tv2; // 10
        let tv5 = tv4.power(&C5.to_be_bytes()); // 11

        // step 12 reads `isQR = tv5 == 1`, which answers false for u = 0 even
        // though 0 is a square: every intermediate above vanishes there, so the
        // test has nothing to go on. Simplified SWU never asks (its numerator
        // would be a point of order 2, and neither curve has one), but the
        // answer should not depend on that.
        let is_qr = CtEqual::ct_eq(&tv5, &<Fp2 as Field>::ONE) | u.is_zero_ct(); // 12

        let tv2 = &tv3 * &C7; // 13
        let tv5 = &tv4 * &tv1; // 14
        tv3 = Fp2::ct_select(is_qr, &tv3, &tv2); // 15
        tv4 = Fp2::ct_select(is_qr, &tv4, &tv5); // 16

        // 17-26: c1 is a constant, so the trip count is fixed and so is the
        // number of squarings each round
        for i in (2..=C1).rev() {
            let mut tv5 = tv4.clone(); // 18, 19
            for _ in 0..i - 2 {
                tv5 = tv5.square(); // 20, tv4^(2^(i-2))
            }
            let e1 = CtEqual::ct_eq(&tv5, &<Fp2 as Field>::ONE); // 21
            let tv2 = &tv3 * &tv1; // 22
            tv1 = tv1.square(); // 23
            let tv5 = &tv4 * &tv1; // 24
            tv3 = Fp2::ct_select(e1, &tv3, &tv2); // 25
            tv4 = Fp2::ct_select(e1, &tv4, &tv5); // 26
        }
        (is_qr, tv3) // 27
    }

    fn is_zero_ct(&self) -> Choice {
        self.ct_zero()
    }
    fn select(cond: Choice, a: &Self, b: &Self) -> Self {
        Fp2::ct_select(cond, a, b)
    }
    fn sgn0(&self) -> Choice {
        // sign_0 OR (zero_0 AND sign_1), the m = 2 case of section 4.1
        self.c0.sgn0() | (self.c0.ct_zero() & self.c1.sgn0())
    }
}

/// The Simplified SWU map of appendix F.2, onto the curve
/// `y^2 = x^3 + a·x + b` — here always an isogenous `E'`, since the map needs
/// `a·b != 0` and both BLS12-381 curves have `a = 0`.
///
/// This is the straight-line form of the specification, step for step, so it is
/// free of branches on the field values.
fn map_to_curve_sswu<F: H2cField>(u: &F, a: &F, b: &F) -> (F, F) {
    let z = &F::Z;
    let tv1 = z.mul(&u.mul(u)); // 1, 2
    let tv2 = tv1.mul(&tv1).add(&tv1); // 3, 4
    let tv3 = b.mul(&tv2.add(&F::ONE)); // 5, 6

    // tv4 is a·Z when tv2 vanishes and a·(-tv2) otherwise, so it is never zero
    // and the division of step 25 is always defined
    let tv4 = a.mul(&F::select(tv2.is_zero_ct(), z, &tv2.neg())); // 7, 8
    let tv6 = tv4.mul(&tv4); // 10
    let tv2 = tv3.mul(&tv3).add(&a.mul(&tv6)).mul(&tv3); // 9, 11, 12, 13
    let tv6 = tv6.mul(&tv4); // 14
    let tv2 = tv2.add(&b.mul(&tv6)); // 15, 16

    let x = tv1.mul(&tv3); // 17
    let (is_gx1_square, y1) = F::sqrt_ratio(&tv2, &tv6); // 18
    let y = tv1.mul(u).mul(&y1); // 19, 20
    let x = F::select(is_gx1_square, &tv3, &x); // 21
    let y = F::select(is_gx1_square, &y1, &y); // 22

    // the sign of y follows the sign of u, which is what makes the map
    // injective on the sign
    let y = F::select(u.sgn0().ct_eq(&y.sgn0()), &y, &y.neg()); // 23, 24

    (x.mul(&tv4.inverse()), y) // 25
}

/// Evaluate `coeffs[0] + coeffs[1]·x + ...` by Horner's rule, with an implicit
/// leading coefficient of 1 when `monic` — the shape the isogeny denominators
/// are given in.
fn poly_eval<F: H2cField>(coeffs: &[F], x: &F, monic: bool) -> F {
    let (mut acc, rest) = if monic {
        (F::ONE, coeffs)
    } else {
        let (last, rest) = coeffs.split_last().expect("non-empty polynomial");
        (last.clone(), rest)
    };
    for c in rest.iter().rev() {
        acc = acc.mul(x).add(c);
    }
    acc
}

/// The rational maps of an isogeny `E' -> E`, as coefficients ascending in `x'`.
///
/// `x = x_num/x_den` and `y = y'·y_num/y_den`, with `x_den` and `y_den` monic
/// (their leading coefficient is not stored). This is the form of appendix E.
struct IsoMap<F: 'static> {
    x_num: &'static [F],
    x_den: &'static [F],
    y_num: &'static [F],
    y_den: &'static [F],
}

/// Push a point of `E'` through the isogeny, in projective coordinates.
///
/// Keeping the result projective avoids an inversion, and takes care of the
/// exceptional case of section 6.6.3 on its own: an input where either
/// denominator vanishes gives `Z = 0`, i.e. exactly the identity the
/// specification asks for.
fn iso_map<F: H2cField>(iso: &IsoMap<F>, x: &F, y: &F) -> projective::Point<F> {
    let x_num = poly_eval(iso.x_num, x, false);
    let x_den = poly_eval(iso.x_den, x, true);
    let y_num = poly_eval(iso.y_num, x, false);
    let y_den = poly_eval(iso.y_den, x, true);

    let z = x_den.mul(&y_den);
    let point = projective::Point {
        x: x_num.mul(&y_den),
        y: y.mul(&y_num).mul(&x_den),
        z: z.clone(),
    };
    // a vanishing denominator would leave (0 : 0 : 0), which is not a valid
    // representation of anything; hand back the canonical identity instead
    let infinity = projective::Point {
        x: F::ZERO,
        y: F::ONE,
        z: F::ZERO,
    };
    projective::Point {
        x: F::select(z.is_zero_ct(), &infinity.x, &point.x),
        y: F::select(z.is_zero_ct(), &infinity.y, &point.y),
        z: F::select(z.is_zero_ct(), &infinity.z, &point.z),
    }
}

/// Build a table of field elements from the byte constants of the parameters.
macro_rules! fe_table {
    ($fe:ty, $bytes:expr; $($i:literal)+) => {
        [$(<$fe>::from_bytes_unchecked(&$bytes[$i])),+]
    };
}

// --- G1: the 11-isogeny of appendix E.2 -------------------------------------

const G1_ISO_A: Fp = Fp::from_bytes_unchecked(&h2c::g1::ISO_A_BYTES);
const G1_ISO_B: Fp = Fp::from_bytes_unchecked(&h2c::g1::ISO_B_BYTES);

const G1_X_NUM: [Fp; 12] = fe_table!(Fp, h2c::g1::ISO_K1_BYTES; 0 1 2 3 4 5 6 7 8 9 10 11);
const G1_X_DEN: [Fp; 10] = fe_table!(Fp, h2c::g1::ISO_K2_BYTES; 0 1 2 3 4 5 6 7 8 9);
const G1_Y_NUM: [Fp; 16] =
    fe_table!(Fp, h2c::g1::ISO_K3_BYTES; 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15);
const G1_Y_DEN: [Fp; 15] = fe_table!(Fp, h2c::g1::ISO_K4_BYTES; 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14);

const G1_ISO: IsoMap<Fp> = IsoMap {
    x_num: &G1_X_NUM,
    x_den: &G1_X_DEN,
    y_num: &G1_Y_NUM,
    y_den: &G1_Y_DEN,
};

// --- G2: the 3-isogeny of appendix E.3 --------------------------------------

const G2_ISO_A: Fp2 = Fp2::from_bytes_unchecked(&h2c::g2::ISO_A_BYTES);
const G2_ISO_B: Fp2 = Fp2::from_bytes_unchecked(&h2c::g2::ISO_B_BYTES);

const G2_X_NUM: [Fp2; 4] = fe_table!(Fp2, h2c::g2::ISO_K1_BYTES; 0 1 2 3);
const G2_X_DEN: [Fp2; 2] = fe_table!(Fp2, h2c::g2::ISO_K2_BYTES; 0 1);
const G2_Y_NUM: [Fp2; 4] = fe_table!(Fp2, h2c::g2::ISO_K3_BYTES; 0 1 2 3);
const G2_Y_DEN: [Fp2; 3] = fe_table!(Fp2, h2c::g2::ISO_K4_BYTES; 0 1 2);

const G2_ISO: IsoMap<Fp2> = IsoMap {
    x_num: &G2_X_NUM,
    x_den: &G2_X_DEN,
    y_num: &G2_Y_NUM,
    y_den: &G2_Y_DEN,
};

/// One `Fp` element from `L` bytes of uniform randomness, reduced modulo `p`.
fn fp_from_uniform(bytes: &[u8]) -> Fp {
    // the reduction takes a double-width buffer, so pad on the left
    let mut wide = [0u8; 2 * Fp::SIZE_BYTES];
    wide[2 * Fp::SIZE_BYTES - L..].copy_from_slice(bytes);
    Fp::init_from_wide_bytes_be(wide)
}

/// `hash_to_field(msg, count)` of section 5.2, for `count` in `1..=2`.
///
/// Each element takes `m` blocks of `L` bytes, one per component of the field.
fn hash_to_field<const COUNT: usize, const M: usize, F, G>(
    msg: &[u8],
    dst: &[u8],
    build: G,
) -> [F; COUNT]
where
    G: Fn(&[u8]) -> F,
{
    assert!(COUNT * M <= 4);
    let mut uniform = [0u8; 4 * L];
    let len = COUNT * M * L;
    expand_message_xmd(msg, dst, &mut uniform[..len]);

    core::array::from_fn(|i| build(&uniform[i * M * L..(i + 1) * M * L]))
}

/// `count` elements of `Fp`.
fn hash_to_field_fp<const COUNT: usize>(msg: &[u8], dst: &[u8]) -> [Fp; COUNT] {
    hash_to_field::<COUNT, 1, _, _>(msg, dst, fp_from_uniform)
}

/// `count` elements of `Fp2`, each `e_0 + e_1·u` from two blocks in that order.
fn hash_to_field_fp2<const COUNT: usize>(msg: &[u8], dst: &[u8]) -> [Fp2; COUNT] {
    hash_to_field::<COUNT, 2, _, _>(msg, dst, |raw| {
        Fp2::new(fp_from_uniform(&raw[..L]), fp_from_uniform(&raw[L..]))
    })
}

/// `map_to_curve` for G1: Simplified SWU onto `E'`, then the 11-isogeny to `E`.
///
/// The result is a point of the curve, *not* of G1: the caller still owes it a
/// cofactor clearing.
pub(super) fn map_to_curve_g1(u: &Fp) -> projective::Point<Fp> {
    let (x, y) = map_to_curve_sswu(u, &G1_ISO_A, &G1_ISO_B);
    iso_map(&G1_ISO, &x, &y)
}

/// `map_to_curve` for G2: Simplified SWU onto `E'`, then the 3-isogeny to `E'`
/// of the twist.
///
/// The result is a point of the twist, *not* of G2; see [`map_to_curve_g1`].
pub(super) fn map_to_curve_g2(u: &Fp2) -> projective::Point<Fp2> {
    let (x, y) = map_to_curve_sswu(u, &G2_ISO_A, &G2_ISO_B);
    iso_map(&G2_ISO, &x, &y)
}

/// The two field elements the `RO_` suites hash a message to, for G1.
pub(super) fn hash_to_field_g1(msg: &[u8], dst: &[u8]) -> [Fp; 2] {
    hash_to_field_fp(msg, dst)
}

/// The single field element the `NU_` suites hash a message to, for G1.
pub(super) fn encode_to_field_g1(msg: &[u8], dst: &[u8]) -> Fp {
    let [u] = hash_to_field_fp::<1>(msg, dst);
    u
}

/// The two field elements the `RO_` suites hash a message to, for G2.
pub(super) fn hash_to_field_g2(msg: &[u8], dst: &[u8]) -> [Fp2; 2] {
    hash_to_field_fp2(msg, dst)
}

/// The single field element the `NU_` suites hash a message to, for G2.
pub(super) fn encode_to_field_g2(msg: &[u8], dst: &[u8]) -> Fp2 {
    let [u] = hash_to_field_fp2::<1>(msg, dst);
    u
}

#[cfg(test)]
#[path = "hash_to_curve_vectors.rs"]
mod vectors;

#[cfg(test)]
mod tests {
    use super::vectors::{
        ExpandVector, Vector, G1_NU, G1_NU_DST, G1_RO, G1_RO_DST, G2_NU, G2_NU_DST, G2_RO,
        G2_RO_DST, XMD, XMD_DST, XMD_LONG, XMD_LONG_DST,
    };
    use super::*;
    use crate::curve::bls12_381::{g1, g2};
    use core::convert::TryInto;

    fn hex(s: &str) -> Vec<u8> {
        assert_eq!(s.len() % 2, 0, "hex needs an even length");
        (0..s.len() / 2)
            .map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).expect("hex"))
            .collect()
    }

    fn fp(s: &str) -> Fp {
        Fp::from_bytes(&hex(s).try_into().expect("48 bytes")).expect("canonical Fp")
    }

    fn fp2(s: &str) -> Fp2 {
        Fp2::from_bytes(&hex(s).try_into().expect("96 bytes")).expect("canonical Fp2")
    }

    /// `expand_message_xmd` against appendix K.1, over every requested output
    /// length in the vectors (32 and 128 bytes, i.e. one and four blocks).
    #[test]
    fn expand_message_xmd_kat() {
        for ExpandVector { msg, uniform } in XMD {
            let expected = hex(uniform);
            let mut out = vec![0u8; expected.len()];
            expand_message_xmd(msg.as_bytes(), XMD_DST.as_bytes(), &mut out);
            assert_eq!(out, expected, "expand_message_xmd(msg = {:?})", msg);
        }
    }

    /// Appendix K.2: a tag longer than 255 bytes is hashed down first, which
    /// the vectors pin — the same messages give different output than K.1.
    #[test]
    fn expand_message_xmd_oversize_dst_kat() {
        assert!(
            XMD_LONG_DST.len() > 255,
            "the long tag must exceed the limit"
        );
        for ExpandVector { msg, uniform } in XMD_LONG {
            let expected = hex(uniform);
            let mut out = vec![0u8; expected.len()];
            expand_message_xmd(msg.as_bytes(), XMD_LONG_DST.as_bytes(), &mut out);
            assert_eq!(out, expected, "oversize dst, msg = {:?}", msg);
        }
    }

    /// The whole G1 pipeline against appendix J.9: the field elements the
    /// message expands to, the points they map to before cofactor clearing,
    /// and the final group element.
    #[test]
    fn g1_kat() {
        for (dst, vectors, ro) in [(G1_RO_DST, G1_RO, true), (G1_NU_DST, G1_NU, false)] {
            for Vector { msg, u, q, p } in vectors {
                let (msg, dst) = (msg.as_bytes(), dst.as_bytes());
                let expected_u: Vec<Fp> = u.iter().map(|s| fp(s)).collect();
                let got_u: Vec<Fp> = if ro {
                    hash_to_field_g1(msg, dst).to_vec()
                } else {
                    vec![encode_to_field_g1(msg, dst)]
                };
                assert_eq!(got_u, expected_u, "hash_to_field, msg = {:?}", msg);

                for (u, (qx, qy)) in got_u.iter().zip(q.iter()) {
                    let mapped = g1::Point::from_projective(map_to_curve_g1(u));
                    let affine = mapped
                        .to_affine()
                        .expect("map_to_curve is not the identity");
                    let (x, y) = affine.to_coordinate();
                    assert_eq!((x, y), (&fp(qx), &fp(qy)), "map_to_curve, msg = {:?}", msg);
                }

                let point = if ro {
                    g1::Point::hash_to_curve(msg, dst)
                } else {
                    g1::Point::encode_to_curve(msg, dst)
                };
                let affine = point.to_affine().expect("the vectors are not the identity");
                let (x, y) = affine.to_coordinate();
                assert_eq!(
                    (x, y),
                    (&fp(p.0), &fp(p.1)),
                    "hash_to_curve, msg = {:?}",
                    msg
                );
                assert!(
                    point.is_in_subgroup().is_true(),
                    "not in G1, msg = {:?}",
                    msg
                );
            }
        }
    }

    /// The same against appendix J.10, for G2.
    #[test]
    fn g2_kat() {
        for (dst, vectors, ro) in [(G2_RO_DST, G2_RO, true), (G2_NU_DST, G2_NU, false)] {
            for Vector { msg, u, q, p } in vectors {
                let (msg, dst) = (msg.as_bytes(), dst.as_bytes());
                let expected_u: Vec<Fp2> = u.iter().map(|s| fp2(s)).collect();
                let got_u: Vec<Fp2> = if ro {
                    hash_to_field_g2(msg, dst).to_vec()
                } else {
                    vec![encode_to_field_g2(msg, dst)]
                };
                assert_eq!(got_u, expected_u, "hash_to_field, msg = {:?}", msg);

                for (u, (qx, qy)) in got_u.iter().zip(q.iter()) {
                    let mapped = g2::Point::from_projective(map_to_curve_g2(u));
                    let affine = mapped
                        .to_affine()
                        .expect("map_to_curve is not the identity");
                    let (x, y) = affine.to_coordinate();
                    assert_eq!(
                        (x, y),
                        (&fp2(qx), &fp2(qy)),
                        "map_to_curve, msg = {:?}",
                        msg
                    );
                }

                let point = if ro {
                    g2::Point::hash_to_curve(msg, dst)
                } else {
                    g2::Point::encode_to_curve(msg, dst)
                };
                let affine = point.to_affine().expect("the vectors are not the identity");
                let (x, y) = affine.to_coordinate();
                assert_eq!(
                    (x, y),
                    (&fp2(p.0), &fp2(p.1)),
                    "hash_to_curve, msg = {:?}",
                    msg
                );
                assert!(
                    point.is_in_subgroup().is_true(),
                    "not in G2, msg = {:?}",
                    msg
                );
            }
        }
    }

    /// Hashing is deterministic, sensitive to both inputs, and lands in the
    /// prime-order group for messages beyond the vectors.
    #[test]
    fn hashing_is_deterministic_and_separated() {
        const DST: &[u8] = b"eccoxide-test-BLS12381G1_XMD:SHA-256_SSWU_RO_";
        let a = g1::Point::hash_to_curve(b"message", DST);
        assert_eq!(a, g1::Point::hash_to_curve(b"message", DST));
        assert_ne!(a, g1::Point::hash_to_curve(b"messagf", DST));
        assert_ne!(a, g1::Point::hash_to_curve(b"message", b"other-dst"));
        assert_ne!(a, g1::Point::encode_to_curve(b"message", DST));
        assert!(a.is_in_subgroup().is_true());

        let b = g2::Point::hash_to_curve(b"message", DST);
        assert_eq!(b, g2::Point::hash_to_curve(b"message", DST));
        assert_ne!(b, g2::Point::hash_to_curve(b"messagf", DST));
        assert!(b.is_in_subgroup().is_true());
    }

    /// The optimized `sqrt_ratio` of each field answers the definition it
    /// replaces: the verdict is the quadratic character of `u/v`, and the root
    /// squares back to `u/v` when that holds and to `Z·u/v` when it does not.
    ///
    /// Worth pinning on its own — the vectors above only see `sqrt_ratio`
    /// through the map, which hides the verdict behind a select and normalizes
    /// the sign of the root away.
    #[test]
    fn sqrt_ratio_matches_the_definition() {
        let (mut fp_squares, mut fp2_squares) = (0, 0);
        for i in 1u64..=40 {
            let (u, v) = (Fp::from_u64(i), Fp::from_u64(1000 + i));
            let (is_qr, y) = <Fp as H2cField>::sqrt_ratio(&u, &v);
            let ratio = &u * &v.inverse();
            assert_eq!(
                is_qr.is_true(),
                ratio.sqrt().is_present().is_true(),
                "Fp: wrong verdict for {}/{}",
                i,
                1000 + i
            );
            let want = if is_qr.is_true() {
                ratio
            } else {
                &ratio * &<Fp as H2cField>::Z
            };
            assert_eq!(y.square(), want, "Fp: wrong root for {}/{}", i, 1000 + i);
            fp_squares += is_qr.is_true() as u32;

            let u = Fp2::new(Fp::from_u64(i), Fp::from_u64(2 * i + 1));
            let v = Fp2::new(Fp::from_u64(3 * i + 1), Fp::from_u64(i + 7));
            let (is_qr, y) = <Fp2 as H2cField>::sqrt_ratio(&u, &v);
            let ratio = &u * &v.inverse();
            assert_eq!(
                is_qr.is_true(),
                ratio.sqrt().is_present().is_true(),
                "Fp2: wrong verdict at i = {}",
                i
            );
            let want = if is_qr.is_true() {
                ratio
            } else {
                &ratio * &<Fp2 as H2cField>::Z
            };
            assert_eq!(y.square(), want, "Fp2: wrong root at i = {}", i);
            fp2_squares += is_qr.is_true() as u32;
        }
        // both branches have to have been taken, or the loop above proves
        // nothing about the one it skipped
        assert!(0 < fp_squares && fp_squares < 40, "Fp: one-sided sample");
        assert!(0 < fp2_squares && fp2_squares < 40, "Fp2: one-sided sample");

        // zero is a square, which the any-field variant needs told
        for v in [Fp::from_u64(7), Fp::one()] {
            let (is_qr, y) = <Fp as H2cField>::sqrt_ratio(&Fp::zero(), &v);
            assert!(is_qr.is_true() && y == Fp::zero());
        }
        for v in [
            Fp2::new(Fp::from_u64(7), Fp::from_u64(3)),
            <Fp2 as Field>::ONE,
        ] {
            let (is_qr, y) = <Fp2 as H2cField>::sqrt_ratio(&<Fp2 as Field>::ZERO, &v);
            assert!(is_qr.is_true() && y == <Fp2 as Field>::ZERO);
        }
    }

    /// Every message lands in the group, including ones whose field elements
    /// are extreme (zero maps to the SSWU exceptional branch).
    #[test]
    fn map_to_curve_handles_edge_inputs() {
        for u in [Fp::zero(), Fp::one(), -Fp::one(), Fp::from_u64(11)] {
            let p = g1::Point::from_projective(map_to_curve_g1(&u));
            let affine = p.to_affine().expect("the map does not give the identity");
            let (x, y) = affine.to_coordinate();
            assert!(
                g1::PointAffine::from_coordinate(x, y).is_some(),
                "map_to_curve left the curve"
            );
            assert!(p.clear_cofactor().is_in_subgroup().is_true());
        }
        for u in [
            <Fp2 as Field>::ZERO,
            <Fp2 as Field>::ONE,
            -<Fp2 as Field>::ONE,
        ] {
            let p = g2::Point::from_projective(map_to_curve_g2(&u));
            let affine = p.to_affine().expect("the map does not give the identity");
            let (x, y) = affine.to_coordinate();
            assert!(
                g2::PointAffine::from_coordinate(x, y).is_some(),
                "map_to_curve left the twist"
            );
            assert!(p.clear_cofactor().is_in_subgroup().is_true());
        }
    }
}
