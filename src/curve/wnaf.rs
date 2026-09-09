//! wNAF module

/// Compute the width-`w` non-adjacent form (wNAF) of the big-endian integer
/// `n`, returned as signed digits, least-significant first.
///
/// Every non-zero digit is odd and lies in `[-(2^(w-1)-1), 2^(w-1)-1]`, and no
/// two consecutive digits are both non-zero (so on average only `1/(w+1)` of
/// the digits are non-zero).
///
/// This is used solely by the *variable-time* scalar multiplications and
/// is deliberately not constant-time.
pub(crate) fn wnaf(n: &[u8], w: u32) -> Vec<i8> {
    debug_assert!((2..=8).contains(&w));

    // Little-endian working copy of `n` with one extra byte of headroom: the
    // `k -= digit` step transiently *increases* k when the digit is negative,
    // and the spare byte guarantees the carry never runs off the end.
    let mut k = Vec::with_capacity(n.len() + 1);
    k.extend(n.iter().rev().copied());
    k.push(0);

    let width = 1i32 << w; // 2^w
    let half = 1i32 << (w - 1); // 2^(w-1)

    let mut naf = Vec::with_capacity(n.len() * 8 + 1);
    while k.iter().any(|&b| b != 0) {
        let mut digit = 0i32;
        if k[0] & 1 == 1 {
            // k mod 2^w (w <= 8, so only the low byte contributes)
            let m = (k[0] as i32) & (width - 1);
            digit = if m >= half { m - width } else { m };

            // k -= digit (signed: adds |digit| when digit < 0). The arithmetic
            // right shift turns the per-byte overflow into a sign-extended
            // carry/borrow that propagates to the next byte.
            let mut carry = -digit;
            let mut i = 0;
            while carry != 0 && i < k.len() {
                let v = k[i] as i32 + carry;
                k[i] = (v & 0xff) as u8;
                carry = v >> 8;
                i += 1;
            }
        }
        naf.push(digit as i8);

        // k >>= 1
        let mut prev = 0u8;
        for b in k.iter_mut().rev() {
            let cur = *b;
            *b = (cur >> 1) | (prev << 7);
            prev = cur & 1;
        }
    }
    naf
}
