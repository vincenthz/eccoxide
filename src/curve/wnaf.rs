//! wNAF module

mod sealed {
    pub trait Sealed {}
    impl<const N: usize> Sealed for [u8; N] {}
}

/// Big-endian integer accepted by the variable-time scalar multiplications:
/// a byte array `[u8; N]`, whose length sizes the wNAF digit buffer on the
/// stack.
///
/// This trait is sealed and only implemented for `[u8; N]`.
pub trait ScalarBytes: sealed::Sealed {
    /// The width-`w` non-adjacent form (wNAF) of the integer, `2 <= w <= 8`,
    /// as signed digits, most significant first, starting from the most
    /// significant non-zero digit.
    ///
    /// Every non-zero digit is odd and lies in `[-(2^(w-1)-1), 2^(w-1)-1]`,
    /// and two non-zero digits are at least `w` positions apart. A `B`-bit
    /// integer has at most `B + 1` digits.
    ///
    /// Not constant-time.
    fn wnaf(&self, w: u32) -> impl Iterator<Item = i8>;
}

impl<const N: usize> ScalarBytes for [u8; N] {
    fn wnaf(&self, w: u32) -> impl Iterator<Item = i8> {
        let naf = Wnaf::new(self, w);
        (0..naf.len).rev().map(move |i| naf.digit(i))
    }
}

/// Width-`w` non-adjacent form (wNAF) of an `N`-byte big-endian integer.
///
/// Every non-zero digit is odd and lies in `[-(2^(w-1)-1), 2^(w-1)-1]`, and two
/// non-zero digits are at least `w` positions apart (so on average only
/// `1/(w+1)` of the digits are non-zero).
///
/// A `B`-bit integer has at most `B + 1` digits: the `8N` digits of the bits,
/// plus one past the top bit for a final carry.
///
/// This is used solely by the *variable-time* scalar multiplications and
/// is deliberately not constant-time.
pub(crate) struct Wnaf<const N: usize> {
    /// `low[i / 8][i % 8]` is the digit of weight `2^i`, for `i < 8N`
    low: [[i8; 8]; N],
    /// the digit of weight `2^(8N)`
    top: i8,
    /// number of significant digits
    len: usize,
}

impl<const N: usize> Wnaf<N> {
    pub(crate) fn new(n: &[u8; N], w: u32) -> Self {
        debug_assert!((2..=8).contains(&w));

        let width = 1u32 << w; // 2^w
        let half = 1u32 << (w - 1); // 2^(w-1)

        let mut naf = Wnaf {
            low: [[0; 8]; N],
            top: 0,
            len: 0,
        };

        // Scan the bits from the least significant, with `carry` the pending
        // +1 at `pos` left by a negative digit. An even window gives a zero
        // digit and moves one bit up; an odd one gives the digit `window mods
        // 2^w`, which clears the next w-1 positions, so the scan jumps w bits.
        let mut carry = 0u32;
        let mut pos = 0;
        while pos <= 8 * N {
            let window = (bits_at(n, pos) & (width - 1)) + carry;
            if window & 1 == 0 {
                pos += 1;
                continue;
            }
            let digit = if window < half {
                carry = 0;
                window as i32
            } else {
                carry = 1;
                window as i32 - width as i32
            };
            naf.set(pos, digit as i8);
            naf.len = pos + 1;
            pos += w as usize;
        }
        // a B-bit integer has at most B+1 digits, so the carry is always
        // absorbed by position 8N
        debug_assert_eq!(carry, 0);
        naf
    }

    /// The digit of weight `2^i`, for `i <= 8N`.
    #[inline]
    fn digit(&self, i: usize) -> i8 {
        if i < 8 * N {
            self.low[i / 8][i % 8]
        } else {
            self.top
        }
    }

    #[inline]
    fn set(&mut self, i: usize, d: i8) {
        if i < 8 * N {
            self.low[i / 8][i % 8] = d
        } else {
            self.top = d
        }
    }
}

/// Bits `pos..pos + 9` of the big-endian integer `n` (at least, as the low bits
/// of the result), reading zero above its top bit.
#[inline]
fn bits_at<const N: usize>(n: &[u8; N], pos: usize) -> u32 {
    let byte = |i: usize| if i < N { n[N - 1 - i] as u32 } else { 0 };
    let b = pos / 8;
    (byte(b) | (byte(b + 1) << 8)) >> (pos % 8)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Evaluate `Σ digits[i]·2^i` back into a big-endian integer of `N`
    /// bytes, by Horner's rule from the most significant digit.
    ///
    /// The sum is taken modulo `2^(8(N+2))`, which its magnitude stays well
    /// under, so two zero high bytes at the end mean it really fits `N` bytes.
    fn eval<const N: usize>(naf: &Wnaf<N>) -> [u8; N] {
        let mut acc = vec![0u8; N + 2]; // little-endian
        for i in (0..naf.len).rev() {
            // acc = 2 * acc + d
            let mut carry = naf.digit(i) as i32;
            for b in acc.iter_mut() {
                let v = ((*b as i32) << 1) + carry;
                *b = (v & 0xff) as u8;
                carry = v >> 8;
            }
        }
        assert_eq!(acc[N..], [0, 0], "value does not fit {} bytes", N);
        let mut out = [0u8; N];
        for (o, a) in out.iter_mut().rev().zip(acc.iter()) {
            *o = *a;
        }
        out
    }

    fn check<const N: usize>(n: &[u8; N], w: u32) {
        let naf = Wnaf::new(n, w);

        assert_eq!(&eval(&naf), n, "w={} n={:x?}", w, n);
        assert!(naf.len <= 8 * N + 1);
        if naf.len > 0 {
            assert_ne!(naf.digit(naf.len - 1), 0, "most significant digit is zero");
        }
        let bound = (1i32 << (w - 1)) - 1;
        for i in 0..naf.len {
            let d = naf.digit(i);
            if d != 0 {
                assert_eq!(d & 1, 1, "even digit {} at {}", d, i);
                assert!((d as i32).abs() <= bound, "digit {} out of range", d);
                // the next w-1 digits are zero
                for j in i + 1..(i + w as usize).min(naf.len) {
                    assert_eq!(naf.digit(j), 0, "non-zero digits closer than {} apart", w);
                }
            }
        }

        // the iterator yields the same digits, most significant first
        let mut it = n.wnaf(w);
        for i in (0..naf.len).rev() {
            assert_eq!(it.next(), Some(naf.digit(i)));
        }
        assert_eq!(it.next(), None);
    }

    #[test]
    fn small_values() {
        for w in 2..=8 {
            for v in 0u16..=0x3ff {
                check(&v.to_be_bytes(), w);
            }
            for v in 0u8..=0xff {
                check(&[v], w);
            }
        }
    }

    #[test]
    fn zero_is_empty() {
        assert_eq!([0u8; 0].wnaf(5).count(), 0);
        assert_eq!([0u8; 300].wnaf(5).count(), 0);
    }

    #[test]
    fn wide_values() {
        // all ones (maximal carry chains), a top bit alone, and a byte pattern,
        // each also zero-padded
        let mut x = [0u8; 80];
        for (i, b) in x.iter_mut().enumerate() {
            *b = (i as u8).wrapping_mul(37).wrapping_add(11);
        }
        let mut top = [0u8; 80];
        top[0] = 0x80;
        for n in [[0xffu8; 80], top, x] {
            for w in 2..=8 {
                check(&n, w);
                let mut padded = [0u8; 87];
                padded[7..].copy_from_slice(&n);
                check(&padded, w);
                assert!(n.wnaf(w).eq(padded.wnaf(w)));
            }
        }
        check(&[0xff; 66], 5);
        check(&[0xff; 32], 5);
    }
}
