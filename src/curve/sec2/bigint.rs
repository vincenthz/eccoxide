//! macros to generate curve definition using num_trait/num_bigint dependencies

use crate::bigint_prime_curve;

#[cfg(feature = "p112r1")]
bigint_prime_curve!(p112r1, 112);
#[cfg(feature = "p112r2")]
bigint_prime_curve!(p112r2, 112);
#[cfg(feature = "p128r1")]
bigint_prime_curve!(p128r1, 128);
#[cfg(feature = "p128r2")]
bigint_prime_curve!(p128r2, 128);
#[cfg(feature = "p160k1")]
bigint_prime_curve!(p160k1, 160);
#[cfg(feature = "p160r1")]
bigint_prime_curve!(p160r1, 160);
#[cfg(feature = "p160r2")]
bigint_prime_curve!(p160r2, 160);
#[cfg(feature = "p192k1")]
bigint_prime_curve!(p192k1, 192);
#[cfg(feature = "p192r1")]
bigint_prime_curve!(p192r1, 192);
#[cfg(feature = "p224k1")]
bigint_prime_curve!(p224k1, 224);
#[cfg(feature = "p224r1")]
bigint_prime_curve!(p224r1, 224);
//#[cfg(feature = "p256k1")]
//bigint_prime_curve!(p256k1, 256);
//#[cfg(feature = "p256r1")]
//bigint_prime_curve!(p256r1, 256);
//#[cfg(feature = "p384r1")]
//bigint_prime_curve!(p384r1, 384);
//#[cfg(feature = "p521r1")]
//bigint_prime_curve!(p521r1, 521);

#[cfg(test)]
mod tests {
    /*
    #[cfg(feature = "p256r1")]
    mod p256r1 {
        use super::super::p256r1::*;

        #[test]
        fn point_add() {
            let k2_x = b"\x7C\xF2\x7B\x18\x8D\x03\x4F\x7E\x8A\x52\x38\x03\x04\xB5\x1A\xC3\xC0\x89\x69\xE2\x77\xF2\x1B\x35\xA6\x0B\x48\xFC\x47\x66\x99\x78";
            let k2_y = b"\x07\x77\x55\x10\xDB\x8E\xD0\x40\x29\x3D\x9A\xC6\x9F\x74\x30\xDB\xBA\x7D\xAD\xE6\x3C\xE9\x82\x29\x9E\x04\xB7\x9D\x22\x78\x73\xD1";
            let x = FieldElement::from_bytes(k2_x).unwrap();
            let y = FieldElement::from_bytes(k2_y).unwrap();

            let p_expected = PointAffine::from_coordinate(&x, &y).unwrap();
            let p_got = &Point::generator().double();
            let p_got_affine = p_got.to_affine().unwrap();

            assert_eq!(p_expected, p_got_affine);
        }
    }
    */

    #[cfg(feature = "p192r1")]
    mod p192r1 {
        use super::super::p192r1::*;

        #[test]
        fn point_add() {
            const KAT0: (&[u8], &[u8]) = (
                &[
                    0x18, 0x8d, 0xa8, 0x0e, 0xb0, 0x30, 0x90, 0xf6, 0x7c, 0xbf, 0x20, 0xeb, 0x43,
                    0xa1, 0x88, 0x00, 0xf4, 0xff, 0x0a, 0xfd, 0x82, 0xff, 0x10, 0x12,
                ],
                &[
                    0x07, 0x19, 0x2b, 0x95, 0xff, 0xc8, 0xda, 0x78, 0x63, 0x10, 0x11, 0xed, 0x6b,
                    0x24, 0xcd, 0xd5, 0x73, 0xf9, 0x77, 0xa1, 0x1e, 0x79, 0x48, 0x11,
                ],
            );

            const KAT1: (&[u8], &[u8]) = (
                &[
                    0xda, 0xfe, 0xbf, 0x58, 0x28, 0x78, 0x3f, 0x2a, 0xd3, 0x55, 0x34, 0x63, 0x15,
                    0x88, 0xa3, 0xf6, 0x29, 0xa7, 0x0f, 0xb1, 0x69, 0x82, 0xa8, 0x88,
                ],
                &[
                    0xdd, 0x6b, 0xda, 0x0d, 0x99, 0x3d, 0xa0, 0xfa, 0x46, 0xb2, 0x7b, 0xbc, 0x14,
                    0x1b, 0x86, 0x8f, 0x59, 0x33, 0x1a, 0xfa, 0x5c, 0x7e, 0x93, 0xab,
                ],
            );
            const KAT2: (&[u8], &[u8]) = (
                &[
                    0x76, 0xe3, 0x2a, 0x25, 0x57, 0x59, 0x9e, 0x6e, 0xdc, 0xd2, 0x83, 0x20, 0x1f,
                    0xb2, 0xb9, 0xaa, 0xdf, 0xd0, 0xd3, 0x59, 0xcb, 0xb2, 0x63, 0xda,
                ],
                &[
                    0x78, 0x2c, 0x37, 0xe3, 0x72, 0xba, 0x45, 0x20, 0xaa, 0x62, 0xe0, 0xfe, 0xd1,
                    0x21, 0xd4, 0x9e, 0xf3, 0xb5, 0x43, 0x66, 0x0c, 0xfd, 0x05, 0xfd,
                ],
            );

            let kats = [KAT0, KAT1, KAT2];

            let mut points = Vec::new();
            for k in &kats {
                let x = FieldElement::from_slice(k.0).unwrap();
                let y = FieldElement::from_slice(k.1).unwrap();

                let p_expected = PointAffine::from_coordinate(&x, &y)
                    .expect("cannot convert from affine: test vector value failed");
                points.push(p_expected);
            }

            let p = PointAffine::generator();
            assert_eq!(p, points[0]);

            let p2 = p.double();
            assert_eq!(p2, points[1]);

            let p3 = &p + &p2;
            assert_eq!(p3, points[2]);

            let p3x = Point::from_affine(&p) + Point::from_affine(&p2);
            assert_eq!(p3x, Point::from_affine(&points[2]));

            //let p = Point::generator() + Point::generator();
            //assert_eq!(p.to_affine().unwrap(), points[1]);
        }
    }
}
