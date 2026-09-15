macro_rules! test_kats_mul {
    ($curve: ident) => {
        #[test]
        fn $curve() {
            use super::kats_data::$curve::KATS;
            use crate::curve::field::Field;
            use crate::curve::sec2::$curve::{FieldElement, Point, PointAffine, Scalar};

            for kv in KATS.iter() {
                let x = FieldElement::from_bytes(&kv.x).expect("x fits");
                let y = FieldElement::from_bytes(&kv.y).expect("y fits");

                let k = if kv.n >= 0 {
                    Scalar::from_u64(kv.n as u64)
                } else {
                    Scalar::ZERO - Scalar::from_u64(kv.n.abs() as u64)
                };

                let expected_affine = PointAffine::from_coordinate(&x, &y).unwrap();
                let expected = Point::from_affine(&expected_affine);
                let got = &Point::GENERATOR * &k;
                let got_affine = got.to_affine().expect("KAT has affine");
                assert_eq!(expected_affine, got_affine, "KAT {}", kv.n);
                assert_eq!(expected, got, "KAT {}", kv.n);
            }
        }
    };
}

#[cfg(feature = "p112r1")]
test_kats_mul!(p112r1);
#[cfg(feature = "p128r1")]
test_kats_mul!(p128r1);
#[cfg(feature = "p160k1")]
test_kats_mul!(p160k1);
#[cfg(feature = "p160r1")]
test_kats_mul!(p160r1);
#[cfg(feature = "p160r2")]
test_kats_mul!(p160r2);

#[cfg(feature = "p192r1")]
test_kats_mul!(p192r1);
#[cfg(feature = "p192k1")]
test_kats_mul!(p192k1);
#[cfg(feature = "p224r1")]
test_kats_mul!(p224r1);
#[cfg(feature = "p224k1")]
test_kats_mul!(p224k1);
#[cfg(feature = "p256r1")]
test_kats_mul!(p256r1);
#[cfg(feature = "p256k1")]
test_kats_mul!(p256k1);
#[cfg(feature = "p384r1")]
test_kats_mul!(p384r1);
#[cfg(feature = "p521r1")]
test_kats_mul!(p521r1);
