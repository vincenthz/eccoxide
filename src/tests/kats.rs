macro_rules! test_kats_mul {
    ($curve: ident) => {
        #[test]
        fn $curve() {
            use super::kats_data::$curve::KATS;
            use crate::curve::sec2::$curve::{FieldElement, Point, PointAffine, Scalar};

            //let kats: &[KV] = &KATS[$start..$end];
            for kv in KATS.iter() {
                /*
                let mut xraw = [0u8; FieldElement::SIZE_BYTES];
                let mut yraw = [0u8; FieldElement::SIZE_BYTES];
                let mut kraw = [0u8; Scalar::SIZE_BYTES];

                xraw[FieldElement::SIZE_BYTES - kv.x.len()..].copy_from_slice(&kv.x);
                yraw[FieldElement::SIZE_BYTES - kv.y.len()..].copy_from_slice(&kv.y);
                kraw[Scalar::SIZE_BYTES - kv.k.len()..].copy_from_slice(&kv.k);

                let x = FieldElement::from_bytes(&xraw).expect("x fits");
                let y = FieldElement::from_bytes(&yraw).expect("y fits");
                let k = Scalar::from_bytes(&kraw).unwrap();
                */

                let x = FieldElement::from_bytes(&kv.x).expect("x fits");
                let y = FieldElement::from_bytes(&kv.y).expect("y fits");
                let k = Scalar::from_u64(kv.n);

                let paffine = PointAffine::from_coordinate(&x, &y).unwrap();
                let expected = Point::from_affine(&paffine);
                let got = &Point::GENERATOR * &k;
                assert_eq!(expected, got);
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
