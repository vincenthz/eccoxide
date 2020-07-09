use super::kats_data::{KATS, KV};

macro_rules! test_kats_mul {
    ($curve: ident, $start: literal, $end: literal) => {
        #[test]
        fn $curve() {
            use crate::curve::sec2::$curve::{FieldElement, Point, PointAffine, Scalar};

            let kats: &[KV] = &KATS[$start..$end];
            for kv in kats.iter() {
                let mut xraw = [0u8; FieldElement::SIZE_BYTES];
                let mut yraw = [0u8; FieldElement::SIZE_BYTES];
                let mut kraw = [0u8; Scalar::SIZE_BYTES];

                xraw[FieldElement::SIZE_BYTES - kv.x.len()..].copy_from_slice(&kv.x);
                yraw[FieldElement::SIZE_BYTES - kv.y.len()..].copy_from_slice(&kv.y);
                kraw[Scalar::SIZE_BYTES - kv.k.len()..].copy_from_slice(&kv.k);

                let x = FieldElement::from_bytes(&xraw).unwrap();
                let y = FieldElement::from_bytes(&yraw).unwrap();
                let k = Scalar::from_bytes(&kraw).unwrap();
                let paffine = PointAffine::from_coordinate(&x, &y).unwrap();
                let expected = Point::from_affine(&paffine);
                let got = &Point::generator() * &k;
                assert_eq!(expected, got);
            }
        }
    };
}

test_kats_mul!(p192r1, 0, 52);
test_kats_mul!(p224r1, 52, 104);
test_kats_mul!(p256r1, 104, 156);
test_kats_mul!(p384r1, 156, 208);
test_kats_mul!(p521r1, 208, 260);
