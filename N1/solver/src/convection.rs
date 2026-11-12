//! Film-coefficient models for the exposed faces (D6 correction).
//!
//! The published type-2 solver silently replaced the input `convective_heat_transfer_film_coefficient` by a
//! height-dependent expression h(z) = 3 + 12 z^2 (z = element centroid height in metres), so the input had no effect.
//! * `Constant(h)` (model 0, **default**): the input film coefficient is applied to every exposed face.
//! * `HeightDependent` (model 1): the legacy h(z) = 3 + 12 z^2, kept as an explicit, documented option.

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Convection {
    Constant(f64),
    HeightDependent,
}

impl Convection {
    /// `model`: optional input key `convection_model` (absent => 0); `h_input`: `convective_heat_transfer_film_coefficient`.
    pub fn from_input(model: Option<usize>, h_input: f64) -> Convection {
        match model.unwrap_or(0) {
            0 => {
                assert!(h_input.is_finite() && h_input >= 0.0, "convective_heat_transfer_film_coefficient must be >= 0");
                Convection::Constant(h_input)
            }
            1 => Convection::HeightDependent,
            other => panic!("unknown convection_model {} (use 0 = constant input h, 1 = legacy h(z)=3+12z^2)", other),
        }
    }
    /// film coefficient [W/m2K] for an element whose centroid is at height `z` [m]
    pub fn h(&self, z: f64) -> f64 {
        match *self {
            Convection::Constant(h) => h,
            Convection::HeightDependent => 3.0 + 12.0 * z * z,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn constant_uses_the_input_everywhere() {
        let c = Convection::from_input(None, 37.5);
        assert_eq!(c.h(0.0), 37.5); assert_eq!(c.h(0.3), 37.5);
    }
    #[test]
    fn legacy_matches_the_published_expression() {
        let c = Convection::from_input(Some(1), 999.0);
        assert!((c.h(0.1) - 3.12).abs() < 1e-12); assert!((c.h(0.0) - 3.0).abs() < 1e-15);
        assert_eq!(c.h(0.2), 3.0 + 12.0 * 0.04); // input value (999) is ignored in legacy mode by design
    }
}
