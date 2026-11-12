//! Bed-contact boundary models for the face of a bed-contact element that rests on the print bed.
//!
//! The bed is held at a fixed temperature T_bed. For an element of height dz and plan area A (top_area),
//! `coefficient()` returns the heat-flow conductance G [W/K] such that the heat leaving the element through
//! its bottom face is  Q = G (T - T_bed)  [W].
//!
//! * `Original` (M0): Bhandari's formulation, G = k_bed A / dz. The "distance" is the element height, so the
//!   conductance per area k_bed/dz changes whenever the mesh is refined: the boundary condition is not
//!   mesh-consistent.
//! * `Perfect` (M1): perfect thermal contact with the bed at T_bed. Heat flows from the element centroid to the
//!   bottom face through half an element of the printed material: G = k_z A / (dz/2). As dz -> 0 this tends to
//!   the Dirichlet condition T(z=0) = T_bed, independent of the mesh.
//! * `Resistance(Rc)` (M2): as M1 plus a finite interface resistance Rc [m^2 K/W] in series:
//!   G = A / (dz/(2 k_z) + Rc). As dz -> 0 this tends to the Robin condition -k dT/dz = (T - T_bed)/Rc,
//!   independent of the mesh. Rc -> 0 recovers M1.
//!
//! The solvers pass `ta_d_td` = A/dz (already available in their stencil), from which dz = A / ta_d_td.

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BedContact {
    Original,
    Perfect,
    Resistance(f64),
}

impl BedContact {
    /// Build from the optional trailing input keys `bed_contact_model` and `bed_contact_resistance`.
    /// A missing model key means `Original`, so existing input files behave exactly as before.
    pub fn from_input(model: Option<usize>, resistance: Option<f64>) -> BedContact {
        match model.unwrap_or(0) {
            0 => BedContact::Original,
            1 => BedContact::Perfect,
            2 => {
                let rc = resistance.expect("bed_contact_model 2 needs bed_contact_resistance [m2.K/W]");
                assert!(rc.is_finite() && rc >= 0.0, "bed_contact_resistance must be >= 0, got {}", rc);
                BedContact::Resistance(rc)
            }
            other => panic!("unknown bed_contact_model {} (use 0, 1 or 2)", other),
        }
    }

    /// Conductance G [W/K] of the bed face of one element (see module docs).
    /// `bed_k`: input `conductivity_at_bed_surface_interface` (used by M0 only);
    /// `kz`: through-thickness conductivity of the printed material;
    /// `top_area` = A; `ta_d_td` = A/dz.
    pub fn coefficient(&self, bed_k: f64, kz: f64, top_area: f64, ta_d_td: f64) -> f64 {
        match *self {
            BedContact::Original => bed_k * ta_d_td,
            BedContact::Perfect => 2.0 * kz * ta_d_td,
            BedContact::Resistance(rc) => top_area / (top_area / (2.0 * kz * ta_d_td) + rc),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coefficients_match_closed_forms() {
        let (a, dz, kz, kb) = (2.52e-4, 5.08e-3, 0.35, 0.17);
        let g = a / dz;
        assert_eq!(BedContact::Original.coefficient(kb, kz, a, g), kb * g);
        assert!((BedContact::Perfect.coefficient(kb, kz, a, g) - kz * a / (dz / 2.0)).abs() < 1e-15);
        let rc = 5e-3;
        let m2 = BedContact::Resistance(rc).coefficient(kb, kz, a, g);
        assert!((m2 - a / (dz / (2.0 * kz) + rc)).abs() < 1e-15);
        // Rc -> 0 recovers perfect contact; Rc -> infinity isolates the bed
        assert!((BedContact::Resistance(0.0).coefficient(kb, kz, a, g) - BedContact::Perfect.coefficient(kb, kz, a, g)).abs() < 1e-12);
        assert!(BedContact::Resistance(1e12).coefficient(kb, kz, a, g) < 1e-9);
    }

    #[test]
    fn m1_and_m2_are_mesh_consistent_m0_is_not() {
        // conductance PER AREA as the element is refined vertically
        let (a, kz, kb, rc) = (1.0, 0.35, 0.17, 5e-3);
        let per_area = |bc: BedContact, dz: f64| bc.coefficient(kb, kz, a, a / dz);
        let (c_m0, c_m1, c_m2) = (
            [1, 2, 4, 8].map(|n| per_area(BedContact::Original, 5.08e-3 / n as f64)),
            [1, 2, 4, 8].map(|n| per_area(BedContact::Perfect, 5.08e-3 / n as f64)),
            [1, 2, 4, 8].map(|n| per_area(BedContact::Resistance(rc), 5.08e-3 / n as f64)),
        );
        assert!((c_m0[3] / c_m0[0] - 8.0).abs() < 1e-9, "M0 conductance per area scales with 1/dz");
        // M2 saturates at 1/Rc as dz -> 0 (finite, mesh-independent limit)
        assert!((c_m2[3] - 1.0 / rc).abs() / (1.0 / rc) < 0.2 && c_m2[3] < 1.0 / rc);
        // M1 grows like 1/dz BUT its meaning is a Dirichlet limit (infinite conductance), which is the correct limit for perfect contact
        assert!(c_m1[3] > c_m1[0]);
    }

    #[test]
    fn input_parsing_defaults_to_original() {
        assert_eq!(BedContact::from_input(None, None), BedContact::Original);
        assert_eq!(BedContact::from_input(Some(1), None), BedContact::Perfect);
        assert_eq!(BedContact::from_input(Some(2), Some(0.005)), BedContact::Resistance(0.005));
    }
}
