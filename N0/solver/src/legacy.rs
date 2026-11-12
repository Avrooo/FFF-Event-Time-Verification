//! R8.4 ablation switches: restore ONE original (Bhandari 2020) implementation behaviour at a time in the corrected solver.
//! All default to false (= corrected behaviour, bit-identical to R6/R8.2). They exist only for the controlled error-attribution study.
//! * `legacy_cp_step`  D5: use the DECLARED specific-heat table step instead of the table's real spacing
//! * `legacy_timing`   D1-D3: partial deposition step integrates a full dt; cooldown integrates the stale first interval; trailing partial step is not committed
//! * `legacy_bed_flag` D4: bed contact flag is `layer == 0` instead of geometric (z = zmin)
//! * `legacy_mesh`     D8: node count rounded and domain divided into it (spacing != element size)
//! (original bed-contact law M0 and legacy h(z) are the existing `bed_contact_model 0` / `convection_model 1` keys)
use std::sync::atomic::{AtomicBool, Ordering};
pub static LEGACY_CP_STEP: AtomicBool = AtomicBool::new(false);
pub fn cp_step() -> bool { LEGACY_CP_STEP.load(Ordering::Relaxed) }
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Legacy { pub cp_step: bool, pub timing: bool, pub bed_flag: bool, pub mesh: bool }
impl Legacy {
    pub fn from_keys(get: impl Fn(&str) -> Option<usize>) -> Legacy {
        let f = |k: &str| get(k).unwrap_or(0) != 0;
        Legacy { cp_step: f("legacy_cp_step"), timing: f("legacy_timing"), bed_flag: f("legacy_bed_flag"), mesh: f("legacy_mesh") }
    }
}
