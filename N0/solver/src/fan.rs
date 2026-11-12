//! Fan-aware forced convection (R8.2). Kept separate from `convection.rs` so the frozen constant-h baseline is untouched.
//!
//! Model (documented assumptions):
//! * The part-cooling fan state is a step function of solver time, read from a schedule file (`fan_schedule_file`,
//!   lines `time_s,speed` with speed = the G-code `S` value 0..255; M106 S<n> -> n, M107 -> 0). The schedule is produced
//!   from the executed G-code by `fan_schedule.py`; the solver does not parse M106/M107 itself.
//! * Forced coefficient: h_forced(S) = h_full * (S/255)^n, with `fan_h_forced_full` (= h at S255) and exponent
//!   `fan_speed_exponent` (default 0.6, between laminar 0.5 and turbulent 0.8 flat-plate scaling).
//! * Mixed convection: h_eff = (h_nat^3 + h_forced^3)^(1/3) (Churchill 1977 superposition for mixed convection). When the fan is
//!   off (or the fan model is disabled) h_eff == h_nat EXACTLY (no arithmetic), so the baseline is bit-identical.
//! * The same h_eff is applied to every exposed face of every element (no directional duct / shadowing model).
//! * Radiation is unaffected by the fan.

#[derive(Clone, Debug, PartialEq)]
pub struct FanModel {
    times: Vec<f64>,
    speeds: Vec<f64>,
    h_full: f64,
    exponent: f64,
}

impl FanModel {
    pub fn off() -> FanModel { FanModel { times: vec![], speeds: vec![], h_full: 0.0, exponent: 0.6 } }
    pub fn new(mut events: Vec<(f64, f64)>, h_full: f64, exponent: f64) -> FanModel {
        assert!(h_full.is_finite() && h_full >= 0.0, "fan_h_forced_full must be >= 0");
        assert!(exponent > 0.0 && exponent <= 1.0, "fan_speed_exponent must be in (0,1]");
        events.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for e in &events { assert!(e.1 >= 0.0 && e.1 <= 255.0, "fan speed must be within 0..255"); }
        FanModel { times: events.iter().map(|e| e.0).collect(), speeds: events.iter().map(|e| e.1).collect(), h_full, exponent }
    }
    pub fn from_schedule_text(text: &str, h_full: f64, exponent: f64) -> FanModel {
        let mut ev = vec![];
        for l in text.lines() {
            let l = l.trim();
            if l.is_empty() || l.starts_with('#') || l.starts_with("time") { continue; }
            let mut p = l.split(',');
            let t: f64 = p.next().unwrap().trim().parse().expect("fan schedule time");
            let s: f64 = p.next().expect("fan schedule speed").trim().parse().expect("fan schedule speed");
            ev.push((t, s));
        }
        FanModel::new(ev, h_full, exponent)
    }
    /// fan speed (0..255) at solver time `t` (step function; 0 before the first event)
    pub fn speed_at(&self, t: f64) -> f64 {
        match self.times.iter().rposition(|&x| x <= t) { Some(i) => self.speeds[i], None => 0.0 }
    }
    /// forced coefficient [W/m2K] at solver time `t`
    pub fn h_forced(&self, t: f64) -> f64 {
        let s = self.speed_at(t) / 255.0;
        if s <= 0.0 || self.h_full == 0.0 { 0.0 } else { self.h_full * s.powf(self.exponent) }
    }
}

/// effective film coefficient for natural coefficient `h_nat` and forced coefficient `h_forced` (exactly `h_nat` when 0)
#[inline]
pub fn combine(h_nat: f64, h_forced: f64) -> f64 {
    if h_forced == 0.0 { h_nat } else { (h_nat * h_nat * h_nat + h_forced * h_forced * h_forced).cbrt() }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn off_is_identity_bitwise() {
        let f = FanModel::off();
        assert_eq!(f.h_forced(123.0), 0.0);
        for h in [0.0, 3.0, 12.7, 1e3] { assert_eq!(combine(h, f.h_forced(5.0)).to_bits(), h.to_bits()); }
    }
    #[test]
    fn schedule_is_a_step_function() {
        let f = FanModel::from_schedule_text("time_s,speed\n100,128\n500,0\n", 40.0, 0.6);
        assert_eq!(f.speed_at(99.9), 0.0); assert_eq!(f.speed_at(100.0), 128.0); assert_eq!(f.speed_at(499.0), 128.0); assert_eq!(f.speed_at(500.0), 0.0);
        assert!((f.h_forced(200.0) - 40.0 * (128.0f64 / 255.0).powf(0.6)).abs() < 1e-12);
    }
    #[test]
    fn mixed_convection_limits() {
        assert!((combine(3.0, 1e-3) - 3.0).abs() < 1e-6);               // weak forcing -> natural
        assert!((combine(1e-3, 30.0) - 30.0).abs() < 1e-6);             // strong forcing -> forced
        let h = combine(3.0, 30.0); assert!(h > 30.0 && h < 33.0);      // between max and sum
    }
}
