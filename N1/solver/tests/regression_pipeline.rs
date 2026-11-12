//! Pipeline-level regression tests: run the built binary on tiny G-code cases whose answer is analytic.
//! They exercise `update_model` (deposition clock, partial steps, cooldown), i.e. defects D1 and D2.
//!
//! Isolated elements, constant c, no radiation, adiabatic bed => lumped cooling
//!   T(t) = T_env + (T_dep - T_env) exp(-t/tau),  tau = rho V c / (h A)
//! with h = 3 + 12 z^2 (the solver's built-in film coefficient, z = element centroid height) and
//! A = 4 w ht + w^2 (four sides + top; bottom face is adiabatic on the bed with bed_k = 0).
//! The output clock starts at the first deposition.
use std::fs;
use std::path::PathBuf;
use std::process::Command;

const W: f64 = 15.875e-3;
const HT: f64 = 5.08e-3;
const RHO: f64 = 1271.1;
const T_DEP: f64 = 473.15;
const T_ENV: f64 = 313.15;

fn case_dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("des_regress_{}_{}", name, std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(d.join("inputfiles")).unwrap();
    fs::create_dir_all(d.join("outputfiles")).unwrap();
    d
}

/// Legacy-compatible helper used by the R4 tests: constant c, legacy film law h(z)=3+12z^2 (convection_model 1).
fn write_input(dir: &PathBuf, gcode: &str, c: f64, dt: f64, cooldown: f64) {
    let cp: String = (0..500).map(|i| format!("{},{}\n", c, 288.15 + i as f64)).collect();
    write_input_ext(dir, gcode, &cp, dt, cooldown, &[("convection_model", "1".to_string())]);
}

/// General helper: `cp_rows` is the text of the c(T) table (rows "c,T", 1 K spacing); `extra` keys replace an
/// existing input line or are appended (optional keys such as convection_model / bed_contact_model).
fn write_input_ext(dir: &PathBuf, gcode: &str, cp_rows: &str, dt: f64, cooldown: f64, extra: &[(&str, String)]) {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let base = fs::read_to_string(format!("{}/inputfiles/Input_file.txt", manifest)).unwrap();
    fs::write(dir.join("inputfiles/cp.csv"), cp_rows).unwrap();
    fs::write(dir.join("inputfiles/test.gcode"), gcode).unwrap();
    let mut sets: Vec<(&str, String)> = vec![
        ("GCode_file", "inputfiles/test.gcode".to_string()),
        ("specific_heat_capacity_file", "inputfiles/cp.csv".to_string()),
        ("specific_heat_capacity_file_temp_step", "1".to_string()),
        ("emissivity", "0.0".to_string()),
        ("conductivity_at_bed_surface_interface", "0.0".to_string()),
        ("time_step", format!("{}", dt)),
        ("cooldown_period", format!("{}", cooldown)),
        ("element_min_temp_diff_store", "0.00001".to_string()),
        ("Number_of_threads", "1".to_string()),
    ];
    for (k, v) in extra { sets.retain(|(k2, _)| k2 != k); sets.push((k, v.clone())); }
    let mut out = String::new();
    let mut used = vec![false; sets.len()];
    for line in base.lines() {
        let mut replaced = false;
        for (i, (k, v)) in sets.iter().enumerate() {
            if line.starts_with(&format!("{}:", k)) {
                out.push_str(&format!("{}: {} //\n", k, v)); // the input parser splits on spaces and needs a trailing token
                used[i] = true; replaced = true;
                break;
            }
        }
        if !replaced { out.push_str(line); out.push('\n'); }
    }
    for (i, (k, v)) in sets.iter().enumerate() { if !used[i] { out.push_str(&format!("{}: {} //\n", k, v)); } }
    fs::write(dir.join("inputfiles/Input_file.txt"), out).unwrap();
}

fn run(dir: &PathBuf) -> (Vec<f64>, Vec<f64>) {
    let st = Command::new(env!("CARGO_BIN_EXE_DES_thermal_simulation")).current_dir(dir).output().unwrap();
    assert!(st.status.success(), "solver failed: {}", String::from_utf8_lossy(&st.stderr));
    let txt = fs::read_to_string(dir.join("outputfiles/elem_temps.csv")).unwrap();
    let first = txt.lines().next().unwrap();
    let v: Vec<f64> = first.split(',').filter(|s| !s.trim().is_empty()).map(|s| s.parse().unwrap()).collect();
    (v.iter().step_by(2).cloned().collect(), v.iter().skip(1).step_by(2).cloned().collect())
}

fn interp(t: &[f64], y: &[f64], x: f64) -> f64 {
    let i = t.partition_point(|v| *v <= x).clamp(1, t.len() - 1);
    y[i - 1] + (y[i] - y[i - 1]) * (x - t[i - 1]) / (t[i] - t[i - 1])
}

fn tau(c: f64) -> f64 {
    let z = HT / 2.0;
    let h = 3.0 + 12.0 * z * z;
    RHO * W * W * HT * c / (h * (4.0 * W * HT + W * W))
}

const TWO_BEADS: &str = ";Printing starts here\n; layer 1\nG1 Z5.080 F2000\nG1 X5 Y10 F3240\nG1 X10 Y10 F3240\nG1 X10.5 Y10 E1 F3240\nG1 X73.5 Y10 F12000\nG1 X74 Y10 E1 F3240\n; layer end\n";

fn ten_beads() -> String {
    let mut g = String::from(";Printing starts here\n; layer 1\nG1 Z5.080 F2000\nG1 X5 Y10 F3240\nG1 X10 Y10 F3240\n");
    let mut x = 10.0;
    for i in 0..10 {
        g.push_str(&format!("G1 X{:.3} Y10 E1 F3240\n", x + 0.5));
        if i < 9 { g.push_str(&format!("G1 X{:.3} Y10 F12000\n", x + 63.5)); x += 63.5; }
    }
    g.push_str("; layer end\n");
    g
}

/// D2 (cooldown) + overall clock consistency: 400 s of cooldown after two isolated depositions.
#[test]
fn isolated_element_cooldown_matches_analytic_for_every_dt() {
    let c = 2000.0;
    for dt in [1.0, 0.5, 0.25, 0.1, 0.01] {
        let d = case_dir(&format!("cooldown_{}", dt));
        write_input(&d, TWO_BEADS, c, dt, 400.0);
        let (t, y) = run(&d);
        let at = 300.0;
        let exact = T_ENV + (T_DEP - T_ENV) * (-at / tau(c)).exp();
        let err = interp(&t, &y, at) - exact;
        println!("cooldown test dt={:>5}: error at t=300 s = {:+.4} K", dt, err);
        assert!(err.abs() < 0.05, "dt={}: {:+.3} K from analytic (published code: +20 K .. -129 K)", dt, err);
    }
}

/// D1 (partial step at the end of each deposition interval): ten isolated beads, 0.32 s apart, no cooldown.
#[test]
fn deposition_phase_physics_time_equals_clock_time() {
    let c = 200.0;
    for dt in [1.0, 0.5, 0.3, 0.1, 0.01] {
        let d = case_dir(&format!("deposit_{}", dt));
        write_input(&d, &ten_beads(), c, dt, 0.0);
        let (t, y) = run(&d);
        let tl = *t.last().unwrap();
        let exact = T_ENV + (T_DEP - T_ENV) * (-tl / tau(c)).exp();
        let err = y.last().unwrap() - exact;
        println!("deposition test dt={:>5}: error at last stored clock ({:.3} s) = {:+.4} K", dt, tl, err);
        assert!(err.abs() < 0.03, "dt={}: {:+.3} K from analytic (published code: -4.5 .. -0.1 K)", dt, err);
    }
}

/// Result must not depend on dt beyond first-order Euler error.
#[test]
fn result_is_nearly_dt_independent() {
    let c = 200.0;
    let mut vals = vec![];
    for dt in [1.0, 0.1] {
        let d = case_dir(&format!("dtind_{}", dt));
        write_input(&d, &ten_beads(), c, dt, 0.0);
        let (t, y) = run(&d);
        vals.push(interp(&t, &y, 2.5));
    }
    assert!((vals[0] - vals[1]).abs() < 0.03, "dt=1 vs dt=0.1 differ by {} K", vals[0] - vals[1]);
}

fn tau_h(c: f64, h: f64) -> f64 { RHO * W * W * HT * c / (h * (4.0 * W * HT + W * W)) }

/// D6: the input film coefficient must control the cooling (it was silently ignored by the published type-2 solver).
#[test]
fn input_film_coefficient_is_used_and_matches_analytic() {
    let c = 2000.0;
    let cp: String = (0..500).map(|i| format!("{},{}\n", c, 288.15 + i as f64)).collect();
    let mut at300 = vec![];
    for h in [3.0, 25.0, 100.0] {
        let d = case_dir(&format!("conv_{}", h));
        write_input_ext(&d, TWO_BEADS, &cp, 0.25, 400.0, &[("convective_heat_transfer_film_coefficient", format!("{}", h))]);
        let (t, y) = run(&d);
        let exact = T_ENV + (T_DEP - T_ENV) * (-300.0 / tau_h(c, h)).exp();
        let got = interp(&t, &y, 300.0);
        println!("constant h={:>5}: T(300 s) = {:.3}, analytic {:.3}", h, got, exact);
        assert!((got - exact).abs() < 0.05, "h={}: {:+.3} K from analytic", h, got - exact);
        at300.push(got);
    }
    assert!(at300[0] - at300[2] > 20.0, "h=3 vs h=100 must differ strongly, got {:?}", at300);
}

/// D6: the legacy height-dependent law stays available only by explicit selection and ignores the input h.
#[test]
fn legacy_height_dependent_convection_is_explicit_option() {
    let c = 2000.0;
    let cp: String = (0..500).map(|i| format!("{},{}\n", c, 288.15 + i as f64)).collect();
    let d = case_dir("legacy_conv");
    write_input_ext(&d, TWO_BEADS, &cp, 0.25, 400.0, &[("convective_heat_transfer_film_coefficient", "500".to_string()), ("convection_model", "1".to_string())]);
    let (t, y) = run(&d);
    let exact = T_ENV + (T_DEP - T_ENV) * (-300.0 / tau(c)).exp();
    assert!((interp(&t, &y, 300.0) - exact).abs() < 0.05);
}

/// D5: a wrongly declared table spacing must not change the curve that is read. c(T) is a ramp; the result is
/// compared with a fine explicit integration of the same ODE using the exact table.
#[test]
fn declared_table_step_cannot_corrupt_cp_curve() {
    let cfun = |t: f64| 1500.0 + 5.0 * (t - 288.15);          // 1500 .. 3500 J/kgK over the table range
    let cp: String = (0..500).map(|i| format!("{},{}\n", cfun(288.15 + i as f64), 288.15 + i as f64)).collect();
    let h = 25.0; let a = 4.0 * W * HT + W * W; let m = RHO * W * W * HT;
    let (mut tt, mut time, dt) = (T_DEP, 0.0, 0.001);
    while time < 300.0 { tt -= dt * h * a * (tt - T_ENV) / (m * cfun(tt)); time += dt; }
    let mut got = vec![];
    for declared in ["1", "5", "0.25"] {
        let d = case_dir(&format!("cpstep_{}", declared));
        write_input_ext(&d, TWO_BEADS, &cp, 0.1, 400.0, &[("convective_heat_transfer_film_coefficient", format!("{}", h)), ("specific_heat_capacity_file_temp_step", declared.to_string())]);
        let (t, y) = run(&d);
        got.push(interp(&t, &y, 300.0));
    }
    println!("T(300) declared 1/5/0.25: {:?}, reference {:.3}", got, tt);
    assert!(got.iter().all(|g| (g - tt).abs() < 0.15), "{:?} vs {}", got, tt);
    assert!(got[0] == got[1] && got[1] == got[2], "declared step must not matter: {:?}", got);
}

fn read_mesh(dir: &PathBuf) -> (Vec<[f64; 3]>, Vec<Vec<usize>>) {
    let nodes: Vec<[f64; 3]> = fs::read_to_string(dir.join("outputfiles/nodefile.csv")).unwrap().lines().map(|l| {
        let v: Vec<f64> = l.split(',').map(|s| s.trim().parse().unwrap()).collect(); [v[1], v[2], v[3]] }).collect();
    let els: Vec<Vec<usize>> = fs::read_to_string(dir.join("outputfiles/elementfile.csv")).unwrap().lines().map(|l| {
        l.split(',').skip(1).map(|s| s.trim().parse::<usize>().unwrap()).collect() }).collect();
    (nodes, els)
}

/// D8: generated cell sizes must equal beadwidth/div and beadheight/div exactly, for extents that are NOT whole
/// numbers of elements (published code: 15.097 .. 15.969 mm for 15.875 mm elements, up to 4.9 %).
#[test]
fn generated_mesh_spacing_equals_element_size() {
    let cp: String = (0..500).map(|i| format!("{},{}\n", 1900.0, 288.15 + i as f64)).collect();
    for (len, dw, dh) in [(120.0, 1, 1), (60.0, 2, 2), (48.0, 1, 3), (77.7, 3, 1)] {
        let g = format!(";Printing starts here\n; layer 1\nG1 Z5.080 F2000\nG1 X7.935 Y7.935 F3240\nG1 X{} Y7.935 E1 F3240\n; layer 2\nG1 Z10.160 F2000\nG1 X7.935 Y7.935 E1 F3240\n", 7.935 + len);
        let d = case_dir(&format!("mesh_{}_{}_{}", len, dw, dh));
        write_input_ext(&d, &g, &cp, 0.5, 1.0, &[("divisions_per_bead_width", format!("{}", dw)), ("divisions_per_bead_height", format!("{}", dh))]);
        let _ = run(&d);
        let (nodes, els) = read_mesh(&d);
        let (ew, eh) = (W / dw as f64, HT / dh as f64);
        for e in &els {
            let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
            for &n in e { for k in 0..3 { lo[k] = lo[k].min(nodes[n - 1][k]); hi[k] = hi[k].max(nodes[n - 1][k]); } }
            assert!(((hi[0] - lo[0]) - ew).abs() < 1e-9 * ew, "len {} dx {} vs {}", len, hi[0] - lo[0], ew);
            assert!(((hi[1] - lo[1]) - ew).abs() < 1e-9 * ew, "dy {} vs {}", hi[1] - lo[1], ew);
            assert!(((hi[2] - lo[2]) - eh).abs() < 1e-9 * eh, "dz {} vs {}", hi[2] - lo[2], eh);
        }
    }
}

/// Bed-contact models M0/M1/M2 against their closed forms: isolated bed element, h = 0, no radiation, so
///   T(t) = T_bed + (T_dep - T_bed) exp(-t G / (rho V c)),   G = bed-face conductance of the model.
#[test]
fn bed_contact_models_match_closed_form_cooling() {
    let c = 2000.0; let (kz, kb, rc) = (0.35, 0.17, 5e-3);
    let cp: String = (0..500).map(|i| format!("{},{}\n", c, 288.15 + i as f64)).collect();
    let a = W * W; let m = RHO * a * HT; let t_bed = 363.15;
    let cases = [
        ("M0", vec![("bed_contact_model", "0".to_string())], 0.17 * a / HT),
        ("M1", vec![("bed_contact_model", "1".to_string())], 2.0 * kz * a / HT),
        ("M2", vec![("bed_contact_model", "2".to_string()), ("bed_contact_resistance", format!("{}", rc))], a / (HT / (2.0 * kz) + rc)),
    ];
    for (name, mut extra, g) in cases {
        extra.push(("conductivity_at_bed_surface_interface", format!("{}", kb)));
        extra.push(("convective_heat_transfer_film_coefficient", "0.0".to_string()));
        extra.push(("kz", format!("{}", kz)));
        let d = case_dir(&format!("bed_{}", name));
        write_input_ext(&d, TWO_BEADS, &cp, 0.1, 300.0, &extra);
        let (t, y) = run(&d);
        let exact = t_bed + (T_DEP - t_bed) * (-100.0 * g / (m * c)).exp();
        let got = interp(&t, &y, 100.0);
        println!("{}: T(100 s) = {:.3}, analytic {:.3}", name, got, exact);
        assert!((got - exact).abs() < 0.08, "{}: {:+.3} K from analytic", name, got - exact);
    }
}

// ---------------------------------------------------------------- R8.2: fan-aware convection
fn fan_keys(d: &PathBuf, sched: &str, hfull: f64) -> Vec<(&'static str, String)> {
    let p = d.join("inputfiles/fan.csv");
    fs::write(&p, sched).unwrap();
    vec![("fan_model", "1".to_string()), ("fan_schedule_file", "inputfiles/fan.csv".to_string()),
         ("fan_h_forced_full", format!("{}", hfull)), ("fan_speed_exponent", "0.6".to_string())]
}

/// Fan off (speed 0 for the whole run) must reproduce the no-fan run BIT-FOR-BIT; fan always on must follow the analytic
/// lumped cooling with h_eff = (h_nat^3 + h_forced^3)^(1/3); a mid-run switch must follow the piecewise analytic solution.
#[test]
fn fan_model_matches_analytic_and_off_is_bit_identical() {
    let c = 2000.0;
    let cp: String = (0..500).map(|i| format!("{},{}\n", c, 288.15 + i as f64)).collect();
    // reference: no fan keys
    let d0 = case_dir("fan_none"); write_input_ext(&d0, TWO_BEADS, &cp, 0.25, 400.0, &[]); let (t0, y0) = run(&d0);
    // fan model active but speed 0 everywhere
    let d1 = case_dir("fan_zero"); let k1 = fan_keys(&d1, "time_s,speed\n0,0\n", 40.0);
    write_input_ext(&d1, TWO_BEADS, &cp, 0.25, 400.0, &k1); let (t1, y1) = run(&d1);
    assert_eq!(t0.len(), t1.len());
    for i in 0..y0.len() { assert_eq!(y0[i].to_bits(), y1[i].to_bits(), "fan speed 0 changed the result at sample {}", i); }
    // fan always on at S255: forced h = 40
    let hf: f64 = 40.0; let heff = (27.0f64 + hf.powi(3)).cbrt();
    let d2 = case_dir("fan_on"); let k2 = fan_keys(&d2, "time_s,speed\n0,255\n", hf);
    write_input_ext(&d2, TWO_BEADS, &cp, 0.25, 400.0, &k2); let (t2, y2) = run(&d2);
    let exact = T_ENV + (T_DEP - T_ENV) * (-300.0 / tau_h(c, heff)).exp(); let got = interp(&t2, &y2, 300.0);
    println!("fan on: h_eff={:.3} T(300)={:.3} analytic {:.3}", heff, got, exact);
    assert!((got - exact).abs() < 0.05, "fan-on: {:+.3} K from analytic", got - exact);
    // switch: natural until 150 s, forced afterwards
    let d3 = case_dir("fan_switch"); let k3 = fan_keys(&d3, "time_s,speed\n0,0\n150,255\n", hf);
    write_input_ext(&d3, TWO_BEADS, &cp, 0.25, 400.0, &k3); let (t3, y3) = run(&d3);
    let t_150 = T_ENV + (T_DEP - T_ENV) * (-150.0 / tau_h(c, 3.0)).exp();
    let exact3 = T_ENV + (t_150 - T_ENV) * (-150.0 / tau_h(c, heff)).exp(); let got3 = interp(&t3, &y3, 300.0);
    println!("switch: T(300)={:.3} analytic {:.3}", got3, exact3);
    assert!((got3 - exact3).abs() < 0.1, "fan switch: {:+.3} K from analytic", got3 - exact3);
    let _ = t1;
}
