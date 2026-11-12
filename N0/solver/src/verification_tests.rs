//! R4 verification tests for `model_orthotropic_td_shc_variable_h::Model` (solver for `model_type: 2`).
//!
//! These tests drive the solver's own stencil (`addCell` + `run_model`) on a fixed domain. They do NOT
//! change any solver physics: the Dirichlet boundaries of the manufactured-solution benchmark are
//! imposed through the solver's existing exposed-face film condition (h_eff = 2k/dx, T_env = T_boundary),
//! which is the standard half-cell Dirichlet condition. Run with `cargo test --release -- --nocapture`.
//! Set VERIF_OUT=<dir> to also write CSV result tables.
//!
//! Methodology: Salari & Knupp (SAND2000-1444) sec. 4.2 (observed order p = log(E1/E2)/log r, with spatial
//! and temporal errors separated); Oberkampf & Trucano (2002) for code-vs-numerical-error terminology.

use super::*;
use crate::interpolator::Interpolator;
use std::fs::File;
use std::io::{BufWriter, Write};

const PI: f64 = std::f64::consts::PI;

fn sink() -> BufWriter<File> {
    let p = std::env::temp_dir().join(format!("des_verif_sink_{}.tmp", std::process::id()));
    BufWriter::new(File::create(p).unwrap())
}
fn pool() -> ThreadPool {
    rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap()
}
fn const_c(c: f64) -> Interpolator {
    let x: Vec<f64> = (0..600).map(|i| 200.0 + i as f64).collect();
    Interpolator::new(x, vec![c; 600], 1.0, 1)
}
fn pla_table() -> Interpolator {
    let p = format!("{}/inputfiles/sp_heat_cap_data.csv", env!("CARGO_MANIFEST_DIR"));
    // the shipped table has 1 K spacing (the shipped Input_file.txt wrongly says 5; see report D5)
    Interpolator::read_data_from_file(&p, 1.0, 1)
}
fn out_dir() -> Option<std::path::PathBuf> {
    std::env::var("VERIF_OUT").ok().map(|d| {
        let p = std::path::PathBuf::from(d);
        std::fs::create_dir_all(&p).unwrap();
        p
    })
}

// ---------------------------------------------------------------------------------------------
// Manufactured solution (Roache/Salari-Knupp style):
//   T = T0 + A exp(-beta t) sin(pi x/L) sin(pi y/L) sin(pi z/L),  beta = alpha pi^2 (3/L^2)
// exact solution of rho c T_t = k lap T with T = T0 on all six faces of the cube [0,L]^3.
// ---------------------------------------------------------------------------------------------
#[derive(Clone, Copy)]
struct Mms { k: f64, rho: f64, c: f64, l: f64, t0: f64, amp: f64 }
impl Mms {
    fn pla_like() -> Mms { Mms { k: 0.2, rho: 1240.0, c: 1800.0, l: 0.04, t0: 300.0, amp: 50.0 } }
    fn alpha(&self) -> f64 { self.k / (self.rho * self.c) }
    fn beta(&self) -> f64 { self.alpha() * PI * PI * 3.0 / (self.l * self.l) }
    fn exact(&self, x: f64, y: f64, z: f64, t: f64) -> f64 {
        self.t0 + self.amp * (-self.beta() * t).exp()
            * (PI * x / self.l).sin() * (PI * y / self.l).sin() * (PI * z / self.l).sin()
    }
}

struct MmsRun { linf: f64, l2: f64, energy_err_rel: f64, field: Vec<f64> }

fn run_mms(n: usize, nsteps: usize, t_end: f64, m: &Mms) -> MmsRun {
    let d = m.l / n as f64;
    let ncell = n * n * n;
    let mut nodes = Vec::with_capacity(8 * ncell);
    for k in 0..n { for j in 0..n { for i in 0..n {
        for dz in 0..2 { for dy in 0..2 { for dx in 0..2 {
            nodes.push([(i + dx) as f64 * d, (j + dy) as f64 * d, (k + dz) as f64 * d]);
        }}}
    }}}
    // zmin = -1 => no cell touches the bed (is_bot_cell false). All cells are deliberately given layer 0:
    // the pre-D4 rule `layer == 0` would have turned every cell into a bed cell and failed this test.
    let mut mdl = Model::new(nodes, ncell, const_c(m.c), d, d, 1000, -1.0, m.rho, 0.0, 0.0, m.t0, 1);
    let mut tmp = Vec::with_capacity(ncell);
    let mut cell_indices = Vec::with_capacity(ncell);
    let mut storer: Vec<Vec<f64>> = Vec::new();
    for k in 0..n { for j in 0..n { for i in 0..n {
        let idx = i + n * (j + n * k);
        let mut isn = [false; 6];
        let mut nb = [0isize; 6];
        if i > 0 { isn[1] = true; nb[1] = (idx - 1) as isize; }
        if j > 0 { isn[3] = true; nb[3] = (idx - n) as isize; }
        if k > 0 { isn[5] = true; nb[5] = (idx - n * n) as isize; }
        let nn = [8 * idx, 8 * idx + 1, 8 * idx + 2, 8 * idx + 3, 8 * idx + 4, 8 * idx + 5, 8 * idx + 6, 8 * idx + 7];
        mdl.addCell(nn, m.c, m.rho, m.k, m.k, m.k, m.t0, isn, nb, idx, &mut tmp, [1.0, 0.0], 0, &mut cell_indices, &mut storer);
    }}}
    let h_eff = 2.0 * m.k / d; // half-cell Dirichlet: flux = h_eff (T_c - T0) = k (T_c - T0)/(d/2)
    let centre = |i: usize| (i as f64 + 0.5) * d;
    for k in 0..n { for j in 0..n { for i in 0..n {
        let idx = i + n * (j + n * k);
        mdl.templist[idx] = m.exact(centre(i), centre(j), centre(k), 0.0);
        mdl.h[idx] = h_eff;
    }}}
    let areas = [d * d, d * d, d, d, d * d * d];
    let dt = t_end / nsteps as f64;
    let p = pool();
    let mut bw = sink();
    let mut gt = 0.0;
    let e0: f64 = mdl.templist.iter().map(|t| m.rho * m.c * d * d * d * t).sum();
    let mut loss = 0.0;
    for _ in 0..nsteps {
        // independent (test-side) boundary-loss accounting from the pre-step state
        for c in 0..ncell {
            let nexp = (0..6).filter(|f| !mdl.has_neighbor[c][*f]).count() as f64;
            loss += nexp * h_eff * (mdl.templist[c] - m.t0) * d * d * dt;
        }
        mdl.run_model(dt, dt, 0.0, m.t0, &mut bw, &mut gt, &p, areas, &mut tmp, 1, &cell_indices);
    }
    let e1: f64 = mdl.templist.iter().map(|t| m.rho * m.c * d * d * d * t).sum();
    let (mut linf, mut s2) = (0.0f64, 0.0f64);
    for k in 0..n { for j in 0..n { for i in 0..n {
        let idx = i + n * (j + n * k);
        let e = (mdl.templist[idx] - m.exact(centre(i), centre(j), centre(k), t_end)).abs();
        linf = linf.max(e);
        s2 += e * e;
    }}}
    MmsRun { linf, l2: (s2 / ncell as f64).sqrt(), energy_err_rel: ((e1 - e0) + loss).abs() / (e1 - e0).abs(), field: mdl.templist.clone() }
}

fn order(e1: f64, e2: f64, r: f64) -> f64 { (e1 / e2).ln() / r.ln() }

/// (Linf, L2) error of a cell-centred field against the manufactured solution at time t.
fn err_norms(field: &[f64], n: usize, t: f64, m: &Mms) -> (f64, f64) {
    let d = m.l / n as f64;
    let (mut linf, mut s2) = (0.0f64, 0.0f64);
    for k in 0..n { for j in 0..n { for i in 0..n {
        let idx = i + n * (j + n * k);
        let e = (field[idx] - m.exact((i as f64 + 0.5) * d, (j as f64 + 0.5) * d, (k as f64 + 0.5) * d, t)).abs();
        linf = linf.max(e);
        s2 += e * e;
    }}}
    (linf, (s2 / field.len() as f64).sqrt())
}

#[test]
fn mms_spatial_convergence_is_second_order() {
    let m = Mms::pla_like();
    let t_end = 300.0;
    // Forward Euler leaves an O(dt) error that would contaminate the spatial order. It is removed by
    // first-order Richardson extrapolation in time, T* = 2 T(dt/2) - T(dt), leaving O(dt^2) ~ 1e-5 K.
    let (dt0, ns0) = (0.2, 1500usize);
    let sizes = [4usize, 8, 16, 32];
    let mut rows = vec![];
    println!("\nMMS spatial convergence (time-extrapolated from dt = {} and {} s; t_end = {} s)", dt0, dt0 / 2.0, t_end);
    println!("{:>4} {:>10} {:>12} {:>12} {:>8} {:>8} {:>12}", "N", "dx [mm]", "Linf [K]", "L2 [K]", "p(Linf)", "p(L2)", "E-bal err");
    let mut prev: Option<(f64, f64)> = None;
    let mut last_p = (0.0, 0.0);
    for &n in &sizes {
        let a = run_mms(n, ns0, t_end, &m);
        let b = run_mms(n, 2 * ns0, t_end, &m);
        let ext: Vec<f64> = a.field.iter().zip(&b.field).map(|(x, y)| 2.0 * y - x).collect();
        let (linf, l2) = err_norms(&ext, n, t_end, &m);
        let (pi, p2) = match prev { Some((x, y)) => (order(x, linf, 2.0), order(y, l2, 2.0)), None => (f64::NAN, f64::NAN) };
        println!("{:>4} {:>10.4} {:>12.3e} {:>12.3e} {:>8.3} {:>8.3} {:>12.2e}", n, m.l / n as f64 * 1e3, linf, l2, pi, p2, b.energy_err_rel);
        assert!(b.energy_err_rel < 1e-9, "energy balance must close to round-off for constant c");
        rows.push(format!("{},{},{},{},{},{},{}", n, m.l / n as f64, linf, l2, pi, p2, b.energy_err_rel));
        prev = Some((linf, l2));
        last_p = (pi, p2);
    }
    let _ = dt0;
    if let Some(dir) = out_dir() {
        let mut f = File::create(dir.join("mms_spatial.csv")).unwrap();
        writeln!(f, "N,dx_m,linf_K,l2_K,p_linf,p_l2,energy_err_rel").unwrap();
        for r in rows { writeln!(f, "{}", r).unwrap(); }
    }
    assert!((last_p.0 - 2.0).abs() < 0.15 && (last_p.1 - 2.0).abs() < 0.15, "observed spatial order {:?} is not ~2", last_p);
}

#[test]
fn mms_temporal_convergence_is_first_order() {
    let m = Mms::pla_like();
    let n = 16;
    let t_end = 240.0;
    // stability for N=16: dt < rho c dx^2 /(6 k) ~ 11.6 s
    let steps = [30usize, 60, 120, 240, 480, 960]; // dt = 8, 4, 2, 1, 0.5, 0.25 s
    let runs: Vec<MmsRun> = steps.iter().map(|&s| run_mms(n, s, t_end, &m)).collect();
    println!("\nMMS temporal convergence (N = {}, t_end = {} s). Self-convergence |T(dt) - T(dt/2)| removes the spatial error.", n, t_end);
    println!("{:>8} {:>12} {:>14} {:>9} {:>12}", "dt [s]", "Linf vs exact", "|T(dt)-T(dt/2)|", "p_t", "E-bal err");
    let mut diffs = vec![];
    let mut rows = vec![];
    for i in 0..runs.len() {
        let dt = t_end / steps[i] as f64;
        let diff = if i + 1 < runs.len() {
            runs[i].field.iter().zip(&runs[i + 1].field).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)
        } else { f64::NAN };
        let p = if i >= 1 && diff.is_finite() { order(diffs[i - 1], diff, 2.0) } else { f64::NAN };
        println!("{:>8.3} {:>12.3e} {:>14.3e} {:>9.3} {:>12.2e}", dt, runs[i].linf, diff, p, runs[i].energy_err_rel);
        assert!(runs[i].energy_err_rel < 1e-9);
        rows.push(format!("{},{},{},{},{}", dt, runs[i].linf, diff, p, runs[i].energy_err_rel));
        diffs.push(diff);
    }
    if let Some(dir) = out_dir() {
        let mut f = File::create(dir.join("mms_temporal.csv")).unwrap();
        writeln!(f, "dt_s,linf_vs_exact_K,selfdiff_K,p_t,energy_err_rel").unwrap();
        for r in rows { writeln!(f, "{}", r).unwrap(); }
    }
    let p_last = order(diffs[2], diffs[3], 2.0);
    assert!((p_last - 1.0).abs() < 0.15, "observed temporal order {} is not ~1 (forward Euler)", p_last);
}

#[test]
#[ignore] // ~1 min in release; adds the N = 64 point
fn mms_spatial_convergence_n64() {
    let m = Mms::pla_like();
    let ext = |n: usize| {
        let a = run_mms(n, 1500, 300.0, &m);
        let b = run_mms(n, 3000, 300.0, &m);
        let f: Vec<f64> = a.field.iter().zip(&b.field).map(|(x, y)| 2.0 * y - x).collect();
        err_norms(&f, n, 300.0, &m).0
    };
    let (e32, e64) = (ext(32), ext(64));
    println!("N=64: Linf {:.3e}  p(32->64) = {:.3}", e64, order(e32, e64, 2.0));
    assert!((order(e32, e64, 2.0) - 2.0).abs() < 0.15);
}

// ---------------------------------------------------------------------------------------------
// run_model must advance the state by exactly `time`, including a trailing partial step (D3).
// ---------------------------------------------------------------------------------------------
fn single_cell(c: f64, h: f64, t_init: f64) -> (Model, Vec<f64>, Vec<usize>) {
    let (w, ht) = (15.875e-3, 5.08e-3);
    let nodes: Vec<[f64; 3]> = (0..8).map(|i| [(i & 1) as f64 * w, ((i >> 1) & 1) as f64 * w, ((i >> 2) & 1) as f64 * ht]).collect();
    let mut mdl = Model::new(nodes, 1, const_c(c), w, ht, 1000, -1.0, 1271.1, 0.0, 0.0, 300.0, 1);
    let mut tmp = vec![];
    let mut ci = vec![];
    let mut st: Vec<Vec<f64>> = vec![];
    mdl.addCell([0, 1, 2, 3, 4, 5, 6, 7], c, 1271.1, 0.5, 0.5, 0.5, t_init, [false; 6], [0; 6], 0, &mut tmp, [1.0, 0.0], 0, &mut ci, &mut st);
    mdl.h[0] = h;
    (mdl, tmp, ci)
}

#[test]
fn run_model_advances_state_by_exactly_the_requested_time() {
    let (w, ht, rho, c, h, t0, tenv) = (15.875e-3, 5.08e-3, 1271.1, 2000.0, 20.0, 473.15, 313.15);
    let a = 6.0 * 0.0 + 4.0 * w * ht + 2.0 * w * w; // all six faces exposed (zmin = -1: no bed)
    let tau = rho * w * w * ht * c / (h * a);
    let areas = [w * ht, w * w, ht, w * w / ht, w * w * ht];
    let p = pool();
    for &(time, dt) in &[(0.37, 0.1), (0.05, 0.1), (1.0, 0.3), (2.0, 2.0), (7.0, 2.0)] {
        let (mut mdl, mut tmp, ci) = single_cell(c, h, t0);
        let mut bw = sink();
        let mut gt = 0.0;
        mdl.run_model(time, dt, 0.0, tenv, &mut bw, &mut gt, &p, areas, &mut tmp, 1, &ci);
        // exact composition of the explicit scheme: floor(time/dt) full steps plus one partial step.
        let n = (time / dt).floor() as i32;
        let rem = time - n as f64 * dt;
        let mut f = (1.0 - dt / tau).powi(n);
        if rem > 1e-12 { f *= 1.0 - rem / tau; }
        let expected = tenv + (t0 - tenv) * f;
        let err = (mdl.templist[0] - expected).abs();
        let trunc = (mdl.templist[0] - (tenv + (t0 - tenv) * (-time / tau).exp())).abs();
        println!("run_model(time={}, dt={}): bookkeeping error {:.2e} K  (scheme truncation vs exponential {:.2e} K)", time, dt, err, trunc);
        assert!(err < 1e-9, "state not advanced by exactly `time` (err {} K): trailing partial step lost?", err);
    }
}

// ---------------------------------------------------------------------------------------------
// D4: only elements whose bottom face lies on the bed are bed-contact cells.
// ---------------------------------------------------------------------------------------------
#[test]
fn bed_flag_is_geometric_not_layer_based() {
    let (w, ht) = (4.0e-3, 1.0e-3);
    // two vertically stacked sub-elements of the same bead layer (layer 0): z in [0, ht] and [ht, 2ht]
    let mut nodes: Vec<[f64; 3]> = vec![];
    for zb in [0.0, ht] {
        for i in 0..8 { nodes.push([(i & 1) as f64 * w, ((i >> 1) & 1) as f64 * w, zb + ((i >> 2) & 1) as f64 * ht]); }
    }
    let mut mdl = Model::new(nodes, 2, const_c(1800.0), w, ht, 1000, 0.0, 1240.0, 0.0, 0.17, 363.15, 1);
    let (mut tmp, mut ci, mut st) = (vec![], vec![], vec![]);
    mdl.addCell([0, 1, 2, 3, 4, 5, 6, 7], 1800.0, 1240.0, 0.2, 0.2, 0.2, 473.15, [false; 6], [0; 6], 0, &mut tmp, [1.0, 0.0], 0, &mut ci, &mut st);
    let mut isn = [false; 6];
    isn[5] = true;
    let mut nb = [0isize; 6];
    nb[5] = 0;
    mdl.addCell([8, 9, 10, 11, 12, 13, 14, 15], 1800.0, 1240.0, 0.2, 0.2, 0.2, 473.15, isn, nb, 1, &mut tmp, [1.0, 0.0], 0, &mut ci, &mut st);
    assert_eq!(mdl.is_bot_cell, vec![true, false], "only the sub-element on the bed may carry the bed condition");
}

// ---------------------------------------------------------------------------------------------
// Energy balance with deposition, convection + radiation, bed conduction and T-dependent PLA c(T).
//   E_stored(end) - sum(E_deposited) + sum(boundary losses) = 0
// Boundary losses are computed by the test from the pre-step state with its own formulas, using the
// exact piecewise-linear enthalpy H(T) = int c dT of the same table the solver interpolates.
// ---------------------------------------------------------------------------------------------
fn enthalpy_fn(interp: &Interpolator) -> impl Fn(f64) -> f64 {
    let xs = interp.xvalues.clone();
    let ys = interp.yvalues.clone();
    let mut cum = vec![0.0; xs.len()];
    for i in 1..xs.len() { cum[i] = cum[i - 1] + 0.5 * (ys[i] + ys[i - 1]) * (xs[i] - xs[i - 1]); }
    move |t: f64| {
        let i = (((t - xs[0]) / (xs[1] - xs[0])) as usize).min(xs.len() - 2);
        let f = t - xs[i];
        let slope = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]);
        cum[i] + ys[i] * f + 0.5 * slope * f * f
    }
}

struct EnergyResult { heat_in: f64, e_loss: f64, e_end: f64, err_over_loss: f64, err_over_heat_in: f64 }

fn deposition_energy_case(dt: f64, eps: f64, bed_k: f64, use_pla: bool) -> EnergyResult {
    let (w, ht) = (15.875e-3, 5.08e-3);
    let (rho, t_dep, t_env, t_bed) = (1271.1, 473.15, 313.15, 363.15);
    let (nx, nl) = (4usize, 3usize);
    let ncell = nx * nl;
    let sig = 5.67e-8;
    let mut nodes = vec![];
    for l in 0..nl { for i in 0..nx { for q in 0..8 {
        nodes.push([(i + (q & 1)) as f64 * w, ((q >> 1) & 1) as f64 * w, (l + (q >> 2)) as f64 * ht]);
    }}}
    let table = if use_pla { pla_table() } else { const_c(1929.0) };
    let hfun = enthalpy_fn(&table);
    let mut mdl = Model::new(nodes, ncell, if use_pla { pla_table() } else { const_c(1929.0) }, w, ht, 1000, 0.0, rho, eps, bed_k, t_bed, 1);
    let (mut tmp, mut ci, mut st) = (vec![], vec![], vec![]);
    let areas = [w * ht, w * w, ht, w * w / ht, w * w * ht];
    let vol = w * w * ht;
    let p = pool();
    let mut bw = sink();
    let mut gt = 0.0;
    let (mut e_dep, mut loss) = (0.0, 0.0);
    let advance = |mdl: &mut Model, tmp: &mut Vec<f64>, ci: &Vec<usize>, bw: &mut BufWriter<File>, gt: &mut f64, span: f64, loss: &mut f64| {
        let n = (span / dt).floor() as usize;
        let mut run = |mdl: &mut Model, tmp: &mut Vec<f64>, s: f64, loss: &mut f64| {
            for c in 0..mdl.templist.len() {
                let t = mdl.templist[c];
                let q_surf = mdl.h[c] * (t - t_env) + eps * sig * (t.powi(4) - t_env.powi(4));
                for f in 0..6 {
                    let area = if f < 4 { w * ht } else { w * w };
                    if f == 5 && mdl.is_bot_cell[c] { *loss += bed_k * (t - t_bed) * (w * w / ht) * s; }
                    else if !mdl.has_neighbor[c][f] { *loss += q_surf * area * s; }
                }
            }
            mdl.run_model(s, s, 0.0, t_env, bw, gt, &p, areas, tmp, 1, ci);
        };
        for _ in 0..n { run(mdl, tmp, dt, loss); }
        let rem = span - n as f64 * dt;
        if rem > 1e-12 { run(mdl, tmp, rem, loss); }
    };
    let mut idx = 0usize;
    for l in 0..nl {
        for i in 0..nx {
            let mut isn = [false; 6];
            let mut nb = [0isize; 6];
            if i > 0 { isn[1] = true; nb[1] = (idx - 1) as isize; }
            if l > 0 { isn[5] = true; nb[5] = (idx - nx) as isize; }
            let nn: [usize; 8] = core::array::from_fn(|q| 8 * idx + q);
            mdl.addCell(nn, 1929.0, rho, 0.59, 0.48, 0.35, t_dep, isn, nb, idx, &mut tmp, [1.0, 0.0], l, &mut ci, &mut st);
            e_dep += rho * vol * hfun(t_dep);
            idx += 1;
            let span = if i + 1 < nx { 0.2924 } else { 30.0 };
            advance(&mut mdl, &mut tmp, &ci, &mut bw, &mut gt, span, &mut loss);
        }
    }
    advance(&mut mdl, &mut tmp, &ci, &mut bw, &mut gt, 200.0, &mut loss);
    let e_end: f64 = mdl.templist.iter().map(|t| rho * vol * hfun(*t)).sum();
    let err = (e_end - e_dep + loss).abs();
    let heat_in = e_dep - ncell as f64 * rho * vol * hfun(t_env);
    EnergyResult { heat_in,e_loss: loss, e_end, err_over_loss: err / loss, err_over_heat_in: err / heat_in }
}

#[test]
fn energy_balance_with_deposition_bed_radiation_and_pla_cp() {
    println!("\nEnergy balance: 4 beads x 3 layers deposited at 0.29 s spacing, 30 s layer dwell, 200 s cooldown");
    println!("{:<34} {:>7} {:>13} {:>13} {:>11} {:>11}", "case", "dt [s]", "heat above T_env [J]", "lost [J]", "err/lost", "err/heat-in");
    let mut rows = vec![];
    let cases: [(&str, f64, f64, bool); 3] = [
        ("const c, conv only, adiabatic bed", 0.0, 0.0, false),
        ("const c, conv+rad, bed 0.17", 0.92, 0.17, false),
        ("PLA c(T), conv+rad, bed 0.17", 0.92, 0.17, true),
    ];
    for (name, eps, bed, pla) in cases {
        for dt in [1.0, 0.25, 0.05] {
            let r = deposition_energy_case(dt, eps, bed, pla);
            println!("{:<34} {:>7} {:>13.3} {:>13.3} {:>11.2e} {:>11.2e}", name, dt, r.heat_in, r.e_loss, r.err_over_loss, r.err_over_heat_in);
            rows.push(format!("{},{},{},{},{},{}", name, dt, r.e_end, r.e_loss, r.err_over_loss, r.err_over_heat_in));
            assert!(r.err_over_heat_in < 0.01, "{} dt={}: energy balance error {:.3}% exceeds 1%", name, dt, 100.0 * r.err_over_heat_in);
            if !pla { assert!(r.err_over_heat_in < 1e-9, "constant c must conserve to round-off"); }
        }
    }
    if let Some(dir) = out_dir() {
        let mut f = File::create(dir.join("energy_deposition.csv")).unwrap();
        writeln!(f, "case,dt_s,e_end_J,e_loss_J,err_over_loss,err_over_heat_in").unwrap();
        for r in rows { writeln!(f, "{}", r).unwrap(); }
    }
}
