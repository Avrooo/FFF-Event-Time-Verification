use rayon::ThreadPool;
use std::fs::File;
use std::io::{BufReader, BufRead};

/// A linear interpolator that interpolates between evenly spaced x-values and gets an equivalent
/// value based on x-value. The x-values need to be uniformly spaced. X-values and Y-values need to
/// have equal length
pub struct Interpolator{
    pub (crate) xvalues:Vec<f64>,
    pub (crate) yvalues:Vec<f64>,
    pub (crate) xstep:f64,
    pub (crate) maxthreads: usize,
}
/// Implementation of methods for interpolator
impl Interpolator{
    pub fn new(xvalues:Vec<f64>, yvalues:Vec<f64>, xstep:f64, maxthreads: usize)->Interpolator{
        return Interpolator {
            xvalues,
            yvalues,
            xstep,
            maxthreads
        }
    }

    pub fn interpolate(&self, input_values: & Vec<f64>, output_values: &mut Vec<f64>, pool: &ThreadPool){
        self.interpolate_inner(input_values.as_slice(), output_values.as_mut_slice(), 1,self.maxthreads, pool);
    }

    pub fn interpolate_inner(&self, input_values: &[f64], output_values:&mut[f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        if numthreads <= maxthreads{
            let split_pt = input_values.len()/2;
            let(iv1, iv2) = input_values.split_at(split_pt);
            let(ov1, ov2) = output_values.split_at_mut(split_pt);
            pool.install(||rayon::join(
                ||self.interpolate_inner(iv1,ov1,numthreads*2, maxthreads, pool),
                ||self.interpolate_inner(iv2,ov2,numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..input_values.len(){
                let current_value = input_values[i];
                unsafe {
                    let location = ((current_value - self.xvalues.get_unchecked(0)) / self.xstep) as usize;
                    if location >= self.xvalues.len() {panic!("The temperature is outside the range of input interpolation values. Input table does not have the ranges suitable for simulation or the simulation is unstable due to large time step")}
                    //println!("temperature {} index {}", current_value, i);
                    let slope = (self.yvalues.get_unchecked(location + 1) - self.yvalues.get_unchecked(location)) / self.xstep;
                    output_values[i] = self.yvalues.get_unchecked(location) + slope * (input_values.get_unchecked(i) - self.xvalues.get_unchecked(location));
                }
            }
        }
    }


    pub fn interpolate1(&self, input_values: & Vec<f64>,  output_values: &mut Vec<f64>, pool:&ThreadPool){
        let input_values_size = input_values.len();
        let (ipval1, ipval2) = input_values.split_at(input_values_size/2);
        let (ipval11, ipval12) = ipval1.split_at(input_values_size/4);
        let (ipval21,ipval22) = ipval2.split_at(input_values_size/4);
        let (ipval111,ipval112) = ipval11.split_at(input_values_size/8);
        let (ipval121,ipval122) = ipval12.split_at(input_values_size/8);
        let (ipval211,ipval212) = ipval21.split_at(input_values_size/8);
        let (ipval221,ipval222) = ipval22.split_at(input_values_size/8);
        let (opval1, opval2) = output_values.split_at_mut(input_values_size/2);
        let (opval11, opval12) = opval1.split_at_mut(input_values_size/4);
        let (opval21,opval22) = opval2.split_at_mut(input_values_size/4);
        let (opval111,opval112) = opval11.split_at_mut(input_values_size/8);
        let (opval121,opval122) = opval12.split_at_mut(input_values_size/8);
        let (opval211,opval212) = opval21.split_at_mut(input_values_size/8);
        let (opval221,opval222) = opval22.split_at_mut(input_values_size/8);
        pool.install(|| rayon::join(
            ||rayon::join(
                ||rayon::join(
                    ||self.interpolate_split(ipval111,opval111),
                    ||self.interpolate_split(ipval112,opval112)),
                ||rayon::join(
                    ||self.interpolate_split(ipval121,opval121),
                    ||self.interpolate_split(ipval122,opval122)
                )
            ),||rayon::join(
                ||rayon::join(
                    ||self.interpolate_split(ipval211,opval211),
                    ||self.interpolate_split(ipval212,opval212)),
                ||rayon::join(
                    ||self.interpolate_split(ipval221,opval221),
                    ||self.interpolate_split(ipval222,opval222)
                ))));


    }



    pub fn interpolate_slice(&self, input_values: & [f64],  output_values: &mut [f64], pool:&ThreadPool){
        let input_values_size = input_values.len();
        let (ipval1, ipval2) = input_values.split_at(input_values_size/2);
        let (ipval11, ipval12) = ipval1.split_at(input_values_size/4);
        let (ipval21,ipval22) = ipval2.split_at(input_values_size/4);
        let (ipval111,ipval112) = ipval11.split_at(input_values_size/8);
        let (ipval121,ipval122) = ipval12.split_at(input_values_size/8);
        let (ipval211,ipval212) = ipval21.split_at(input_values_size/8);
        let (ipval221,ipval222) = ipval22.split_at(input_values_size/8);
        let (opval1, opval2) = output_values.split_at_mut(input_values_size/2);
        let (opval11, opval12) = opval1.split_at_mut(input_values_size/4);
        let (opval21,opval22) = opval2.split_at_mut(input_values_size/4);
        let (opval111,opval112) = opval11.split_at_mut(input_values_size/8);
        let (opval121,opval122) = opval12.split_at_mut(input_values_size/8);
        let (opval211,opval212) = opval21.split_at_mut(input_values_size/8);
        let (opval221,opval222) = opval22.split_at_mut(input_values_size/8);
        pool.install(|| rayon::join(
            ||rayon::join(
                ||rayon::join(
                    ||self.interpolate_split(ipval111,opval111),
                    ||self.interpolate_split(ipval112,opval112)),
                ||rayon::join(
                    ||self.interpolate_split(ipval121,opval121),
                    ||self.interpolate_split(ipval122,opval122)
                )
            ),||rayon::join(
                ||rayon::join(
                    ||self.interpolate_split(ipval211,opval211),
                    ||self.interpolate_split(ipval212,opval212)),
                ||rayon::join(
                    ||self.interpolate_split(ipval221,opval221),
                    ||self.interpolate_split(ipval222,opval222)
                ))));


    }

    fn interpolate_split(&self, input_values:&[f64],output_values: &mut [f64]){
        for i in 0..input_values.len(){
            let current_value = input_values[i];
            let location = ((current_value-self.xvalues[0]) / self.xstep) as usize;
            //println!("temperature {} index {}", current_value, i);
            let slope = (self.yvalues[location+1]-self.yvalues[location]) / self.xstep;
            output_values[i] = self.yvalues[location] + slope * (input_values[i] - self.xvalues[location]);
        }
    }

    pub(crate) fn read_data_from_file(filename: &str, xstep:f64, maxthreads: usize) ->Interpolator{
        let mut temps =  Vec::new();
        let mut sp_ht_caps = Vec::new();
        let spheat_data_file = File::open(filename).expect("can't find file");
        let spht_buf = BufReader::with_capacity(10000,spheat_data_file);
        for line in spht_buf.lines(){
            let line2 = line.expect("cant read the line").clone();
            let mut dat = line2.split(",");
            let sphtcap = dat.next().expect("row does not have first column");
            let temp = dat.next().expect("row does not have second column");
            temps.push(temp.parse::<f64>().expect("can't parse temp to f64"));
            sp_ht_caps.push(sphtcap.parse::<f64>().expect("can't parse sp heat cap to f64"));
        }
        // R8.4 ablation: legacy_cp_step restores the published behaviour (declared step trusted)
        let step = if crate::legacy::cp_step() { xstep } else { Interpolator::validated_step(filename, &temps, xstep) };
        return Interpolator::new(temps,sp_ht_caps,step, maxthreads);
    }

    /// D5 fix. The lookup index is (T - T0)/step, so a declared step that differs from the real spacing of the
    /// table silently reads a different (stretched or compressed) curve - the published example declares 5 K for a
    /// 1 K table and therefore only ever used 288-333 K of it (c effectively constant). The spacing is now taken
    /// from the table itself: it must be uniform and ascending (otherwise the run stops), and a declared step that
    /// disagrees is reported and overridden.
    pub(crate) fn validated_step(filename: &str, temps: &[f64], declared: f64) -> f64 {
        assert!(temps.len() >= 2, "table {} needs at least two rows", filename);
        let step = (temps[temps.len() - 1] - temps[0]) / (temps.len() - 1) as f64;
        assert!(step > 0.0, "table {}: temperatures must be ascending", filename);
        for w in temps.windows(2) {
            assert!(((w[1] - w[0]) - step).abs() <= 1e-6 * step.max(1.0),
                    "table {}: temperature spacing is not uniform ({} then {}); the lookup needs equally spaced rows", filename, step, w[1] - w[0]);
        }
        if (declared - step).abs() > 1e-6 * step {
            eprintln!("WARNING: {} declares a temperature step of {} but its rows are spaced {} apart; using {} (D5 correction)", filename, declared, step, step);
        }
        step
    }

}



#[cfg(test)]
mod tests {
    use super::*;
    fn table(step: f64, n: usize) -> Interpolator {
        // c(T) = 1000 + 10 (T-300) + 0.5 (T-300)^2 is piecewise-linearised on the grid
        let t: Vec<f64> = (0..n).map(|i| 300.0 + step * i as f64).collect();
        let c: Vec<f64> = t.iter().map(|x| 1000.0 + 10.0 * (x - 300.0) + 0.5 * (x - 300.0).powi(2)).collect();
        Interpolator::new(t, c, step, 1)
    }
    fn at(ip: &Interpolator, x: f64) -> f64 {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let mut out = vec![0.0]; ip.interpolate_inner(&[x], &mut out, 2, 1, &pool); out[0]
    }
    #[test]
    fn validated_step_is_taken_from_the_table() {
        let t: Vec<f64> = (0..50).map(|i| 288.15 + i as f64).collect();
        assert!((Interpolator::validated_step("x", &t, 5.0) - 1.0).abs() < 1e-12, "declared 5 must be overridden by the real 1 K spacing");
        assert!((Interpolator::validated_step("x", &t, 1.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    #[should_panic(expected = "not uniform")]
    fn non_uniform_table_is_rejected() {
        Interpolator::validated_step("x", &[300.0, 301.0, 303.0, 304.0], 1.0);
    }
    #[test]
    fn lookup_reproduces_the_table_at_nodes_and_is_linear_between() {
        let ip = table(1.0, 100);
        assert!((at(&ip, 310.0) - (1000.0 + 100.0 + 50.0)).abs() < 1e-9);
        let mid = 0.5 * ((1000.0 + 10.0 * 10.0 + 0.5 * 100.0) + (1000.0 + 10.0 * 11.0 + 0.5 * 121.0));
        assert!((at(&ip, 310.5) - mid).abs() < 1e-9);
    }
    #[test]
    fn wrong_declared_step_reads_the_wrong_curve_negative_control() {
        // the OLD behaviour: step declared 5 for a 1 K table -> value at 330 K is read from row (330-300)/5=6 (i.e. 306 K)
        let t: Vec<f64> = (0..100).map(|i| 300.0 + i as f64).collect();
        let c: Vec<f64> = t.iter().map(|x| 1000.0 + 10.0 * (x - 300.0)).collect();
        let wrong = Interpolator::new(t.clone(), c.clone(), 5.0, 1);
        let right = Interpolator::new(t, c, 1.0, 1);
        assert!((at(&right, 330.0) - 1300.0).abs() < 1e-9);
        assert!((at(&wrong, 330.0) - 1300.0).abs() > 50.0, "declared-step mismatch must visibly corrupt the lookup");
    }
}
