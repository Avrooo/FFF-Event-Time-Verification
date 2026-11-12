#![allow(warnings)]
mod bed_contact;
mod convection;
mod fan;
mod legacy;
mod gcode_reader;
mod interpolator;
mod model;
mod model_generator;
mod model_iso_td_shc_td_con;
mod model_orthotropic_td_shc;
mod model_updater;
mod primitives;
mod simpleModel;
mod model_iso_td_shc_td;
mod model_orthotropic_td_shc_variable_h;
//mod model_with_structural;

use num_cpus;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::borrow::BorrowMut;
use std::cmp::{max, Ordering};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter, Read};
use std::iter::FromIterator;
use std::ops::Index;
use std::thread::current;
use std::{fmt, fs};

use crate::bed_contact::BedContact;
use crate::convection::Convection;
use crate::fan::FanModel;
use crate::legacy::Legacy;
use crate::gcode_reader::GCodeReader;
use crate::interpolator::*;
use crate::model::Model;
use crate::model_generator::ModelGenerator;
use crate::model_updater::ModelUpdater;
use crate::primitives::{Node, Point};
use std::str::FromStr;

use std::time::Instant;

fn read_next_item<'a>(inputfilebufreader: &'a mut BufReader<File>, line: &'a mut String) -> Option<&'a str>{
    line.clear();
    inputfilebufreader.read_line(line);
    let mut linesplit = line.split(" ");
    let v1 = linesplit.next();
    //println!("v1 {:?}", v1?);
    //linesplit.next();
    let value =  linesplit.next();
    // let v2 = value.clone();
    //println!("{:?}", v2.unwrap());
    return value
}

/// Optional trailing `key: value` lines (blank lines skipped, anything after the value ignored).
fn read_optional_keys(inputfilebufreader: &mut BufReader<File>) -> HashMap<String, String> {
    let mut keys = HashMap::new();
    let mut line = String::new();
    loop {
        line.clear();
        if inputfilebufreader.read_line(&mut line).unwrap_or(0) == 0 { break; }
        let mut parts = line.trim().splitn(2, ':');
        if let (Some(k), Some(v)) = (parts.next(), parts.next()) {
            if let Some(val) = v.trim().split_whitespace().next() { keys.insert(k.trim().to_string(), val.to_string()); }
        }
    }
    keys
}

/// Number of nodes and far-boundary coordinate of a mesh axis starting at `lo` with elements of exactly `h`
/// that covers at least [lo, hi] (a tolerance avoids adding a whole extra element for round-off).
pub(crate) fn mesh_axis(lo: f64, hi: f64, h: f64) -> (usize, f64) {
    let cells = ((hi - lo) / h - 1e-6).ceil().max(1.0) as usize;
    (cells + 1, lo + cells as f64 * h)
}

fn main() {
    let inputfile = File::open("inputfiles/Input_file.txt").expect("InputFile not found");
    let mut inputfilebufreader = BufReader::with_capacity(10000, inputfile);
    let mut line = String::with_capacity(100);
    inputfilebufreader.read_line(&mut line).expect("cant read line 1");
    line.clear();
    inputfilebufreader.read_line(&mut line).expect("cant read line 2");
    let mut linesplit = line.split(" "); linesplit.next();
    let mut num_cpus = usize::from_str(linesplit.next().expect("cant read number of cpus")).expect("cant parse number of cpus provided to integer");
    if num_cpus == 0 {
        num_cpus = num_cpus::get_physical();
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_cpus+2)
        .build()
        .expect("can't create threads");
    println!("using {} cpu threads", num_cpus);
    let maxthreads = num_cpus*2;
    line.clear();
    inputfilebufreader.read_line(&mut line).expect("cant read line 3");
    let mut linesplit = line.split(" "); linesplit.next();
    let mut gcodefile = linesplit.next().expect("cant read the input gcode filename");
    println!("reading gcode file |{}|", gcodefile);
    let gcr = GCodeReader::new(gcodefile);
    line.clear();
    inputfilebufreader.read_line(&mut line).expect("cant read line 3");
    let mut linesplit = line.split(" "); linesplit.next();
    let divs_per_bead= usize::from_str(linesplit.next().expect("cant read divisions per beadwidth")).expect("cant read divisions per beadwidth as positive integer");
    line.clear();
    inputfilebufreader.read_line(&mut line).expect("cant read line 3");
    let mut linesplit = line.split(" "); linesplit.next();
    let divs_per_bead_z =  usize::from_str(linesplit.next().expect("cant read divisions per beadheight")).expect("cant read divisions per beadheight as positive integer");
    let beadwidth = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read beadwidth")).expect("cant parse beadwidth to float");
    let beadheight = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read beadheight")).expect("cant parse beadheight to float");
    
    let model_type = usize::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read model type")).expect("cant parse model type to integer");
    if model_type == 0 {
        
    }
    else if model_type == 1 {
        let sp_ht_cap_filename = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read specific heat capacity filename"));
        let sp_heat_cap_step = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read specific heat capacity file temperature spacing")).expect("cant parse specific heat capacity file temperature spacing to float");
        let conductivity_filename = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read conductivity filename"));
        let conductivity_step = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read conductivity file temperature spacing")).expect("cant parse conductivity file temperature spacing to float");
        let density = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read density")).expect("cant parse density to float");
        let h = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read convective heat film transfer coefficient")).expect("cant parse convective heat transfer film coefficient to float");
        let temp_bed = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read bed temperature")).expect("cant parse bed temperature to float");
        let bed_k = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read bed interface conductivity")).expect("cant parse bed interface conductivity to float");
        let ambient_temp = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read ambient temperature")).expect("cant parse ambient temperature to float");

        let extrusion_temperature = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read extrusion temperature")).expect("cant parse extrusion temperature to float");
        let emissivity = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read extrusion temperature")).expect("cant parse extrusion temperature to float");
        let time_step = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read bed time step")).expect("cant parse time step to float");
        let cooldown_period = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read cooldown period")).expect("cant parse cooldown period to float");

        let node_file_output_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read node file output name"));

        let element_file_output_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read element file output name"));
        let activation_times_file_output_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read activation times file output name"));
        let abaqus_input_file_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read abaqus input file output name"));
        let elemental_temperature_output_file_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read elemental temperature output filename"));
        let min_temp_change_store = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read min temp change to store file output name")).expect("cant parse min temp change to store to float");
        let nodal_temperature_output_file_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read nodal temperature output filename"));
        let turn_off_layers_at = usize::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read turn off layer at")).expect("cant parse turn off layers at as integer");
        // optional trailing keys (absent in the published input files => defaults)
        let optional_keys = read_optional_keys(&mut inputfilebufreader);
        let bed_contact = BedContact::from_input(
            optional_keys.get("bed_contact_model").and_then(|s| usize::from_str(s).ok()),
            optional_keys.get("bed_contact_resistance").and_then(|s| f64::from_str(s).ok()));
        let convection = Convection::from_input(optional_keys.get("convection_model").and_then(|s| usize::from_str(s).ok()), h);
        println!("bed contact model: {:?}", bed_contact);
        println!("convection model: {:?}", convection);

        // all inputs done. calculations start here
        let element_width = beadwidth / divs_per_bead as f64;

        let element_height = beadheight / divs_per_bead_z as f64;
        // D8 fix: the mesh spacing must EQUAL the element size the solver uses for volumes, areas and conduction
        // distances. The domain is anchored at its minimum corner, covers the whole toolpath, and its far boundary
        // is moved outward to the next whole number of elements (previously the number of nodes was rounded and the
        // domain was divided into it, so the generated spacing differed from `element_width` by up to ~0.5 %).
        let (xdiv, xmax_dom) = mesh_axis(gcr.xmin - beadwidth / 2.0, gcr.xmax + beadwidth / 2.0, element_width);
        let (ydiv, ymax_dom) = mesh_axis(gcr.ymin - beadwidth / 2.0, gcr.ymax + beadwidth / 2.0, element_width);
        let (zdiv, zmax_dom) = mesh_axis(gcr.zmin - beadheight, gcr.zmax, element_height);

        let zmin = gcr.zmin.clone();
        println!(
            "xmin {} xmax {} ymin {} ymax {} zmin {} zmax {} xdiv {} ydiv {} zdiv {}",
            gcr.xmin, gcr.xmax, gcr.ymin, gcr.ymax, gcr.zmin, gcr.zmax, xdiv, ydiv, zdiv
        );
        let init_temp = extrusion_temperature;
        let tic = Instant::now();
        let mut m = ModelGenerator::new(
            gcr.xmin - beadwidth / 2.0,
            xmax_dom,
            gcr.ymin - beadwidth / 2.0,
            ymax_dom,
            gcr.zmin - beadheight,
            zmax_dom,
            xdiv,
            ydiv,
            zdiv,
            init_temp,
            element_width,
            divs_per_bead_z,
        );
        println!("time taken to generate model {:?}", tic.elapsed());

        let mut activation_times_all: Vec<(f64,usize,[f64;2],usize)> =  m.generate_activation_times_all_layers(
            gcr.segment,
            gcr.is_extrusion_on,
            gcr.speed,
            beadwidth,
            beadheight,
            &pool,
            maxthreads,
        );//Vec::with_capacity(10000000);




        let mut activation_times_file = File::create(activation_times_file_output_name.as_str()).expect("can't create file");
        // let mut activation_times_file = File::create("/mnt/c/rustFiles/activation_times_store_file.csv").expect("can't create file");
        let mut write_buffer = BufWriter::with_capacity(100000,activation_times_file);
        for i in activation_times_all.clone(){
            writeln!(write_buffer,"{},{},{},{},{}",i.0,i.1,i.2[0],i.2[1],i.3);
        }
        println!("writing node and element files");
        let mut activated_elements = Vec::with_capacity(activation_times_all.len());
        for i in 0..activation_times_all.len(){
            activated_elements.push(activation_times_all[i].1);
        }
        //println!("{}",node_file_output_name);
        let (nodeveclen, elemveclen, nodesnew, elementsupdated,nodenum_old_to_new) = m.write_nodes_and_elements_to_file(node_file_output_name.as_str(),element_file_output_name.as_str(), activated_elements.clone());
        //println!("{}", abaqus_input_file_name);
        //m.write_abaqus_input_mesh_file(abaqus_input_file_name.as_str(),activated_elements.clone());

        let sp_heat_cap_interpolation_table = Interpolator::read_data_from_file(sp_ht_cap_filename.as_str(),sp_heat_cap_step, maxthreads);
        let conductivity_interpolation_table = Interpolator::read_data_from_file(conductivity_filename.as_str(),conductivity_step, maxthreads);

        println!("specific heat capacity and conductivity data files read");
        let mut mu = ModelUpdater::new(activation_times_all, m, nodenum_old_to_new);

        //default values
        let input_data = ([0.205, 0.205, 0.205], density, 1500.0, h, [init_temp, ambient_temp]);

        //convert old nodes to new nodes

        let mut mdl = model_iso_td_shc_td::Model::new(
            nodesnew.clone(),
            mu.activation_times.len(),
            sp_heat_cap_interpolation_table,
            element_width,
            element_height,
            turn_off_layers_at,
            zmin-beadheight,input_data.1,emissivity,bed_k,temp_bed,maxthreads);
        mdl.set_bed_contact(bed_contact);
        assert!(convection == Convection::Constant(h), "model_type 1 always uses the constant input film coefficient; convection_model 1 applies to model_type 2 only");


        let areas_and_dists = [
            element_width * element_height,
            element_width * element_width,
            element_height * element_width / element_width,
            element_width * element_width / element_height,
            element_width * element_width * element_height,
        ];
        let mut tmpfile = File::create("tempfile.csv").expect("cant create temporary scratch file");
        let mut bw= BufWriter::with_capacity(10000,tmpfile);

        let mut datastorer = Vec::with_capacity(activated_elements.len());
        for i in 0..activated_elements.len(){
            datastorer.push(Vec::with_capacity(1000));
        }

        let mut datastorernode = Vec::with_capacity(nodeveclen);
        // for i in 0..nodeveclen{
        //   datastorernode.push(Vec::with_capacity(1000));
        // }

        //create node to element map
        let mut nd_to_elem = Vec::with_capacity(nodeveclen);
        for i in 0..nodeveclen{
            nd_to_elem.push(Vec::with_capacity(6));
        }

        for i in 0..elementsupdated.len(){
            for j in 0..8{
                nd_to_elem[elementsupdated[i][j]].push(i);
            }
        }

        println!("starting simulation timestep {}", time_step);
        let tic = Instant::now();
        // assuming inputfile has model type 1
        model_iso_td_shc_td::Model::update_model(
            &mut mu,
            &mut mdl,
            &mut bw,
            areas_and_dists,
            time_step,
            &conductivity_interpolation_table,
            maxthreads,
            input_data,
            min_temp_change_store,
            &mut datastorer,
            &mut datastorernode,
            nd_to_elem.as_slice(),
            cooldown_period,
            &pool,
        );
        println!("simulation ended in {:?} . start writing elemental temperature output file", tic.elapsed());
        let tic = Instant::now();
        let mut outfile = File::create(elemental_temperature_output_file_name).expect("cant create elemental temperature output file");
        let mut outfilebw = BufWriter::with_capacity(100000,outfile);
        for i in 0..datastorer.len(){
            for j in 0..datastorer[i].len(){
                write!(outfilebw,"{},",datastorer[i][j]);
            }
            write!(outfilebw,"\n");
        }
        println!("writing elemental temperature output file took {:?}", tic.elapsed());

        let tic = Instant::now();
        let mut outfilenode = File::create(nodal_temperature_output_file_name).expect("cant create nodal temperature output file");
        let mut outfilenodebw = BufWriter::with_capacity(10000, outfilenode);
        //println!("{}",datastorernode[0][0]);
        for i in 0..datastorernode.len(){
            for j in 0..datastorernode[i].len(){
                write!(outfilenodebw,"{},",datastorernode[i][j]);
            }
            write!(outfilenodebw,"\n");
        }
        println!("writing nodal temperature output file took {:?}", tic.elapsed());
    }
    
    else if model_type == 2 {
        //println!("entered model 2");
        let sp_ht_cap_filename = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read specific heat capacity filename"));
        let sp_heat_cap_step = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read specific heat capacity file temperature spacing")).expect("cant parse specific heat capacity file temperature spacing to float");
        let kx = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read kx")).expect("cant parse kx to float");
        let ky = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read ky")).expect("cant parse ky to float");
        let kz = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read kz")).expect("cant parse kz to float");     
        
        println!("kx ky kz {},{},{}", kx, ky, kz);
        
        let density = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read density")).expect("cant parse density to float");
        let h = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read convective heat film transfer coefficient")).expect("cant parse convective heat transfer film coefficient to float");
        let temp_bed = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read bed temperature")).expect("cant parse bed temperature to float");
        let bed_k = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read bed interface conductivity")).expect("cant parse bed interface conductivity to float");
        let ambient_temp = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read ambient temperature")).expect("cant parse ambient temperature to float");

        let extrusion_temperature = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read extrusion temperature")).expect("cant parse extrusion temperature to float");
        let emissivity = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read extrusion temperature")).expect("cant parse extrusion temperature to float");
        let time_step = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read bed time step")).expect("cant parse time step to float");
        let cooldown_period = f64::from_str(read_next_item(&mut inputfilebufreader,&mut line).expect("cant read cooldown period")).expect("cant parse cooldown period to float");

        let node_file_output_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read node file output name"));

        let element_file_output_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read element file output name"));
        let activation_times_file_output_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read activation times file output name"));
        let abaqus_input_file_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read abaqus input file output name"));
        let elemental_temperature_output_file_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read elemental temperature output filename"));
        let min_temp_change_store = f64::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read min temp change to store file output name")).expect("cant parse min temp change to store to float");
        let nodal_temperature_output_file_name = String::from(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read nodal temperature output filename"));
        let turn_off_layers_at = usize::from_str(read_next_item(&mut inputfilebufreader, &mut line).expect("cant read turn off layer at")).expect("cant parse turn off layers at as integer");
        // optional trailing keys (absent in the published input files => defaults)
        let optional_keys = read_optional_keys(&mut inputfilebufreader);
        let bed_contact = BedContact::from_input(
            optional_keys.get("bed_contact_model").and_then(|s| usize::from_str(s).ok()),
            optional_keys.get("bed_contact_resistance").and_then(|s| f64::from_str(s).ok()));
        let convection = Convection::from_input(optional_keys.get("convection_model").and_then(|s| usize::from_str(s).ok()), h);
        println!("bed contact model: {:?}", bed_contact);
        println!("convection model: {:?}", convection);
        // R8.2 optional fan model (absent => off => identical to the frozen baseline)
        let fan = match optional_keys.get("fan_model").and_then(|s| usize::from_str(s).ok()).unwrap_or(0) {
            0 => FanModel::off(),
            1 => {
                let path = optional_keys.get("fan_schedule_file").expect("fan_model 1 needs fan_schedule_file");
                let hf = f64::from_str(optional_keys.get("fan_h_forced_full").expect("fan_model 1 needs fan_h_forced_full")).expect("fan_h_forced_full");
                let ex = optional_keys.get("fan_speed_exponent").and_then(|s| f64::from_str(s).ok()).unwrap_or(0.6);
                let text = std::fs::read_to_string(path).expect("cannot read fan schedule file");
                FanModel::from_schedule_text(&text, hf, ex)
            }
            o => panic!("unknown fan_model {} (0 = off, 1 = schedule file)", o),
        };
        println!("fan model: {:?}", fan);
        let legacy = Legacy::from_keys(|k| optional_keys.get(k).and_then(|s| usize::from_str(s).ok()));
        legacy::LEGACY_CP_STEP.store(legacy.cp_step, std::sync::atomic::Ordering::Relaxed);
        println!("legacy switches: {:?}", legacy);

        // all inputs done. calculations start here
        let element_width = beadwidth / divs_per_bead as f64;

        let element_height = beadheight / divs_per_bead_z as f64;
        // D8 fix: the mesh spacing must EQUAL the element size the solver uses for volumes, areas and conduction
        // distances. The domain is anchored at its minimum corner, covers the whole toolpath, and its far boundary
        // is moved outward to the next whole number of elements (previously the number of nodes was rounded and the
        // domain was divided into it, so the generated spacing differed from `element_width` by up to ~0.5 %).
        let (mut xdiv, mut xmax_dom) = mesh_axis(gcr.xmin - beadwidth / 2.0, gcr.xmax + beadwidth / 2.0, element_width);
        let (mut ydiv, mut ymax_dom) = mesh_axis(gcr.ymin - beadwidth / 2.0, gcr.ymax + beadwidth / 2.0, element_width);
        let (mut zdiv, mut zmax_dom) = mesh_axis(gcr.zmin - beadheight, gcr.zmax, element_height);
        if legacy.mesh {
            // D8 (original): node counts rounded, domain divided into them (spacing differs from the element size)
            xdiv = ((gcr.xmax - gcr.xmin + beadwidth) / element_width + 1.0).round() as usize; xmax_dom = gcr.xmax + beadwidth / 2.0;
            ydiv = ((gcr.ymax - gcr.ymin + beadwidth) / element_width + 1.0).round() as usize; ymax_dom = gcr.ymax + beadwidth / 2.0;
            zdiv = ((gcr.zmax - gcr.zmin + beadheight) / element_height + 1.0).round() as usize; zmax_dom = gcr.zmax;
        }

        let zmin = gcr.zmin.clone();
        println!(
            "xmin {} xmax {} ymin {} ymax {} zmin {} zmax {} xdiv {} ydiv {} zdiv {}",
            gcr.xmin, gcr.xmax, gcr.ymin, gcr.ymax, gcr.zmin, gcr.zmax, xdiv, ydiv, zdiv
        );
        let init_temp = extrusion_temperature;
        let tic = Instant::now();
        let mut m = ModelGenerator::new(
            gcr.xmin - beadwidth / 2.0,
            xmax_dom,
            gcr.ymin - beadwidth / 2.0,
            ymax_dom,
            gcr.zmin - beadheight,
            zmax_dom,
            xdiv,
            ydiv,
            zdiv,
            init_temp,
            element_width,
            divs_per_bead_z,
        );
        println!("time taken to generate model {:?}", tic.elapsed());

        let mut activation_times_all: Vec<(f64,usize,[f64;2],usize)> =  m.generate_activation_times_all_layers(
            gcr.segment,
            gcr.is_extrusion_on,
            gcr.speed,
            beadwidth,
            beadheight,
            &pool,
            maxthreads,
        );//Vec::with_capacity(10000000);




        let mut activation_times_file = File::create(activation_times_file_output_name.as_str()).expect("can't create file");
        // let mut activation_times_file = File::create("/mnt/c/rustFiles/activation_times_store_file.csv").expect("can't create file");
        let mut write_buffer = BufWriter::with_capacity(100000,activation_times_file);
        for i in activation_times_all.clone(){
            writeln!(write_buffer,"{},{},{},{},{}",i.0,i.1,i.2[0],i.2[1],i.3);
        }
        println!("writing node and element files");
        let mut activated_elements = Vec::with_capacity(activation_times_all.len());
        for i in 0..activation_times_all.len(){
            activated_elements.push(activation_times_all[i].1);
        }
        //println!("{}",node_file_output_name);
        let (nodeveclen, elemveclen, nodesnew, elementsupdated,nodenum_old_to_new) = m.write_nodes_and_elements_to_file(node_file_output_name.as_str(),element_file_output_name.as_str(), activated_elements.clone());
        //println!("{}", abaqus_input_file_name);
        //m.write_abaqus_input_mesh_file(abaqus_input_file_name.as_str(),activated_elements.clone());

        let sp_heat_cap_interpolation_table = Interpolator::read_data_from_file(sp_ht_cap_filename.as_str(),sp_heat_cap_step, maxthreads);
        //let conductivity_interpolation_table = Interpolator::read_data_from_file(conductivity_filename.as_str(),conductivity_step, maxthreads);

        println!("specific heat capacity and conductivity data files read");
        let mut mu = ModelUpdater::new(activation_times_all, m, nodenum_old_to_new);

        //default values
        let input_data = ([kx, ky, kz], density, 1500.0, h, [init_temp, ambient_temp]);

        //convert old nodes to new nodes

        let mut mdl = model_orthotropic_td_shc_variable_h::Model::new(
            nodesnew.clone(),
            mu.activation_times.len(),
            sp_heat_cap_interpolation_table,
            element_width,
            element_height,
            turn_off_layers_at,
            zmin-beadheight,input_data.1,emissivity,bed_k,temp_bed,maxthreads);
        mdl.set_bed_contact(bed_contact);
        mdl.set_convection(convection);
        mdl.set_fan(fan);
        mdl.set_legacy(legacy);


        let areas_and_dists = [
            element_width * element_height,
            element_width * element_width,
            element_height * element_width / element_width,
            element_width * element_width / element_height,
            element_width * element_width * element_height,
        ];
        let mut tmpfile = File::create("tempfile.csv").expect("cant create temporary scratch file");
        let mut bw= BufWriter::with_capacity(10000,tmpfile);

        let mut datastorer = Vec::with_capacity(activated_elements.len());
        for i in 0..activated_elements.len(){
            datastorer.push(Vec::with_capacity(1000));
        }

        let mut datastorernode = Vec::with_capacity(nodeveclen);
        // for i in 0..nodeveclen{
        //   datastorernode.push(Vec::with_capacity(1000));
        // }

        //create node to element map
        let mut nd_to_elem = Vec::with_capacity(nodeveclen);
        for i in 0..nodeveclen{
            nd_to_elem.push(Vec::with_capacity(6));
        }

        for i in 0..elementsupdated.len(){
            for j in 0..8{
                nd_to_elem[elementsupdated[i][j]].push(i);
            }
        }

        println!("starting simulation timestep {}", time_step);
        let tic = Instant::now();
        // assuming inputfile has model type 1
        model_orthotropic_td_shc_variable_h::Model::update_model(
            &mut mu,
            &mut mdl,
            &mut bw,
            areas_and_dists,
            time_step,
            maxthreads,
            input_data,
            min_temp_change_store,
            &mut datastorer,
            &mut datastorernode,
            nd_to_elem.as_slice(),
            cooldown_period,
            &pool,
        );
        println!("simulation ended in {:?} . start writing elemental temperature output file", tic.elapsed());
        let tic = Instant::now();
        let mut outfile = File::create(elemental_temperature_output_file_name).expect("cant create elemental temperature output file");
        let mut outfilebw = BufWriter::with_capacity(100000,outfile);
        for i in 0..datastorer.len(){
            for j in 0..datastorer[i].len(){
                write!(outfilebw,"{},",datastorer[i][j]);
            }
            write!(outfilebw,"\n");
        }
        println!("writing elemental temperature output file took {:?}", tic.elapsed());

        let tic = Instant::now();
        let mut outfilenode = File::create(nodal_temperature_output_file_name).expect("cant create nodal temperature output file");
        let mut outfilenodebw = BufWriter::with_capacity(10000, outfilenode);
        //println!("{}",datastorernode[0][0]);
        for i in 0..datastorernode.len(){
            for j in 0..datastorernode[i].len(){
                write!(outfilenodebw,"{},",datastorernode[i][j]);
            }
            write!(outfilenodebw,"\n");
        }
        println!("writing nodal temperature output file took {:?}", tic.elapsed());
    }
    
}

#[cfg(test)]
mod mesh_axis_tests {
    use super::mesh_axis;
    #[test]
    fn spacing_is_exactly_the_element_size_and_covers_the_range() {
        for (lo, hi, h) in [(0.0, 0.135875, 0.015875), (0.0, 0.12, 0.015875), (0.0, 0.015875, 0.015875), (0.01, 0.0763, 0.0079375), (0.0, 0.0254, 0.00127)] {
            let (n, top) = mesh_axis(lo, hi, h);
            let spacing = (top - lo) / (n - 1) as f64;
            assert!((spacing - h).abs() <= 1e-12 * h, "spacing {} vs {}", spacing, h);
            assert!(top >= hi - 1e-9 && top - hi < h, "must cover the range by less than one element: {} {}", top, hi);
        }
    }
    #[test]
    fn exact_multiples_do_not_gain_an_element() {
        assert_eq!(mesh_axis(0.0, 0.3175, 0.015875).0, 21); // 20 cells -> 21 nodes
    }
}
