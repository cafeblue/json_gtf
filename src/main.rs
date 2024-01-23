use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::collections::HashMap;
use bzip2::read::BzDecoder;
use serde_json::{Result, Value};

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn main() {

    let mut exon_posset: HashMap<String, HashMap<u32, bool>> = HashMap::new();
  // Read in gtf values
  if let Ok(lines) = read_lines("FOXG1.gtf".to_string()) {
	for line in lines {
		if let Ok(ip) = line {
			let spl: Vec<&str> = ip.split("\t").collect();
            if spl[2] != "exon" {
                continue;
            }
			let startp = spl[3].parse::<u32>().unwrap();
			let stopp = spl[4].parse::<u32>().unwrap();
            let mut currexon: HashMap<u32, bool> = HashMap::new();
            for i in startp..(stopp+1) {
                currexon.insert(i, true);
            }
            let desc = spl[8];
            let desc_parts: Vec<&str> = desc.split("; ").collect();
            for part in desc_parts {
                let part_parts: Vec<&str> = part.split(" ").collect();
                if part_parts[0] == "gene_name" {
                    let gene_id = part_parts[1].replace("\"", "");
                    if !exon_posset.contains_key(&gene_id) {
                        exon_posset.insert(gene_id, currexon.clone());
                    }
                    else {
                        let tmp = exon_posset.get(&gene_id).unwrap();
                        currexon.extend(tmp.clone().into_iter());
                        exon_posset.insert(gene_id, currexon.clone());
                    }
                }
            }
		}
	}
  }



    let path = "FOXG1.json.bz2";
    let file = File::open(path).unwrap();
    let reader = BufReader::new(BzDecoder::new(file));
    for line in reader.lines() {
        let l = line.unwrap();
        let v: Value = serde_json::from_str(&l).unwrap();

        let mut a_count: u32 = 0;
        let mut c_count: u32 = 0;
        let mut freq_idx: usize = 0;
        let mut tmp_pos: u32 = 0;
        let mut on_exon: bool = false;

        'marker: for idx in 0..v["primary_snapshot_data"]["placements_with_allele"].as_array().unwrap().len() {
            if v["primary_snapshot_data"]["placements_with_allele"][idx]["is_ptlp"].as_bool().unwrap() {
                for idx1 in 0..v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"].as_array().unwrap().len() {
                    if v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"][idx1]["allele"]["spdi"]["deleted_sequence"].as_str().unwrap() == v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"][idx1]["allele"]["spdi"]["inserted_sequence"].as_str().unwrap() {
                        tmp_pos = v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"][idx1]["allele"]["spdi"]["position"].as_u64().unwrap() as u32;
                        freq_idx = idx1;
                        break 'marker;
                    }
                }
            }
        }

        
        let gene_symbol = v["primary_snapshot_data"]["allele_annotations"][freq_idx]["assembly_annotation"][0]["genes"][0]["locus"].as_str().unwrap();
        if exon_posset.contains_key(&gene_symbol.to_string()) {
            let tmp = exon_posset.get(&gene_symbol.to_string()).unwrap();
            if tmp.contains_key(&tmp_pos) {
                on_exon = true;
            }
        }         
        let vt1 = v["primary_snapshot_data"]["allele_annotations"][freq_idx]["assembly_annotation"][0]["genes"][0]["rnas"][0]["sequence_ontology"][0]["name"].to_string();
        for idx in 0..v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"].as_array().unwrap().len() {
            a_count += v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"][idx]["allele_count"].as_u64().unwrap() as u32;
            c_count += v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"][idx]["total_count"].as_u64().unwrap() as u32;
        }
        let vt2: String = v["primary_snapshot_data"]["variant_type"].to_string();
        println!("{}:\tidx: {}\tposition: {}\tvt1: {}\tvt2: {}\ta_count: {}\tt_count: {}\ton_exon: {}", gene_symbol, freq_idx, tmp_pos, vt1, vt2, a_count, c_count, on_exon );
    }
    for key in exon_posset.keys() {
        println!("{}", key);
    }
}