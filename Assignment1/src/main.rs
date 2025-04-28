use serde::Deserialize;
use std::env::current_exe;
use std::fs;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::time::{Duration, Instant};

fn main(){
    // let start = Instant::now();
    let ns = [1000];
    
    // Write the JSON string to a file
    for i  in ns{
        println!("At value: {}", i as usize);
        let data = avarage_values(100, i);
        let json = serde_json::to_string_pretty(&data).expect("Failed to convert Hasmap to jsonstring");
        let new_path = format!("JSON_files/Square_Grid_simulation_{}.json", i);
        let mut file = File::create(new_path ).expect("failed to create file");
        let e = file.write_all(json.as_bytes());
        match e{
        Ok(())=>{
            println!("Data has been written to json");


        }
        Err(err)=>{
            println!("The issue happend at {}",i);
            print!("Its so over it failed to write, due to {}", err)
        }
    }
    }
    // let duration = start.elapsed();
    // println!("Completed simulation after: {:?}", duration);
}


// fn find_root_node(list_index: usize, sites: &Vec<f64>) -> usize {
//     let mut current = list_index; 
//     // Traverse to find the root node
//     while sites[current] >= 0.0 {
//         current = sites[current] as usize;
//     }
//     current // Return the root node
// }

fn find_root_node(list_index: usize, sites: &Vec<f64>) -> usize {
    let mut current = list_index;

    // Ensure the initial index is within bounds
    if current >= sites.len() {
        current -= 1;
    }

    // Traverse to find the root node
    while current < sites.len() && sites[current] >= 0.0 {
        current = sites[current] as usize;

        // Check if the new current index is out of bounds
        if current >= sites.len() {
            panic!("Traversal led to an out-of-bounds index: {}", current);
        }
    }

    current // Return the root node
}

fn swapping_bonds(bonds_list: &mut Vec<Vec<usize>>){
    for i in 1..(bonds_list.len()-2){
        let candidate_num = quad_rand::gen_range(i as f64,bonds_list.len() as f64) as usize;
        bonds_list.swap(i,candidate_num);
    }
}

#[derive(Deserialize, Debug)]
struct BondsList {
    #[serde(rename = "Bonds")]
    bonds: Vec<Vec<usize>>,
    #[serde(rename = "Number of bonds")]
    num_bonds: usize,
}


fn simulate_bonds_rust(num:i64, iteration:usize)-> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>){
    let path = format!("JSON_files/Square_grid_bonds_{}.json", num);
    let data = fs::read_to_string(path).expect("didnt find json file");
    let bonds:BondsList = serde_json::from_str(&data).expect("couldn't convert to rust type");
    // let n = bonds.bonds.last().expect("no last value")[0];
    let n = num.pow(2);
    let mut bonds_list: Vec<Vec<usize>> = bonds.bonds;
    for _ in 0..iteration{
        swapping_bonds(&mut bonds_list);
    }
    let num_bonds: f64 = bonds.num_bonds as f64;
    let mut sites = vec![-1.0; n as usize];
    let mut largest_clust= [1.0, -1.0]; // Initial largest cluster
    let mut number_of_activated_bonds: f64 = 0.0; // Counter for activated bonds
    println!("{:?}",bonds_list.len());
    let mut p_inf: Vec<f64> = vec![0.0; bonds_list.len()];
    let mut p_0: Vec<f64> = vec![0.0; bonds_list.len()];
    let mut p_inf_2: Vec<f64> = vec![(1.0/n as f64).powf(2.0); bonds_list.len()];
    let mut s: Vec<f64> = vec![0.0; bonds_list.len()];;
    let mut susept: Vec<f64> = vec![(n as f64) * (p_inf_2[0].powf(2.0) - p_inf[0].powf(2.0)).powf(1.0 / 2.0); bonds_list.len()];
    let mut avarage_s: f64 = n as f64; // Average size
    let mut s_step: f64 = (avarage_s - ((n as f64) * p_inf[0]).powf(2.0)) / ((n as f64) * (1.0 - p_inf[0]));
    for i in 0..(bonds_list.len()){
        number_of_activated_bonds += 1.0;
        let current_bonds =  &bonds_list[i];
        let node1= current_bonds[0];
        let node2 = current_bonds[1];
        let root1 = find_root_node(node1, &mut sites);
        let root2 = find_root_node(node2, &mut sites);
        let new_root: usize;
        if root1 != root2{
            avarage_s += -(sites[root1 as usize] as f64).powf(2.0) - (sites[root2] as f64).powf(2.0);
            if sites[root1] < sites[root2]{
                sites[root1] += sites[root2];
                sites[root2] = root1 as f64;
                new_root = root1;
            }
            else{
                sites[root2] += sites[root1];
                sites[root1] = root2 as f64;
                new_root = root2;
            }
            if sites[new_root] < largest_clust[1]{
                largest_clust[0] = new_root as f64;
                largest_clust[1] = sites[new_root];
            }
            else {
                s_step = (avarage_s - ((n as f64) * p_inf[i-1]).powi(2)) / ((n as f64) * (1.0 - p_inf[i-1]));
                s_step = if s_step.is_finite() { s_step } else { 0.0 };
            }
            avarage_s += sites[new_root].powf(2.0);
        }
        s[i] = s_step;
        p_0[i] = number_of_activated_bonds/num_bonds;
        p_inf_2[i] =(largest_clust[1]/(n as f64)).powf(2.0);
        let p_inf_2_mean = p_inf_2.iter().sum::<f64>() / p_inf_2.len() as f64;
        let p_inf_mean = p_inf.iter().sum::<f64>() / p_inf.len() as f64;
        susept[i] =  (n as f64) * ((p_inf_2_mean - p_inf_mean.powi(2)).max(0.0)).sqrt();
    }
    (p_0, p_inf, p_inf_2, susept, s)
}


fn avarage_values(iteration:usize, n:i64)-> HashMap<String, Vec<f64>>{
    let mut p_inf_list: Vec<Vec<f64>> = Vec::new();
    let mut susept_list: Vec<Vec<f64>> = Vec::new();
    let mut s_list: Vec<Vec<f64>> = Vec::new();
    let mut p_list: Vec<Vec<f64>> = Vec::new();
    let mut p_inf_2_list: Vec<Vec<f64>> = Vec::new();
    // let mut result: (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);
    let mut scores = HashMap::new();
    for i in 0..iteration{
        let result = simulate_bonds_rust(n, i);
        let p = result.0;
        let p_inf = result.1;
        let p_inf_2 = result.2;
        let susept = result.3;
        let s = result.4;
        p_list.push(p);
        p_inf_list.push(p_inf);
        p_inf_2_list.push(p_inf_2);
        susept_list.push(susept);
        s_list.push(s);
        println!("Done with {} iteration", i)
    }
    let p_inf_avarage =  calculate_average(&p_inf_list);
    scores.insert("p_inf".to_owned(), p_inf_avarage);
    let p_inf_2_avarage = calculate_average(&p_inf_2_list);
    scores.insert("p_inf_2".to_owned(), p_inf_2_avarage);
    let p_avarage = calculate_average(&p_list);
    scores.insert("p".to_owned(), p_avarage);
    let s_avarage = calculate_average(&s_list);
    scores.insert("s".to_owned(),s_avarage);
    let susept_avarage = calculate_average(&susept_list);
    scores.insert("susept".to_owned(), susept_avarage);
    return scores
}

fn calculate_average(data: &Vec<Vec<f64>>) -> Vec<f64> {
    if data.is_empty() {
        return vec![]; // Return an empty vector if there's no data
    }

    let num_columns = data[0].len();
    let mut averages = vec![0.0; num_columns];

    for row in data {
        for (i, &value) in row.iter().enumerate() {
            averages[i] += value;
        }
    }

    let num_rows = data.len() as f64;
    for avg in &mut averages {
        *avg /= num_rows; // Calculate the average
    }

    averages
}