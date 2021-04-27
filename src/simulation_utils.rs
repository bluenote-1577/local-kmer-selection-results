use rand::distributions::{Bernoulli, Distribution};
use rand::{thread_rng, Rng};
use fxhash::{FxHashSet,FxHashMap};

pub fn gen_rand_string(n: usize) -> Vec<u8> {
    const CHARSET: &[u8] = b"ATCG";
    let mut rng = thread_rng();
    let return_string: Vec<u8> = (0..n)
        .map(|_| {
            let idx = rng.gen_range(0..CHARSET.len());
            CHARSET[idx] as u8
        })
        .collect();

    //println!("{:?}", return_string);
    return return_string;
}

pub fn gen_mutated_string(sequence: &[u8], theta: f64) -> (Vec<u8>,Vec<bool>) {
    let mut return_vec = vec![];
    let mut return_vec_bool = vec![];
    let d = Bernoulli::new(theta).unwrap();
    let mut rng = thread_rng();
    //Turn A to T, and C to G with probability theta. 
    for base in sequence {
        let x = d.sample(&mut rng);
        let mut new_base = b' ';
        if x {
            if *base == b'A' {
                new_base = b'C';
            }
            if *base == b'T' {
                new_base = b'G';
            }
            if *base == b'C' {
                new_base = b'T';
            }
            if *base == b'G' {
                new_base = b'A';
            }

            return_vec.push(new_base);
            return_vec_bool.push(false);
        } else {
            return_vec.push(*base);
            return_vec_bool.push(true);
        }
    }
    return (return_vec,return_vec_bool);
}

pub fn check_context_dependent_mutation(
    seeds_orig: &FxHashMap<&[u8], FxHashSet<usize>>,
    seeds_mut: &FxHashMap<&[u8], FxHashSet<usize>>,
    string_mut: &[u8],
    string_orig: &[u8],
    k: usize){

    let mut num_context_mutations = 0;
    let set1: FxHashSet<_> = seeds_orig.keys().cloned().collect();
    let set2: FxHashSet<_> = seeds_mut.keys().cloned().collect();
    let difference: FxHashSet<_> = set1.difference(&set2).collect();

    //We only take a subset of context mutated k-mers,
    //if somehow a mismatch k-mer was generated through a mutation,
    //this won't capture it.
    for kmer in difference{
        let positions_s = seeds_orig.get(kmer).unwrap();
        for pos in positions_s{
            //dbg!(String::from_utf8(string_mut[*pos..*pos+k].to_vec()),String::from_utf8(kmer.to_vec()));
            if string_mut[*pos..*pos+k] == **kmer{
                num_context_mutations += 1;
            }
            if string_orig[*pos..*pos+k] != **kmer{
                panic!("something went wrong");
            }
        }
    }

    dbg!(num_context_mutations);

}

pub fn get_conservation(
    seeds_orig: FxHashMap<&[u8], FxHashSet<usize>>,
    seeds_mut: FxHashMap<&[u8], FxHashSet<usize>>,
    k: usize,
    str_len: usize,
) -> (f64,f64) {
    let mut num_spurious_matches = 0;
    let set1: FxHashSet<_> = seeds_orig.keys().cloned().collect();
    let set2: FxHashSet<_> = seeds_mut.keys().cloned().collect();
    let intersection: FxHashSet<_> = set1.intersection(&set2).collect();
    let mut homologous_positions: FxHashSet<usize> = FxHashSet::default();
    let mut spurious_positions: FxHashSet<usize> = FxHashSet::default();
    //If we get the same k-mer on both strings, make sure that
    //they are at the same positions.
    //dbg!(intersection.len());
    for kmer in intersection {
        let positions_s = seeds_orig.get(kmer).unwrap();
        let positions_sp = seeds_mut.get(kmer).unwrap();

        let int_pos: FxHashSet<_> = positions_s.intersection(&positions_sp).collect();
        for pos in int_pos {
            for i in 0..k {
                homologous_positions.insert(pos + i);
            }
        }

        let bad_matches = positions_sp.difference(&positions_s).collect::<Vec<_>>();
        num_spurious_matches += bad_matches.len();
        for pos in bad_matches{
            for i in 0..k{
                spurious_positions.insert(pos+i);
            }
        }
    }

    //dbg!(num_spurious_matches);
    return (homologous_positions.len() as f64 / str_len as f64,
    spurious_positions.len() as f64 / str_len as f64);

}

pub fn get_schem_prob_minimizers(
    seeds_orig: &FxHashMap<&[u8], FxHashSet<usize>>,
    seeds_mut: &FxHashMap<&[u8], FxHashSet<usize>>,
    positions_changed: &Vec<bool>,
    k: usize
) -> Vec<f64>{

    let set1: FxHashSet<_> = seeds_orig.keys().cloned().collect();
    let set2: FxHashSet<_> = seeds_mut.keys().cloned().collect();
    let intersection: FxHashSet<_> = set1.intersection(&set2).collect();
    let mut union = FxHashSet::default();

    for key in intersection{
        let positions = seeds_orig.get(key).unwrap();
        for pos in positions{
            union.insert(pos);
        }
    }

    let mut successes = vec![0;k];
    let mut total_counts = vec![0;k];
    let mut running_length = 0;
    for (i,pos) in positions_changed.iter().enumerate(){
        if *pos == true && running_length < 2*k - 1{
            running_length += 1;
        }
        else{
            if running_length >= k{
                for j in i-running_length..i-k+1{
                    if union.contains(&j){
                        successes[running_length - k] += 1;
                        break;
                    }
                }
                total_counts[running_length - k] += 1;
            }
            if *pos == false{
                running_length = 0;
            }
        }
    }

    let mut probabilities = vec![0.0;k];
    for i in 0..k{
        probabilities[i] = successes[i] as f64/ total_counts[i] as f64;
    }

    return probabilities;

    //dbg!(probabilities,total_counts);
}


