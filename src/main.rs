use bio::alignment::sparse;
use fxhash::{FxHashMap, FxHashSet};
use simulations::seeding_methods;
use simulations::simulation_utils;
use statistical;

fn main() {
    let n = 2;
    let k = 24;
    let w = 15;
    let s = k - w / 2 + 1 / 2 as usize;
    //let s = 10;
    let t = (k - s + 1) / 2 + 1 as usize;
    let k_0 = k - w+2;
    let add_w = 2;
    let num_iters = 100;
    //let t = 3;
    let thetas = vec![
        0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14,0.15
    ];
    let string_length = 5 * 10_usize.pow(5);
    println!("k = {}, s = {}, k-s+1 = {}, t = {}", k, s, k - s + 1, t);


    let mut os_cons_over_theta = vec![];
    let mut cept_cons_over_theta = vec![];
    let mut cust_cons_over_theta = vec![];
    let mut min_cons_over_theta = vec![];

    let mut min_cons_std_over_theta = vec![];
    let mut cept_cons_std_over_theta = vec![];
    let mut cust_cons_std_over_theta = vec![];

    let mut min_dens_trials = vec![];
    let mut cept_dens_trials = vec![];
    for theta in thetas {
        let mut os_cons_trials = vec![];
        let mut min_cons_trials = vec![];
        let mut cept_cons_trials = vec![];
        let mut min_pr_vec_trials = vec![];
        let mut cust_cons_trials = vec![];

        for i in 0..num_iters {
            println!("iteration {},{}", i, theta);
            let s_orig = simulation_utils::gen_rand_string(string_length);
            let (s_prime, positions_changed) =
                simulation_utils::gen_mutated_string(s_orig.as_slice(), theta);
            //dbg!(String::from_utf8(s.clone()),String::from_utf8(s_prime.clone()));

            //MINIMIZERS
            //println!("Seeding minimizers");
            let (minimizers_s, minimizers_s_num) =
                seeding_methods::minimizer_seeds(s_orig.as_slice(), w, k);
            let (minimizers_sprime, _minimizers_sp_num) =
                seeding_methods::minimizer_seeds(s_prime.as_slice(), w, k);
            let dens_min = minimizers_s_num as f64 / (string_length - k + 1) as f64;
            min_dens_trials.push(dens_min);
            let min_pr_vec = simulation_utils::get_schem_prob_minimizers(
                &minimizers_s,
                &minimizers_sprime,
                &positions_changed,
                k,
            );
            min_pr_vec_trials.push(min_pr_vec);
            //    simulation_utils::check_context_dependent_mutation(&minimizers_s,&minimizers_sprime,s_prime.as_slice(),s_orig.as_slice(),k);

            //OPEN SYNCMERS
            //println!("Seeding syncmers");
            let (syncmer_s, syncmer_s_num) =
                seeding_methods::open_sync_seeds(s_orig.as_slice(), k, s, t);
            let (syncmer_sp, _syncmer_sp_num) =
                seeding_methods::open_sync_seeds(s_prime.as_slice(), k, s, t);
            let dens_sync = syncmer_s_num as f64 / (string_length - k + 1) as f64;
            simulation_utils::get_schem_prob_minimizers(
                &syncmer_s,
                &syncmer_sp,
                &positions_changed,
                k,
            );

            //MINICEPTION
            //println!("Seeding miniception");
            let (miniception_s, miniception_s_num) =
                seeding_methods::miniception_seeds(s_orig.as_slice(), w - add_w, k, k_0);
            let (miniception_sprime, _miniception_sprime_num) =
                seeding_methods::miniception_seeds(s_prime.as_slice(), w - add_w, k, k_0);
            let dens_minicept = miniception_s_num as f64 / (string_length - k + 1) as f64;
            cept_dens_trials.push(dens_minicept);
            //println!("Done seeding");
            simulation_utils::get_schem_prob_minimizers(
                &miniception_s,
                &miniception_sprime,
                &positions_changed,
                k,
            );

            //WORDS
            //println!("Seeding wordss");
            let (words_s, words_s_num) = seeding_methods::words_seeds(s_orig.as_slice(), n, k);
            let (words_sp, _words_sp_num) = seeding_methods::words_seeds(s_prime.as_slice(), n, k);
            let dens_words = words_s_num as f64 / (string_length - k + 1) as f64;
            simulation_utils::get_schem_prob_minimizers(&words_s, &words_sp, &positions_changed, k);

            println!("Density for minimizer/syncmer = {},{}",dens_min,dens_sync);
            println!("Density for miniception = {}",dens_minicept);
            //println!("Density for words = {}",dens_words);
            
            //CUST WORDS_5_4
            let (cust_words_s, cust_words_s_num) = seeding_methods::custom_words_seeds_8_8(s_orig.as_slice(), n, k);
            let (cust_words_sp, _cust_words_sp_num) = seeding_methods::custom_words_seeds_8_8(s_prime.as_slice(), n, k);
            let dens_cust_words = cust_words_s_num as f64 / (string_length - k + 1) as f64;
            //simulation_utils::get_schem_prob_minimizers(&words_s, &words_sp, &positions_changed, k);

            let (cons1, spur1) = simulation_utils::get_conservation(
                minimizers_s,
                minimizers_sprime,
                k,
                string_length,
            );
            min_cons_trials.push(cons1);
            let (cons2, spur2) =
                simulation_utils::get_conservation(syncmer_s, syncmer_sp, k, string_length);
            os_cons_trials.push(cons2);
            let (cons3, spur3) = simulation_utils::get_conservation(
                miniception_s,
                miniception_sprime,
                k,
                string_length,
            );
            cept_cons_trials.push(cons3);
            let (cons4, spur4) =
                simulation_utils::get_conservation(words_s, words_sp, k, string_length);

            let (cons5, spur5) =
                simulation_utils::get_conservation(cust_words_s, cust_words_sp, k, string_length);
            cust_cons_trials.push(cons5);

//            println!(
//                "Cons | Spur | C-S for Minmizers = {},{},{}",
//                cons1,
//                spur1,
//                cons1 - spur1
//            );
            //        println!("Cons | Spur | C-S for OS = {},{},{}",cons2,spur2,cons2-spur2);
            //        println!("Cons | Spur | C-S for Minicept = {},{},{}",cons3,spur3,cons3-spur3);
            //        println!("Cons | Spur | C-S for Words = {},{},{}",cons4,spur4,cons4-spur4);
        }

        let min_cons_mean = statistical::mean(&min_cons_trials);
        let os_cons_mean = statistical::mean(&os_cons_trials);
        let cept_cons_mean = statistical::mean(&cept_cons_trials);
        let cust_cons_mean = statistical::mean(&cust_cons_trials);

        let min_cons_std = statistical::standard_deviation(&min_cons_trials, Some(min_cons_mean));
        let cept_cons_std = statistical::standard_deviation(&cept_cons_trials, Some(cept_cons_mean));
        let cust_cons_std = statistical::standard_deviation(&cust_cons_trials, Some(cust_cons_mean));

        min_cons_over_theta.push(min_cons_mean);
        os_cons_over_theta.push(os_cons_mean);
        cept_cons_over_theta.push(cept_cons_mean);
        cust_cons_over_theta.push(cust_cons_mean);

        min_cons_std_over_theta.push(min_cons_std);
        cept_cons_std_over_theta.push(cept_cons_std);
        cust_cons_std_over_theta.push(cust_cons_std);
    }

    let cept_dens_mean = statistical::mean(&cept_dens_trials);
    dbg!(os_cons_over_theta);
    dbg!(cept_dens_mean);

    dbg!(min_cons_over_theta);
    dbg!(min_cons_std_over_theta);

    dbg!(cept_cons_over_theta);
    dbg!(cept_cons_std_over_theta);

    dbg!(cust_cons_over_theta);
    dbg!(cust_cons_std_over_theta);
}
