# Results directory for paper

Simulation code, plotting code, and miscellaneous supplementary files/materials for the paper "Theory of local k-mer selection with application to long-read alignment" are found in this folder. 

# Simulation software

Computing conservation via simulations is done in rust. Files are in src. We compute empirical values for 
1. Words based method with W_8
2. Miniception 
3. Random minimizers
4. Open syncmers

Parameters can be changed in `src/main.rs`. Rust needs to be installed to run simulations.

To run simulations:
```
cargo build
cargo run
```

outputs a set of lists of values, which can be copy and pasted with some slight modifications into a cell in ``plots_kmer_select.ipynb``. 

# Comparing .sam outputs for minimizers and open syncmers

Given two sam files `output_using_syncmers.sam` and `output_using_minimizers.sam` that contain the same reads and same reference, we analyze them using these scripts.  

`analyze_sam_reads.py output_using_syncmers.sam output_using_minimizers.sam` tells the user how many unmapped reads are successfully mapped by each of the methods.

`plot_chain_scores.py output_using_syncmers.sam output_using_minimizers.am` is used for plotting histogram of chaining scores. Title and labels need to be changed. 

# Notebook 

## plots_kmer_select.ipynb

This notebook contains all plots from the paper. For the last plot, empirical results were taken from the simulations outside of the notebook and copy and pasted into the notebook.
