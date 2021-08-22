# Transcriptome mapping experiment reproduction

## Code for running the experiment outlined in Section 5.2.3 of the paper for ONT cDNA mapping. 

### Required software and files:

1. [NanoSim](https://github.com/bcgsc/NanoSim)  Note that you **must unzip the pre-trained models before using this pipeline**. The pre-trained model used is NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy.tar.gz.
    1. Making sure Python3-pip is available
        ```
        sudo apt install python3-pip
        ```
    1. Installing NanoSim
        ```
        mkdir -p $HOME/software
        cd $HOME/software
        git clone https://github.com/bcgsc/NanoSim.git
        cd $HOME/software/NanoSim
        pip3 install -r requirements.txt
        ```
    1. Unzipping pre-trained models
        ```
        cd $HOME/software/NanoSim/pre-trained_models
        for file in *.tar.gz; do tar -xzf $file; done
2. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
3. [os-minimap2](https://github.com/bluenote-1577/os-minimap2) 
```
 git clone https://github.com/bluenote-1577/os-minimap2
 cd os-minimap2 && make
```
4. GRCh38 reference genome (can be found at https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz). This exact reference transcriptome must be used, as the pre-trained model is based on this specific transcriptome. 
5. ensembl GRCh38 reference transcriptome (can be found at http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz)
6. [SciencePlots](https://github.com/garrettj403/SciencePlots) is used for plotting. 
    1. Installing SciencePlots via Pip3
        ```
        pip3 install matplotlib
        pip3 install SciencePlots
        
### Running the experiment

The first step is to modify the `Snakefile.smk` so that the paths are correct.

1. EXP_FILE - expression file used by NanoSim for simulating transcripts. Should be `human_NA12878_cDNA_Bham1_guppy/expression_abundance.tsv` where `human_NA12878_cDNA_Bham1_guppy` is the **unzipped** pre-trained model folder.
2. TRANSCRIPTOME_REF - the ensembl GRCh38 reference transcriptome. 
3. TRAINING_FOLDER - we use the model `human_NA12878_cDNA_Bham1_guppy` from NanoSim's pre-trained model folder. Again, this is the **unzipped** pre-trained model folder. 
4. REF_GENOME - GRCh38 reference genome.
5. OS_MINIMAP_BIN - the `minimap2` executable file compiled from the os-minimap2 directory. Important: this is _not_ the standard minimap2, this is the minimap2 binary with open syncmer support.
6. NANOSIM_BIN - ``NanoSim/src/simulator.py`` is the python file which is in the src folder of NanoSim. 

After modification, run 
0. Index reference files using `samtools faidx` 
1. `./multiple_iteration.sh > experiment.log 2>&1` to align simulated reads over a range of parameters. Output to log is needed to get running time.
2. `scripts/get_times_from_log.py experiment.log` to get the runtimes
3. `scripts/transcriptome_plot (TRANSCRIPTOME_FILE) run_times.pkl aln_SMK*` to generate the plot for the experiment. 
