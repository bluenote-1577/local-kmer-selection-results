RAN = range(11,25+1)
density = 3
it = 3

EXP_FILE = "NanoSim/pre-trained_models/human_NA12878_cDNA_Bham1_guppy/expression_abundance.tsv"
TRANSCRIPTOME_REF = "ref/Homo_sapiens.GRCh38.cdna.all.fa"
TRAINING_FOLDER = "NanoSim/pre-trained_models/human_NA12878_cDNA_Bham1_guppy/"
REF_GENOME = "ref/GCF_000001405.26_GRCh38_genomic.fna"
OS_MINIMAP_BIN = "os-minimap2/minimap2"
NANOSIM_BIN = "NanoSim/src/simulator.py"

rule all:
    input:
        expand("aln_SMK%d/k{ran}_sync_d%d.sam" %(it,density), ran=RAN)

rule generate_reads:
    input:
        exp_file = EXP_FILE,
        transcriptome = TRANSCRIPTOME_REF,
        training_folder = TRAINING_FOLDER,
        ref_genome = REF_GENOME

    output: 
        "reads/transcript/simulated_d%di%d_aligned_reads.fasta" %(density,it)
    shell:
        "%s transcriptome -rt {input.transcriptome} -e {input.exp_file} -c {input.training_folder}/training -rg {input.ref_genome} -o reads/transcript/simulated_d%di%d" %(NANOSIM_BIN, density,it)

rule sync_map:
    input:
        REF_GENOME,
        "reads/transcript/simulated_d%di%d_aligned_reads.fasta" %(density,it)
    output:
        expand("aln_SMK%d/k{ran}_sync_d%d.sam" %(it,density),ran=RAN),
        expand("aln_SMK%d/k{ran}_mini_d%d.sam" %(it,density),ran=RAN)
    shell:
        "for i in {{11..25}}; do %s -t 20 -ax splice -k $i --syncs $((i - %d)) {input} > aln_SMK%d/k${{i}}_sync_d%d.sam; done;for i in {{11..25}}; do %s -t 20 -ax splice -k $i -w %d {input} > aln_SMK%d/k${{i}}_mini_d%d.sam; done" %(OS_MINIMAP_BIN, density - 1, it, density, OS_MINIMAP_BIN, 2 * density - 1, it, density)



        


