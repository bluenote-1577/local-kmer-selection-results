for d in {3,5,7}
do
    for i in {1..2}
    do 
        sed -i "s/it = [0-9]/it = ${i}/g" Snakefile.smk
        sed -i "s/density = [0-9]/density = ${d}/g" Snakefile.smk
        snakemake -s Snakefile.smk -p --cores 20
    done
done

