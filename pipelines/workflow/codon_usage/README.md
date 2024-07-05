# Workflow for codon usage

This pipeline allows for the study of codon usage based on expression levels, tRNA pool, and within genomes given a list of assemblies.

A command to test this pipeline is:

``` bash
snakemake codon_usage_analysis --configfile config.json --use-singularity --singularity-args "--bind /beegfs/banque/gtdrift/:/beegfs/banque/gtdrift/" --cluster "sbatch -J {params.name} -p {params.partition} -N 1 --ntasks={params.ntasks} --mem={params.mem} -t {params.time} -o {params.out} -e {params.err}" --rerun-incomplete --rerun-triggers mtime -j 100 -n
```

To create a dag file :

``` bash
snakemake codon_usage_analysis --configfile config.json --forceall --dag | dot -Tpdf > dag-GTDrift.pdf
```
