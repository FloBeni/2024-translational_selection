# Collecting genome annotation from NCBI

This pipeline allows for the downloading of genomes and their annotations given a list of assemblies.

A command to test this pipeline is:

``` bash
snakemake collect_everything --configfile collecting_genome_annotation/config.json --cluster "sbatch -J {params.name} -p {params.partition} -N 1 --ntasks={params.ntasks} --mem={params.mem} -t {params.time} -o {params.out} -e {params.err}" --rerun-incomplete --rerun-triggers mtime -j 100 -n
```

To create a dag file :

``` bash
snakemake collect_everything --configfilecollecting_genome_annotation/config.json --forceall --dag | dot -Tpdf > dag-GTDrift.pdf
```
