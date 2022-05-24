# BSA Snakemake Pipeline

This is a snakemake pipeline for bulk segregation analysis.

This pipeline was designed to be running in IGDB-HPCC.

## Dependencies

All dependencies have been listed in `bsa-pipeline-env.yaml`, and can be
easily installed using conda.

```bash
mamba env create --file bsa-pipeline-env.yaml
```

## Configuration

To run this pipeline, you need these three type of data.

1. Reference Genome FASTA data (downloaded from EnsemblPlants)
2. Annotation GFF3 data (downloaded from EnsemblPlants)
3. Your bulk sequencing FASTQ data

All file paths need to be filled in the config file `bsa-pipeline-config.yaml`.

## Run Pipeline

After all preparation was down, using these commands to run this pipeline
in IGDB-HPCC.

```bash
conda activate bsa
snakemake \
    --snakefile bsa-pipeline.smk \
	--cluster 'bsub -n {resources.core_num} -R "rusage[mem={resources.mem_gb}]" -o {log}' \
	--default-resources core_num=1 mem_gb=1 \
	--jobs 24
```

A series of files will be generated in the working directory, including results.

