# Genotype phasing pipeline

## 1. Description
Pipeline for genotype phasing using [SHAPEIT4](https://odelaneau.github.io/shapeit4/). The pipeline automatically chunks chromosomes into overlapping regions, distributes chunks for phasing accross multiple compute nodes, and concatinates chunks based on overlapping flanking regions. The pipeline was developed using [Nextflow](https://www.nextflow.io/) and was tested on [SLURM](https://slurm.schedmd.com/documentation.html) job scheduler.

## 2. Prerequisites
The following software is required:
- Nextflow
- bcftools

## 3. Installation
To run this pipline you will need:
1. Download and install (SHAPEIT4](https://odelaneau.github.io/shapeit4/)

## 4. Execution
1. Modify `nextflow.config` configuration file:
* `params.unphased_vcfs` -- path to VCF/BCF files with unphased genotypes. Each VCF/BCF file must have the corresponding tbi/csi index.
* `params.window` --  size of the chromosome chunk in base-pairs for distributed phasing. See [SHAPEIT4 documentation](https://odelaneau.github.io/shapeit4/#documentation) for details.
* `params.flank` -- size of the flanking region around chunk in base-pairs.
* `params.shapeit_exec` --  path to the SHAPEIT4 executable
* `params.shapeit_maps` -- path to the SHAPEIT4 folder with genetic maps. NOTE: For GRCh37, you need to modify `shapeit` command line inside `phase_chunks` process in `Phasing.nf`.
* `process.*` and `executor.*` -- set this arguments according to your compute cluster configuration.
2. Run pipleine. Example of interactive SLURM job:
```
salloc --time=12:00:00 --ntasks=1 --mem-per-cpu=16G
module load nextflow
module load bcftools
nextflow run Phasing.nf
```
