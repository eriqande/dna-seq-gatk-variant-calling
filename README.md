# yukon-chinookomes branch

This is a snakemake workflow.  This branch of it has been
configured to run through the Yukon Chinookomes.

Some modifications are as follows.

## Getting the Genome

Otsh_V2.0 is not up on Ensembl yet, so I can't have it downloaded directly
using the Snakemake workflow.  Instead I am going to just download it directly,
gunzip it and then put it into resources/genome.fasta.

I will download it from: 

https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/296/145/GCA_018296145.1_Otsh_v2.0/GCA_018296145.1_Otsh_v2.0_genomic.fna.gz

So, I got it with:
```{r}
(snakemake) [node34: resources]--% pwd
/home/eanderson/Documents/projects/yukon-chinookomes/resources
(snakemake) [node34: resources]--% wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/296/145/GCA_018296145.1_Otsh_v2.0/GCA_018296145.1_Otsh_v2.0_genomic.fna.gz

# that was ridiculously fast.


```
Then I gunzipped it and renamed it `resources/genome.fasta`

## Indexing the genome.

While waiting to get more quota, etc., I will index the genome.
```sh
snakemake --use-conda  --cores 8  resources/genome.fasta.{fai,bwt} resources/genome.dict
```


# Tweaking the workflow

I have done a fair bit, including adding a place for NMFS_DNA_ID.

Right now, I am working on just getting bams to the dedup stage,
and running some multiqc on all that.  The dag for it looks like
this:
  
  
![DAG on a small ficitious set](images/first-phase-dag.svg)
  

## Testing/Developing

I want to run it all the way through for one of the samples.  This one
is only 98 Mb:  `NVS127B_Columbus_T2088_R1/T200097_T2088_5C_S35_L002_R1_001.fastq.gz`

So, it is probably crappy DNA, but I can at least map it.  Looking at `units.tsv`
I see that we have:
```
sample_id = T200097
sample = s129
```
The final output before multiqc that we would want is of the form:
`qc/samtools-stats/{sample}-{unit}.txt`. So, `qc/samtools-stats/s129-1.txt`.
A quick,
```sh
snakemake --use-conda  -np  qc/samtools-stats/s129-1.txt
```
tells me that the logic of the workflow is correct.  I will now get all the
conda packages installed:
```sh
snakemake --use-conda  --conda-create-envs-only --cores 8  qc/samtools-stats/s129-1.txt
```

while that is running, I will get a slurm-sedna profile together.

I got that, now let us test it:
```sh
snakemake --use-conda  --profile ./slurm_profile  --jobs 100 qc/samtools-stats/s129-1.txt
```
That produced this output, from which we can get a read on how long
each of these steps will take.
```sh
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cluster nodes: 100
Provided resources: cpus=500, mem_mb=2300000
Singularity containers: ignored
Job counts:
        count   jobs
        1       map_reads
        1       mark_duplicates
        1       samtools_stats
        1       trim_reads_pe
        4
Select jobs to execute...

[Thu Jan  6 13:01:27 2022]
rule trim_reads_pe:
    input: /home/eanderson/Documents/projects/yukon-chinookomes/resources/fastq/NVS127B_Columbus_T2088_R1/T200097_T2088_5C_S35_L002_R1_001.fastq.gz, /home/eanderson/Documents/projects/yukon-chinookomes/resources/fastq/NVS127B_Columbus_T2088_R2/T200097_T2088_5C_S35_L002_R2_001.fastq.gz
    output: trimmed/s129-1.1.fastq.gz, trimmed/s129-1.2.fastq.gz, trimmed/s129-1.1.unpaired.fastq.gz, trimmed/s129-1.2.unpaired.fastq.gz, trimmed/s129-1.trimlog.txt
    log: logs/trimmomatic/s129-1.log
    jobid: 3
    wildcards: sample=s129, unit=1
    resources: cpus=1, mem_mb=4600

Submitted job 3 with external jobid '104829'.
Removing temporary output file trimmed/s129-1.1.unpaired.fastq.gz.
Removing temporary output file trimmed/s129-1.2.unpaired.fastq.gz.
[Thu Jan  6 13:03:59 2022]
Finished job 3.
1 of 4 steps (25%) done
Select jobs to execute...

[Thu Jan  6 13:03:59 2022]
rule map_reads:
    input: trimmed/s129-1.1.fastq.gz, trimmed/s129-1.2.fastq.gz, resources/genome.fasta.amb, resources/genome.fasta.ann, resources/genome.fasta.bwt, resources/genome.fasta.pac, resources/genome.fasta.sa
    output: mapped/s129-1.sorted.bam
    log: logs/bwa_mem/s129-1.log
    jobid: 2
    wildcards: sample=s129, unit=1
    threads: 4
    resources: cpus=5, mem_mb=4600

Submitted job 2 with external jobid '104830'.
Removing temporary output file trimmed/s129-1.1.fastq.gz.
Removing temporary output file trimmed/s129-1.2.fastq.gz.
[Thu Jan  6 13:18:53 2022]
Finished job 2.
2 of 4 steps (50%) done
Select jobs to execute...

[Thu Jan  6 13:18:53 2022]
rule mark_duplicates:
    input: mapped/s129-1.sorted.bam
    output: dedup/s129-1.bam, qc/dedup/s129-1.metrics.txt
    log: logs/picard/dedup/s129-1.log
    jobid: 1
    wildcards: sample=s129, unit=1
    resources: cpus=1, mem_mb=4600

Submitted job 1 with external jobid '104831'.
Removing temporary output file mapped/s129-1.sorted.bam.
[Thu Jan  6 13:20:14 2022]
Finished job 1.
3 of 4 steps (75%) done
Select jobs to execute...

[Thu Jan  6 13:20:14 2022]
rule samtools_stats:
    input: dedup/s129-1.bam
    output: qc/samtools-stats/s129-1.txt
    log: logs/samtools-stats/s129-1.log
    jobid: 0
    wildcards: sample=s129, unit=1
    resources: cpus=1, mem_mb=4600

Submitted job 0 with external jobid '104880'.
Removing temporary output file dedup/s129-1.bam.
[Thu Jan  6 13:20:34 2022]
Finished job 0.
4 of 4 steps (100%) done
Complete log: /home/eanderson/Documents/projects/yukon-chinookomes-dna-seq-gatk-variant-calling/.snakemake/log/2022-01-06T130125.375371.snakemake.log

```
So, how long did these take? Time in minutes
```r
trim_time <- 2.366667
map_time <- 14.9
mkdup_time <- 1.35
stats_time <- 0.3333
```
Here are the sizes for that:
```
(base) [sedna: fastq]--% du -h NVS127B_Columbus_T2088_R*/T200097*
98M	NVS127B_Columbus_T2088_R1/T200097_T2088_5C_S35_L002_R1_001.fastq.gz
100M	NVS127B_Columbus_T2088_R2/T200097_T2088_5C_S35_L002_R2_001.fastq.gz
```
So, that was 198M total sequences.  

The average size of a pair of fastqs in Gb is
```sh
(base) [sedna: fastq]--% du */*.fastq.gz | awk '{sum += $1; n++} END {print (2*sum/n) / 1e6 }'
4.55943
```
So, the average fastq pair is 4.55943 / 0.198 = 23.027 times larger than
our test case.

So, the average fastq will take 18.95 * 23.027 = 436.36
minutes, or 7.27 hours to finish.  That is 2791.68 hours.  But,
we can get 100 sets of 4 cpus, so that would be 27 hours.

Rippin!

## Further Tweaks

I must:

1. Move input/output directories and logs into `results`, just to keep it cleaner
in there.
2. Set MarkDuplicates to not discard them.  This means --REMOVE_DUPLICATES = false.
3. Rename the `dedup` directory to `mkdup`.
4. Make sure the `mkdup` outputs are protected rather than temp.  
5. I thought for a while that I might add barcode sequences to the PU, for completeness, but, since I don't know how gatk really uses them, and it is
unlikely they would be used correctly for non-model organisms, I ditched that.

## First run

```sh
# did this:
snakemake --use-conda  --profile ./slurm_profile  --jobs 500  results/qc/multiqc.html

# it mapped 374 of 384.  I think I had 9 trimmomatic failures and
# 1 mkdup failure

# full log is in here:
/home/eanderson/Documents/projects/yukon-chinookomes-dna-seq-gatk-variant-calling/.snakemake/log/2022-01-08T204846.726634.snakemake.log
```
Gonna restart and see if it was just transient read errors...
```sh
(snakemake) [node34: yukon-chinookomes-dna-seq-gatk-variant-calling]--% snakemake --use-conda  --profile ./slurm_profile  --jobs 200   results/qc/multiqc.html

Complete log: /home/eanderson/Documents/projects/yukon-chinookomes-dna-seq-gatk-variant-calling/.snakemake/log/2022-01-10T052721.802307.snakemake.log
```
So, those are bona fide failures!

Let's tally them up:
```sh
# get the output from snakemake
snakemake --use-conda  --profile ./slurm_profile  -np   results/qc/multiqc.html > /tmp/sm.txt

# then process it to summarize
awk '/^rule/ {rule = $2; next} /output:/ {out = $2; print rule, out}' /tmp/sm.txt | sort

```
I copied the output from that and put it into `prepare/failed-steps.txt`.
Then, I analyze it with: `prepare/002-mapping-failures.Rmd`.


# Snakemake workflow: dna-seq-gatk-variant-calling

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.14.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![GitHub actions status](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/workflows/Tests/badge.svg?branch=master)](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/actions?query=branch%3Amaster+workflow%3ATests)

This Snakemake pipeline implements the [GATK best-practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) for calling small germline variants.

## Authors

* Johannes Köster (https://koesterlab.github.io)


## Usage

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).


#### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files `config.yaml`, `samples.tsv` and `units.tsv`.

#### Step 3: Execute workflow

This workflow will automatically download reference genomes and annotation.
In order to save time and space, consider to use [between workflow caching](https://snakemake.readthedocs.io/en/stable/executing/caching.html) by adding the flag `--cache` to any of the commands below.
The workflow already defines which rules are eligible for caching, so no further arguments are required.
When caching is enabled, Snakemake will automatically share those steps between different instances of this workflow.

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details (e.g. cloud execution).

#### Step 4: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen [here](https://cdn.rawgit.com/snakemake-workflows/dna-seq-gatk-variant-calling/master/.test/report.html).

#### Step 5: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

#### Step 6: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/dna-seq-gatk-variant-calling.git` or `git remote add -f upstream https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master rules scripts envs schemas report > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config.yaml samples.tsv units.tsv`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.

#### Step 7: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

Test cases are in the subfolder `.test`. They are automtically executed via continuous integration with Github actions.
