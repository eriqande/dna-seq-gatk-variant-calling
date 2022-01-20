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


### Playing with HaplotypeCaller

I really need to know, roughly, how long it will take to make the individual
gVCFs.

So, I will do a test run with an average size bam file. Let'd do
`results/mkdup/s001-1.bam` which is 5.5 Gb.  Let's try:
```sh
# going to use the gatk that comes with sedna
module load bio/gatk
conda activate /opt/bioinformatics/miniconda3/envs/gatk

gatk --java-options "-Xmx8g" HaplotypeCaller  \
   -R resources/genome.fasta \
   -I results/mkdup/s001-1.bam \
   -O results/gvcf/s001.g.vcf.gz \
   -ERC GVCF
```
Note that we will have to index all these bams before doing this for real.

I started that at 20:18:36 and it seems to be doing about 3 Mb per minute which would suggest about 32 hours for the whole schmeer.

### How much slower is it if we use just one core?

Let's try this out.  I will put this on just one core
```sh
module load bio/gatk
conda activate /opt/bioinformatics/miniconda3/envs/gatk

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R resources/genome.fasta \
   -I results/mkdup/s001-1.bam \
   -O results/gvcf/s001.g.vcf.gz \
   -ERC GVCF > gatk-1core4thread.stdout 2> gatk-1core4thread.stderr
```

### How about if we set threads to one at the same time?

```sh
module load bio/gatk
conda activate /opt/bioinformatics/miniconda3/envs/gatk

FF=s001
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R resources/genome.fasta \
   -I results/mkdup/${FF}-1.bam \
   -O results/gvcf/${FF}.g.vcf.gz \
   -ERC GVCF > gatk-1core1thread.stdout 2> gatk-1core1thread.stderr
```

So, preliminary results are suggesting to me that using 4 cores with 4 threads
is only 1.26X faster than using 1 core with 4 threads (which is, itself,
slightly faster that 1 core with 1 thread.)  So, that is settled.  When I do
this for real, it will be 1 core per sample, and we can expect it to be about 36
to 40 hours for each sample.  But, if no one is on SEDNA, we can do them all at
once.

So, I killed those and I will now start three runs so I have gvcfs to play
with for importing genomic data bases, etc.
```sh
mkdir results/logs/gvcf-play

module load bio/gatk
conda activate /opt/bioinformatics/miniconda3/envs/gatk

FF=s003
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R resources/genome.fasta \
   -I results/mkdup/${FF}-1.bam \
   -O results/gvcf/${FF}.g.vcf.gz \
   -ERC GVCF > results/logs/gvcf-play/$FF.stdout 2> results/logs/gvcf-play/$FF.stderr

```

Aw crap! It looks like I didn't really set threads to 1.  Oh well.  I will just do it for the
FF=s003 job:
```sh
gatk --java-options "-Xmx4g" HaplotypeCaller     -R resources/genome.fasta    -I results/mkdup/${FF}-1.bam    -O results/gvcf/${FF}.g.vcf.gz --native-pair-hmm-threads 1   -ERC GVCF > results/logs/gvcf-play/$FF.stdout 2> results/logs/gvcf-play/$FF.stderr
```

#### Looking at the results of all that

Note that by the time the program has reached CM031232.1, it has processed
1992217379 bases.

So, let's see where the first time they hit that was:
```sh
(gatk) [node21: yukon-chinookomes-dna-seq-gatk-variant-calling]--% for i in results/logs/gvcf-play/*.stderr; do echo =====$i======; awk '/ProgressMeter/ && /CM031232.1/ {print; exit}' $i; done
=====results/logs/gvcf-play/s001.stderr======
08:07:30.817 INFO  ProgressMeter -    CM031232.1:635749            665.4              12690740          19071.8
=====results/logs/gvcf-play/s002.stderr======
09:35:00.226 INFO  ProgressMeter -    CM031232.1:932671            741.3              12418360          16751.9
=====results/logs/gvcf-play/s003.stderr======
05:51:39.382 INFO  ProgressMeter -    CM031232.1:820622            510.8              12131850          23748.4
```
OK, and how big were each of the bams?
```sh
(gatk) [node21: yukon-chinookomes-dna-seq-gatk-variant-calling]--% du -h results/mkdup/s00[1-3]-1.bam
5.6G	results/mkdup/s001-1.bam
6.6G	results/mkdup/s002-1.bam
4.0G	results/mkdup/s003-1.bam
```
So, let's scale all those together and compare:
```r
> c(665, 741, 510) / 510
[1] 1.303922 1.452941 1.000000

> c(5.6, 6.6, 4.0) / 4.0 
[1] 1.40 1.65 1.00
```
So, there might be a small advantage to setting the threads to 1.
Or maybe not.  Either way, these all ripped through all but one of the autosomes
in under 12.5 hours, with a single core.  So this is great news!


## GenomicsDBImport

Playing around with this.  There is no problem getting it to run on my three example
gVCFs above. Using a single core and no modifications to the memory for the Java machine:
```sh
# this is recommended for NFS filesystems
export TILEDB_DISABLE_FILE_LOCKING=1

# Doing just the first chromosome
gatk GenomicsDBImport -V results/gvcf-play/s001.g.vcf.gz -V results/gvcf-play/s002.g.vcf.gz -V results/gvcf-play/s003.g.vcf.gz --genomicsdb-workspace-path sandbox --intervals CM031199.1
```
It took about 8 minutes to run though that.  I don't know how it scales with number of samples,
but let's imagine that it is linear.  If so, then let's estimate 3 minutes per sample on this
average sized (maybe somewhat large-ish) chromosome.  That would be roughly 19 hours for all
384 samples.  OK.  

It doesn't seem that multithreading will improve this, but I should test it.  NIH HPC did some
benchmarking and adding more cores does not seem to help (nor does adding more memory).  For production, I think that running
it on two cores with "-Xmx4g" to the java engine and the other 5.4 Gb to the C++ parts of
GenomicsDBImport should be fine.  

With lots of samples, the program opens and maintains connections to many files (like one per
sample).  In extreme cases, this can start to slow things down and even crash your process.  So,
the recommend reading things in in batches of 50.  Apparently this can be slower than doing everyone
at once, but should be more fault tolerant.  A good discussion of these issues can
be found at: [https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GenomicsDBImport-usage-and-performance-guidelines](https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GenomicsDBImport-usage-and-performance-guidelines).  This is set with the `--batch-size` option.

That page mentions that multiple reader threads might be useful when opening files
in batches and dealing with filesystem latency.  So, `--reader-threads` might be a good
thing to use.  Also, they recommend setting
`--genomicsdb-shared-posixfs-optimizations` to `true` in this case.

Scattering it over chromosomes is definitely what we will want to do.  We would have one class of
jobs, each of which just operated on a single chromosome over all samples.
But then, we have to do something with the 8,000 or so small contigs.  In the Otsh_v2.0 genome
they are organized in the genome from shortest to longest, it appears.  GenomicsDBImport now
provides the `--merge-contigs-into-num-partitions` option, which is just what we want!

I think it would be worthwhile to carve those contigs into sets that are each of about the
same cumulative size---say, maybe about 50% of the average length of a chromosome---and then
we could do a separate job for each of those, merging all those contigs into a single partition.

To set this up, it would be easy to read the fai file and have the user provide a prefix for the
chromosomes (like CM).  Then we have a rule that uses R to set that all up, just writing the intervals
to text files that can be inserted onto the command line.  Killer!


I am continuing with this and now have a rule `genomics_db_import`.  I set it up to do s001--s004, for
whatever chromosome I name.  



## GenotypeGVCFs

I want to play around with this utility a bit, as well.

```sh
mkdir -p  /scratch/eanderson/tmp

gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R resources/genome.fasta \
   -V gendb://sandbox \
   -O CM031199.1.vcf.gz \
   --tmp-dir  /scratch/eanderson/tmp
```
That failed with an error.  But I think it is a bug in GATK.  I am
going to rerun with genomeDB's from the same, latest, version of GATK.

```sh
mkdir -p results/vcf_parts
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R resources/genome.fasta \
   -V gendb://results/genomics_db/chromosomes/CM031200.1 \
   -O results/vcf_parts/CM031200.1.vcf.gz
```

* with CM031199 that fails after writing 29760738 with `java.lang.IllegalStateException: Genotype [T199970 ATATATAT/T GQ 49 DP 4 AD 0,2,0,0,0,2,0,0 {SB=0,0,2,2}] does not contain likelihoods necessary to calculate posteriors.`
* with CM031200.1 that fails after writing 15955241 with: `java.lang.IllegalStateException: Genotype [T199968 AAA/TA GQ 57 DP 5 AD 0,0,0,3,2,0,0,0 {SB=0,0,3,2}] does not contain likelihoods necessary to calculate posteriors.`
* with CM031201.1 it fails after writing 45940161 with: `java.lang.IllegalStateException: Genotype [T199967 ATTT/TT GQ 34 DP 5 AD 0,3,2,0,0,0,0,0 {SB=0,0,4,1}] does not contain likelihoods necessary to calculate posteriors.` And the two warnings immediately before it are
```
Chromosome CM031201.1 position 45940163 (TileDB column 217279462) has too many alleles in the combined VCF record : 7 : current limit : 6. Fields, such as  PL, with length equal to the number of genotypes will NOT be added for this location.

Chromosome CM031201.1 position 45940171 (TileDB column 217279470) has too many alleles in the combined VCF record : 7 : current limit : 6. Fields, such as  PL, with length equal to the number of genotypes will NOT be added for this location.
```
So, that looks like two adjacent positions that both have too many alleles, and
hence won't have genotype likelihoods.

But, I think it might be that it fails when there is exactly one more than the
max number of alleles.  All of them fail when there are 7 alleles, and the max allowed number of alleles is 6.  

CM031199.1 actually doesn't have a problem at 6786063 where there are 8
alleles.  It just charges right through that.  So, I will try an experiment:
set the `--max-alternate-alleles` to 7 and see if it bombs on position 6786063.

It did!  I filed an issue on GitHub.

## Moving forward

It turns out that the bug occurred in GATK 4.2.4.1. As user seboyden pointed out:
```
To elaborate & clarify, in my use case:

GenomicsDBImport 4.2.4.0 plus GenotypeGVCFs 4.2.4.0: works
GenomicsDBImport 4.2.4.0 plus GenotypeGVCFs 4.2.4.1: fails
GenomicsDBImport 4.2.4.1 plus GenotypeGVCFs 4.2.4.0: works
GenomicsDBImport 4.2.4.1 plus GenotypeGVCFs 4.2.4.1: fails
```

So, I think that I should be fine just defining a new GATK environment
for the GenotypeGVCFs steps.

I had to reinstall conda (I used mambaforge) to get that to work, but I am
back now.

### Tossing the 9 corrupted files for now

This gets the files that we are still waiting for:
```sh
(snakemake) [node36: yukon-chinookomes-dna-seq-gatk-variant-calling]--% snakemake --use-conda  --profile ./slurm_profile --jobs 50 -np   results/qc/multiqc.html | awk '/wildcards: sample=/ {print $2}' | sed 's/sample=//g; s/,//g;' | sort | uniq > prepare/failures_1.txt
(snakemake) [node36: yukon-chinookomes-dna-seq-gatk-variant-calling]--% cat prepare/failures_1.txt
s221
s227
s243
s244
s250
s265
s305
s343
s360
```
So, I should be able to merely remove them from the `samples.tsv` and 
`units.tsv`, and then run through everything.  I'll at least get a chance to
see how long everything takes.
```sh
cp samples.tsv full_samples.tsv
cp units.tsv full_units.tsv

# get the non-faily samples
(cat prepare/failures_1.txt; echo xxxxxxxxxxxxxxx; cat full_samples.tsv) | awk '/xxxxxxxxxx/ {go=1; next} go==0 {tosser[$1]++} go==1 {if(!($1 in tosser)) print}'  > samples.tsv

# get the non-failing units
(cat prepare/failures_1.txt; echo xxxxxxxxxxxxxxx; cat full_units.tsv) | awk '/xxxxxxxxxx/ {go=1; next} go==0 {tosser[$1]++} go==1 {if(!($1 in tosser)) print}' > units.tsv

```
Now we can look at the DAG like this:
```sh
snakemake --use-conda  --profile ./slurm_profile --jobs 50 -np  --dag | condense_dag | dot -Tsvg > toss-9-corrupted-fastqs.svg
```
Which looks like this:
![DAG for 375 non-corrupted fastqs](images/toss-9-corrupted-fastqs.svg)


So, let's give it a whirl:
```sh
(snakemake) [node36: yukon-chinookomes-dna-seq-gatk-variant-calling]--% snakemake --use-conda  --profile ./slurm_profile --jobs 60
```

That cruised along happily, but I forgot that SEDNA now imposes an 8 hour default
time limit on jobs.  So, about six of the genomics_db_imports failed, and almost
all of the genotypegvcfs were headed to failure.  So I scancelled everything and 
removed all the vcf_sections.  It looks like the latter were getting through
about 25 megabases in 8 hours.  And the longest chromosomes are about
100 megabases.  So, 24 hours should be sufficient.  I will give each of those
jobs 3 days by default....And the genomicsDBimports should all get 36 hours,
I believe, just to be on the safe side.  (It will all need more with more
individuals).

Also, some of the genotype_gvcfs failed with a runtime memory of 4Gb.  So, I 
really should up that to 8 Gb.  But, if I do that, then I probably ought to up
things to 4 cores...Actually, I will do 2 CPUs but up the memory a bit.  So, I 
will do 8Gb to the Java engine, but 11750 total.  That will let me get 8 jobs on each node.

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
