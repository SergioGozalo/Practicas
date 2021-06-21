# EMOSE Gene Catalogue

###### Based on TARA and Pablo Sánchez pipeline.
###### Written by Lidia Montiel.

## Introduction

The EMOSE gene catalogue includes **50 metagenomes** from one sampling day with different filter size and filtered water volume. This markdown explains the methods used to build the **gene catalogue** and the **gene** abundance tables starting from the predicted genes. Functional abundance analysis is already done by people from Banyuls.

### 1. Gene Catalogue building

Starting from the CDS of each sample, identifiers are modified to simplify them. Initially, identifiers have this format `>EMOSE-GC_ERZXXXXXX_XX` but later they are changed to `>EMOSE-GC_000000001` to simplify them more.

```bash
ls *gz > pre_list
sed 's/_.*$//' pre_list > CDS.list
```


```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="rename"
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/log_rename_%A_%a.out
#SBATCH --error=logs/log_rename_%A_%a.err
#SBATCH --array=1-50%10

# Variables
sample=$(awk "NR==$SLURM_ARRAY_TASK_ID" CDS.list)

# renaming
zless ${sample}_FASTA_CDS_annotated.ffn.gz | awk '/^>/{print ">EMOSE-GC_" ++i; next}{print}' | sed 's/GC_/GC_'${sample}'_/g' > ${sample}_renamed_CDS_annotated.ffn  # | head -n 100

gzip ${sample}_renamed_CDS_annotated.ffn
```

Building a preliminary and the final catalogue with the best identifier format.

```bash
zcat ERZ80*named* > EMOSE.pre-catalogue.fasta
perl -ane 'if(/\>/){$a++;printf ">EMOSE-GC_%09d\n", $a}else{print;}' EMOSE.pre-catalogue.fasta > EMOSE.catalogue.fasta
```

Also, a coordinate file is created to correlate new headers with the old ones.

```bash
grep ">" EMOSE.catalogue.raw.fasta > raw_headers
grep ">" EMOSE.pre-catalogue.fasta > headers_pre
grep ">" EMOSE.catalogue.fasta > headers
paste headers headers_pre old/raw_headers | sed 's/>//g' > coordinates.txt
```

To remove redundancy, we tried `linclust` from [MMseqs2 v.9-d36de](https://github.com/soedinglab/MMseqs2) since [CD-HIT](http://weizhongli-lab.org/cd-hit/) can take months to cluster. We clustered at **95% of identity** and __80% of coverage__.

```bash
#!/bin/bash

#SBATCH --account=emm1
#SBATCH --job-name="clust"
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=logs/create_linclust_v2.log_%j.out
#SBATCH --error=logs/create_linclust_v2.log_%j.err

# modules
module load mmseqs2

# variables
input=EMOSE.catalogue.fasta
output=EMOSE.catalogue_95id_clust
OUT=renamed_emose.create-lincluster
DB=${OUT}/EMOSE.catalogue_noclust_db

# job
mkdir ${OUT}

echo "Creating database before clustering..."
date

# Before clustering, convert your database into the MMseqs2 database format
mmseqs createdb ${input} ${OUT}/EMOSE.catalogue_noclust_db

echo ""
echo "--- createdb IS DONE --- $(date) "
echo ""
echo "Clustering... $(date)  "

# linear time clustering (faster but less sensitive)
mmseqs linclust ${OUT}/EMOSE.catalogue_noclust_db ${OUT}/${output} ${OUT} -c 0.8 --min-seq-id 0.95 # --cluster-mode 1 --threads 3
# --min-seq-id: matches above this sequence identity (0.95 as we usually did with CD-HIT)
# -c: is the coverage
# --cluster-mode: 0 is setcover, 1 is connected component, 2 is gredy clustering by sequence length

echo ""
echo "--- CLUSTERING IS DONE! --- $(date) "
echo ""
echo "Running now the conversion to fasta..."

# editing results in an informative way
mmseqs createtsv ${DB} ${DB} ${OUT}/${output} ${OUT}/resultsDB_clu_mode0.tsv

# representative sequence of the clusters
mmseqs result2repseq ${DB} ${OUT}/${output} ${OUT}/EMOSE.catalogue_95id_repseq_mode0

# results in fasta
mmseqs result2flat ${DB} ${DB} ${OUT}/EMOSE.catalogue_95id_repseq_mode0 ${OUT}/${output}.fasta --use-fasta-header

echo ""
echo "--- ALL THE PROGRAM IS DONE! --- "
date
```

Clustering information is gathered in the following table:

|                       | EMOSE gene catalogue |
| --------------------- |:--------------------:|
| Total predicted genes |      22,192,039      |
| Non-redundant genes   |    __11,339,535__    |

### 2. Gene Annotation

In order to obtain the gene abundance table normalized by single copy genes, an annotation of the catalogue with COG database is needed. First, catalogue is translated using `transeq`, a program inside the [EMBOSS v6.6.0](http://emboss.sourceforge.net/index.html) Open Software Suite, and splitted into 50 files using `Perl` scripts called `fasta_spliter.pl`:

```bash
#!/bin/bash

#SBATCH --job-name="trans_split"
#SBATCH --account=emm1               # to request the queue
#SBATCH --mem=50G
#SBATCH --time=10-00:00:00           # limit wall clock time
#SBATCH -D .                        # the working directory of the job
#SBATCH --error=logs/log_transplit_%j.err        # file to collect the standard error output (stderr)
#SBATCH --output=logs/log_transplit_%j.out       # file to collect the standard output (stdout)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Translation
module load emboss

indir=../gene.catalogue
file=EMOSE.catalogue_95id_last

## Command
date
echo "Translating..."
transeq ${indir}/${file}.fasta ${indir}/${file}.faa
echo "Translation finished!"
date
echo " "

## Splitting
module load perl

# CMD
date
echo "Splitting..."
perl fasta_spliter.pl --n-parts 50 --measure count --part-num-prefix bgc_ --out-dir ../gene.annotation/chunks/ ${indir}/${file}.faa
echo "Splitting finished!"
date
```

Gene annotation using `rpsblast` from [BLAST v2.7.1](https://blast.ncbi.nlm.nih.gov/Blast.cgi):

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH -D .
#SBATCH --job-name=emose_cog
#SBATCH --time=100-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/log_annot_cog_%A_%a.out
#SBATCH --error=logs/log_annot_cog_%A_%a.err
#SBATCH --array=1-50%4

# Load modules
module load blast

# vars
infile=$(awk "NR==$SLURM_ARRAY_TASK_ID" input.list)
indir=chunks
outdir=annotation
#db=~/dbs/databases/cdd/Cog/Cog
cog=/mnt/lustre/repos/bio/databases/public/NCBI/CDD/CDD.NCBI_release_3.14/Cog/Cog
pfam=/mnt/lustre/repos/bio/databases/public/pfam/pfam_release_31.0/Pfam-A.hmm

date
echo "Starting COG annotation..."
# create output dir
mkdir -p ${outdir}/COG

# annotation COG
rpsblast -soft_masking true -query ${indir}/${infile} -db ${cog} -evalue 0.1 -outfmt 6 -out ${outdir}/COG/${infile}.cog

date
echo "COG anno finished!"
```

#### Parsing annotation

When annotation finishes, all the output files are put all together in one single file.

```bash
find *cog | while read i; do cat ${i}; done > EMOSE_COG_results
```

If the IDs are sorted, we start to parse the hits using the Perl script `hmmScanFilter.pl`.

```bash
sort -k1,1 EMOSE_COG_results > sorted_EMOSE_COG_results
perl ../../../scripts/hmmScanFilter.pl -i sorted_EMOSE_COG_results -r
```

Before to continue with the parsing of COG annotation, the `_1` from the gene IDs are removed.

```bash
sed 's/_1\t/\t/' sorted_EMOSE_COG_results.annot > sorted_EMOSE_COG_results.annot.raw
```

Finally, descriptions to the CDD IDs are added to be more informative using the `cog_descriptions.sh` script:

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="cog_postpro"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/log_postpro_%J.out
#SBATCH --error=logs/log_postpro_%J.err

# Variables
infile=annotation/COG/sorted_EMOSE_COG_results.annot.raw
outfile=annotation/COG/sorted_EMOSE_COG_results.annot.raw.description.txt
cddid=/mnt/lustre/repos/bio/databases/public/NCBI/CDD/CDD.NCBI_release_3.16/cddid_cog.tbl

cat ${infile} | bash cog-rpsblast-postprocess.sh - ${cddid} > ${outfile}
```

And COG parsing is done!
**9,535,099 genes** from the EMOSE gene catalogue have been annotated with COG database.

### 3. Abundance Analysis

We built abundance profiles to the EMOSE gene catalogue. To do so, we mapped the cleaned reads to the catalogue using [BWA v0.7.17](http://bio-bwa.sourceforge.net/). We first did `bwa index` and then we mapped with `bwa mem` running [Samtools v1.8](http://www.htslib.org/) in one step. Thus, we skipped building huge SAM files and get filtered BAM files automatically. Unmapped reads, alignments with map quality smaller than 10% and secondary hits were removed from these bam files.

`bwa index`:
```bash
#!/bin/bash

#SBATCH --job-name="bwa_index"
#SBATCH --account=emm1               # to request the queue
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00           # limit wall clock time
#SBATCH -D .                        # the working directory of the job
#SBATCH --error=bwa_index.err        # file to collect the standard error output (stderr)
#SBATCH --output=bwa_index.out       # file to collect the standard output (stdout)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20          # 10 tasks in the same node, not filling all the node (48 tasks is the top)

date
## BWA
module load bwa

## Var
assembly=EMOSE.catalogue_95id_clust.fasta

## CMD
bwa index $assembly
echo ""
echo "finished!"
date
```

`bwa mem`:
```bash
#!/bin/bash

#SBATCH --job-name="emose_bwa_samtools"
#SBATCH --time=30-00:00:00
#SBATCH --account=emm1
#SBATCH --mem=200G
#SBATCH -D .
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
##SBATCH --constraint="largemem"
##SBATCH --partition=batch
#SBATCH --output=logs/log_bwa_samtools_%A_%a.out
#SBATCH --error=logs/log_bwa_samtools_%A_%a.err
#SBATCH --array=1-50%4

# Module
module load bwa gcc samtools

# Variables
r1=$(awk "NR==$SLURM_ARRAY_TASK_ID" R1.list)
r2=$(awk "NR==$SLURM_ARRAY_TASK_ID" R2.list)
contigs_index=/mnt/lustre/bio/shared/EMOSE2017/gene.catalog/mapping/index/EMOSE.catalogue_95id_clust.fasta
DATA_DIR=/mnt/lustre/bio/shared/EMOSE2017/metaG/data/url.txt
SAMPLENAME=$(basename ${r1} | perl -pe 's/_\d.fastq.gz//')

# CMD
echo "Starting mapping with BWA..."
echo ""
# map and filter out UNMAP (4) and Secondary (256): -F 260
bwa mem -t 32 ${contigs_index} ${DATA_DIR}/${r1} ${DATA_DIR}/${r2} | samtools view -h -q 10 -F 260 -b - > ${SAMPLENAME}.bam 2>> logs/${SAMPLENAME}.log

# sort and index bam file
samtools sort -@ 32 ${SAMPLENAME}.bam -o ${SAMPLENAME}.sorted.bam 2>> logs/${SAMPLENAME}.log
samtools index ${SAMPLENAME}.sorted.bam > ${SAMPLENAME}.sorted.bam.bai 2>> logs/${SAMPLENAME}.log

echo ""
echo "Job array finished!"
date
```


The __gene catalogue fasta file__ was converted to __gff file__ using a custom Perl script written by Pablo Sánchez named `fasta2gff.pl`. To start counting genes, we will use [HTSeq v0.10.0](https://htseq.readthedocs.io/en/release_0.11.1/) which requires the gff file.

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="fasta2gff"
#SBATCH --time=01-00:00:00
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/log_fasta2gff_%j.out
#SBATCH --error=logs/log_fasta2gff_%j.err

# Modules

# Variables
indir=/mnt/lustre/bio/shared/EMOSE2017/gene.catalog/clustering/renamed_emose.create-lincluster
infile=EMOSE.catalogue_95id_clust
outdir=../htseq
outfile=${infile}
fasta2gff=fasta2gff.pl

# CMD
perl ${fasta2gff} ${indir}/${infile}.fasta > ${outdir}/${outfile}.gff
```

The **gff file conversion** is done!


Once the __gff file__ and the __mapping__ are done, we start counting how many reads mapped to each gene from the catalogue with `htseq-count`.

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="htseq"
#SBATCH --time=5-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120000
#SBATCH --output=logs/log_htseq_%j.out
#SBATCH --error=logs/log_htseq_%j.err
#SBATCH --array=1-50%10

# Modules
module load htseq

# Variables
#export PATH=/home/x_sanchep/miniconda2/bin:$PATH
indir=../mapping
sample=$(awk "NR==$SLURM_ARRAY_TASK_ID" sample.list)
outdir=../htseq/read_counts
gff=../htseq/EMOSE.catalogue_95id_clust.gff

# CMD
mkdir ${outdir}
htseq-count -f bam -s no -t CDS --idattr ID -r pos --nonunique all ${indir}/${sample}.sorted.bam ${gff} > ${outdir}/${sample}.htseq-count.txt 2>> logs/htseq-count.${sample}.log
```

After `htseq` run, we have counts for each fragment for each gene in the catalogue. Now, we count all fragments in each sample (reads). As reads are in __fastq__ format (each read belongs to 4 lines), we do a basic count `wc -l` of pair 1 file lines and divide by 2 to obtain counts from pair 1 and pair 2.

`count_library.sh`:
```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="count"
#SBATCH --time=30-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --mem=120000
#SBATCH --output=logs/count_log_%A_%a.out
#SBATCH --error=logs/count_log_%A_%a.err
#SBATCH --array=1-50%20

# Variables

indir=/mnt/lustre/bio/shared/EMOSE2017/metaG/data/url.txt
sample=$(awk "NR==$SLURM_ARRAY_TASK_ID" sample.list)
outdir=../htseq/fragment_counts
# count lines in pair1 fastq file (each read belongs to 4 lines)
count1P=$(zcat ${indir}/${sample}_1.fastq.gz | wc -l)
# divide by 2 to get the number of fragments (reads P1 + reads P2)
countPE=$(echo "$count1P/2" | bc)

mkdir -p ${outdir}
echo -e "${sample}\t${countPE}" > ${outdir}/${sample}_read.count.txt
```

With this, a file for each sample will be created. These new files will have a value to normalize the read counts per library size and we need to put together these values to a common file:

```bash
cat *.txt | sort -k1,1 > fragment.count.per.sample.txt
```

#### Raw counts table

Now, we'll start to build the raw counts table...

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="build_tbl"
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/build_raw_log_%j.out
#SBATCH --error=logs/build_raw_log_%j.err

# Variables
indir="../htseq/read_counts"
outdir="../tables"
first=$(ls ${indir}/*htseq-count.txt | head -n 1)
echo "awk '{print \$1}' ${first} > ${outdir}/EMOSE_GC_95id.gene.raw.counts.tbl" | bash -x

while read sample; do find ${indir}/ -name "${sample}.htseq-count.txt" | while read i; do join -t $'\t' ${outdir}/EMOSE_GC_95id.gene.raw.counts.tbl ${i} > ${outdir}/raw.tmp && mv ${outdir}/raw.tmp ${outdir}/EMOSE_GC_95id.gene.raw.counts.tbl; done; done < sample_list
```

The header is added to the raw table and last lines generated in the htseq-count program are deleted:

```bash
echo "gene" > header1; cat header1 ../scripts/sample.list | tr '\n' '\t' > header.txt
cat header.txt EMOSE_GC_95id.gene.raw.counts.tbl > tpm
grep -v '^_' tpm > EMOSE_GC_95id.gene_raw_counts.tbl
rm -rf tpm
```

NB: `header.txt` had a tabulation at the end of the line. It should has been removed to build the following tables properly.

#### Converting raw counts to TPM by gene length and the scaling factor

We'll convert the raw counts to TPM for each sample by first dividing the raw count per gene length in kbp. Then divide that value by the quocient between the sample read count and the 10^6 scaling factor.

In the command:
- `$7` = read length
- `($2/($7*0.001))/(rc/1000000)` = raw count (`$2`) divided by the read length in kbp in the numerator and the sample fragment count (`rc`) divided by the scaling factor in the denominator.

`get_TPM.sh`:
```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="get_TPM"
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/get_TPM_log_%j.out
#SBATCH --error=logs/get_TPM_log_%j.err

read_counts=../htseq/read_counts
fragments=../htseq/fragment_counts

while read sample; do rc=$(awk '{print $2}' ${fragments}/${sample}_read.count.txt); grep -v '^_' ${read_counts}/${sample}.htseq-count.txt | paste - ../htseq/EMOSE.catalogue_95id_last.gff | awk -v rc=${rc} '{print $1"\t"$2"\t"$7"\t"($2/($7*0.001))/(rc/1000000)}' > ${read_counts}/${sample}.TPM.txt; done < sample.list
```

Now we have done the proper calculations, we put together the TPMs and build a single table:

`build_TPM_table.sh`:
```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="build_TPM"
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/build_TPM_log_%j.out
#SBATCH --error=logs/build_TPM_log_%j.err

# Variables
indir="../htseq/read_counts"
outdir="../tables"
first=$(ls ${indir}/*TPM.txt | head -n 1)
echo "awk '{print \$1}' ${first} > ${outdir}/EMOSE-GC_gene_TPM_no_header.tbl" | bash -x

while read sample; do find ${indir}/ -name "${sample}.TPM.txt" | while read i; do awk '{print $1"\t"$4}' ${i} | join -t $'\t' ${outdir}/EMOSE-GC_gene_TPM_no_header.tbl -  > ${outdir}/raw.tmp && mv ${outdir}/raw.tmp ${outdir}/EMOSE-GC_gene_TPM_no_header.tbl; done; done < sample.list
```

Adding the header...

```bash
cat tables/header.txt tables/BBMO_Megahit.v1_gene_TPM_no_header.tbl > BBMO_Megahit.v1_gene_TPM.tbl
```

### Table based on normalized counts

Taken from the **Malaspina catalogue Markdown**, Pablo Sánchez textually says the following:
> *TARA is not using TPMs and it may be a bit confusing (mainly because the T in TPM means transcripts). I think it's a better metric to get the gene length normalized counts by dividing the raw counts by the gene length in kb, and then use any of the normalization strategies to compare between samples (namely the correction by the mean of several universal single-copy genes).*

This means we will build normalized tables by the gene length (normalized within samples) and later by the mean of 10 universal single-copy genes (normalized between samples).

#### Gene-length normalized table

This is the **gene table**, it is a correction of the gene counts by dividing the values with the gene length in kb. The following jobscript uses the bash script called `normalizeTableByGeneLength.sh` where gathers all the calculations done.

`normalize_by_gene_length.sh`:

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="norm_geneLength"
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/norm_genelength_log_%j.out
#SBATCH --error=logs/norm_genelength_log_%j.err

# Variables
rawCounts=../tables/EMOSE_GC_95id.gene_raw_counts.tbl
gffFile=../htseq/EMOSE.catalogue_95id_last.gff
outfile=../tables/EMOSE-GC_gene_lengthNorm_counts.tbl

bash normalizeTableByGeneLength.sh ${rawCounts} ${gffFile} > ${outfile}
```

This left the gene abundance table normalized by gene length.


#### Table based on Single Copy Genes

The recently built functional tables can be used to compare abundances of genes within a sample or metagenome. To compare across samples, we have to add another layer of normalization and this is about to calculating the abundance of the genes that have been annotated with 10 Universal Single Copy genes:

* **COG0012 GTP1**; Ribosome-binding ATPase YchF, GTP1/OBG family
* **COG0016 PheS**; Phenylalanyl-tRNA synthetase alpha subunit
* **COG0018 ArgS**; Arginyl-tRNA synthetase
* **COG0172 SerS**; Seryl-tRNA synthetase
* **COG0215 CysS**; Cysteinyl-tRNA synthetase
* **COG0495 LeuS**; Leucyl-tRNA synthetase
* **COG0525 ValS**; Valyl-tRNA synthetase
* **COG0533 TsaD**; tRNA A37 threonylcarbamoyltransferase TsaD
* **COG0541 Ffh**; Signal recognition particle GTPase
* **COG0552 FtsY**; Signal recognition particle GTPase

We'll take first the genes annotated for these 10 single copy genes annotated in the COG annotation file.

```bash
bash scripts/get_SCG.sh gene.annotation/annotation/COG/sorted_EMOSE_COG_results.annot.raw.description.txt > htseq/EMOSE_single_copy_genes.txt
```

**99,458 genes** are one of the 10 universal single-copy genes!

Another table is built with the sum of each single copy gene per sample. But, since I had sort problems with this step, I tried to fix it adding an `a` in the beginning of the `gene` from the header. Then I proceeded:

```bash
sed 's/gene/agene/' ../tables/EMOSE-GC_gene_lengthNorm_counts.tbl > tpm.tbl
bash scripts/condense_SCG_table.sh htseq/EMOSE_single_copy_genes.txt htseq/tpm.tbl > tables/EMOSE_single_copy_gene_with_counts.txt
rm -rf htseq/tpm.tbl
```

The recently table has 10 rows belonging to each universal single-copy gene and it has columns as number of samples or metagenomes. After this, geometric mean of single-copy gene COG counts per sample are obtained:

```bash
module load R
Rscript scripts/single_copy_gene_Median.R tables/EMOSE_single_copy_gene_with_counts.txt > tables/EMOSE_single_copy_genes_geometric_mean.txt
```

Finally, we build a **gene table normalised by gene length but now normalized by single-copy gene counts**.

`normalize_by_Single_Copy_COGs.sh`:

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="norm_SCG"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/norm_SCG_log_%J.out
#SBATCH --error=logs/norm_SCG_log_%J.err

# Module
module load python/2.7.15

# Variables
countTable=../tables/EMOSE-GC_gene_lengthNorm_counts.tbl
geomMeanTable=../tables/EMOSE_single_copy_genes_geometric_mean.txt
outfile=../tables/EMOSE-GC_gene_lengthNorm_SingleCopyGeneNorm_counts.tbl

python normalizeCountsToSingleCopyCOGs.py --counts ${countTable} --metric ${geomMeanTable} > ${outfile}
```

And now, **the abundance gene table normalized by gene length and single copy genes is finally done!**

### Normalization by Metagenome Size

Since every metagenome has different size (sequencing depth, more or less reads, etc), it will be used to normalize among samples. The metagenome-size normalization layer is added in the gene-length normalization. This normalization is crucial for our eukaryotic gene catalogues, since single-copy COG normalization is not considered yet, and it is relevant even for catalogues from different datasets.

Getting first the number of reads in reverse and forward `fastq` files:

`number.bases.sh`:
```bash
#!/bin/bash

#SBATCH --account=emm1
#SBATCH --job-name="bases_mp"
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/bases.log_%A_%a.out
#SBATCH --error=logs/bases.log_%A_%a.err
#SBATCH --array=1-155%10


sample=$(awk "NR==$SLURM_ARRAY_TASK_ID" sample.list)
dir=/mnt/lustre/bio/shared/ecomics/malaspina/metagenomes/Malaspina.metaG/malaspina.all

number1=$( pigz -dc ${dir}/ST*_${sample}*R1.fastq.gz |
     awk 'NR%4==2{c++; l+=length($0)}
          END{
                #print "Number of reads: "c;
                print "Number of bases in reads: "l
              }')

echo -e "${sample}\t${number1}" >> R1.numbers

number2=$( pigz -dc ${dir}/ST*_${sample}*R2.fastq.gz |
     awk 'NR%4==2{c++; l+=length($0)}
          END{
                #print "Number of reads: "c;
                print "Number of bases in reads: "l
              }')

echo -e "${sample}\t${number2}" >> R2.numbers
```

Summing R1+R2, determining the Gb size and creating the required-file format:

```bash
sort -k1,1 R1.numbers > R1.numbers.sorted
sort -k1,1 R2.numbers > R2.numbers.sorted
paste *sorted | awk '{print $1 "\t" ($7+$14)/1000000000}' > metaG.Gb
echo -e "Sample\tGb" > header
cat header metaG.Gb > EMOSE_metag_size_Gb.txt
```

The python script `normalizeCountsToSingleCopyCOGs.py` for the single-copy COGs normalization will be reused to generate the final table:

`normalize_by_MetaG_size_Gb.sh`:
```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="norm_Gb_EMOSE"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/norm_Gb_log_%J.out
#SBATCH --error=logs/norm_Gb_log_%J.err

# Module
module load python/2.7.15

# Variables
countTable=../tables/EMOSE-GC_gene_lengthNorm_counts.tbl
metaGsizeTable=../number.bases/EMOSE_metag_size_Gb.txt
outfile=../tables/EMOSE-GC_gene.lengthNorm.metaGsizeGbNorm.counts.tbl

python normalizeCountsToSingleCopyCOGs.py --counts ${countTable} --metric ${metaGsizeTable} > ${outfile}
```

## Final tables

Filenames to each gene abundance table:

| Table                     | Filename                                                 |
|:------------------------- | -------------------------------------------------------- |
| Raw counts                | "EMOSE-GC_95id.gene_raw_counts.tbl"                      |
| TPM                       | "EMOSE-GC_gene_TPM.tbl"                                  |
| Normalized by gene length | "EMOSE-GC_gene_lengthNorm_counts.tbl"                    |
| Normalized by SCG         | "EMOSE-GC_gene_lengthNorm_SingleCopyGeneNorm_counts.tbl" |
| Normalized by MetaG Size  | "EMOSE-GC_gene.lengthNorm.metaGsizeGbNorm.counts.tbl"    |

All these tables are found in `/mnt/lustre/bio/shared/EMOSE2017/gene.catalog/tables`.

# EMOSE gene catalogue following our pipeline

### 1. Assembly

_Ramiro' steps here_

### 2. Gene Prediction

**IMPORTANT POINT:** We detected that `linclust` adds two white spaces at the end of the IDs of each gene. To remove these final spaces, you can run the following `awk` command:

```bash
awk '{$1=$1;print}' catalog.fasta > fixed_catalog.fasta
```

#### Filtering minimum length

In others gene catalogues, we also filtered by minimum length to avoid very small genes. We filtered those gene catalogues by >250 bp length. To filter them, we run a _Perl_ script:

`remove_smalls.sh`:
```bash
```

Gathering all the information and comparing it to the Alex's gene catalogue version:

| EMOSE gene catalogue  |  Alex version  |  Our pipeline  |
| --------------------- |:--------------:|:--------------:|
| Total predicted genes |   22,192,039   |  135,199,238   |
| Non-redundant genes   | __11,339,535__ | __78,222,851__ |
| > 250 bp genes        |   9,615,623    |   36,726,233   |

### 3. Gene Annotation

To start to annotate our genes, the gene catalogue **fasta** is converted to a **faa file** using `transeq`, a program inside the [EMBOSS v6.6.0](http://emboss.sourceforge.net/index.html) Open Software Suite. Also, to optimize the annotation, the catalogue is splitted to 50 subfiles using a `Perl` script named `fasta_spliter.pl`. The next jobscript translates the sequences into proteins and split the **faa file**.

`trans.sh`:
```bash
#!/bin/bash

#SBATCH --job-name="trans"
#SBATCH --account=emm1               # to request the queue
#SBATCH --mem=30G
#SBATCH --time=10-00:00:00           # limit wall clock time
#SBATCH -D .                        # the working directory of the job
#SBATCH --error=logs/log_trans_%j.err        # file to collect the standard error output (stderr)
#SBATCH --output=logs/log_trans_%j.out       # file to collect the standard output (stdout)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Translation
module load emboss

# Variables
gc=EMOSE-GC_95id_min250
dir=../gene.catalogue

## Command
date
echo "Translating..."
transeq  ${dir}/${gc}.fasta  ${dir}/${gc}_pre.faa
echo "Translation finished!"
date
```

After translation, an `_1` is added at the end of the gene id. It is removed before proceeding to the split:

```bash
sed 's/_1$//' EMOSE-GC_95id_min250_pre.faa > EMOSE-GC_95id_min250.faa
```

`split.sh`:
```bash
#!/bin/bash

#SBATCH --job-name="split"
#SBATCH --account=emm1               # to request the queue
#SBATCH --mem=30G
#SBATCH --time=10-00:00:00           # limit wall clock time
#SBATCH -D .                        # the working directory of the job
#SBATCH --error=logs/log_split_%j.err        # file to collect the standard error output (stderr)
#SBATCH --output=logs/log_split_%j.out       # file to collect the standard output (stdout)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

## Splitting
module load perl

# Variables
gc=EMOSE-GC_95id_min250
dir=../gene.catalogue
outdir=../gene.annotation

# CMD
date
echo "Splitting..."
perl fasta_spliter.pl --n-parts 200 --measure count --part-num-prefix bk_ --out-dir ${outdir}/chunks/ ${dir}/${gc}.faa
echo "Splitting finished!"
date
```

We will annotate functions using these 4 databases:
- **COG** with `rpsblast`
- **Pfam** with `hmmsearch`
- **KEGG** with `diamond blastp`
- **eggNOG** with `hmmsearch`

#### COG

To annotate with COG:

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH -D .
#SBATCH --job-name=emose_cog
#SBATCH --time=100-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/cog.annot.log_%A_%a.out
#SBATCH --error=logs/cog.annot.log_%A_%a.err
#SBATCH --array=1-200%4

# Load modules
module load blast

# vars
infile=$(awk "NR==$SLURM_ARRAY_TASK_ID" input.list)
indir=chunks
outdir=annotation
#db=~/dbs/databases/cdd/Cog/Cog
cog=/mnt/lustre/repos/bio/databases/public/NCBI/CDD/CDD.NCBI_release_3.14/Cog/Cog

date
echo "Starting COG annotation..."
# create output dir
mkdir -p ${outdir}/COG

# annotation COG
rpsblast -soft_masking true -query ${indir}/${infile} -db ${cog} -evalue 0.1 -outfmt 6 -out ${outdir}/COG/${infile}.cog

date
echo "COG anno finished!"
```


#### Pfam

To annotate with Pfam:

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH -D .
#SBATCH --job-name=emose_pfam
#SBATCH --time=100-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/log_annot_pfam_%A_%a.out
#SBATCH --error=logs/log_annot_pfam_%A_%a.err
#SBATCH --array=1-200%4

# Load modules
module load blast

# vars
infile=$(awk "NR==$SLURM_ARRAY_TASK_ID" input.list)
indir=chunks
outdir=annotation
#db=~/dbs/databases/cdd/Cog/Cog
pfam=/mnt/lustre/repos/bio/databases/public/pfam/pfam_release_31.0/Pfam-A.hmm

date
echo "Pfam annotation is starting..."

# Load modules
module load gcc hmmer

# create output dir
mkdir -p ${outdir}/Pfam

# annotation Pfam
hmmsearch --noali --domtblout ${outdir}/Pfam/${infile}.hmmsearch.pfam.dom -E 0.1 --cpu 8 ${pfam} ${indir}/${infile} # > /dev/null

date
echo "Pfam anno is finished!"
```

It took 3 days to complete the annotation with Pfam. There are 49,047,180 hits.

##### Parsing Pfam

Once Pfam annotation is finished, all the comment lines that start with `#` are deteled while ORF IDs are sorted in chunk files:

`pfam_sort.sh`:
```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="pfam_sort"
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/pfam_sort_%A_%a.out
#SBATCH --error=logs/pfam_sort_%A_%a.err
#SBATCH --array=1-200%10

# Variables
infile=$(awk "NR==$SLURM_ARRAY_TASK_ID" ../gene.annotation/input.list)
dir=../gene.annotation/annotation/Pfam/raw.output

# job
grep -v '^#' ${dir}/${infile}.hmmsearch.pfam.dom | sort -k1,1 > ${dir}/${infile}.hmmsearch.pfam.dom.sorted
```

All these sorted files are concatenated in one single file and removed afterwards:

```bash
cat *sorted > ../EMOSE-GC_95id_min250.pfam.dom.raw
rm -rf *sorted
sort -k1,1 EMOSE-GC_95id_min250.pfam.dom.raw > EMOSE-GC_95id_min250.pfam.dom.raw.sorted
```

Finally, best domain(s) for each ORF is parsed with the same script `hmmScanFilter.pl`:

`parse_pfam.sh`:
```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="parse_pfam"
#SBATCH --time=30-00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/log_parse_pfam_%J.out
#SBATCH --error=logs/log_parse_pfam_%J.err

# Variables
indir="../gene.annotation/annotation/Pfam"
infile="EMOSE-GC_95id_min250.pfam.dom.raw"

perl hmmScanFilter.pl -i ${indir}/${infile}.sorted -t
```

Pfam parsing is done! __23,360,506 genes__ from the EMOSE-ICM Gene Catalogue version have been annotated with Pfam database. There are some peptides with 2 different HMM models for the same peptide, as it is considered that it has 2 or more valid domains.

#### KEGG

To annotate KEGG:

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH -D .
#SBATCH --job-name=emose_kegg
#SBATCH --time=100-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --output=logs/log_annot_kegg_%j.out
#SBATCH --error=logs/log_annot_kegg_%j.err

# Load modules
module load diamond/0.9.22

# vars
keggDB=/mnt/lustre/repos/bio/databases/public/kegg/kegg_release_20170220/diamond-db/kegg20170220-genes.dmnd
outdir=annotation/kegg
infile=../gene.catalogue/EMOSE-GC_95id_min250.faa

# Make temporary directory in scratch
project=$(pwd | xargs -n 1 basename)
tmp=/mnt/lustre/scratch/lmontiel/${project}
mkdir -p ${tmp}
mkdir -p ${tmp}/diamond_emose-tmp

# Blast with diamond
diamond blastp -d ${keggDB} -q ${infile} -o ${tmp}/EMOSE-GC_95id_250bp.faa.kegg.raw -e 0.1 --sensitive --tmpdir ${tmp}/diamond_emose-tmp -p 22

# clean up
mkdir -p ${outdir}
mv ${tmp}/EMOSE-GC_95id_250bp.faa.kegg.raw ${outdir}
rm -rf ${tmp}/diamond_emose-tmp
```

It took 1 days and 18 hours to annotate with KEGG.

##### Parsing KEGG

Annotations have to be sorted as:

```bash
sort -k1,1 EMOSE-GC_95id_250bp.faa.kegg.raw > EMOSE-GC_95id_250bp.faa.kegg.raw.sorted
```

As Pfam databases, KEGG results are parsed using the same script `hmmScanFilter.pl`:

`parse_kegg.sh`:
```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="parse_kegg"
#SBATCH --time=30-00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/log_parse_kegg_%J.out
#SBATCH --error=logs/log_parse_kegg_%J.err

# Variables
indir="../gene.annotation/annotation/kegg"
infile="EMOSE-GC_95id_250bp.faa.kegg.raw.sorted"

perl hmmScanFilter.pl -i ${indir}/${infile} -r
```

It took 3 hours to parse the results and left __22,288,565 parsed annotated genes__ from a total of 498,684,893 annotated genes with KEGG.

Then, **KEGG genes hits** are parsed to **KEGG KOs** with the following jobscript called `kegg-genes-2-kegg-ko.slurm`. This step is going to link the KEGG genes hits with their KEGG KO numbers using the file `ko_genes-sorted-by-gene.list` located in the KEGG public database path in *Marbits*. The resulting file will be a two-column table: `<gene ID> \t <KO>` sorted by gene ID.

```bash
#!/bin/bash

#SBATCH --job-name="kegg2ko"
#SBATCH --account=emm1
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/log_kegg2ko_%j.out
#SBATCH --error=logs/log_kegg2ko_%j.err

# Variables
dir=../gene.annotation/annotation/kegg
infile=EMOSE-GC_95id_250bp.faa.kegg.raw.sorted.annot
outfile=EMOSE-GC_95id_250bp.faa.kegg.annot.ko
keggGene2ko=/mnt/lustre/repos/bio/databases/public/kegg/kegg_release_20170220/genes/ko/ko_genes-sorted-by-gene.list

bash kegg-genes-2-kegg-ko.sh ${dir}/${infile} ${keggGene2ko} > ${dir}/${outfile}
```

This left __15,142,793 genes with KEGG KO numbers__.

As some genes will have more than one KO, we collapsed the KOs using the Python script `condenseRowsByColumn1Key.py` running `kegg-ko-collapsed-by-geneId.sh`. If one gene has two or more KOs, these KO IDs will be separated with commas `,`.

```bash
#!/bin/bash

#SBATCH --job-name="collapse_ko"
#SBATCH --account=emm1
#SBATCH --time=1-30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/collapseKO_%j.out
#SBATCH --error=logs/collapseKO_%j.err

# Variables
dir=../gene.annotation/annotation/kegg
infile=EMOSE-GC_95id_250bp.faa.kegg.annot.ko
outfile=EMOSE-GC_95id_250bp.faa.kegg.annot.ko.out

python condenseRowsByColumn1Key.py -i ${dir}/${infile} -o ${dir}/${outfile}
```

Finally, the `ko:` prefixes are removed which it could be done earlier. It is done with a basic Perl command.

`kegg_clean_ko_prefix.slurm`:
```bash
#!/bin/bash

#SBATCH --job-name="clean_ko"
#SBATCH --account=emm1
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/clean_ko_%j.out
#SBATCH --error=logs/clean_ko_%j.err

# Variables
dir=../gene.annotation/annotation/kegg
infile=EMOSE-GC_95id_250bp.faa.kegg.annot.ko.out

bash kegg_clean_ko_prefix.sh ${dir}/${infile}
```

KEGG parsing is done!

This final file (`EMOSE-GC_95id_250bp.faa.kegg.annot.ko.out`) left __15,108,436 genes__ that may be linked to one or more KEGG KOs.

#### EggNOG

As EggNOG annotation can last 4 months, I splitted the protein file into 200 subfiles for it. Then, the 200 subfiles were annotated as:


### 3. Abundance Analysis

We built abundance profiles to the EMOSE gene catalogue. To do so, we mapped the cleaned reads to the catalogue using [BWA v0.7.17](http://bio-bwa.sourceforge.net/). We first did `bwa index` and then we mapped with `bwa mem` running [Samtools v1.8](http://www.htslib.org/) in one step. Thus, we skipped building huge SAM files and get filtered BAM files automatically. Unmapped reads, alignments with map quality smaller than 10% and secondary hits were removed from these bam files.

`bwa index`:
```bash
#!/bin/bash

#SBATCH --job-name="bwa_index"
#SBATCH --account=emm1               # to request the queue
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00           # limit wall clock time
#SBATCH -D .                        # the working directory of the job
#SBATCH --error=bwa_index.err        # file to collect the standard error output (stderr)
#SBATCH --output=bwa_index.out       # file to collect the standard output (stdout)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

date
## BWA
module load bwa

## Var
assembly=EMOSE-GC_95id_min250.fasta

## CMD
bwa index $assembly

echo "bwa index finished!"
date
```

`bwa mem`:
```bash
#!/bin/bash

#SBATCH --job-name="emose_bwa"
#SBATCH --time=30-00:00:00
#SBATCH --account=emm1
#SBATCH --mem=200G
#SBATCH -D .
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
##SBATCH --constraint="largemem"
##SBATCH --partition=batch
#SBATCH --output=logs/log_bwa_samtools_%A_%a.out
#SBATCH --error=logs/log_bwa_samtools_%A_%a.err
#SBATCH --array=1-50%5

# Module
module load bwa gcc samtools

# Variables
r1=$(awk "NR==$SLURM_ARRAY_TASK_ID" R1.list)
r2=$(awk "NR==$SLURM_ARRAY_TASK_ID" R2.list)
contigs_index=index/EMOSE-GC_95id_min250.fasta
DATA_DIR=/mnt/lustre/bio/shared/EMOSE2017/metaG/data/url.txt
SAMPLENAME=$(basename ${r1} | perl -pe 's/_\d.fastq.gz//')

# CMD
echo "Starting mapping with BWA..."
echo ""
# map and filter out UNMAP (4) and Secondary (256): -F 260
bwa mem -t 32 ${contigs_index} ${DATA_DIR}/${r1} ${DATA_DIR}/${r2} | samtools view -h -q 10 -F 260 -b - > ${SAMPLENAME}.bam 2>> logs/${SAMPLENAME}.log

# sort and index bam file
samtools sort -@ 32 ${SAMPLENAME}.bam -o ${SAMPLENAME}.sorted.bam 2>> logs/${SAMPLENAME}.log
samtools index ${SAMPLENAME}.sorted.bam > ${SAMPLENAME}.sorted.bam.bai 2>> logs/${SAMPLENAME}.log

echo ""
echo "Job array finished!"
date
```

The __gene catalogue fasta file__ was converted to __gff file__ using a custom Perl script written by Pablo Sánchez named `fasta2gff.pl`. To start counting genes, we will use [HTSeq v0.10.0](https://htseq.readthedocs.io/en/release_0.11.1/) which requires the gff file.

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="fasta2gff"
#SBATCH --time=01-00:00:00
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/log_fasta2gff_%j.out
#SBATCH --error=logs/log_fasta2gff_%j.err

# Modules

# Variables
dir=../gene.catalogue
infile=EMOSE-GC_95id_min250
outfile=${infile}
fasta2gff=fasta2gff.pl

# CMD
perl ${fasta2gff} ${dir}/${infile}.fasta > ${dir}/${outfile}.gff
```

The **gff file conversion** is done!


Once the __gff file__ and the __mapping__ are done, we start counting how many reads mapped to each gene from the catalogue with `htseq-count`.

```bash
#!/bin/sh

#SBATCH --account=emm1
#SBATCH --job-name="htseq"
#SBATCH --time=5-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120000
#SBATCH --output=logs/log_htseq_%j.out
#SBATCH --error=logs/log_htseq_%j.err
#SBATCH --array=1-50%10

# Modules
module load htseq

# Variables
#export PATH=/home/x_sanchep/miniconda2/bin:$PATH
indir=../mapping
sample=$(awk "NR==$SLURM_ARRAY_TASK_ID" sample.list)
outdir=../htseq/read_counts
gff=../htseq/EMOSE.catalogue_95id_clust.gff

# CMD
mkdir ${outdir}
htseq-count -f bam -s no -t CDS --idattr ID -r pos --nonunique all ${indir}/${sample}.sorted.bam ${gff} > ${outdir}/${sample}.htseq-count.txt 2>> logs/htseq-count.${sample}.log
```

After `htseq` finished, we have counts for each fragment for each gene in the catalogue. Now, all fragments should be counted for each sample (reads) but it is already done with the previous version of EMOSE gene catalogue. For downstream steps, the file `fragment.count.per.sample.txt` with all values of this library size will be used to normalize tables.

#### Raw counts table
