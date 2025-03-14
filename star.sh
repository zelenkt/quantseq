#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8GB
#SBATCH -t 12:00:00
#SBATCH --job-name=star
#SBATCH --output=star.%j.out
#SBATCH --error=star.%j.err

##### LAST UPDATED 01/24/2023 Tomas Zelenka #####
#####  2021.11.23 Brian Connolly #####

##### ARGUMENTS AND USAGE #####
usage="$(basename $0) [-h] [-d dir] [-n string] [-i string] [-g file] [-f file] [-r file] -- Program to align RNA-seq BAM files to genome.

Where:
    -h  Show this help text
    -d  Set target directory (directory containing RNA-seq fastq files)
    -n  Set sample base name (i.e. RNA_Th2_WT_1452) (default: Name of directory containing files)
    -i	Set star index location (default: /share/lab_avram/HPC_Cluster/annotations_genomes/data/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/star_index/transcript/genome_tran)
    -g	Set reference GTF file (default: /share/lab_avram/HPC_Cluster/annotations_genomes/alias/mm10/ensembl_gtf/default/mm10.gtf.gz)
    -f  Set forward fastq file (default: \$base_1P.fastq.gz)
    -r	Set reverse fastq file (default: \$base_2P.fastq.gz)

Example usage:
    sbatch --workdir=/directory/containing/fastq/files --export=n=RNA_Th2_WT_1452 $HOME/bin/quantseq/$(basename $0)

    or

    $(basename $0) -d /directory/containing/fastq/files -n RNA_Th2_WT_1452"

while getopts ':hd:n:i:g:f:r:' option; do
    case "$option" in
	h) echo "$usage"
	   exit 0
	   ;;
	d) if [[ -d "$OPTARG" ]]; then
	       cd "$OPTARG"
	       printf "Working directory set to %s\n" "$OPTARG"
	   else
	       printf "%s doesn't exist or isn't a valid directory.\n" "$OPTARG"
	       exit 1
	   fi
	   ;;
	n) n="$OPTARG"
	   printf "Sample base name set to %s.\n" "$OPTARG"
	   ;;
       \?) printf "Invalid option: -%s.\n" "$OPTARG"
	   echo "$usage"
	   exit 1
	   ;;
	i) i="$OPTARG"
	   printf "STAR index location set to: %s\n" "$OPTARG"
	   ;;
	g) g="$OPTARG"
	   printf "Reference GTF file set to: %s\n" "$OPTARG"
	   ;;
	f) f="$OPTARG"
	   printf "Forward fastq file set to: %s\n" "$OPTARG"
	   ;;
	r) r="$OPTARG"
	   printf "Reverse fastq file set to: %s\n" "$OPTARG"
	   ;;
	:) printf "Invalid option: -%s requires an argument.\n" "$OPTARG"
	   echo "$usage"
	   exit 1
	   ;;
    esac
done
shift $((OPTIND -1))

date;pwd

##### LOAD MODULES #####
# module load STAR/2.7.7a-GCC-10.2.0 # this version does not support the genome index I have from refgenie
module load STAR/2.7.2b-GCC-8.3.0
module load SAMtools
module load picard

##### SET VARIABLES #####
base="${n-${PWD##*/}}"
mkdir "$PWD"/tmp
TMPDIR="$PWD"/tmp


index="/share/lab_avram/HPC_Cluster/annotations_genomes/data/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/star_index/default" # Verify index files exist at specified location
[[ -f "$index"/SAindex ]] || { printf "Star index not found at %s. Please specify correct index location and rerun %s. Exiting...\n" "$index" "$(basename $0)"; exit 1; }

gtf="${g-/share/lab_avram/HPC_Cluster/annotations_genomes/alias/mm10/ensembl_gtf/default/mm10.gtf.gz}" # Verify GTF file exists at specified location
[[ -f "$gtf" ]] || { printf "%s not found. Please specify correct gtf location and rerun %s. Exiting...\n" "$gtf" "$(basename $0)"; exit 1; }

forward_reads="${f:-${base}_1P.fastq.gz}"
reverse_reads="${r:-${base}_2P.fastq.gz}"

##### SET JAVA&TMPDIR OPTIONS #####
export _JAVA_OPTIONS="-Xmx4g"
#CREATE LOCAL TMP DIR IF IT DOESN'T ALREADY EXIST
[[ -d tmp ]] || mkdir tmp
export TMPDIR="$PWD"/tmp

##### BEGIN SCRIPT #####
#CREATES STAR DIRECTORY IF IT DOESN'T EXIST
[[ -d star ]] || mkdir star

#ALIGN TRIMMED READS
if ! [[ -f star/"$base"_filtered.sam || -f star/"$base"_full.bam || -f star/"$base".bam ]]; then
    printf "Aligning trimmed reads to mm10 via star...\n"
    if [[ -f "$forward_reads" && -f "$reverse_reads" ]]; then
	printf "\tPaired-end mode used.\n"
	printf "\tMODIFY COMMAND BY ADDING THE RIGHT PARAMETERS.\n"
    elif [[ -f "$forward_reads" ]]; then
	printf "\tSingle-end mode used.\n"
	STAR --genomeDir "$index" --clip5pNbases 12 --clip3pAdapterSeq polyA --runThreadN 10 --readFilesIn "$forward_reads" --outFileNamePrefix star/"$base"_ --outSAMtype SAM --readFilesCommand zcat
    #| sed '/chrM/d;/chrY/d;/random/d;/chrUn/d' S13_ctrl_1P_Aligned.out.bam > star/"$base"_filtered.sam
    sed '/chrM/d;/chrY/d;/random/d;/chrUn/d' star/"$base"_Aligned.out.sam > star/"$base"_filtered.sam

# # testing file - RUNNING 6-5-24
# STAR --genomeDir /share/lab_avram/HPC_Cluster/annotations_genomes/data/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/star_index/default --clip5pNbases 12 --clip3pAdapterSeq polyA --runThreadN 10 --readFilesIn S13_ctrl_1P.fastq.gz  --outFileNamePrefix star/S13_ctrl_1P_ --outSAMtype SAM --readFilesCommand zcat; sed '/chrM/d;/chrY/d;/random/d;/chrUn/d' star/S13_ctrl_1P_Aligned.out.sam > star/S13_ctrl_1P_filtered.sam

# # THIS DOES NOT WORK
# STAR --genomeDir /share/lab_avram/HPC_Cluster/annotations_genomes/data/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/star_index/default --clip5pNbases 12 --clip3pAdapterSeq polyA --runThreadN 10 --readFilesIn S13_ctrl_1P.fastq.gz  --outFileNamePrefix star/S13_ctrl_1P_ --outSAMtype BAM Unsorted --quantMode GeneCounts --readFilesCommand zcat --sjdbGTFfile /share/lab_avram/HPC_Cluster/annotations_genomes/data/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/ensembl_gtf/default/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1.gtf


# /share/lab_avram/HPC_Cluster/annotations_genomes/alias/mm10/ensembl_gtf/default/mm10.gtf.gz


# quantMode GeneCounts
# --genomeLoad LoadAndRemove
# https://groups.google.com/g/rna-star/c/oGvChGG2gto/m/UQrP3c8yDQAJ >>> 1. For the stranded RNA-Seq data, I don't need to specify the strandness information for STAR. But I need to specify the strandness information for featureCount by using strandSpecific. Am I right?
# You are right, use -s 0/1/2 option from featureCounts

#  --soloBarcodeReadLength 0
# # all these are probably only for single cell star analysis
# --soloStrand Forward
# FASTQ="fastq/test.fastq.gz"
# INDEX="index"
# CORES=12
# PASSLIST="passlist.txt"
# OUTDIR="results"

# echo "TATA" > "${PASSLIST}"

# https://github.com/alexdobin/STAR/issues/1232
# STAR \
#   --readFilesIn "${FASTQ}" \
#   --genomeDir "${INDEX}" \
#   --runThreadN "${CORES}" \
#   --readFilesCommand "gzip -cd" \
#   --clip3pAdapterSeq polyA \
#   --outFilterType BySJout \
#   --outFilterMultimapNmax 20 \
#   --alignSJoverhangMin 8 \
#   --alignSJDBoverhangMin 1 \
#   --outFilterMismatchNmax 999 \
#   --outFilterMismatchNoverLmax 0.1 \
#   --alignIntronMin 20 \
#   --alignIntronMax 1000000 \
#   --alignMatesGapMax 1000000 \
#   --outSAMattributes NH HI NM MD \
#   --outSAMtype BAM SortedByCoordinate \
#   --soloBarcodeMate 1 \
#   --soloType CB_UMI_Simple \
#   --soloCBstart 7 \
#   --soloCBlen 4 \
#   --soloUMIstart 1 \
#   --soloUMIlen 6 \
#   --soloCellFilter None \
#   --soloCBwhitelist "${PASSLIST}" \
#   --outFileNamePrefix "${OUTDIR}/" \
#   --soloMultiMappers Uniform EM \
#   --soloStrand Forward





	# original version was modified to include strandness and to increase trimming of 5 prime from 8 to 12 nt as recommended in the quantseq guidelines as there was random primer with possible mismatches
	# 	hisat2 -p 10 --no-unal --time --met-stderr --dta-cufflinks --trim5 8 -x "$index" -U "$forward_reads" 2> star/"$base"_star.log | sed '/chrM/d;/chrY/d;/random/d;/chrUn/d' > star/"$base"_filtered.sam
    else
	printf "Did not detect valid trimmed fastq files...exiting.\n"
	exit 1
    fi
    printf "Done.\n"
    date
fi

# samtools sort -l 9 -m 3G -o S13_ctrl_1P_full.bam -@ 10 S13_ctrl_1P_filtered.sam
# samtools index S13_ctrl_1P_full.bam
#SORT FILTERED SAM FILE
if ! [[ -f star/"$base"_full.bam || -f star/"$base".bam ]]; then
    printf "Sorting %s_filtered.sam file via samtools..." "$base"
    samtools sort -l 9 -m 3G -o star/"$base"_full.bam -T star/"$SLURM_JOB_ID" -@ 10 star/"$base"_filtered.sam
    printf "Done.\n"
    date
fi

#DEDUPLICATE FULL BAM FILE
if [[ ! -f star/"$base".bam ]]; then
    printf "Deduplicating %s_full.bam via Picardtools..." "$base"
	# Touch output files due to bug in picard
	touch star/"$base".bam
	touch star/"$base"_picard.metrics
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=star/"$base"_full.bam O=star/"$base".bam REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE m=star/"$base"_picard.metrics TMP_DIR="$PWD"/tmp
    printf "Done.\n"
    date
fi

#REMOVE TMP DIR AND INTERMEDIATE FILES
rm -r tmp
rm star/"$base"_Aligned.out.sam
rm star/"$base"_filtered.sam
# rm star/"$base"_full.bam


# add bam index
cd star
samtools index "$base".bam
cd ../

printf "star.sh has completed for sample %s...Program terminating.\n" "$base"


#BEGIN FASTQC
if [[ ! -f FastQC_star/"$base"_fastqc.html ]]; then
    printf "Analyzing %s.bam via FastQC...\n" "$base"
    module load FastQC
    [[ -d "FastQC_star" ]] || mkdir FastQC_star #Create FastQC directory if it doesn't already exist
    fastqc star/"$base".bam
    mv -v star/"$base"_fastqc.* FastQC_star/
    printf "Done.\n"
else
    printf "%s_fastqc.html already exists...skipping FastQC analysis.\n" "$base"
fi

date

