#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8GB
#SBATCH -t 12:00:00
#SBATCH --job-name=hisat2
#SBATCH --output=hisat2.%j.out
#SBATCH --error=hisat2.%j.err

##### LAST UPDATED 01/24/2023 Tomas Zelenka #####
#####  2021.11.23 Brian Connolly #####

##### ARGUMENTS AND USAGE #####
usage="$(basename $0) [-h] [-d dir] [-n string] [-i string] [-g file] [-f file] [-r file] -- Program to align RNA-seq BAM files to genome.

Where:
    -h  Show this help text
    -d  Set target directory (directory containing RNA-seq fastq files)
    -n  Set sample base name (i.e. RNA_Th2_WT_1452) (default: Name of directory containing files)
    -i	Set hisat2 index location (default: /share/lab_avram/HPC_Cluster/annotations_genomes/data/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/hisat2_index/transcript/genome_tran)
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
	   printf "Hisat2 index location set to: %s\n" "$OPTARG"
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
module load HISAT2
module load SAMtools
module load picard

##### SET VARIABLES #####
base="${n-${PWD##*/}}"
mkdir "$PWD"/tmp
TMPDIR="$PWD"/tmp

index="/share/lab_avram/HPC_Cluster/annotations_genomes/data/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/hisat2_index/transcript/genome_tran" # Verify index files exist at specified location
[[ -f "$index".1.ht2 ]] || { printf "Hisat2 index not found at %s. Please specify correct index location and rerun %s. Exiting...\n" "$index" "$(basename $0)"; exit 1; }

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
#CREATES HISAT2 DIRECTORY IF IT DOESN'T EXIST
[[ -d hisat2 ]] || mkdir hisat2

#ALIGN TRIMMED READS
if ! [[ -f hisat2/"$base"_filtered.sam || -f hisat2/"$base"_full.bam || -f hisat2/"$base".bam ]]; then
    printf "Aligning trimmed reads to mm10 via hisat2...\n"
    if [[ -f "$forward_reads" && -f "$reverse_reads" ]]; then
	printf "\tPaired-end mode used.\n"
	hisat2 -p 10 --no-unal --time --met-stderr --dta-cufflinks --trim5 8 -x "$index" -1 "$forward_reads" -2 "$reverse_reads" 2> hisat2/"$base"_hisat2.log | sed '/chrM/d;/chrY/d;/random/d;/chrUn/d' > hisat2/"$base"_filtered.sam
    elif [[ -f "$forward_reads" ]]; then
	printf "\tSingle-end mode used.\n"
	hisat2 -p 10 --no-unal --time --met-stderr --dta-cufflinks --trim5 12 --rna-strandness F -x "$index" -U "$forward_reads" 2> hisat2/"$base"_hisat2.log | sed '/chrM/d;/chrY/d;/random/d;/chrUn/d' > hisat2/"$base"_filtered.sam
	# original version was modified to include strandness and to increase trimming of 5 prime from 8 to 12 nt as recommended in the quantseq guidelines as there was random primer with possible mismatches
	# 	hisat2 -p 10 --no-unal --time --met-stderr --dta-cufflinks --trim5 8 -x "$index" -U "$forward_reads" 2> hisat2/"$base"_hisat2.log | sed '/chrM/d;/chrY/d;/random/d;/chrUn/d' > hisat2/"$base"_filtered.sam
    else
	printf "Did not detect valid trimmed fastq files...exiting.\n"
	exit 1
    fi
    printf "Done.\n"
    date
fi

#SORT FILTERED SAM FILE
if ! [[ -f hisat2/"$base"_full.bam || -f hisat2/"$base".bam ]]; then
    printf "Sorting %s_filtered.sam file via samtools..." "$base"
    samtools sort -l 9 -m 3G -o hisat2/"$base"_full.bam -T hisat2/"$SLURM_JOB_ID" -@ 10 hisat2/"$base"_filtered.sam
    printf "Done.\n"
    date
fi

#DEDUPLICATE FULL BAM FILE
if [[ ! -f hisat2/"$base".bam ]]; then
    printf "Deduplicating %s_full.bam via Picardtools..." "$base"
	# Touch output files due to bug in picard
	touch hisat2/"$base".bam
	touch hisat2/"$base"_picard.metrics
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=hisat2/"$base"_full.bam O=hisat2/"$base".bam REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE m=hisat2/"$base"_picard.metrics TMP_DIR="$PWD"/tmp
    printf "Done.\n"
    date
fi

#REMOVE TMP DIR AND INTERMEDIATE FILES
rm -r tmp
rm hisat2/"$base"_filtered.sam
# rm hisat2/"$base"_full.bam

# add bam index
cd hisat2
samtools index "$base".bam
cd ../


printf "hisat2.sh has completed for sample %s...Program terminating.\n" "$base"


#BEGIN FASTQC
if [[ ! -f hisat2/"$base"_fastqc.html ]]; then
    printf "Analyzing %s.bam via FastQC...\n" "$base"
    module load FastQC
    [[ -d "FastQC" ]] || mkdir FastQC #Create FastQC directory if it doesn't already exist
    fastqc hisat2/"$base".bam
    mv -v hisat2/"$base"_fastqc.* FastQC/
    printf "Done.\n"
else
    printf "%s_fastqc.html already exists...skipping FastQC analysis.\n" "$base"
fi

date

