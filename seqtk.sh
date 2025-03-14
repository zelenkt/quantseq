#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --job-name=seqtk
#SBATCH --output seqtk.%j.out
#SBATCH --error seqtk.%j.err

##### LAST UPDATED 01/24/2023 Tomas Zelenka #####
#####  2021.11.23 Brian Connolly #####

date;pwd

# module load glibc/2.14.1
module load seqtk
module load parallel

trim_file=$HOME/lib/IlluminaAdaptorsALL.fa  #TruSeqPE.fa or NexteraPE.fa or IlluminaAdaptorsALL.fa
base="${PWD##*/}"

###added by Tomas
#DETERMINE IF SAMPLE IS PAIRED END OR SINGLE END
if [[ -z $is_paired_end ]]; then
    if test -n "$(find -L "$PWD" -maxdepth 1 -name '*_R2_001.fastq.gz' -o -iname '*_Reverse.fastq.gz')"; then
	is_paired_end=true
	printf "%s detected to be paired end.\n" "$base"
    else
	is_paired_end=false
	printf "%s detected to be single end.\n" "$base"
    fi
fi

export _JAVA_OPTIONS="-Xmx4g"



## CHECK IF READS ARE CONCATENATED ##
if [[ ! -f "$base"_forward.fastq.gz && "$(find . -name "*R1_001.fastq.gz" -type f -print -quit)" ]]; then # Checking if (1) concatenated forward read file already exists, and, if not, (2) if the unconcatenated forward read files are present.
	cat *R1_001.fastq.gz > "$base"_forward.fastq.gz # If this doesn't result in the creation of the concatenated file, make sure the unconcatenated fastq files follow standard Illumina naming convention (e.g. "$sample_name"_S4_L001_R1_001.fastq.gz)
fi

if [[ ! -f "$base"_reverse.fastq.gz && "$(find . -name "*R2_001.fastq.gz" -type f -print -quit)" ]]; then
	cat *R2_001.fastq.gz > "$base"_reverse.fastq.gz
fi

## TRIM READS WITH SEQTK AND COMPRESS OUTPUT ##
if [[ -f "$base"_forward.fastq.gz && ! -f "$base"_forward_seqtk.fastq.gz ]]; then # If the concatenated forward read file exists, proceed
    seqtk trimfq "$base"_forward.fastq.gz | parallel --pipe -k --block 10M gzip -9 > "$base"_forward_seqtk.fastq.gz
fi

if [[ -f "$base"_reverse.fastq.gz && ! -f "$base"_reverse_seqtk.fastq.gz ]]; then # If the concatenated reverse read file exists, proceed. Will not run if no reverse file was found (i.e. run was Single End)
    seqtk trimfq "$base"_reverse.fastq.gz | parallel --pipe -k --block 10M gzip -9 > "$base"_reverse_seqtk.fastq.gz
fi

echo "seqtk trimming complete...program terminating"


### added by Tomas
#BEGIN TRIMMING
if [[ ! -f "$base"_1P.fastq.gz ]]; then
    printf "Trimming read via Trimmomatic..."
    module load Trimmomatic/0.39-Java-1.8
    if $is_paired_end; then
	trimmomatic PE -threads 6 -phred33 -baseout "$base".fastq.gz "$base"_forward_seqtk.fastq.gz "$base"_reverse_seqtk.fastq.gz ILLUMINACLIP:"$trim_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> "$base".trim.log
    else
	trimmomatic SE -threads 6 -phred33 "$base"_forward_seqtk.fastq.gz "$base"_1P.fastq.gz ILLUMINACLIP:"$trim_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> "$base".trim.log
    fi
    printf "Done.\n"
else
    printf "Trimmed fastq files already detected...skipping trimmomatic.\n"
fi

[[ -d "FastQC" ]] || mkdir FastQC #Create FastQC directory if it doesn't already exist
#BEGIN FASTQC
if [[ ! -f FastQC/"$base"*_fastqc.html ]]; then
    printf "Analyzing %s_001.fastq.gz via FastQC...\n" "$base"
    module load FastQC
    # fastqc "$base".bam -o FastQC
    fastqc *R1_001.fastq.gz
    fastqc *R2_001.fastq.gz
    mv -v *_fastqc.* FastQC/
    printf "Done.\n"
else
    printf "%s_fastqc.html already exists...skipping FastQC analysis.\n" "$base"
fi
### /////added by Tomas


date

