#!/bin/bash

##### LAST UPDATED 01/24/2023 Tomas Zelenka #####
##### 2021.11.23 Brian Connolly #####

##### ARGUMENTS AND USAGE #####
usage="$(basename $0) [-h] [-d dir] [-r file] [-t string] [-c string] -- Program to initiate DESeq2 analysis pipeline on RNA-Seq samples.

Optional Arguments:
    -h  Show this help text
    -d  Set target directory (directory to write analysis files) (default: $PWD)
    -r	Set RNA-Seq sample list file (default: \$HOME/lib/RNA_samples.txt)
    -t	Set test condition (column from RNA-Seq sample list file) (i.e. Genotype or Cell_subset) (default: Genotype)
    -c	Set control group (i.e. WT if test condition is Genotype, or Th0 if test condition is Cell_subset (default: WT)
	
Example usage:
    sbatch --chdir=/share/lab_avram/HPC_Cluster/user/tomas/12_Nina_Pilon-Thomas/CICPT_4749_Lexogen-RNAseq_5-2024 --export=t=\Genotype\,c=ctrl ~/bin/quantseq/DESeq2_loop.sh

    or

    $(basename $0) -d /share/lab_avram/HPC_Cluster/user/kyle.lorentsen/RNA/CD8_RNA_Seq_Run -c Th0 -t Cell_subset"

while getopts ':hd:r:t:c:' option; do
    case "$option" in
    h)  echo "$usage"
	    exit 0
	    ;;
	d)  if [[ -d "$OPTARG" ]]; then
	    d="$OPTARG"
		printf "Working directory set to %s\n" "$OPTARG"
	    else
		printf "%s doesn't exist or is not a valid directory. Please specify a valid directory and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"
		exit 1
	    fi
	    ;;
	r)  if [[ -s "$OPTARG" ]]; then
		r="$OPTARG"
		printf "RNA_Seq sample list file manually set to: %s\n" "$OPTARG"
	    else
		printf "%s is empty or does not exist. Please specify valid RNA_Seq sample list file and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"
	    fi
	    ;;
	t)  t="$OPTARG"
	    printf "Test condition manually set to: %s\n" "$OPTARG"
	    ;;
	c)  c="$OPTARG"
	    printf "Control group manually set to: %s\n" "$OPTARG"
	    ;;
	\?) printf "Invalid option: -%s.\n" "$OPTARG"
	    echo "$usage"
	    exit 1
	    ;;
	:)  printf "Invalid option: -%s requires an argument.\n" "$OPTARG"
	    echo "$usage"
	    exit 1
	    ;;
    esac
done
shift $((OPTIND -1))

##### SET VARIABLES #####
master_dir="${d:-$PWD}"
#TMPDIR="$PWD"/tmp
current_dir="$PWD"

RNA_samples="${r:-${HOME}/lib/RNA_samples.txt}"
[[ -s "$RNA_samples" ]] || { printf "%s is empty or does not exist. Please specify valid RNA_Seq sample list file and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"; exit 1; } # VERIFY RNA SAMPLE LIST FILE IS VALID

condition="${t:-Genotype}"
read -ra condition_list <<<$(head -n1 "$RNA_samples") # GENERATE ARRAY OF POSSIBLE CONDITIONS
[[ "${condition_list[@]}" =~ "$condition" ]] || { printf "%s is not a valid condition.  Please select a condition from the following list: %s. Exiting...\n" "$condition" "${condition_list[*]}"; exit 1; } #VERIFY SELECTED CONDITION IS VALID

baseline="${c:-WT}"
baseline_check=$(awk -F$'\t' -v condition=$condition -v baseline=$baseline 'BEGIN {s=0} 
	NR==1{for(i=1; i<=NF; i++) if($i==condition) {a[i]++; break} }
	{for (i in a) if ($i==baseline) {s++; exit} }
	END {print s}' "$RNA_samples")
[[ $baseline_check > 0 ]] || { printf "There are no RNA samples of type \"%s\" under condition \"%s\". Please select a valid control group and rerun %s. Exiting...\n" "$baseline" "$condition" "$(basename $0)"; exit 1; } #VERIFY THAT AT LEAST ONE SAMPLE UNDER SELECTED CONDITION MATCHES BASELINE TYPE

folders=($( find "$master_dir" -mindepth 1 -maxdepth 1 -type d | grep -v analysis))
job_list=""
bam_files=""


##### BEGIN SCRIPT #####
for folder in "${folders[@]}"; do

    base="${folder##*/}"

    first=""
    if [[ ! -f "$folder"/"$base"_forward_trimmed.fastq.gz ]]; then
	first=$(sbatch --chdir="$folder" $HOME/bin/quantseq/seqtk.sh)
	first=$(echo afterok:$(echo $first | cut -f 4 -d " "))
    fi

    second=""
    if [[ ! -f "$folder"/hisat2/"$base".bam ]]; then
	second=$(sbatch --chdir="$folder" --dependency=$first $HOME/bin/quantseq/hisat2.sh)
	second=$(echo $second | cut -f 4 -d " ")
	if [[ -z "$job_list" ]]; then
        job_list=$(echo afterok:$second)
	else
	    job_list=$(echo "$job_list",$second)
	fi
    fi

    third=""
    if [[ ! -f "$folder"/star/"$base".bam ]]; then
	third=$(sbatch --chdir="$folder" --dependency=$first $HOME/bin/quantseq/star.sh)
	third=$(echo $third | cut -f 4 -d " ")
	if [[ -z "$job_list" ]]; then
        job_list=$(echo afterok:$third)
	else
	    job_list=$(echo "$job_list",$third)
	fi
    fi


# uncomment depending on whether you want to use hisat2 or star mapped data for differential expression analysis
    # if [[ -z "$bam_files" ]]; then
	# bam_files="$folder"/hisat2/"$base".bam
    # else
	# bam_files=$(echo "$bam_files" "$folder"/hisat2/"$base".bam)
    # fi

    if [[ -z "$bam_files" ]]; then
	bam_files="$folder"/star/"$base".bam
    else
	bam_files=$(echo "$bam_files" "$folder"/star/"$base".bam)
    fi


done

printf "Variables sent to DESeq2.sh:\n"
printf "bam_list = %s\n" "$bam_files"
printf "condition = %s\n" "$condition"
printf "baseline = %s\n" "$baseline"
printf "RNA sample list file = %s\n" "$RNA_samples"

sbatch --chdir="$master_dir" --dependency=$job_list --export=bam_list="$bam_files",t="$condition",c="$baseline",r="$RNA_Samples" $HOME/bin/quantseq/DESeq2.sh

