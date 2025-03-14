#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=8:00:00
#SBATCH --job-name=DESeq2
#SBATCH --output=DESEq2.%j.out
#SBATCH --error=DESEq2.%j.err

##### LAST UPDATED 01/24/2023 Tomas Zelenka #####
#####  2021.11.23 Brian Connolly #####

##### ARGUMENTS AND USAGE #####
usage="$(basename $0) [-h] [-d dir] [-r file] [-g file] [-t string] [-c string] bam_list -- Program to run DESeq2 differential analysis on RNA-Seq bam files.

Where: bam_list is a space-separated list of RNA-Seq bam files (e.g. RNA_WT_Th0_sample1.bam RNA_WT_Th0_sample2.bam RNA_WT_Th2_sample1.bam RNA_WT_Th2_sample2.bam)

Optional Arguments:
    -h  Show this help text
    -d  Set target directory (directory to write analysis files) (default: $PWD/analysis/DESeq2)
    -r  Set RNA-Seq sample list file (default: \$HOME/lib/RNA_samples.txt)
    -g  Set annotated GTF file (default: \$HOME/lib/refgenie_ensembl_mm10.gtf)
    -t  Set test condition (column from RNA-Seq sample list file) (i.e. Genotype or Cell_subset) (default: Genotype)
    -c  Set control group (i.e. WT if test condition is Genotype, or Th0 if test condition is Cell_subset (default: WT)


Example usage:
    sbatch --chdir=/share/lab_avram/HPC_Cluster/user/kyle.lorentsen/RNA/CD8_RNA_Seq_Run --export=t=Cell_subset,c=Th0,bam_list=\"RNA_WT_Th0_sample1.bam RNA_WT_Th0_sample2.bam RNA_WT_Th2_sample1.bam RNA_WT_Th2_sample2.bam\" $HOME/bin/quantseq/$(basename $0)

    or

    $(basename $0) -d /share/lab_avram/HPC_Cluster/user/kyle.lorentsen/RNA/CD8_RNA_Seq_Run -c Th0 -t Cell_subset RNA_WT_Th0_sample1.bam RNA_WT_Th0_sample2.bam RNA_WT_Th2_sample1.bam RNA_WT_Th2_sample2.bam"

while getopts ':hd:r:g:t:c:' option; do
    case "$option" in
    h)  echo "$usage"
        exit 0
        ;;
	d)  if [[ -d "$OPTARG" ]]; then
            d="$OPTARG"
            printf "Target directory manually set to %s\n" "$OPTARG"
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
    g)  g="$OPTARG"
        printf "Annotated GTF file manually set to %s\n" "$OPTARG"
        ;;
    t)  t="$OPTARG"
        printf "Test condition manually set to: %s\n" "$OPTARG"
        ;;
    c)  c="$OPTARG"
        printf "Control group manually set to: %s\n" "$OPTARG"
        ;;
   \?)  printf "Invalid option: -%s.\n" "$OPTARG"
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

date;pwd

##### LOAD MODULES #####
module load R/4.1.0-foss-2020b

##### Use Avram Lab's shared R library #####
export R_LIBS="/share/lab_avram/HPC_Cluster/share/R/x86_64-pc-linux-gnu-library/4.1"

##### SET VARIABLES #####
target_dir="${d:-${PWD}/analysis/DESeq2}"
if [[ "$target_dir" == "$PWD/analysis/DESeq2" && ! -d "$target_dir" ]]; then # IF TARGET DIR DOES NOT EXIST
    mkdir -p "$target_dir"
elif [[ ! -d "$target_dir" ]]; then # THIS CAN ONLY BE TRUE IF SOMEONE CALLED THIS FUNCTION THROUGH SBATCH AND DIRECTORY EXPORTED AN INVALID "d" variable. (e.g. sbatch --export=d=/some/invalid/dir DESeq2.sh)
    printf "%s doesn't exist or is not a valid directory. Please specify a valid directory and rerun %s. Exiting...\n" "$target_dir" "$(basename $0)"
    exit 1
fi

gtf="${g-${HOME}/lib/refgenie_ensembl_mm10.gtf}"
[[ -f "$gtf" ]] || { printf "%s does not exist.  Please specify valid annotated GTF file and rerun %s. Exiting...\n" "$gtf" "$(basename $0)"; exit 1; } # VERIFY ANNOTATED GTF FILE IS VALID

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

#ADD TEST TO CHECK IF AT LEAST ONE SAMPLE (BUT NOT ALL SAMPLES) HAVE BASELINE AS VALUE UNDER CONDITION
bam_list="${bam_list:-$@}"
if [[ -z "$bam_list" ]]; then
    printf "No list of BAM files detected. Please provide valid list of bam files and rerun %s. Exiting...\n" "$(basename $0)"
    exit 1
else
    IFS=' ' read -ra bam_array<<<"$bam_list"
    [[ "${#bam_array[@]}" > 1 ]] || { printf "At least two BAM files must be provided. Please update BAM file list and rerun %s. Exiting...\n" "$(basename $0)"; exit 1; } # VERIFY THAT AT LEAST TWO BAM FILES HAVE BEEN PROVIDED
    for bam in "${bam_array[@]}"; do # VERIFY THAT ALL THE GIVEN BAM FILES EXIST
        [[ -f "$bam" ]] || { printf "%s can not be found.  Please give valid bam file list and rerun %s. Exiting...\n" "$bam" "$(basename $0)"; exit 1; } 
    done

    uniqueBams=$(printf "%s\n" "${bam_array[@]}" | awk '!($0 in seen){seen[$0];c++} END {print c}')
    [[ $uniqueBams == "${#bam_array[@]}" ]] || { printf "Duplicate BAM files detected in given list.  Each BAM file can only be given once, please remove duplicates and rerun $0. Exiting...\n" "$(basename $0)"; exit 1; } # VERIFY THAT THERE ARE NO DUPLICATE BAM FILES IN GIVEN BAM LIST

fi

##### BEGIN SCRIPT #####
printf "Condition = %s\n" "$condition"
printf "Baseline = %s\n" "$baseline"
printf "RNA Sample list file = %s\n" "$RNA_samples"
printf "Annotated GTF file = %s\n" "$gtf"
printf "Bam_list = %s\n" "$bam_list"

( cd "$target_dir" && Rscript "$HOME"/bin/quantseq/DESeq2.R -c "$condition" -b "$baseline" -r "$RNA_samples" -g "$gtf" -l "$bam_list" ) # RUN DESEQ2.R WITHIN SUBSHELL TO CHANGE DIRECTORY TO ANALYSIS DIR

# ( cd "$target_dir" && Rscript "$HOME"/bin/quantseq/annotate_genes.R -i DESeq2_analysis.txt -o DESeq2_analysis_annotated.txt ) # RUN ANNOTATE_GENES_DESE2.R WITHIN SUBSHELL TO CHANGE DIRECTORY TO ANALYSIS DIR




# added by Tomas to gene annotate the resulting files
if [[ ! -f ${HOME}/lib/geneNames_refgenie_ensembl_mm10.gtf ]]; then
    printf "Extracting gene names from GTF file...\n" "$base"
        
    awk '$3=="gene"{print $0}' "$gtf" | awk -F'gene_id "|gene_name "|"; ' -v OFS="\t" '{print $2,$5}' | sort | uniq > ${HOME}/lib/geneNames_refgenie_ensembl_mm10.gtf

    printf "Done.\n"
else
    printf "Gene names extracted from the GTF annotation file already exist\n" "$base"
fi


# prepare annotation files
# awk -F'\t' -v OFS="\t" 'NR>1{print $2,$3,$4}' humanALLEnsemblIDannotated.txt > humanALLEnsemblIDannotated_forHPC.txt 
# awk -F'\t' -v OFS="\t" 'NR>1{print $2,$3,$4}' mouseALLEnsemblIDannotated.txt > mouseALLEnsemblIDannotated_forHPC.txt
# awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' humanALLEnsemblIDannotated_forHPC.txt geneNames_refgenie_ensembl_hg38.gtf | awk -F'\t' -v OFS="\t" '{print $1,$2,$5}' > geneNames_refgenie_ensembl_hg38_description.gtf
# awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' mouseALLEnsemblIDannotated_forHPC.txt geneNames_refgenie_ensembl_mm10.gtf | awk -F'\t' -v OFS="\t" '{print $1,$2,$5}' > geneNames_refgenie_ensembl_mm10_description.gtf


# # to normalize to gene length I add to the pipeline the lines to extract the longest isoform for each gene (sum of exons)
# awk -v OFS="\t" '$3=="exon" {print $14,  $10, $20,  $5-$4 }' refgenie_ensembl_mm10.gtf | tr -d '"|;' > tmpIDs449848423; awk -F'\t' -v OFS="\t" '{ seen[$1] += $4 } END { for (i in seen) print i, seen[i] }' tmpIDs449848423 > tmpLengths198489

# # note that this needs to be sorted from smallest to largest length because the next "vlookup" code takes the last lane
# awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' tmpIDs449848423 tmpLengths198489 | sort -k4,4 -k2n,2 | awk -F'\t' -v OFS="\t" '{print $4, $5, $1, $2}'  > tmp_sorted_iso_lengths

# awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' tmp_sorted_iso_lengths geneNames_refgenie_ensembl_mm10_description.gtf | awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $7}' > geneNames_refgenie_ensembl_mm10_description_longestIsoLength.gtf; rm tmp_sorted_iso_lengths tmpLengths198489 tmpIDs449848423

# # this is for human data
# awk -v OFS="\t" '$3=="exon" {print $14,  $10, $20,  $5-$4 }' refgenie_ensembl_hg38.gtf | tr -d '"|;' > tmpIDs449848423; awk -F'\t' -v OFS="\t" '{ seen[$1] += $4 } END { for (i in seen) print i, seen[i] }' tmpIDs449848423 > tmpLengths198489

# # note that this needs to be sorted from smallest to largest length because the next "vlookup" code takes the last lane
# awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' tmpIDs449848423 tmpLengths198489 | sort -k4,4 -k2n,2 | awk -F'\t' -v OFS="\t" '{print $4, $5, $1, $2}'  > tmp_sorted_iso_lengths

# awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' tmp_sorted_iso_lengths geneNames_refgenie_ensembl_hg38_description.gtf | awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $7}' > geneNames_refgenie_ensembl_hg38_description_longestIsoLength.gtf; rm tmp_sorted_iso_lengths tmpLengths198489 tmpIDs449848423





geneNames="${g-${HOME}/lib/geneNames_refgenie_ensembl_mm10_description_longestIsoLength.gtf}"

#Add gene name annotation to DESeq2 files
if [[ ! -f analysis/DESeq2/gene_annotated* ]]; then
    printf "Gene annotating DESeq2 output files...\n" "$base"

    cat analysis/DESeq2/DESeq2_analysis.txt | tr -d '"'| awk 'NR>1{print $0}' > analysis/DESeq2/tmp1315; awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $geneNames analysis/DESeq2/tmp1315 | awk -F'\t' -v OFS="\t" '{print $9,$1,$2,$3,$4,$5,$6,$7,$10}' > analysis/DESeq2/gene_annotated_DESeq2_analysis.txt; rm analysis/DESeq2/tmp1315


    # cat analysis/DESeq2/DESeq2_analysis.txt | tr -d '"'| awk 'NR>1{print $0}' > analysis/DESeq2/tmp1315; awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' /Volumes/lab_avram/HPC_Cluster/user/tomas/lib/geneNames_refgenie_ensembl_mm10_description_longestIsoLength.gtf analysis/DESeq2/tmp1315 | awk -F'\t' -v OFS="\t" '{print $9,$1,$2,$3,$4,$5,$6,$7,$10}' > analysis/DESeq2/gene_annotated_DESeq2_analysis.txt; rm analysis/DESeq2/tmp1315

    # cat analysis/DESeq2/DESeq2_analysis.txt | tr -d '"'| awk 'NR>1{print $0}' > analysis/DESeq2/tmp1315; awk 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $geneNames analysis/DESeq2/tmp1315 | awk -v OFS="\t" '{print $9,$1,$2,$3,$4,$5,$6,$7}' > analysis/DESeq2/gene_annotated_DESeq2_analysis.txt; rm analysis/DESeq2/tmp1315

    cat analysis/DESeq2/DESeq2_normalized_counts.txt | tr -d '"'| awk 'NR>1{print $0}' > analysis/DESeq2/tmp131662

    colCount=$(awk '{print NF; exit}' analysis/DESeq2/tmp131662)

    awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $geneNames analysis/DESeq2/tmp131662 | awk -F'\t' -v OFS="\t" -v colCount="$colCount" '{print $(colCount+2),$0}' > analysis/DESeq2/gene_annotated_DESeq2_normalized_counts.txt; rm analysis/DESeq2/tmp131662

    # awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' /Volumes/lab_avram/HPC_Cluster/user/tomas/lib/geneNames_refgenie_ensembl_mm10_description_longestIsoLength.gtf analysis/DESeq2/tmp131662 | awk -F'\t' -v OFS="\t" -v colCount="5" '{print $(colCount+2),$0}' > analysis/DESeq2/gene_annotated_DESeq2_normalized_counts.txt; rm analysis/DESeq2/tmp131662

    # awk 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $geneNames analysis/DESeq2/tmp131662 | awk -v OFS="\t" -v colCount="$colCount" '{print $(colCount+2),$0}' | awk 'NF{NF-=2}1' > analysis/DESeq2/gene_annotated_DESeq2_normalized_counts2.txt; rm analysis/DESeq2/tmp131662

    cat analysis/DESeq2/Raw_counts.txt | tr -d '"'| awk -F'\t' -v OFS="\t" 'NR>1{print $0}' > analysis/DESeq2/tmp13199; awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $geneNames analysis/DESeq2/tmp13199 | awk -F'\t' -v OFS="\t" -v colCount="$colCount" '{print $(colCount+2),$0}' > analysis/DESeq2/gene_annotated_Raw_counts.txt; rm analysis/DESeq2/tmp13199


    # cat analysis/DESeq2/Raw_counts.txt | tr -d '"'| awk -F'\t' -v OFS="\t" 'NR>1{print $0}' > analysis/DESeq2/tmp13199; awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' /Volumes/lab_avram/HPC_Cluster/user/tomas/lib/geneNames_refgenie_ensembl_mm10_description_longestIsoLength.gtf analysis/DESeq2/tmp13199 | awk -F'\t' -v OFS="\t" -v colCount="5" '{print $(colCount+2),$0}' > analysis/DESeq2/gene_annotated_Raw_counts.txt; rm analysis/DESeq2/tmp13199


    # cat analysis/DESeq2/Raw_counts.txt | tr -d '"'| awk 'NR>1{print $0}' > analysis/DESeq2/tmp13199; awk 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' $geneNames analysis/DESeq2/tmp13199 | awk -v OFS="\t"  -v colCount="$colCount" '{print $(colCount+2),$0}' | awk 'NF{NF-=2}1' > analysis/DESeq2/gene_annotated_Raw_counts.txt; rm analysis/DESeq2/tmp13199

    # for local testing
    #awk -F'"|\t' -v OFS="\t" 'NR>1{print $2,$4,$5,$6,$7,$8,$9}' DESeq2_analysis.txt > tmp1315; awk 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"zzzzzzzzzzzz#NA"} }' ../../../lib/geneNames_gencode.vM25.annotation_GRCm38p6.gtf tmp1315 | awk -v OFS="\t" '{print $9,$1,$2,$3,$4,$5,$6,$7}' | head

    #awk -F'"|\t' -v OFS="\t" 'NR>1{print $2,$4,$5,$6,$7}' DESeq2_normalized_counts.txt > tmp1315; awk 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"zzzzzzzzzzzz#NA"} }' ../../../lib/geneNames_gencode.vM25.annotation_GRCm38p6.gtf tmp1315 | awk -v OFS="\t" '{print $7,$1,$2,$3,$4,$5}' | head





    # # this is to normalize the DESeq2 normalized counts using also gene-length - using the sum of the exons for the longest isoform
    # # if we have 4 replicates
    # awk -F'\t' -v OFS="\t" '{print $1,$2,$3/$10*1000, $4/$10*1000, $5/$10*1000, $6/$10*1000, $9, $10}' analysis/DESeq2/gene_annotated_DESeq2_normalized_counts.txt > analysis/DESeq2/gene_annotated_gene-length-normalized-DESeq2_normalized_counts.txt

    # # if we have 6 replicates
    awk -F'\t' -v OFS="\t" '{print $1,$2,$3/$12*1000, $4/$12*1000, $5/$12*1000, $6/$12*1000,$7/$12*1000,$8/$12*1000, $11, $12}' analysis/DESeq2/gene_annotated_DESeq2_normalized_counts.txt > analysis/DESeq2/gene_annotated_gene-length-normalized-DESeq2_normalized_counts.txt
    

    printf "Done.\n"
else
    printf "Gene annotated DESeq2 files already exist\n" "$base"
fi
### /////added by Tomas


date

