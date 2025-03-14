#!/usr/bin/Rscript

##### LAST UPDATED 01/24/2023 Tomas Zelenka #####
##### 2019.12.16 EYH #####

##### LOAD LIBRARIES #####
suppressMessages(
if  (!require("DESeq2")) {
    install.packages("DESeq2", dependencies = TRUE)
    library(DESeq2)
    }
)

suppressMessages(
if  (!require("Rsubread")) {
    install.packages("Rsubread", dependencies = TRUE)
    library(Rsubread)
    }
)

suppressMessages(
if  (!require("getopt")) {
    install.packages("getopt", dependencies = TRUE)
    library(getopt)
    }
)

##### GET ARGUMENTS AND SET VARIABLES #####

spec = matrix(c(
    'help'      , 'h', 0, "logical"  , "Prints usage text",
    'condition' , 'c', 1, "character", "Test condition (default: 'Genotype')",
    'baseline'  , 'b', 1, "character", "Control group identifier (default: 'WT')",
    'rna'       , 'r', 1, "character", "RNA-Seq sample list file (default: '$HOME/lib/RNA_samples.txt')",
    'gtf'       , 'g', 1, "character", "Annotated GTF file (default: '$HOME/lib/refgenie_ensembl_mm10.gtf')",
	'bams'      , 'l', 1, "character", "Space-separated list of bam files (e.g. 'Bam1.bam Bam2.bam') [REQUIRED]"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

# CALL HELP MESSAGE IF NEEDED
if (!is.null(opt$help) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# REQUIRED ARGS
if (is.null(opt$bams) ) {
   cat("bams argument is required. Please rerun with bam list provided.\n")
   q(status=1)
}

# OPTIONAL ARGS
if ( is.null(opt$condition ) ) { opt$condition = "Genotype"                                             }
if ( is.null(opt$baseline  ) ) { opt$baseline  = "WT"                                                   }
if ( is.null(opt$rna       ) ) { opt$rna       = paste0(path.expand("~"),"/lib/RNA_samples.txt")        }
if ( is.null(opt$gtf       ) ) { opt$gtf       = paste0(path.expand("~"),"/lib/refgenie_ensembl_mm10.gtf") }

# SPLIT BAM-LIST INTO CHARACTER VECTORres <
bam_list<-unlist(strsplit(opt$bams, ' '))

##### BEGIN SCRIPT #####
# GENERATE RAW COUNTS FILE IF IT DOESN'T ALREADY EXIST
# TZ included also strand-specific parameter
if (!file.exists("Raw_counts.txt")) {
    myFeatureCounts <- featureCounts(bam_list,
                       annot.ext=opt$gtf,
                       isGTFAnnotationFile=TRUE,
                       isPairedEnd=FALSE,
                       strandSpecific=1,
                       nthreads=4)

    myCounts <- myFeatureCounts$counts

    ## GET BASE NAMES
        ## featureCounts converts "/ufrc/user/sample1.bam" to "X.ufrc.user.sample1.bam".  This command extracts "sample1"
        colnames(myCounts) <- gsub('\\.', '_', x=colnames(myCounts))
        colnames(myCounts) <- gsub('_bam', '', x=colnames(myCounts))

        ## WRITE COUNT DATA TO FILE		
        write.table(myCounts, file = "Raw_counts.txt", sep = "\t")

} else {
	cat("Raw_counts.txt file already detected. Skipping featureCounts.\n")
    myCounts <- as.matrix(read.table("Raw_counts.txt", sep = "\t", row.names=1, header=1))
}

## IMPORT SAMPLE INFO TABLE
    sampleInfo <- read.table(opt$rna, sep = "\t", row.names=1, header=1)

## REMOVE UNNECCESARY SAMPLEINFO SAMPLE DATA AND REORDER TO MATCH "myCounts"
    sortedIndices <- match(colnames(myCounts), rownames(sampleInfo))
    sampleInfo <- sampleInfo[sortedIndices,,drop=FALSE]

## INPUT DATA INTO DESeq2
    dds <- DESeqDataSetFromMatrix(countData = myCounts,
                                  colData = sampleInfo,
                                  design = as.formula(sprintf("~ %s",opt$condition)))

## SPECIFY REFERENCE LEVEL(GENOTYPE)
    dds[[opt$condition]] <- relevel(dds[[opt$condition]], ref=opt$baseline)

## MINIMAL DATA FILTERING
    dds <- dds[ rowSums(counts(dds)) > 1, ]

## PERFORM DIFFERENTIAL ANALYSIS
    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered <- res[order(res$padj),]
    normalizedCounts <- counts(dds, normalized=TRUE)

## EXPORT FILES
    write.table(resOrdered, file = "DESeq2_analysis.txt", sep = "\t") # Differential analysis file
    write.table(normalizedCounts, file = "DESeq2_normalized_counts.txt", sep = "\t") # Normalized counts file

