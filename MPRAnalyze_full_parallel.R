#install.packages('devtools')

# Load packages

#library("devtools")

# install the correct version of MPRAnalyze
#devtools::install_github("YosefLab/MPRAnalyze", force=TRUE)

library("rSubmitter")
library("MPRAnalyze")
sessionInfo()

nchunks <- 300
chunksize <- 50

# Read in data
 
annot = read.csv(
	"/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/barcode_counts/NPC_col_annot.txt", 
	sep = '\t', 
	row.names = 1
	)
dcounts_table = read.table(
	"/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/barcode_counts/NPC_DNA_MPRAnalyze.txt", 
	header=T, 
	sep='\t', 
	stringsAsFactors=F, 
	quote="", 
	row.names=1
	)
rcounts_table = read.table(
	"/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/barcode_counts/NPC_RNA_MPRAnalyze.txt", 
	header=T, 
	sep='\t', 
	stringsAsFactors=F, 
	quote="", 
	row.names=1
	)
#dcounts = as.matrix(read.csv("/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/barcode_counts/Hob_DNA_MPRAnalyze.txt", sep='\t', row.names = 1))
#rcounts = as.matrix(read.csv("/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/barcode_counts/Hob_RNA_MPRAnalyze.txt", sep='\t', row.names = 1))


# Coerce to matrices for later

dcounts = as.matrix(as.data.frame(dcounts_table))
rcounts = as.matrix(as.data.frame(rcounts_table))


# Estimate and extract depth factors on entire dataset

obj <- MpraObject(dnaCounts = dcounts, rnaCounts = rcounts, colAnnot = annot)
obj <- estimateDepthFactors(obj, lib.factor = "replicate")

dnadepth <- obj@dnaDepth
rnadepth <- obj@rnaDepth


# Split original tables, keep colnames and rownames

dcounts_table_split <- split(dcounts_table, rep(1:nchunks, chunksize))
rcounts_table_split <- split(rcounts_table, rep(1:nchunks, chunksize))

data_item <- list()
for (i in 1:nchunks) {
	df1 <- dcounts_table_split[i]
	df2 <- rcounts_table_split[i]
	tmp <- c(df1, df2)
	data_item[[i]] <- tmp
}


# Define a function to analyze each chunk of data

scalable_function <- function(stuff, a, ddepth, rdepth)  {
d <- as.matrix(stuff[[1]])
r <- as.matrix(stuff[[2]])

obj2 <- MpraObject(dnaCounts = d, rnaCounts = r, colAnnot = a)
obj2 <- setDepthFactors(obj2, dnaDepth = ddepth, rnaDepth = rdepth)
obj2 <- analyzeComparative(
	obj2,
	dnaDesign = ~ replicate + barcode_allele,
	rnaDesign = ~ replicate + allele,
	correctControls = FALSE,
	reducedDesign = ~ replicate
	)

res <- testLrt(obj2)

return(res)
}


# superApply to submit many Sherlock jobs in parallel

sapOut <- superApply(
		data_item, 
		scalable_function, 
		a=annot, 
		ddepth=dnadepth, 
		rdepth=rnadepth, 
		tasks=nchunks, 
		extraBashLines="module load R", 
		partition="hns", 
		time="1-23:59:00", 
		mem="20GB"
		)


# Merge results of sapOut into one df

final_df <- do.call(rbind, sapOut)


# Fix FDRs to correct values
pvals <- final_df[,'pval']
FDR <- p.adjust(pvals, method='BH')
final_df$fdr <- FDR


# Write results to out

write.table(
			final_df, 
            file="/scratch/users/cweiss19/MPRA/DATA/second_RNA_DNA_sequencing/barcode_counts/NPC_results_full_parallel_test.txt", 
            sep="\t", 
            quote = F, 
            na = "", 
            row.names = T, 
            col.names = T
            )

