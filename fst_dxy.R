#fasta_from_tags_tsv <- function(catalog_file,sumstats_file) {

catalog_file <- "/home/a499a400/ceyx_stacks_output/ceyx.catalog.tags.tsv"
structure_file <- "/home/a499a400/ceyx_stacks_output/ceyx_60_m5.structure.tsv" 

catalog <- as.matrix(read.table(catalog_file))
structure <- readLines(structure_file)
headerline <- grep("# Stacks",structure,fixed=TRUE)
if(length(headerline)==0) {
  headerline <- 0
}

struct <- structure[(headerline+1):(length(structure))]
nrow_struct <- length(struct)
struct <- unlist(strsplit(struct,"\t",fixed=TRUE))
struct <- matrix(struct,nrow=nrow_struct,byrow=TRUE)

if(struct[1,2]=="") {
  struct <- struct[,-2]
}

unique_loci <- unique(as.numeric(struct[1,2:(dim(struct)[2])]))
catalog <- catalog[(which(as.numeric(catalog[,3]) %in% as.numeric(unique_loci))),]

for (i in 1:(dim(catalog)[2])) {
  locus_name <- as.numeric(catalog[i,3])
  base_changes <- cbind(struct[,1],struct[,(which(as.numeric(struct[1,])==locus_name))])
  locus_length <- nchar(catalog[i,9])
