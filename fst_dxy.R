#fasta_from_tags_tsv <- function(catalog_file,sumstats_file, pop_map_file) {

catalog_file <- "/home/a499a400/ceyx_stacks_output/ceyx.catalog.tags.tsv"
structure_file <- "/home/a499a400/ceyx_stacks_output/ceyx_60_m5.structure.tsv"
pop_map_file <- "pop_map"

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

pop_map <- as.matrix(read.table(pop_map_file))
no_pops <- unique(pop_map[which(!(pop_map[,2] %in% c("outgroup","exclude"))),2])

comparisons <- combn(no_pops,2)
outputname <- NULL
for (i in 1:(dim(comparisons)[2])) {
  outputname <- c(outputname,(paste(comparisons[1,i],"&",comparisons[2,i],"fst",sep="")))
  outputname <- c(outputname,(paste(comparisons[1,i],"&",comparisons[2,i],"dxy",sep="")))
  outputname <- c(outputname,(paste(comparisons[1,i],"&",comparisons[2,i],"dstd",sep="")))
}

outputname <- c("locus_name","locus_length",outputname)

for (i in 1:(dim(catalog)[1])) {
  locus_name <- as.numeric(catalog[i,3])
  base_changes <- cbind(struct[,1],struct[,(which(as.numeric(struct[1,])==locus_name))])
  locus_length <- nchar(catalog[i,9])
  outgroup_haplotypes <- base_changes[(which(base_changes[,1] %in% pop_map[(which(pop_map[,2]=="outgroup")),1])),2:(dim(base_changes)[2])]
  if (dim(base_changes)[2]==2) {
    outgroup_haplotypes <- matrix(outgroup_haplotypes,ncol=1)
   }
   remove <- (unique(which(outgroup_haplotypes=="0",arr.ind=TRUE)[,1]))
   if(length(remove)>0) {
     outgroup_haplotypes <- outgroup_haplotypes[-remove,]
   }
   if (dim(base_changes)[2]==2) {
     outgroup_haplotypes <- matrix(outgroup_haplotypes,ncol=1)
   }
  locustemp <- matrix(c(locus_name,locus_length),nrow=1)
  for (j in 1:(length(no_pops)-1)) {
    fst <- NA
    dxy <- NA
    dstd <- NA
    j_haplotypes <- base_changes[(which(base_changes[,1] %in% pop_map[(which(pop_map[,2]==no_pops[j])),1])),2:(dim(base_changes)[2])]
    if (dim(base_changes)[2]==2) {
      j_haplotypes <- matrix(j_haplotypes,ncol=1)
    }
    remove <- (unique(which(j_haplotypes=="0",arr.ind=TRUE)[,1]))
    if(length(remove)>0) {
      j_haplotypes <- j_haplotypes[-remove,]
    }
    if (dim(base_changes)[2]==2) {
      j_haplotypes <- matrix(j_haplotypes,ncol=1)
    }
    j_nos <- dim(j_haplotypes)[1]
    if(j_nos>0) {
      j_uniques <- unique(j_haplotypes)
      j_Hs <- 0
      for (k in 1:(dim(j_uniques)[1])) {
        j_Hs <- j_Hs + (sum(paste(j_uniques[k,], collapse="")==apply(j_haplotypes,1,paste,collapse=""))/j_nos)^2
      }
      j_Hs <- 1-j_Hs
      for (m in (j+1):(length(no_pops))) {
        m_haplotypes <- base_changes[(which(base_changes[,1] %in% pop_map[(which(pop_map[,2]==no_pops[m])),1])),2:(dim(base_changes)[2])]
        if (dim(base_changes)[2]==2) {
          m_haplotypes <- matrix(m_haplotypes,ncol=1)
        }
        remove <- (unique(which(m_haplotypes=="0",arr.ind=TRUE)[,1]))
        if(length(remove)>0) {
          m_haplotypes <- m_haplotypes[-remove,]
        }
        if (dim(base_changes)[2]==2) {
          m_haplotypes <- matrix(m_haplotypes,ncol=1)
        }
        m_nos <- dim(m_haplotypes)[1]
        if(m_nos>0) {
          m_uniques <- unique(m_haplotypes)
          m_Hs <- 0
          for (k in 1:(dim(m_uniques)[1])) {
            m_Hs <- m_Hs + (sum(paste(m_uniques[k,], collapse="")==apply(m_haplotypes,1,paste,collapse=""))/m_nos)^2
          }
          m_Hs <- 1-m_Hs
          all_nos <- j_nos + m_nos
          all_haplotypes <- rbind(j_haplotypes,m_haplotypes)
          all_uniques <- unique(all_haplotypes)
          all_Ht <- 0
          for (k in 1:(dim(all_uniques)[1])) {
            all_Ht <- all_Ht + (sum(paste(all_uniques[k,], collapse="")==apply(all_haplotypes,1,paste,collapse=""))/all_nos)^2
          }
          all_Ht <- 1-all_Ht
          fst <- (all_Ht-((j_Hs+m_Hs)/2))/all_Ht
          dxy <- 0
          for (n in 1:(dim(j_haplotypes)[1])) {
            for (o in 1:(dim(m_haplotypes)[1])) {
              dxy <- dxy + (dim(j_haplotypes)[2])-sum(j_haplotypes[n,]==m_haplotypes[o,])
            }
          }
          dxy <- dxy/(j_nos*m_nos)
          if(length(outgroup_haplotypes)>0){
            dxy_outgroup <- 0
            for (n in 1:(dim(all_haplotypes)[1])) {
              for (o in 1:(dim(outgroup_haplotypes)[1])) {
                dxy_outgroup <- dxy_outgroup + (dim(all_haplotypes)[2])-sum(all_haplotypes[n,]==outgroup_haplotypes[o,])
              }
            }
          dxy_outgroup <- dxy_outgroup/(all_nos*(dim(outgroup_haplotypes)[1]))
          dstd <- dxy/dxy_outgroup  
          } else {
          dstd <- NA
          }
        }
      }
     } 
     temp <- matrix(c(fst,dxy,dstd),nrow=1)
     locustemp <- cbind(locustemp,temp)
     }
   outputname <- rbind(outputname,locustemp)
}

write.table(outputname,"dxy_fst.table",quote=FALSE,row.names=FALSE,col.names=FALSE)
