fst_dxy <- function(catalog_file,structure_file, pop_map_file) {

# Error/help message if any of the inputs are missing  
if(missing(catalog_file) | missing(structure_file) | missing(pop_map_file)) {
  print(paste("This script needs you to define the location and name of your catalog.tags.tsv file, your structure file (which it expects to be tab-delimited and to have 0 == missing data), and the pop_map file which defines which populations are being used in the comparisons to generate FST and Dxy. Sample names as they appear in your structure file should be in the left column of the pop_map file, the population they map to should be in the right column. Use 'exclude' for any samples you do not want to include, and 'outgroup' for your outgroup samples used to rescale dxy. Example of pop_map file:"))
  print(paste("Ceyx.erit.KU.23190","1"))
  print(paste("Ceyx.erit.KU.23309","1"))
  print(paste("Ceyx.erit.KU.12622","2"))
  print(paste("Ceyx.erit.KU.12774","2"))
  print(paste("Ceyx.erit.UW-73854","exclude"))
  print(paste("Ceyx.lepi.14384","outgroup"))
  print(paste("Ceyx.lepi.19259","outgroup"))
  print(paste("Ceyx.mela.19297","exclude"))
  print(paste(""))
  print(paste("Example of calling fst_dxy:"))
  cat('fst_dxy("/home/a499a400/ceyx_stacks_output/ceyx.catalog.tags.tsv","/home/a499a400/ceyx_stacks_output/ceyx_60_m5.structure.tsv","pop_map")\n\n')
  stop("Please fill in the missing info in your function call")
}

# Reading in the catalog
catalog <- as.matrix(read.table(catalog_file))

# Reading in the structure file and removing the header if there is one
structure <- readLines(structure_file)
headerline <- grep("# Stacks",structure,fixed=TRUE)
if(length(headerline)==0) {
  headerline <- 0
}

# After ditching any header lines, splitting the file out into columns based on the tab-separator 
struct <- structure[(headerline+1):(length(structure))]
nrow_struct <- length(struct)
struct <- unlist(strsplit(struct,"\t",fixed=TRUE))
struct <- matrix(struct,nrow=nrow_struct,byrow=TRUE)

# If there is an existing populations column in the file, getting rid of that  
if(struct[1,2]=="") {
  struct <- struct[,-2]
}

# Getting the unique loci from the structure locus name row  
unique_loci <- unique(as.numeric(struct[1,2:(dim(struct)[2])]))

# Paring the catalog down to only include theseloci  
catalog <- catalog[(which(as.numeric(catalog[,3]) %in% as.numeric(unique_loci))),]

# Getting an array of the population names which are not "outgroup" or "exclude" for our comparisons  
pop_map <- as.matrix(read.table(pop_map_file)) 
no_pops <- unique(pop_map[which(!(pop_map[,2] %in% c("outgroup","exclude"))),2])

# Figuring out the dimensions of our output file by the number of between population comparisons we are carrying out
# Making the headers for our output file  
comparisons <- combn(no_pops,2)
outputname <- NULL
for (i in 1:(dim(comparisons)[2])) {
  outputname <- c(outputname,(paste(comparisons[1,i],"&",comparisons[2,i],"fst",sep="")))
  outputname <- c(outputname,(paste(comparisons[1,i],"&",comparisons[2,i],"dxy",sep="")))
  outputname <- c(outputname,(paste(comparisons[1,i],"&",comparisons[2,i],"dstd",sep="")))
}

outputname <- c("locus_name","locus_length",outputname)

# Stepping through each of the loci remaining in the catalog  
for (i in 1:(dim(catalog)[1])) {
  # Getting locus name from catalog
  locus_name <- as.numeric(catalog[i,3])
  # Extracting all SNPs from the structure file corresponding to that locus name
  base_changes <- cbind(struct[,1],struct[,(which(as.numeric(struct[1,])==locus_name))])
  # Extracting locus length from the catalog
  locus_length <- nchar(catalog[i,9])
  # Pulling out the outgroup haplotypes across these SNPs
  outgroup_haplotypes <- base_changes[(which(base_changes[,1] %in% pop_map[(which(pop_map[,2]=="outgroup")),1])),2:(dim(base_changes)[2])]
  # Massaging the haplotypes into a matrix format if there is only one SNP and/or missing data
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
  # Creating the temporary output array for this locus
  locustemp <- matrix(c(locus_name,locus_length),nrow=1)
  # Now stepping through the populations
  for (j in 1:(length(no_pops)-1)) {
    # These variables, by default, are NA
    fst <- NA
    dxy <- NA
    dstd <- NA
    # Pulling out the haplotypes for population j across locus i's SNPs
    j_haplotypes <- base_changes[(which(base_changes[,1] %in% pop_map[(which(pop_map[,2]==no_pops[j])),1])),2:(dim(base_changes)[2])]
    # Massaging the haplotypes into a matrix format if there is only one SNP and/or missing data
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
    # Counting the number of haplotypes present
    j_nos <- dim(j_haplotypes)[1]
    # Only perform the following calculations if we actually have data for individuals from population j
    if(j_nos>0) {
      # Getting the number of unique haplotypes in population
      j_uniques <- unique(j_haplotypes)
      j_Hs <- 0
      # For each unique haplotype, summing its frequency^2
      for (k in 1:(dim(j_uniques)[1])) {
        j_Hs <- j_Hs + (sum(paste(j_uniques[k,], collapse="")==apply(j_haplotypes,1,paste,collapse=""))/j_nos)^2
      }
      # Heterozygosity equals 1 - sum of (haplotype frequencies ^2)
      j_Hs <- 1-j_Hs
      # Comparing population j to population m
      for (m in (j+1):(length(no_pops))) {
        # Pulling out the haplotypes for population m across locus i's SNPs
        m_haplotypes <- base_changes[(which(base_changes[,1] %in% pop_map[(which(pop_map[,2]==no_pops[m])),1])),2:(dim(base_changes)[2])]
        # Massaging the haplotypes into a matrix format if there is only one SNP and/or missing dat
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
        # Counting the number of haplotypes present
        m_nos <- dim(m_haplotypes)[1]
        # Only perform the following calculations if we actually have data for individuals from population m
        if(m_nos>0) {
          # Getting the number of unique haplotypes in population
          m_uniques <- unique(m_haplotypes)
          m_Hs <- 0
          # For each unique haplotype, summing its frequency^2
          for (k in 1:(dim(m_uniques)[1])) {
            m_Hs <- m_Hs + (sum(paste(m_uniques[k,], collapse="")==apply(m_haplotypes,1,paste,collapse=""))/m_nos)^2
          }
          # Heterozygosity equals 1 - sum of (haplotype frequencies ^2)
          m_Hs <- 1-m_Hs
          # Now calculating total heterozygosity across populations j and m combined
          all_nos <- j_nos + m_nos
          all_haplotypes <- rbind(j_haplotypes,m_haplotypes)
          all_uniques <- unique(all_haplotypes)
          all_Ht <- 0
          for (k in 1:(dim(all_uniques)[1])) {
            all_Ht <- all_Ht + (sum(paste(all_uniques[k,], collapse="")==apply(all_haplotypes,1,paste,collapse=""))/all_nos)^2
          }
          all_Ht <- 1-all_Ht
          # Using subpopulation heterozygosity and total heterozygosity to calculate fst
          fst <- (all_Ht-((j_Hs+m_Hs)/2))/all_Ht
          # Calculating dxy
          dxy <- 0
          for (n in 1:(dim(j_haplotypes)[1])) {
            for (o in 1:(dim(m_haplotypes)[1])) {
              dxy <- dxy + (dim(j_haplotypes)[2])-sum(j_haplotypes[n,]==m_haplotypes[o,])
            }
          }
          # Averaging the count of nucleotide differences over the number of comparisons made
          dxy <- dxy/(j_nos*m_nos)
          # If outgroup data is available calculating dxy between outgroups and populations j+m combined
          if(length(outgroup_haplotypes)>0){
            dxy_outgroup <- 0
            for (n in 1:(dim(all_haplotypes)[1])) {
              for (o in 1:(dim(outgroup_haplotypes)[1])) {
                dxy_outgroup <- dxy_outgroup + (dim(all_haplotypes)[2])-sum(all_haplotypes[n,]==outgroup_haplotypes[o,])
              }
            }
          # Using outgroup data to "standardize" ingroup dxy to control for variation in substitution rate etc.  
          dxy_outgroup <- dxy_outgroup/(all_nos*(dim(outgroup_haplotypes)[1]))
          dstd <- dxy/dxy_outgroup  
          } else {
          dstd <- NA
          }
        }
      }
     } 
     # Binding these output variables to the locus name and locus length
     temp <- matrix(c(fst,dxy,dstd),nrow=1)
     locustemp <- cbind(locustemp,temp)
     }
   # Binding the output for variables for this locus across all population comparisons to those of the previous loci
   outputname <- rbind(outputname,locustemp)
}

write.table(outputname,"dxy_fst.table",quote=FALSE,row.names=FALSE,col.names=FALSE)

}
