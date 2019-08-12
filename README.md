# fst_dxy v0.0
I inherited a STACKs dataset consising of a structure (containing all SNPs for a given locus, not one SNP/locus) and catalog.tags.tsv files. I used the R-code in this repository to generate dxy and fst across the loci for a given comparison between populations defined in a separate pop_map file, where sample name as it appears in the structure file is in the left-hand column, and the population each sample is assigned to is in the right column. Samples to not be included in the comparisons should be labelled with "exclude". Samples that are outgroups (for standardizing dxy to the dxy calculated between the ingroup comparison populations and the outgroup) should be labelled with "outgroup".

```
Ceyx.erit.KU.23190	1
Ceyx.erit.KU.23309	1
Ceyx.erit.KU.12622	2
Ceyx.erit.KU.12774	2
Ceyx.erit.UW-73854	exclude
Ceyx.lepi.14384	outgroup
Ceyx.lepi.19259	outgroup
Ceyx.mela.19297	exclude
```

Example of invoking R-script:
```
fst_dxy("/home/a499a400/ceyx_stacks_output/ceyx.catalog.tags.tsv","/home/a499a400/ceyx_stacks_output/ceyx_60_m5.structure.tsv","pop_map")
```

Output:
A tab-delimited file with locus name in the first column, locus length in the next, and then repeating rows of Fst, Dxy, Dstd (Dxy scaled to outgroup) for each population comparison carried out. Note - Dxy/Dstd are not divided by locus length in this output, so you'll need to do this manually if this is what you are after.

# Suggested citation
This code was first published in TBD Ceyx.  

If you could cite the pub, and the progam as below, that would be lovely:  
Alexander, A. 2017. fst_dxy v0.0 Available from https://github.com/laninsky/fst_dxy

This pipeline also wouldn't be possible without:  
R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/

# Version history
v0.0: ready to rock 'n' roll
