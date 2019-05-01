
##-----------------------------------------------
## To do the analysis with GREAT:
##-----------------------------------------------

setwd('~/Documents/499data/test/')
library(GenomicRanges)
library(rGREAT)
library(rtracklayer)
library(ratelimitr)

##-----------------------------------------------
## Accessory Functions
##-----------------------------------------------

## Make a BED file of the genomic positions of the target sites
make.bed <- function(targets){
  # initialize BED format dataframe
  bedFrame <- data.frame(chr=character(),
                         start=numeric(),
                         end=numeric(),
                         name=character(),
                         score=numeric(),
                         strand=integer(),
                         stringsAsFactors = FALSE)
  # populate dataframe
  for (i in seq(1, nrow(targets))){
    # Chromosome
    bedFrame[i,1] = paste0('chr', targets[i,5])
    # Start
    bedFrame[i,2] = as.numeric(targets[i,7]) + as.numeric(targets[i,2])
    # End
    bedFrame[i,3] = as.numeric(targets[i,7]) + as.numeric(targets[i,3])
    # Name
    bedFrame[i,4] = paste0(targets[i,1], "_", targets[i,4])
    # Score
    bedFrame[i,5] = 0
    # Strand
    if (targets[i,6] == '1'){
      strand = '+'
    } else if (targets[i,6] == '-1'){
      strand = '-'
    } else {
      stop("error: strand information missing")
    }
    bedFrame[i,6] = strand
  }
  #return(bedFrame)
  # Convert BED format data frame to GRanges object
  #bed <- makeGRangesFromDataFrame(bedFrame, seqnames.field = 'chr', end.field = 'end')
  return(bedFrame)
}

## Limit GREAT requests to every 5 minutes
great_lim <- limit_rate(submitGreatJob,rate(n=1, period=300))


##-----------------------------------------------
## Main Funtion
##-----------------------------------------------

## Import and process 3'UTR background
bg.in <- read.delim("bg.bed6", header=FALSE)
for (i in 1:nrow(bg.in)){
  # add chr prefix
  bg.in[i,1] <- paste0("chr", bg.in[i,1])
  # change strand format to +/-
  if (bg.in[i,6]==1 | bg.in[i,6]=='1'){
    bg.in[i,6] <- '+'
  } else if (bg.in[i,6]==-1 | bg.in[i,6]=='-1'){
    bg.in[i,6] <- '-'
  }
  #print(paste(i, 'of', nrow(bg.in)))
}
names(bg.in) <- c('chr', 'start', 'end', 'name', 'score', 'strand')

# for miRNA in LOM, create a new object with all lines containing miRNA
allmiRNA <- read.delim("miR_input2.txt", header=FALSE)
allTargets <- read.delim("enrichment_input2.txt", header=FALSE)


## TargetScan to GREAT
ts.to.great <- function(LOM, Targets){
  # start timer for runtime updates
  ptm <- proc.time()
  
  chain <- import.chain('~/Documents/499data/R/danRer11ToDanRer7.over.chain')
  
  #setwd("~/Documents/499data/test/output/")
  setwd("~/Documents/499data/test/output-all/")
  
  # For each miRNA, perform functional enrichment
  for (miRNA in 1:nrow(LOM)){
    print(paste0("Processing: ", LOM[miRNA,1], " (", miRNA, " of ", nrow(LOM), ") [elapsed: ", round(proc.time()[3] - ptm[3]), "s]"))
    
    # initialize variables
    miR_name <- LOM[miRNA,1]
    tmp <- data.frame(name=character(0),
                      miR=character(0),
                      UTRstart=integer(0),
                      UTRend=integer(0),
                      type=character(0),
                      chr=integer(0),
                      strand=integer(0),
                      start=integer(0),
                      end=integer(0),
                      stringsAsFactors = FALSE)
    
    ## x. Create miRNA Target BED files
    #print(paste0("Compiling Targets [elapsed: ", round(proc.time()[3] - ptm[3]), "s]"))
    
    for (line in 1:nrow(Targets)){
      if (grepl(gsub("^",miR_name,"$"), Targets[line,2], perl=TRUE)){
        tmp <- rbind(tmp, Targets[line,c(1,3:9)], stringsAsFactors = FALSE)
      }
    }
    names(tmp) <- c('name', 'UTRstart', 'UTRend', 'type', 'chr', 'strand', 'start', 'end')
    #print(head(tmp))
    # Make target BED file (GRanges Object)
    targetBed <- make.bed(tmp)
    targetGR <- makeGRangesFromDataFrame(targetBed, seqnames.field = 'chr', end.field = 'end')
    
    # Make background BED file (GRanges Object)
    bgBed <- rbind(bg.in, targetBed)
    bgGR <- makeGRangesFromDataFrame(bgBed, seqnames.field = 'chr', end.field = 'end')
    
    # UCSC LiftOver from danRer11 to danRer7 for GREAT
    # Foreground
    convertedGRList <- liftOver(targetGR, chain)
    convertedGR <- unlist(convertedGRList)
    # Background
    convertedBGList <- liftOver(bgGR, chain)
    convertedBG <- unlist(convertedBGList)
    
    ## x. Submit Files to GREAT
    print(paste0("Submitting to GREAT [elapsed: ", round(proc.time()[3] - ptm[3]), "s]"))
    GREAT <- great_lim(convertedGR,
                       bg=convertedBG,
                       species='danRer7',
                       rule='oneClosest')
    annotations <- getEnrichmentTables(GREAT, category=c('GO', 
                                                         #'Phenotype Data',
                                                         'Pathway Data', 
                                                         'Gene Expression', 
                                                         'Gene Families'))
    ## x. Write everything!
    # create directory for files
    #print(paste0("Writing results [elapsed: ", round(proc.time()[3] - ptm[3]), "s]"))
    dir.create(paste0(miR_name))
    # write targets
    write.table(tmp, 
                file = paste0(miR_name, "/", miR_name, "-Targets.txt"), 
                sep='\t', quote=FALSE, row.names=FALSE, col.names = FALSE)
    
    ## All Terms
    for (category in 1:length(annotations)){
      annotationTable <- annotations[[category]]
      categoryName <- gsub(" ", "", names(annotations)[category])
      write.table(annotationTable,
                  file = paste0(miR_name, "/", miR_name, "-",  categoryName, ".txt"),
                  sep='\t', quote=FALSE, row.names=FALSE)
    }
    
    ## Only Significant Terms (raw p-val < 0.05)
    # for (category in 1:length(annotations)){
    #   annotationTable <- annotations[[category]]
    #   categoryName <- gsub(" ", "", names(annotations)[category])
    #   output <- data.frame(ID=character(),
    #                        name=character(),
    #                        Hyper_Total_Regions=numeric(),
    #                        Hyper_Expected=numeric(),
    #                        Hyper_Foreground_Region_Hits=numeric(),
    #                        Hyper_Fold_Enrichment=numeric(),
    #                        Hyper_Region_Set_Coverage=numeric(),
    #                        Hyper_Term_Region_Coverage=numeric(),
    #                        Hyper_Foreground_Gene_Hits=numeric(),
    #                        Hyper_Background_Gene_Hits=numeric(),
    #                        Total_Genes_Annotated=numeric(),
    #                        Hyper_Raw_PValue=numeric(),
    #                        Hyper_Adjp_BH=numeric())
    #   for (term in 1:nrow(annotationTable))
    #     if (annotationTable[term,12] < 0.05){
    #       output <- rbind(output, annotationTable[term,])
    #     }
    #   write.table(output,
    #               file = paste0(miR_name, "/", miR_name, "-", categoryName, ".txt"),
    #               sep='\t', quote=FALSE, row.names=FALSE)
    # }
  }
  
}

ts.to.great(LOM=allmiRNA, Targets=allTargets)

#testing
#targetTest <- head(allTargets)
#miRNATest <- as.data.frame(allmiRNA[c(1,50,149,200,215),1])
#ts.to.great(LOM=miRNATest, Targets=allTargets)

## Single microRNA 
##-----------------------------------------------
## FOREGROUND
#miR29a_targets <- read.delim("biomart_data/miR-29a_targets-fixed.txt", header=FALSE)
miR29a_targets <- read.delim("~/Documents/499data/test/miR-29a_targets.txt", header=FALSE)
miR29a_bed <- make.bed(miR29a_targets)
miR29a_GRange <- makeGRangesFromDataFrame(miR29a_bed, seqnames.field = 'chr', end.field = 'end')

## BACKGROUND
## Process 3'UTR background
bg.in <- read.delim("bg.bed6", header=FALSE)
for (i in 1:nrow(bg.in)){
  # add chr prefix
  bg.in[i,1] <- paste0("chr", bg.in[i,1])
  # change strand format to +/-
  if (bg.in[i,6]==1 | bg.in[i,6]=='1'){
    bg.in[i,6] <- '+'
  } else if (bg.in[i,6]==-1 | bg.in[i,6]=='-1'){
    bg.in[i,6] <- '-'
  }
  #print(paste(i, 'of', nrow(bg.in)))
}
names(bg.in) <- c('chr', 'start', 'end', 'name', 'score', 'strand')
bg2 <- rbind(bg.in, miR29a_bed)
bgBed <- makeGRangesFromDataFrame(bg2, seqnames.field = 'chr', end.field = 'end')

## LIFTOVER
# UCSC LiftOver from danRer11 to danRer7 for GREAT
chain <- import.chain('../R/danRer11ToDanRer7.over.chain')
# Foreground
convertedGRangeList <- liftOver(miR29a_GRange, chain)
convertedGrange <- unlist(convertedGRangeList)
# Background
converted.bgList <- liftOver(bgBed, chain)
converted.bg <- unlist(converted.bgList)

## GREAT
# Retreive enrichment lists from GREAT
GREATanno.nobg <- submitGreatJob(convertedGrange, 
                            #bg=converted.bg, 
                            species='danRer7', 
                            rule='oneClosest')
GREATanno <- submitGreatJob(convertedGrange, 
                            bg=converted.bg, 
                            species='danRer7', 
                            rule='oneClosest')
#ontologies <- as.vector(availableOntologies(GREATanno))
categories <- as.vector(availableCategories(GREATanno))
miR29a_anno <- getEnrichmentTables(GREATanno, category=c('GO', 
                                                         #'Phenotype Data',
                                                         'Pathway Data', 
                                                         'Gene Expression', 
                                                         'Gene Families'))
miR29a_anno.nobg <- getEnrichmentTables(GREATanno.nobg, category=c('GO', 
                                                         #'Phenotype Data',
                                                         'Pathway Data', 
                                                         'Gene Expression', 
                                                         'Gene Families'))

# ## OUTPUT
# # GO
# miR29a_anno[[1]][1,c(2,12,13)] # GO Molecular Function
# miR29a_anno[[2]][1,c(2,12,13)] # GO Biological Process
# miR29a_anno[[3]][1,c(2,12,13)] # GO Cellular Component
# # Pathway Data
# miR29a_anno[[4]][1,c(2,12,13)] # Wiki Pathways
# # Gene Expression
# miR29a_anno[[5]][1,c(2,12,13)] # Zebrafish WT Expression
# # Gene Families
# miR29a_anno[[6]][1,c(2,12,13)] # InterPro (protein families)
# miR29a_anno[[7]][1,c(2,12,13)] # TreeFam (gene families)

# write.table(miR29a_anno[[1]][,c(2,12,13)], file="miR-29a_MolecularFunction-2.txt", sep='\t', quote=FALSE, row.names=FALSE)
# write.table(miR29a_anno[[2]][,c(2,12,13)], file="miR-29a_BiologicalProcess-2.txt", sep='\t', quote=FALSE, row.names=FALSE)
# write.table(miR29a_anno[[3]][,c(2,12,13)], file="miR-29a_CellularComponent-2.txt", sep='\t', quote=FALSE, row.names=FALSE)
# write.table(miR29a_anno[[4]][,c(2,12,13)], file="miR-29a_Pathway-2.txt", sep='\t', quote=FALSE, row.names=FALSE)
# write.table(miR29a_anno[[5]][,c(2,12,13)], file="miR-29a_Expression-2.txt", sep='\t', quote=FALSE, row.names=FALSE)
# write.table(miR29a_anno[[6]][,c(2,12,13)], file="miR-29a_ProteinFam-2.txt", sep='\t', quote=FALSE, row.names=FALSE)
# write.table(miR29a_anno[[7]][,c(2,12,13)], file="miR-29a_GeneFam-2.txt", sep='\t', quote=FALSE, row.names=FALSE)


##-----------------------------------------------
## To do the analysis with Enrichr (old)
#library(gProfileR)
#library(enrichR)

## Test gconvert with miR-29a target gene list
#ensembl.29a <- scan("miR-29a.txt", character(), quote="") # 2251 Entries
#gconvert.29a <- gconvert(ensembl.29a, organism="drerio", target="ENTREZGENE")
#entrez.29a <- as.vector(gconvert.29a$target) # 1816 Entries ...435 gene names are N/A
#len <- length(entrez.29a)
#cat("Identified", len, "genes")

## Test enrichR with miR-29a target gene list
#fish.dbs <- c("GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018")
#out.29a <- enrichr(entrez.29a, fish.dbs)
#printEnrich(out.29a, "enrichR-29a.txt", columns = c(1, 2, 4, 9))
