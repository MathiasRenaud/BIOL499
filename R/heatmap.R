library(pheatmap)
library(RColorBrewer)
library(gplots)

setwd("~/Documents/499data/test/output-all/")
dirs<- list.dirs()
dirs <- dirs[2:length(dirs)]

## Todo: Remove entirely blue rows (ontology terms with no significant association)
# for each row, check if all rows == 1, if not then output the row - collect in new matrix.

## Binary heatmaps (less than 0.05 == red, greater == blue)
#------------------------------------------------------------------------------
# GO Cellular Component (files[2]) - binary, only significant rows!
#------------------------------------------------------------------------------
example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[2]))
IDs <- example[order(example$ID),][,1]
miRs <- gsub("./", "", dirs)
matr <- matrix(nrow=length(IDs), ncol=length(miRs))
colnames(matr) <- miRs
row.names(matr) <- IDs

d=0
for (d in 1:length(dirs)){
  files <- list.files(dirs[d])
  MolFun <- read.delim(paste0(dirs[d], "/", files[2]))[c(1,13)]
  # order by term ID
  orderedMolFun <- MolFun[order(MolFun$ID),]
  #print(head(orderedMolFun))
  matr[,d] <- as.matrix(orderedMolFun[2])
}

# Remove terms with no significant assocaitions.
mat2 <- matrix(ncol=length(miRs))
r=0
ontTerms <- c()
for (r in 1:nrow(matr)){
 if (!(all(matr[r,]>0.05))){
   mat2 <- rbind(mat2, matr[r,])
   ontTerms <- c(ontTerms, names(matr[,1][r]))
 }
}
mat2 <- mat2[2:nrow(mat2),]
row.names(mat2) <- ontTerms

pheatmap(mat2, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 10,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "~/Desktop/test.pdf",
         #filename = "CellularComponent-heat-bin.pdf",
         border_color = NA
)

#------------------------------------------------------------------------------
# GO Biological Process (files[1]) - binary
#------------------------------------------------------------------------------
example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[1]))
IDs <- example[order(example$ID),][,1]
miRs <- gsub("./", "", dirs)
matr <- matrix(nrow=length(IDs), ncol=length(miRs))
colnames(matr) <- miRs
row.names(matr) <- IDs

d=0
for (d in 1:length(dirs)){
  files <- list.files(dirs[d])
  MolFun <- read.delim(paste0(dirs[d], "/", files[1]))[c(1,13)]
  # order by term ID
  orderedMolFun <- MolFun[order(MolFun$ID),]
  #print(head(orderedMolFun))
  matr[,d] <- as.matrix(orderedMolFun[2])
}

# Remove terms with no significant assocaitions.
mat2 <- matrix(ncol=length(miRs))
r=0
ontTerms <- c()
for (r in 1:nrow(matr)){
  if (!(all(matr[r,]>0.05))){
    mat2 <- rbind(mat2, matr[r,])
    ontTerms <- c(ontTerms, names(matr[,1][r]))
  }
}
mat2 <- mat2[2:nrow(mat2),]
row.names(mat2) <- ontTerms

pheatmap(mat2, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 90,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "BiologicalProcess-heat-bin.pdf",
         border_color = NA
)

# Remove miRNAs with no significant assocaitions.
mat3 <- matrix(nrow=nrow(mat2))
c=0
miRNAs <- c()
for (c in 1:ncol(mat2)){
  if (!(all(mat2[,c]>0.05))){
    mat3 <- cbind(mat3, mat2[,c])
    miRNAs <- c(miRNAs, names(mat2[1,][c]))
  }
}
mat3 <- mat3[,2:ncol(mat3)]
colnames(mat3) <- miRNAs

pheatmap(mat3, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 100,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "BiologicalProcess-heat-bin-2.pdf",
         border_color = NA
)

#------------------------------------------------------------------------------
# GO Cellular Component (files[2]) - binary
#------------------------------------------------------------------------------
example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[2]))
IDs <- example[order(example$ID),][,1]
miRs <- gsub("./", "", dirs)
matr <- matrix(nrow=length(IDs), ncol=length(miRs))
colnames(matr) <- miRs
row.names(matr) <- IDs

d=0
for (d in 1:length(dirs)){
  files <- list.files(dirs[d])
  MolFun <- read.delim(paste0(dirs[d], "/", files[2]))[c(1,13)]
  # order by term ID
  orderedMolFun <- MolFun[order(MolFun$ID),]
  #print(head(orderedMolFun))
  matr[,d] <- as.matrix(orderedMolFun[2])
}

# Remove terms with no significant assocaitions.
mat2 <- matrix(ncol=length(miRs))
r=0
ontTerms <- c()
for (r in 1:nrow(matr)){
  if (!(all(matr[r,]>0.05))){
    mat2 <- rbind(mat2, matr[r,])
    ontTerms <- c(ontTerms, names(matr[,1][r]))
  }
}
mat2 <- mat2[2:nrow(mat2),]
row.names(mat2) <- ontTerms

pheatmap(mat2,  
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 10,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "CellularComponent-heat-bin.pdf",
         border_color = NA
)

# Remove miRNAs with no significant assocaitions.
mat3 <- matrix(nrow=nrow(mat2))
c=0
miRNAs <- c()
for (c in 1:ncol(mat2)){
  if (!(all(mat2[,c]>0.05))){
    mat3 <- cbind(mat3, mat2[,c])
    miRNAs <- c(miRNAs, names(mat2[1,][c]))
  }
}
mat3 <- mat3[,2:ncol(mat3)]
colnames(mat3) <- miRNAs

pheatmap(mat3, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 20,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "CellularComponent-heat-bin-2.pdf",
         border_color = NA
)

#------------------------------------------------------------------------------
# GO Molecular Function (files[3]) - binary
#------------------------------------------------------------------------------
example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[3]))
IDs <- example[order(example$ID),][,1]
miRs <- gsub("./", "", dirs)
matr <- matrix(nrow=length(IDs), ncol=length(miRs))
colnames(matr) <- miRs
row.names(matr) <- IDs

d=0
for (d in 1:length(dirs)){
  files <- list.files(dirs[d])
  MolFun <- read.delim(paste0(dirs[d], "/", files[3]))[c(1,13)]
  # order by term ID
  orderedMolFun <- MolFun[order(MolFun$ID),]
  #print(head(orderedMolFun))
  matr[,d] <- as.matrix(orderedMolFun[2])
}

# Remove terms with no significant assocaitions.
mat2 <- matrix(ncol=length(miRs))
r=0
ontTerms <- c()
for (r in 1:nrow(matr)){
  if (!(all(matr[r,]>0.05))){
    mat2 <- rbind(mat2, matr[r,])
    ontTerms <- c(ontTerms, names(matr[,1][r]))
  }
}
mat2 <- mat2[2:nrow(mat2),]
row.names(mat2) <- ontTerms

pheatmap(mat2,  
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 10,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "MolecularFunction-heat-bin.pdf",
         border_color = NA
)

# Remove miRNAs with no significant assocaitions.
mat3 <- matrix(nrow=nrow(mat2))
c=0
miRNAs <- c()
for (c in 1:ncol(mat2)){
  if (!(all(mat2[,c]>0.05))){
    mat3 <- cbind(mat3, mat2[,c])
    miRNAs <- c(miRNAs, names(mat2[1,][c]))
  }
}
mat3 <- mat3[,2:ncol(mat3)]
colnames(mat3) <- miRNAs

pheatmap(mat3, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 40,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "MolecularFunction-heat-bin-2.pdf",
         border_color = NA
)

#------------------------------------------------------------------------------
# InterPro (files[4]) - binary
#------------------------------------------------------------------------------
example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[4]))
IDs <- example[order(example$ID),][,1]
miRs <- gsub("./", "", dirs)
matr <- matrix(nrow=length(IDs), ncol=length(miRs))
colnames(matr) <- miRs
row.names(matr) <- IDs

d=0
for (d in 1:length(dirs)){
  files <- list.files(dirs[d])
  MolFun <- read.delim(paste0(dirs[d], "/", files[4]))[c(1,13)]
  # order by term ID
  orderedMolFun <- MolFun[order(MolFun$ID),]
  #print(head(orderedMolFun))
  matr[,d] <- as.matrix(orderedMolFun[2])
}

# Remove terms with no significant assocaitions.
mat2 <- matrix(ncol=length(miRs))
r=0
ontTerms <- c()
for (r in 1:nrow(matr)){
  if (!(all(matr[r,]>0.05))){
    mat2 <- rbind(mat2, matr[r,])
    ontTerms <- c(ontTerms, names(matr[,1][r]))
  }
}
mat2 <- mat2[2:nrow(mat2),]
row.names(mat2) <- ontTerms

pheatmap(mat2,
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 40,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "InterPro-heat-bin.pdf",
         border_color = NA
)

# Remove miRNAs with no significant assocaitions.
mat3 <- matrix(nrow=nrow(mat2))
c=0
miRNAs <- c()
for (c in 1:ncol(mat2)){
  if (!(all(mat2[,c]>0.05))){
    mat3 <- cbind(mat3, mat2[,c])
    miRNAs <- c(miRNAs, names(mat2[1,][c]))
  }
}
mat3 <- mat3[,2:ncol(mat3)]
colnames(mat3) <- miRNAs

pheatmap(mat3, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 50,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "InterPro-heat-bin-2.pdf",
         border_color = NA
)

#------------------------------------------------------------------------------
# TreeFam (files[6]) - binary
#------------------------------------------------------------------------------
example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[6]))
IDs <- example[order(example$ID),][,1]
miRs <- gsub("./", "", dirs)
matr <- matrix(nrow=length(IDs), ncol=length(miRs))
colnames(matr) <- miRs
row.names(matr) <- IDs

d=0
for (d in 1:length(dirs)){
  files <- list.files(dirs[d])
  MolFun <- read.delim(paste0(dirs[d], "/", files[6]))[c(1,13)]
  # order by term ID
  orderedMolFun <- MolFun[order(MolFun$ID),]
  #print(head(orderedMolFun))
  matr[,d] <- as.matrix(orderedMolFun[2])
}

# Remove terms with no significant assocaitions.
mat2 <- matrix(ncol=length(miRs))
r=0
ontTerms <- c()
for (r in 1:nrow(matr)){
  if (!(all(matr[r,]>0.05))){
    mat2 <- rbind(mat2, matr[r,])
    ontTerms <- c(ontTerms, names(matr[,1][r]))
  }
}
mat2 <- mat2[2:nrow(mat2),]
row.names(mat2) <- ontTerms

pheatmap(mat2, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 20,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "TreeFam-heat-bin.pdf",
         border_color = NA
)

# Remove miRNAs with no significant assocaitions.
mat3 <- matrix(nrow=nrow(mat2))
c=0
miRNAs <- c()
for (c in 1:ncol(mat2)){
  if (!(all(mat2[,c]>0.05))){
    mat3 <- cbind(mat3, mat2[,c])
    miRNAs <- c(miRNAs, names(mat2[1,][c]))
  }
}
mat3 <- mat3[,2:ncol(mat3)]
colnames(mat3) <- miRNAs

pheatmap(mat3, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 20,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "TreeFam-heat-bin-2.pdf",
         border_color = NA
)

#------------------------------------------------------------------------------
# WikiPathways (files[7]) - binary
#------------------------------------------------------------------------------
example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[7]))
IDs <- example[order(example$ID),][,1]
miRs <- gsub("./", "", dirs)
matr <- matrix(nrow=length(IDs), ncol=length(miRs))
colnames(matr) <- miRs
row.names(matr) <- IDs

d=0
for (d in 1:length(dirs)){
  files <- list.files(dirs[d])
  MolFun <- read.delim(paste0(dirs[d], "/", files[7]))[c(1,13)]
  # order by term ID
  orderedMolFun <- MolFun[order(MolFun$ID),]
  #print(head(orderedMolFun))
  matr[,d] <- as.matrix(orderedMolFun[2])
}


# Remove terms with no significant assocaitions.
mat2 <- matrix(ncol=length(miRs))
r=0
ontTerms <- c()
for (r in 1:nrow(matr)){
  if (!(all(matr[r,]>0.05))){
    mat2 <- rbind(mat2, matr[r,])
    ontTerms <- c(ontTerms, names(matr[,1][r]))
  }
}
mat2 <- mat2[2:nrow(mat2),]
row.names(mat2) <- ontTerms

pheatmap(mat2,  
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 10,
         cellheight = 20,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "WikiPathways-heat-bin.pdf",
         border_color = NA
)
# Remove miRNAs with no significant assocaitions.
mat3 <- matrix(nrow=nrow(mat2))
c=0
miRNAs <- c()
for (c in 1:ncol(mat2)){
  if (!(all(mat2[,c]>0.05))){
    mat3 <- cbind(mat3, mat2[,c])
    miRNAs <- c(miRNAs, names(mat2[1,][c]))
  }
}
mat3 <- mat3[,2:ncol(mat3)]
colnames(mat3) <- miRNAs

pheatmap(mat3, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 10,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "WikiPathways-heat-bin-2.pdf",
         border_color = NA
)

#------------------------------------------------------------------------------
# Expression (files[8]) - binary
#------------------------------------------------------------------------------
example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[1]))
IDs <- example[order(example$ID),][,1]
miRs <- gsub("./", "", dirs)
matr <- matrix(nrow=length(IDs), ncol=length(miRs))
colnames(matr) <- miRs
row.names(matr) <- IDs

d=0
for (d in 1:length(dirs)){
  files <- list.files(dirs[d])
  MolFun <- read.delim(paste0(dirs[d], "/", files[1]))[c(1,13)]
  # order by term ID
  orderedMolFun <- MolFun[order(MolFun$ID),]
  #print(head(orderedMolFun))
  matr[,d] <- as.matrix(orderedMolFun[2])
}

# Remove terms with no significant assocaitions.
mat2 <- matrix(ncol=length(miRs))
r=0
ontTerms <- c()
for (r in 1:nrow(matr)){
  if (!(all(matr[r,]>0.05))){
    mat2 <- rbind(mat2, matr[r,])
    ontTerms <- c(ontTerms, names(matr[,1][r]))
  }
}
mat2 <- mat2[2:nrow(mat2),]
row.names(mat2) <- ontTerms

pheatmap(mat2,  
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 80,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "Expression-heat-bin.pdf",
         border_color = NA
)
# Remove miRNAs with no significant assocaitions.
mat3 <- matrix(nrow=nrow(mat2))
c=0
miRNAs <- c()
for (c in 1:ncol(mat2)){
  if (!(all(mat2[,c]>0.05))){
    mat3 <- cbind(mat3, mat2[,c])
    miRNAs <- c(miRNAs, names(mat2[1,][c]))
  }
}
mat3 <- mat3[,2:ncol(mat3)]
colnames(mat3) <- miRNAs

pheatmap(mat3, 
         #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
         color = colorRampPalette(c('red3', 'steelblue'))(2),
         breaks = c(0, 0.05, 1),
         legend_breaks = c(0.05, 1),
         annotation_names_col = FALSE, 
         show_colnames = T,
         show_rownames = T,
         cellwidth = 100,
         cellheight = 10,
         cluster_cols = T,
         cluster_rows = T,
         width = NA,
         height = NA,
         filename = "Expression-heat-bin-2.pdf",
         border_color = NA
)

#------------------------------------------------------------------------------

# ## 3-colour heatmaps:
# #------------------------------------------------------------------------------
# # GO Biological Process (files[1])
# #------------------------------------------------------------------------------
# example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[1]))
# IDs <- example[order(example$ID),][,1]
# miRs <- gsub("./", "", dirs)
# matr <- matrix(nrow=length(IDs), ncol=length(miRs))
# colnames(matr) <- miRs
# row.names(matr) <- IDs
# 
# d=0
# for (d in 1:length(dirs)){
#   files <- list.files(dirs[d])
#   MolFun <- read.delim(paste0(dirs[d], "/", files[1]))[c(1,13)]
#   # order by term ID
#   orderedMolFun <- MolFun[order(MolFun$ID),]
#   #print(head(orderedMolFun))
#   matr[,d] <- as.matrix(orderedMolFun[2])
# }
# 
# pheatmap(matr, 
#          #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
#          color = colorRampPalette(c('red3', 'orangered', 'steelblue'))(3),
#          breaks = c(0, 0.01, 0.05, 1),
#          legend_breaks = c(0.01, 0.05, 1),
#          annotation_names_col = FALSE, 
#          show_colnames = T,
#          show_rownames = T,
#          cellwidth = 200,
#          cellheight = 10,
#          cluster_cols = T,
#          cluster_rows = T,
#          width = NA,
#          height = NA,
#          filename = "BiologicalProcess-heat.pdf",
#          border_color = NA
# )
# 
# #------------------------------------------------------------------------------
# # GO Cellular Component (files[2])
# #------------------------------------------------------------------------------
# example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[2]))
# IDs <- example[order(example$ID),][,1]
# miRs <- gsub("./", "", dirs)
# matr <- matrix(nrow=length(IDs), ncol=length(miRs))
# colnames(matr) <- miRs
# row.names(matr) <- IDs
# 
# d=0
# for (d in 1:length(dirs)){
#   files <- list.files(dirs[d])
#   MolFun <- read.delim(paste0(dirs[d], "/", files[2]))[c(1,13)]
#   # order by term ID
#   orderedMolFun <- MolFun[order(MolFun$ID),]
#   #print(head(orderedMolFun))
#   matr[,d] <- as.matrix(orderedMolFun[2])
# }
# 
# pheatmap(matr, 
#          #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
#          color = colorRampPalette(c('red3', 'orangered', 'steelblue'))(3),
#          breaks = c(0, 0.01, 0.05, 1),
#          legend_breaks = c(0.01, 0.05, 1),
#          annotation_names_col = FALSE, 
#          show_colnames = T,
#          show_rownames = T,
#          cellwidth = 50,
#          cellheight = 10,
#          cluster_cols = T,
#          cluster_rows = T,
#          width = NA,
#          height = NA,
#          filename = "CellularComponent-heat.pdf",
#          border_color = NA
# )
# 
# #------------------------------------------------------------------------------
# # GO Molecular Function (files[3])
# #------------------------------------------------------------------------------
# example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[3]))
# IDs <- example[order(example$ID),][,1]
# miRs <- gsub("./", "", dirs)
# matr <- matrix(nrow=length(IDs), ncol=length(miRs))
# colnames(matr) <- miRs
# row.names(matr) <- IDs
# 
# d=0
# for (d in 1:length(dirs)){
#   files <- list.files(dirs[d])
#   MolFun <- read.delim(paste0(dirs[d], "/", files[3]))[c(1,13)]
#   # order by term ID
#   orderedMolFun <- MolFun[order(MolFun$ID),]
#   #print(head(orderedMolFun))
#   matr[,d] <- as.matrix(orderedMolFun[2])
# }
# 
# pheatmap(matr, 
#          #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
#          color = colorRampPalette(c('red3', 'orangered', 'steelblue'))(3),
#          breaks = c(0, 0.01, 0.05, 1),
#          legend_breaks = c(0.01, 0.05, 1),
#          annotation_names_col = FALSE, 
#          show_colnames = T,
#          show_rownames = T,
#          cellwidth = 100,
#          cellheight = 10,
#          cluster_cols = T,
#          cluster_rows = T,
#          width = NA,
#          height = NA,
#          filename = "MolecularFunction-heat.pdf",
#          border_color = NA
# )
# 
# #------------------------------------------------------------------------------
# # InterPro (files[4])
# #------------------------------------------------------------------------------
# example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[4]))
# IDs <- example[order(example$ID),][,1]
# miRs <- gsub("./", "", dirs)
# matr <- matrix(nrow=length(IDs), ncol=length(miRs))
# colnames(matr) <- miRs
# row.names(matr) <- IDs
# 
# d=0
# for (d in 1:length(dirs)){
#   files <- list.files(dirs[d])
#   MolFun <- read.delim(paste0(dirs[d], "/", files[4]))[c(1,13)]
#   # order by term ID
#   orderedMolFun <- MolFun[order(MolFun$ID),]
#   #print(head(orderedMolFun))
#   matr[,d] <- as.matrix(orderedMolFun[2])
# }
# 
# pheatmap(matr, 
#          #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
#          color = colorRampPalette(c('red3', 'orangered', 'steelblue'))(3),
#          breaks = c(0, 0.01, 0.05, 1),
#          legend_breaks = c(0.01, 0.05, 1),
#          annotation_names_col = FALSE, 
#          show_colnames = T,
#          show_rownames = T,
#          cellwidth = 800,
#          cellheight = 10,
#          cluster_cols = T,
#          cluster_rows = T,
#          width = NA,
#          height = NA,
#          filename = "InterPro-heat.pdf",
#          border_color = NA
# )
# 
# #------------------------------------------------------------------------------
# # TreeFam (files[6])
# #------------------------------------------------------------------------------
# example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[6]))
# IDs <- example[order(example$ID),][,1]
# miRs <- gsub("./", "", dirs)
# matr <- matrix(nrow=length(IDs), ncol=length(miRs))
# colnames(matr) <- miRs
# row.names(matr) <- IDs
# 
# d=0
# for (d in 1:length(dirs)){
#   files <- list.files(dirs[d])
#   MolFun <- read.delim(paste0(dirs[d], "/", files[6]))[c(1,13)]
#   # order by term ID
#   orderedMolFun <- MolFun[order(MolFun$ID),]
#   #print(head(orderedMolFun))
#   matr[,d] <- as.matrix(orderedMolFun[2])
# }
# 
# pheatmap(matr, 
#          #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
#          color = colorRampPalette(c('red3', 'orangered', 'steelblue'))(3),
#          breaks = c(0, 0.01, 0.05, 1),
#          legend_breaks = c(0.01, 0.05, 1),
#          annotation_names_col = FALSE, 
#          show_colnames = T,
#          show_rownames = T,
#          cellwidth = 600,
#          cellheight = 10,
#          cluster_cols = T,
#          cluster_rows = T,
#          width = NA,
#          height = NA,
#          filename = "TreeFam-heat.pdf",
#          border_color = NA
# )
# #------------------------------------------------------------------------------
# # WikiPathways (files[7])
# #------------------------------------------------------------------------------
# example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[7]))
# IDs <- example[order(example$ID),][,1]
# miRs <- gsub("./", "", dirs)
# matr <- matrix(nrow=length(IDs), ncol=length(miRs))
# colnames(matr) <- miRs
# row.names(matr) <- IDs
# 
# d=0
# for (d in 1:length(dirs)){
#   files <- list.files(dirs[d])
#   MolFun <- read.delim(paste0(dirs[d], "/", files[7]))[c(1,13)]
#   # order by term ID
#   orderedMolFun <- MolFun[order(MolFun$ID),]
#   #print(head(orderedMolFun))
#   matr[,d] <- as.matrix(orderedMolFun[2])
# }
# 
# 
# pheatmap(matr, 
#          #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
#          color = colorRampPalette(c('red3', 'orangered', 'steelblue'))(3),
#          breaks = c(0, 0.01, 0.05, 1),
#          legend_breaks = c(0.01, 0.05, 1),
#          annotation_names_col = FALSE, 
#          show_colnames = T,
#          show_rownames = T,
#          cellwidth = 10,
#          cellheight = 10,
#          cluster_cols = T,
#          cluster_rows = T,
#          width = NA,
#          height = NA,
#          filename = "WikiPathways-heat.pdf",
#          border_color = NA
# )
# #------------------------------------------------------------------------------
# # Expression (files[8])
# #------------------------------------------------------------------------------
# example <- read.delim(paste0(dirs[1], "/", list.files(dirs[1])[1]))
# IDs <- example[order(example$ID),][,1]
# miRs <- gsub("./", "", dirs)
# matr <- matrix(nrow=length(IDs), ncol=length(miRs))
# colnames(matr) <- miRs
# row.names(matr) <- IDs
# 
# d=0
# for (d in 1:length(dirs)){
#   files <- list.files(dirs[d])
#   MolFun <- read.delim(paste0(dirs[d], "/", files[1]))[c(1,13)]
#   # order by term ID
#   orderedMolFun <- MolFun[order(MolFun$ID),]
#   #print(head(orderedMolFun))
#   matr[,d] <- as.matrix(orderedMolFun[2])
# }
# 
# pheatmap(matr, 
#          #color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(3),
#          color = colorRampPalette(c('red3', 'orangered', 'steelblue'))(3),
#          breaks = c(0, 0.01, 0.05, 1),
#          legend_breaks = c(0.01, 0.05, 1),
#          annotation_names_col = FALSE, 
#          show_colnames = T,
#          show_rownames = T,
#          cellwidth = 400,
#          cellheight = 10,
#          cluster_cols = T,
#          cluster_rows = T,
#          width = NA,
#          height = NA,
#          filename = "Expression-heat.pdf",
#          border_color = NA
# )
# 
# 
# #------------------------------------------------------------------------------
# # pHeatMap Default parameters
# #------------------------------------------------------------------------------
# pheatmap(mat,
#          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
#          kmeans_k = NA, 
#          breaks = NA, 
#          border_color = "grey60",
#          cellwidth = NA, 
#          cellheight = NA, 
#          scale = "none", 
#          cluster_rows = TRUE,
#          cluster_cols = TRUE, 
#          clustering_distance_rows = "euclidean",
#          clustering_distance_cols = "euclidean", 
#          clustering_method = "complete",
#          clustering_callback = identity2, 
#          cutree_rows = NA, 
#          cutree_cols = NA,
#          treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows, 50, 0), 
#          treeheight_col = ifelse((class(cluster_cols) == "hclust") || cluster_cols, 50, 0), 
#          legend = TRUE, 
#          legend_breaks = NA,
#          legend_labels = NA, 
#          annotation_row = NA, 
#          annotation_col = NA,
#          annotation = NA, 
#          annotation_colors = NA, 
#          annotation_legend = TRUE,
#          annotation_names_row = TRUE, 
#          annotation_names_col = TRUE,
#          drop_levels = TRUE, 
#          show_rownames = T, 
#          show_colnames = T, 
#          main = NA,
#          fontsize = 10, 
#          fontsize_row = fontsize, 
#          fontsize_col = fontsize,
#          display_numbers = F, 
#          number_format = "%.2f", 
#          number_color = "grey30",
#          fontsize_number = 0.8 * fontsize, 
#          gaps_row = NULL, 
#          gaps_col = NULL,
#          labels_row = NULL, 
#          labels_col = NULL, 
#          filename = NA, 
#          width = NA,
#          height = NA, 
#          silent = FALSE, 
#          na_col = "#DDDDDD",
#          ...)
# 
