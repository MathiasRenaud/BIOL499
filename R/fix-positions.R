setwd('~/Documents/499data/mar16')
UTR_input <- read.delim('fix_pos-input.txt', 
                   header = FALSE, 
                   stringsAsFactors = FALSE,
                   col.names = c('name', 'chr', 'strand', 'start', 'end', 'seq'))

UTRs <- UTR_input[c(1,31,41,304),] #for testing

#------------------------------------------------------------------------------

## METHOD 1: using rbind() to append outputs
## fix.positions(x) takes a dataframe of UTR data (processed biomart download)
##    and splits spliced UTRs into separate entries for each exon.
## The output can then be used in the downstream analyses for this project.
## x: dataframe with 6 columns (name, chr, strand, stard, end, seq)

fix.positions <- function(x, verbose=TRUE){
  
  # Store number of rows in input
  nrows=nrow(x)
  
  # Initalize output dataframe
  out.df <- data.frame(name=character(0),
                       chr=integer(0),
                       strand=integer(0),
                       start=integer(0),
                       end=integer(0),
                       sq=character(0),
                       stringsAsFactors = FALSE)
  
  for (i in 1:nrows) {
    
    if (verbose==TRUE){print(paste(i, "of", nrows))}
    
    name <- as.character(x[i,1])
    chr <- as.integer(x[i,2])
    strand <- as.character(x[i,3])
    start <- as.character(x[i,4])
    end <- as.character(x[i,5])
    sq <- as.character(x[i,6])
    
    #if (strand=='1'){
    #  strd <- '+'
    #} else if (strand=='-1'){
    #  strd <- '-'
    #}
    
    # Check if Start and End Fields contain semicolon
    # If 'start' and 'end' contain at least 1 semicolon
    if (grepl(";", start) && grepl(";", end)){
      # vector of start indices
      starts <- strsplit(as.character(x[i,4]), ";")
      # vector of end indices
      ends <- strsplit(as.character(x[i,5]), ";")
      # number of subsequences
      numSeqs <- length(starts[[1]])
      # initialize vector for lengths of subsequnces
      seqLengths <- c()
      exon = 1
      
      # Populate vector of subsequence lengths
      for (n in seq(1, numSeqs)){
        len <- (as.integer(ends[[1]][n])-as.integer(starts[[1]][n])) + 1
        seqLengths <- c(seqLengths, len)
      }

      #
      for (l in seq(1, length(seqLengths))){
        
        if (l==1){
          out.df <- rbind(out.df, 
                          setNames(as.list(c(name=as.character(paste0(name, ".", exon)),
                                             chr=as.integer(chr),
                                             strand=as.character(strand),
                                             start=as.integer(starts[[1]][l]),
                                             end=as.integer(ends[[1]][l]),
                                             sq=as.character(substr(sq, 1, seqLengths[l])))),
                                   names(out.df)),
                          stringsAsFactors = FALSE)
          index = seqLengths[l] + 1
          exon = exon + 1
        } else {
          out.df <- rbind(out.df, 
                          setNames(as.list(c(name=as.character(paste0(name, ".", exon)),
                                             chr=as.integer(chr),
                                             strand=as.character(strand),
                                             start=as.integer(starts[[1]][l]),
                                             end=as.integer(ends[[1]][l]),
                                             sq=as.character(substr(sq, index, (index + seqLengths[l] - 1))))),
                                   names(out.df)),
                          stringsAsFactors = FALSE)
          index = index + seqLengths[l]
          exon = exon + 1
        }
      }
    
    } else if (!(grepl(";", start) && grepl(";", end))) {
      out.df <- rbind(out.df, 
                      setNames(as.list(c(name=as.character(name),
                                         chr=as.integer(chr),
                                         strand=as.character(strand),
                                         start=as.integer(start),
                                         end=as.integer(end),
                                         sq=as.character(sq))), 
                               names(out.df)),
                      stringsAsFactors = FALSE)
    }
  }
  return(out.df)
}

new_UTRs <- fix.positions(UTR_input)
write.table(new_UTRs, file='fix_pos-output.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

#------------------------------------------------------------------------------

## METHOD 2: using vectors to store outputs
fix.positions2 <- function(input){
  # Initialize number of rows
  nrows=nrow(UTRs)

  A <- c() #names
  B <- c() #chrs
  C <- c() #strands
  D <- c() #starts
  E <- c() #ends
  G <- c() #sqs

  for (i in 1:nrows) {
    # Check if Start and End Fields contain semicolon
    name <- as.character(UTRs[i,1])
    chr <- as.integer(UTRs[i,2])
    strand <- as.integer(UTRs[i,3])
    start <- as.character(UTRs[i,4])
    end <- as.character(UTRs[i,5])
    sq <- as.character(UTRs[i,6])
    
    # If 'start' and 'end' contain at least 1 semicolon
    if (grepl(";", start) && grepl(";", end)){
      # vector of start indices
      starts <- strsplit(as.character(UTRs[i,4]), ";")
      # vector of end indices
      ends <- strsplit(as.character(UTRs[i,5]), ";")
      # number of subsequences
      numSeqs <- length(starts[[1]])
      # initialize vector for lengths of subsequnces
      seqLengths <- c()
      
      # Populate vector of subsequence lengths
      for (n in seq(1, numSeqs)){
        len <- (as.integer(ends[[1]][n])-as.integer(starts[[1]][n])) + 1
        seqLengths <- c(seqLengths, len)
      }
      
      #
      for (l in seq(1, length(seqLengths))){
        
        if (l==1){
          A <- c(A, name) #names
          B <- c(B, chr) #chrs
          C <- c(C, strand) #strands
          D <- c(D, starts[[1]][l]) #starts
          E <- c(E, ends[[1]][l]) #ends
          G <- c(G, substr(sq, 1, seqLengths[l])) #seqs
          
          index = seqLengths[l] + 1
          
        } else {
          A <- c(A, name) #names
          B <- c(B, chr) #chrs
          C <- c(C, strand) #strands
          D <- c(D, starts[[1]][l]) #starts
          E <- c(E, ends[[1]][l]) #ends
          G <- c(G, substr(sq, index, (index + seqLengths[l] - 1))) #seqs
          
          index = index + seqLengths[l]
        }
      }
      
    } else if (!(grepl(";", start) && grepl(";", end))) {
      A <- c(A, name) #names
      B <- c(B, chr) #chrs
      C <- c(C, strand) #strands
      D <- c(D, start) #starts
      E <- c(E, end) #ends
      G <- c(G, sq) #seqs
    }
  }
  df <- data.frame(names=A,
                   chr=B,
                   strand=C,
                   start=D,
                   end=E,
                   seq=G,
                   stringsAsFactors = FALSE)
  return(df)
}
