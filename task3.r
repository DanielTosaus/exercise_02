library("Biostrings")
library("readxl")
fasta_f <- "fishes.fna.gz"
MIDs_f <- "fishes_MIDs.xls"

fmids_list <- c("ACGAGTGCGT", "ACGCTCGACA","AGACGCACTC")
rmids_list <- c("ACGAGTGCGT", "ACGCTCGACA","AGACGCACTC")
rcf <- c("ACGCACTCGT", "TGTCGAGCGT","GAGTGCGTCT")
rcr <- c("ACGCACTCGT", "TGTCGAGCGT","GAGTGCGTCT")

sample_names <- c("Paracheirodon_axelrodi", "Poecilia_sphenops", "Hypostomus_plecostomus")

demultiplex <- function(fasta_f, fmids_list, rmids_list, sample_names){
  seqs2 <- readDNAStringSet(fasta_f)
  #mids <- read_excel(MIDs_f)
  
  # We find reverse complement of MID and search for elements between MID and its reverse complement
  
  
  
  stats <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(stats) <- c('name', 'n_seq')
  
  for(i in 1: 3){
    print(i)
    print(fmids_list[i])
    print(rcf[i])
    res <- grep(paste(fmids_list[i],".+",rcf[i], sep = ""), seqs2, perl=TRUE)
    write.fasta(sequences = trimLRPatterns(Lpattern = fmids_list[i], Rpattern = rcf[i], seqs2[res]), names = res, nbchar = 80, file.out = paste(sample_names[i],".gz", sep =""))
    data <- data.frame(name=sample_names[i], n_seq = toString(length(res)))
    stats <- rbind(stats, data)
  }
  
  
  write.table(stats, file = "report.txt")
  return(stats)
}
try <- demultiplex(fasta_f, fmids_list, rmids_list, sample_names)

seqs2 <- readDNAStringSet(fasta_f)
#mids <- read_excel(MIDs_f)

# We find reverse complement of MID and search for elements between MID and its reverse complement

# grep("ACGAGTGCGT.+ACGCACTCGT", seqs2, perl=TRUE)
  
stats <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(stats) <- c('name', 'n_seq')

for(i in 1: 3){
  print(i)
  print(fmids_list[i])
  print(rcf[i])
  res <- grep(paste(fmids_list[i],".+",rcf[i], sep = ""), seqs2, perl=TRUE)
  write.fasta(sequences = trimLRPatterns(Lpattern = fmids_list[i], Rpattern = rcf[i], seqs2[res]), names = res, nbchar = 80, file.out = paste(sample_names[i],".gz", sep =""))
  data <- data.frame(name=sample_names[i], n_seq = toString(length(res)))
  stats <- rbind(stats, data)
}


write.table(stats, file = "report.txt")

