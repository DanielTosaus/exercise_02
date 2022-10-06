library("seqinr")

seqs <- read.fasta("fishes.fna.gz")

library("Biostrings")

seqs2 <- readDNAStringSet("fishes.fna.gz")

# We find reverse complement of MID and search for elements between MID and its reverse complement

#grep("ACGAGTGCGT.+ACGCACTCGT", seqs2, perl=TRUE)

grep(paste("ACGAGTGCGT",".+","ACGCACTCGT", sep = ""), seqs2, perl=TRUE)
