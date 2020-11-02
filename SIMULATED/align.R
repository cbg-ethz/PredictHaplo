

library("Biostrings")

library("seqinr")


HIV <-  readDNAStringSet("HIV1.fas", format="fasta")
 

Xset.haplo <- readDNAStringSet("clones.fas", format="fasta")

mat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
Align <- pairwiseAlignment(Xset.haplo,  toString(HIV[[1]]),type = "local", substitutionMatrix = mat, gapOpening = -3, gapExtension = -1)


A <- list()
for(i in 1:10){
  A[[i]] <- toString(aligned(Align)[[i]])
}
write.fasta(A,names = as.character(1:10),file = "haplosAligned.fasta")



