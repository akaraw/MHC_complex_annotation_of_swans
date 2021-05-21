#!/usr/bin/env Rscript
#usage: Rscript --vanilla consensus.R aln.fa consensus.fa
args = commandArgs(trailingOnly=TRUE)
aln = args[1]
outtxt = args[2]

library(DECIPHER) 
aln0 <- readDNAStringSet(aln, format='fasta')
(cons0 <- ConsensusSequence(aln0, threshold=0.05, ambiguity=T, ignoreNonBases=T, includeTerminalGaps=F))
writeXStringSet(cons0, filepath=outtxt)
