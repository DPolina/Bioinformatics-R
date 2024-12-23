rm(list=ls())

# deletion penalty
del.penalty <- -2

# nucleotide substitution matrix (diagonal is a reword, all the others cells are penalties)
S <- matrix( -1,  nrow=4, ncol=4 )
diag(S) <- 1
rownames(S) <- c("A", "T", "G", "C")
colnames(S) <- c("A", "T", "G", "C")
print(S)

# input data to align, where dna1 will be in the column header, dna2 - in the raw header.  
dna1 <- unlist(strsplit("AAACGTAA", "")) #Creating the character vector for the first sequence.
dna2 <- unlist(strsplit("ACGTA", "")) #Creating the character vector for the second sequence.

# initialise F matrix
F.row <- length(dna2)+1
F.col <- length(dna1)+1
F <- matrix( nrow=F.row, ncol=F.col )

############### Algorithm starts  ###############
#### PART 1: fill in the first deletion penalty row/column
for (i in 1:F.row) {
  ### YOUR CODE GOES HERE ####
}

for (j in 1:F.col) {
  ### YOUR CODE GOES HERE ####
}


#### PART 2: fill in the matrix F
for (i in 2:F.row){
  for (j in 2:F.col)
  {
    dna1.letter <- dna1[j-1]
    dna2.letter <- dna2[i-1]
    
    ### YOUR CODE GOES HERE ####
    
  }
}
print(F)

#### PART 3: trace back
AlignmentA <- ""
AlignmentB <- ""
i <- F.row
j <- F.col

# reconstruct the back way and construct the actual alignments
while (i > 1 && j > 1)
{
  # retrieve from F the current score (i,j), back-diagonal score as well as Up and Left scores
  Score <- F[i,j]
  ScoreDiag <- F[i-1, j-1]
  ScoreUp <- F[i-1, j]
  ScoreLeft <- F[i, j-1]
  
  dna1.letter <- dna1[j-1]
  dna2.letter <- dna2[i-1]
  
  # go back from down-right corner to up left and retrieve the route we went on a previous step
  if (Score == ScoreDiag + S[dna1.letter, dna2.letter])
  {
    AlignmentA <- paste(dna1.letter, AlignmentA, sep = "")
    AlignmentB <- paste(dna2.letter, AlignmentB, sep = "")
    ### YOUR CODE GOES HERE ####
  }
  else if (Score == ScoreLeft + del.penalty)
  {
    ### YOUR CODE GOES HERE ####
  }
  else if (Score == ScoreUp + del.penalty)
  {
    ### YOUR CODE GOES HERE ####
  }
}

# add head / tail to both sequences
while (i > 1)
{
  ### YOUR CODE GOES HERE ####
}

while (j > 1)
{
  ### YOUR CODE GOES HERE ####
}

# show the result
AlignmentA
AlignmentB
