rm(list=ls())

# Deletion penalty
del.penalty <- -1

# Nucleotide substitution matrix (diagonal is a reward, all the others cells are penalties)
S <- matrix( -1, nrow=4, ncol=4 )
diag(S) <- 1
rownames(S) <- c("A", "T", "G", "C")
colnames(S) <- c("A", "T", "G", "C")
print(S)

# Input data to align, where dna1 will be in the column header, dna2 - in the row header.  
# Creating the character vector for the first sequence.
dna1 <- unlist(strsplit("ACTGGAGAA", ""))
# Creating the character vector for the second sequence.
dna2 <- unlist(strsplit("CATGCAGTCGCCACTTTGAAAACTAGTCGATGC", ""))

# Alignment should look as follows:
# ------------AC-TGGAGAA------------
# CATGCAGTCGCCACTTTGA-AAACTAGTCGATGC

# Initialise F matrix
F.row <- length(dna2)+1
F.col <- length(dna1)+1
F <- matrix( nrow=F.row, ncol=F.col )

############### Algorithm starts  ###############
#### PART 1: fill in the first row and column
for (i in 1:F.row) {s
  ### YOUR CODE GOES HERE ####
}

for (j in 1:F.col) {
  ### YOUR CODE GOES HERE ####
}

# Initialize variables for finding the cell with the highest score
max.F <- 0
max.i <- 0
max.j <- 0

#### PART 2: fill in the matrix F
for (i in 2:F.row){
  for (j in 2:F.col)
  {
    dna1.letter <- dna1[j-1]
    dna2.letter <- dna2[i-1]
    
    ### YOUR CODE GOES HERE ####
    
  }
}
rownames(F) <- c('', dna2)
colnames(F) <- c('', dna1)
print(F)

#### PART 3: trace back
AlignmentA <- ""
AlignmentB <- ""

i <- F.row
j <- max.j

# Initialize variables for calculating identity and gaps
ident <- 0
gaps <- 0

# Add tail to both sequences
#  (see trailing sequence of minuses in "Alignment should look as follows")
while (i > max.i) {
  ### YOUR CODE GOES HERE ###
}

# Trace back from the cell with the highest score
# Reconstruct the back way and construct the actual alignments
while (F[i, j] != 0)
{
  # Retrieve from F the current score (i,j), back-diagonal score as well as Up and Left scores
  Score <- F[i, j]
  ScoreDiag <- F[i-1, j-1]
  ScoreUp <- F[i-1, j]
  ScoreLeft <- F[i, j-1]
  
  dna1.letter <- dna1[j-1]
  dna2.letter <- dna2[i-1]
  
  # Go back from down-right corner to up left and retrieve the route we went on a previous step
  
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

# Alignment spans can be calculated now
query.start <- ### YOUR CODE GOES HERE ###
query.end <- ### YOUR CODE GOES HERE ###
reference.start <- ### YOUR CODE GOES HERE ###
reference.end <- ### YOUR CODE GOES HERE ###

# Add head to both sequences
#  (see leading sequence of minuses in "Alignment should look as follows")
while (i > 1) {
  ### YOUR CODE GOES HERE ###
}

### Calculate alignment statistics ###

# Calculate alignment length
align.len <- ### YOUR CODE GOES HERE ###
# Length of query sequence
query.len <- ### YOUR CODE GOES HERE ###
# Calculate coverage
coverage <- ### YOUR CODE GOES HERE ###
# Calculate percent of identity
p.ident <- ### YOUR CODE GOES HERE ###
# Calculate percent of gaps
p.gaps <- ### YOUR CODE GOES HERE ###

  
# Print alignment
cat(sprintf("  %s\n  %s\n\n",
            AlignmentA, AlignmentB))
