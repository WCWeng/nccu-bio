# hw2_105753039

######################################
# initial
######################################
# read parameters
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("USAGE: Rscript hw2_105753039.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call. = FALSE)
}

# parse parameters
i<-1
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<-args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-args[i+1]
    i<-i+1
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

print("PARAMETERS")
print(paste("input file         :", i_f))
print(paste("output file        :", o_f))
print(paste("score file         :", s_f))
print(paste("gap open penalty   :", g_o))
print(paste("gap extend penalty :", g_e))

######################################
# main
######################################

#pseudocode:

#initialize
pam250 <- read.table("pam250.txt")
pam250 <- as.matrix(pam250)

#library("seqinr")
#library("Biostrings")
s1 <- read first seq. from test.fasta
s2 <- read sceond seq. from test.fasta

#initialize score_matrix
score_matrix <- matrix(row = s1+1, cul = s2+1)
set values of first row = 0 -1 -2 -3 ... -length(s1)
set values of first culumn = 0 -1 -2 -3 ... -length(s2)

tag <- 0  #i=0 : no gap yet
while(s1!=NULL || s2!=NULL){
  if(tag == 0){
    #gap_open penalty
    gap_penalty <- -10
    tag <- 1
  }
  #gap_extend penalty
  gap_penalty <- -2
  for(int i = 1; i <= length(s1); i++)
    for(int j = 1; j <= length(s2); j++)
      score_matrix[i,j] <- MAX(score_matrix[i, j-1] + gap_penalty,
                               score_matrix[i-1, j-1] + pam250[i,j],
                               score_matrix[i-1, j] + gap_penalty)
}
global <- trace-back score_matrix along with hightest score
write.fasta(global, "result.fasta")