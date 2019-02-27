#max print has to be changed for this to work, since to file contains many rows. I just made it arbitrarily huge
options(max.print=999999999)
library(tidyverse)
library(plyr)

#reads data (list of all CHH positions) into a df
CHH <- read.table("CHH_pos.txt", header = FALSE) 
#extracts a vector from said df
CHHvec <- CHH[[1]]

#br is just a vector containing a sequence of numbers from 1 to the largest position rounded up spaced 100 apart. 
#its components will serve as breaks within the hist() function below
br <- seq(1,round_any(max(CHHvec), 100)+101,by=100)

#The hist() function finds the number of CHHs within each of the 100 bp tiles
freq <- hist(CHHvec, breaks = br, include.lowest = TRUE, plot = FALSE)
imp <- cbind(freq[[1]],freq[[2]])

mCHH <- read.table("mCHH_pos.txt", header = FALSE)
mCHHvec <- mCHH[[1]]
mfreq <- hist(mCHHvec, breaks = br, include.lowest = TRUE, plot = FALSE)
mimp <- cbind(mfreq[[1]],mfreq[[2]])

percentages <- mfreq[[2]]/freq[[2]]

fin <- cbind(freq[[1]],freq[[2]], mfreq[[2]], percentages)
fin <- as.data.frame(fin)

#if a percentage for the bin is over 25 and has >5 Cs, it will be 1, otherwise 0
fin$over25 <- ifelse(fin$percentages > 0.25, ifelse(fin[[2]] > 5, 1, 0), 0)

#convert NAs to 0s
fin[is.na(fin)] <- 0
colnames(fin) <- c("startpos", "Cs", "mCs", "percentages", "over25")

write_tsv(fin, "CHHislands_over25.txt")

perc_nozeros <- percentages[ percentages != 0 ]
hist(perc_nozeros)

#create simple df which has one column with CHH context cytosines and one column with methylated: true/false
matched <- as.data.frame(CHHvec)
matched$mCHH <- ifelse(matched$CHHvec %in% mCHHvec, 1, 0)

CHHislands <- fin[fin$over25 == 1,]

#Read in modified gff file to know gene positions. Gene positions given by bash script
genes <- read.table("genepositions.txt", header = FALSE)
colnames(genes) <- c("feature", "start", "stop")
genes$genenumber <- c(1:length(genes$start))

#Add 2 kb ranges before and after genes
genes$beforewindow <- genes$start - 2000
genes$afterwindows <- genes$stop + 2000

#make a function to determine whether each island is proximal to gene
getValue <- function(x, data) {
  tmp <- data %>%
    dplyr::filter(beforewindow <= x, x <= start) %>%
    filter(row_number() == 1)
  return(tmp$genenumber)
}

testydata <- c(40000, 44280, 48900)
test <- sapply(testydata, getValue, data=genes)

#beforegene <- sapply(CHHislands$startpos, getValue, data=genes)
CHHislands$beforegene <- sapply(CHHislands$startpos, getValue, data=genes) 
CHHislands$beforegene <- gsub("integer(0)", "0", CHHislands$beforegene)

#CHHislands1 <- grep("[[:digit:]]", CHHislands$beforegene)

geneswithislands <- unique(CHHislands$beforegene) %>% grep("1","2","3","4","5","6","7","8","9")


#create 1000 random permutations of methylated positions
permutations <- replicate(1000, sample(matched$mCHH))


#
# Code graveyard (nothing beyond here matters)

#https://stackoverflow.com/questions/24766104/checking-if-value-in-vector-is-in-range-of-values-in-different-length-vector#
#getValue <- function(x, data) {
#  tmp <- data %>%
#    filter(V2 <= x, x <= V3)
#  return(tmp$V4)
#}

#x <- c(107599440, 150769180, 155204690)
#sapply(x, getValue, data=df)

#for (i in genes$start){
#  for (j in CHHislands$startpos){
#    for (k in genes$beforewindow){
#      for (l in CHHislands$beforegene){
#        if (j < i & j > k){
#          l = 1
#        }
#      }
#    }
#  }
#}

#for (i in genes$start){
#  for (k in genes$beforewindow){
#    CHHislands$beforegene <- ifelse(CHHislands$startpos < i & CHHislands$startpos > k, 1, 0)
#  }
#}

#for (i in CHHislands$startpos){
#  if (i < genes$start & i > genes$beforewindow){
#    beforegene <- print(1)
#    else 
#      beforegene <- print(0)
#    }
#  }
#}

#for ( i in CHHislands$startpos) {
#  if (i < genes$start){
#    print(genes$start)
#  }
#}

#for (i in genes$start){
#  for (j in CHHislands$startpos){
#    if (i > j){
#      beforegene <- print(1)
#    }
#  }
#}

#CHHislands$beforegene <- ifelse(CHHislands$startpos > genes$beforewindow, 
  #                            ifelse(CHHislands$startpos < genes$start, 1, 0), 
   #                           0)
#CHHislands$beforegene <- ifelse(any(CHHislands$startpos > genes$beforewindow), 1, 0)
#ifelse()
#beforegene <- sapply(CHHislands$startpos, " 1")
#  for (i in CHHislands$startpos) {
#    if (i < genes$start)
#  }


#create 1000 random permutations of methylated positions
#permutations <- replicate(1000, sample(matched$mCHH))

#beforegene <- sapply(CHHislands$startpos, getValue, data=genes)


