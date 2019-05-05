#max print has to be changed for this to work, since to file contains many rows. I just made it arbitrarily huge
options(max.print=999999999)
library(plyr)
library(dplyr)
library(readr)
library(snow)
library(snowfall)
#datatable package in Aline's script
#Append = true


##FOR HPC
##module load R/3.4.1 <- readr included

#reads data (list of all CHH positions) into a df
CHH <- read.table("CHH_pos.txt", header = FALSE) 
#extracts a vector from said df
CHHvec <- CHH[[1]]

#br is just a vector containing a sequence of numbers from 1 to the largest position rounded up spaced 100 apart. 
#its components will serve as breaks within the hist() function below
br <- seq(1,round_any(max(CHHvec), 100)+101,by=100)

#The hist() function finds the number of CHHs within each of the 100 bp tiles
freq <- hist(CHHvec, breaks = br, include.lowest = TRUE, plot = FALSE)
imp <- cbind(head(freq[[1]], n = length(freq[[1]])-1),freq[[2]])

mCHH <- read.table("mCHH_pos.txt", header = FALSE)
mCHHvec <- mCHH[[1]]
mfreq <- hist(mCHHvec, breaks = br, include.lowest = TRUE, plot = FALSE)
mimp <- cbind(head(mfreq[[1]], n = length(mfreq[[1]])-1),mfreq[[2]])

#Cs <- read.table("C_pos.txt", header = FALSE)
#Cvec <- Cs[[1]]
#Cfreq <- hist(Cvec, breaks = br1, include.lowest = TRUE, plot = FALSE)
#Cimp <- cbind(head(Cfreq[[1]], n = length(Cfreq[[1]])-1),Cfreq[[2]])
#Cimp1 <- cbind(Cimp, imp[,2])
#percentages1 <- Cimp1[,3]/Cimp1[,2]

percentages <- mfreq[[2]]/freq[[2]]

fin <- cbind(head(freq[[1]], n = length(freq[[1]])-1),freq[[2]], mfreq[[2]], percentages)
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
colnames(genes) <- c("feature", "start", "stop", "ID")

#Need to go back in here and ad ID (maize chr1)
#colnames(genes) <- c("feature", "start", "stop")

genes$genenumber <- c(1:length(genes$start))

#Add 2 kb ranges before and after genes
genes$beforewindow <- genes$start - 2000
genes$afterwindows <- genes$stop + 2000

#make a function to determine whether each island is proximal to gene
getValue1 <- function(x, data) {
  tmp <- data %>%
    dplyr::filter(beforewindow <= x, x <= start) %>%
    filter(row_number() == 1)
  return(tmp$genenumber)
}

getValue2 <- function(x, data) {
  tmp <- data %>%
    dplyr::filter(afterwindows >= x, x >= stop) %>%
    filter(row_number() == 1)
  return(tmp$genenumber)
}

#testydata <- c(40000, 44280, 48900)
#test <- sapply(testydata, getValue, data=genes)

#beforegene <- sapply(CHHislands$startpos, getValue, data=genes)
CHHislands$beforegene <- sapply(CHHislands$startpos, getValue1, data=genes) 
CHHislands$aftergene <- sapply(CHHislands$startpos, getValue2, data=genes)

genes$beforewindow2 <- genes$start - 1000
genes$afterwindows2 <- genes$stop + 1000

getValue3 <- function(x, data) {
  tmp <- data %>%
    dplyr::filter(beforewindow2 <= x, x <= start) %>%
    filter(row_number() == 1)
  return(tmp$genenumber)
}

getValue4 <- function(x, data) {
  tmp <- data %>%
    dplyr::filter(afterwindows2 >= x, x >= stop) %>%
    filter(row_number() == 1)
  return(tmp$genenumber)
}

#CvecDF <- as.data.frame(Cvec)

#CvecDF$beforegene <- c(1:length(Cvec))
#CvecDF$beforegene <- sapply(CvecDF$Cvec, getValue3, data=genes)
#CvecDF$beforegene <- as.numeric(CvecDF$beforegene)

#testCvecDF <- head(CvecDF, n = 1000000)
#testCvecDF$beforegene <- sapply(testCvecDF$Cvec, getValue3, data=genes)
#testCvecDF$beforegene <- as.numeric(testCvecDF$beforegene)
#testCvecDF$aftergene <- sapply(testCvecDF$Cvec, getValue4, data=genes)
#testCvecDF$aftergene <- as.numeric(testCvecDF$beforegene)
#testCvecDF1 <- unlist(testCvecDF)

#write(testCvecDF1, file = "testCvecDF.txt")

#CvecDF$aftergene <- c(1:length(Cvec))
#CvecDF$aftergene <- sapply(CvecDF$Cvec, getValue4, data=genes)
#CvecDF$aftergene <- as.numeric(CvecDF$aftergene)

#CHHvecDF <- as.data.frame(CHHvec)

#testCHHvecDF <- filter(CHHvecDF, CHHvecDF$CHHvec<max(testCvecDF$Cvec))

#testCHHvecDF$beforegene <- sapply(testCHHvecDF$CHHvec, getValue3, data=genes)
#testCHHvecDF$beforegene <- as.numeric(testCHHvecDF$beforegene)
#testCHHvecDF$aftergene <- sapply(testCHHvecDF$CHHvec, getValue4, data=genes)
#testCHHvecDF$aftergene <- as.numeric(testCHHvecDF$beforegene)
#testCHHvecDF1 <- unlist(testCHHvecDF)

#write(testCHHvecDF1, file = "testCHHvecDF.txt")

#Cisbeforegene <- filter(testCvecDF, testCvecDF$beforegene >= 1)
#Cisaftergene <- filter(testCvecDF, testCvecDF$aftergene >= 1)
#CHHisbeforegene <- filter(testCHHvecDF, testCHHvecDF$beforegene >= 1)
#CHHisaftergene <- filter(testCHHvecDF, testCHHvecDF$aftergene >= 1)

#Cnotbefore <- filter(testCvecDF, is.na(testCvecDF$beforegene))
#Cnotafter <- filter(testCvecDF, is.na(testCvecDF$aftergene))
#CHHnotbefore <- filter(testCHHvecDF, is.na(testCHHvecDF$beforegene))
#CHHnotafter <- filter(testCHHvecDF, is.na(testCHHvecDF$aftergene))

#percentCHHbefore <- length(CHHisbeforegene$CHHvec)/length(Cisbeforegene$Cvec)
#percentCHHnotbefore <- length(CHHnotbefore$CHHvec)/length(Cnotbefore$Cvec)

#percentCHHafter <- length(CHHisaftergene$CHHvec)/length(Cisaftergene$Cvec)
#percentCHHnotafter <- length(CHHnotafter$CHHvec)/length(Cnotafter$Cvec)

# It's 53.49% (maize chr1)
#genomewidepercentCHH <- 231858368/433461618

#fin$beforegene <- sapply(fin$startpos, getValue1, data=genes)
#fin$aftergene <- sapply(CHHislands$startpos, getValue2, data=genes)
#CHHislands$beforegene <- gsub("integer(0)", "0", CHHislands$beforegene)

#CHHislands1 <- grep("[[:digit:]]", CHHislands$beforegene)

CHHislands$beforegene <- as.numeric(CHHislands$beforegene)
CHHislands$aftergene <- as.numeric(CHHislands$aftergene)
#geneswithislands <- unique(CHHislands$beforegene) %>% grep("1","2","3","4","5","6","7","8","9")
#geneswithislands <- grep()

#create 1000 random permutations of methylated positions

CHHislands$actual <- ifelse(is.na(CHHislands$beforegene), ifelse(is.na(CHHislands$aftergene), 0, 1), 1)
percentislands <- sum(CHHislands$actual)/length(freq[[1]])
write(percentislands, file = "percentislands.txt")

geneswithislandsBEFORE <- unique(CHHislands$beforegene) %>% length()
geneswithislandsAFTER <- unique(CHHislands$aftergene) %>% length()
percentGWIbefore <- geneswithislandsBEFORE/length(genes$feature)
percentGWIafter <- geneswithislandsAFTER/length(genes$feature)

write(percentGWIafter, file = "percentGWIafter.txt")
write(percentGWIbefore, file = "percentGWIbefore.txt")
percentGWIbefore
percentGWIafter

permutationstest <- replicate(10, sample(matched$mCHH), simplify = FALSE)
#permutations <- cbind(CHHvec, permutations)
#permutations1 <- t(permutations)
#test <- apply(1, )

#test <- numeric(length = length(permutationstest))
#test <- as.list(test)
#for ( i in seq_along(permutationstest)) {
#  test[i] <- tail(permutationstest[i])
#}
#head(test)


#permutationstest1 <- Map(cbind, list(CHHvec), permutationstest)
permutationstest1 <- sfLapply(permutationstest, cbind, CHHvec)
permutationstest1 <- sfLapply(permutationstest1, as.data.frame)
#permutationstest2 <- Map(as.numeric, )
PMTmCHH <- sfLapply(permutationstest1, subset, V1 == 1)
#PMTmCHHvec <- sfLapply(PMTmCHH, subset, select = V1)
PMTmCHHvec <- sfLapply(PMTmCHH, dplyr::pull, var = CHHvec)
#PMTmCHHvec <- sfLapply(PMTmCHHvec, as.vector, mode = "any")
#PMTmCHHvec <- sfLapply(PMTmCHHvec, as.numeric)
PMTmfreq <- sfLapply(PMTmCHHvec, hist, breaks = br, include.lowest = TRUE, plot = FALSE)
#PMTmfreq1 <- sfLapply(PMTmfreq, getElement, counts)
PMTcountslist <- list()
for ( i in c(1:length(PMTmfreq))){
  PMTcountslist[[i]] <- PMTmfreq[[i]]$counts
}

#PMTmatched <- Map(cbind, list(freq$counts), PMTcountslist)
PMTmatched <- sfLapply(PMTcountslist, cbind, freq$counts)
#PMTmatched <- Map(as.data.frame, PMTmatched)
PMTmatched <- sfLapply(PMTmatched, as.data.frame)
PMTmatched1 <- sfLapply(PMTmatched, subset, V2 > 5)
#PMTpercent <- 

PMTpercent <- list()
for ( i in c(1:length(PMTmatched1))){
  PMTpercent[[i]] <- PMTmatched1[[i]]$V1/PMTmatched1[[i]]$V2
}

PMTmfreq.combined <- unlist(PMTpercent)

PMTmfreq.sorted <- sort(PMTmfreq.combined, decreasing = TRUE) 

write_tsv(as.data.frame(PMTmfreq.sorted), "PMTfreqs.txt")
#Bernie <- head(P
#MTmfreq.sorted, n = 


#https://stackoverflow.com/questions/1563961/how-to-find-top-n-of-records-in-a-column-of-a-dataframe-using-r
#subset(data, V2 > quantile(V2, prob = 1 - n/100))
Bernie <- quantile(PMTmfreq.sorted, prob = .99)
write.table(Bernie, file = "onepercent.txt", sep = ,)
PMTsummary <- summary(permutationstest)
write_tsv(as.data.frame(PMTsummary), "PMTsummary.txt")

#intervals <- sapply(genes, seq, beforewindow, start)


## Gent 2kb region test ##

intervals1 <- list()
for ( i in c(1:(length(genes$start)))){
  intervals1[[i]] <- cbind(genes$beforewindow[[i]], genes$start[[i]])
}

#intervals1 <- apply(intervals, hist, breaks = cbind(min(test), max(test)), plot = FALSE)
intervals2 <- list()
#length(intervals2) <- length(intervals1)
#for (i in c(1:length(intervals))){
#  intervals1[[i]] <- hist(mCHHvec, breaks = c(1, min(intervals[[i]])-1, intervals[[i]], max(intervals[[i]])+1, max(mCHHvec)), plot = FALSE)
#}


#intervals1 <- lapply(intervals1, as.vector)

for (i in c(1:length(intervals1))){
  intervals2[[i]] <- cut(mCHHvec, breaks = c(1, intervals1[[i]], max(mCHHvec)))
}

for (i in c(1:length(intervals2))){
  levels(intervals2[[i]]) <- c("no1", "yes", "no2")
}

#intervals3 <- lapply(intervals2, levels, except = intervals2[2])
intervals3 <- lapply(intervals2, droplevels.factor, exclude=c("no1", "no2"))

intervals3 <- lapply(intervals3, na.omit)

for (i in c(1:length(intervals3))){
  genes[i,8] <- length(intervals3[[i]])
}

##
rm(intervals1)
rm(intervals2)
rm(intervals3)
rm(freq)
rm(permutationstest)
rm(permutationstest1)
rm(PMTmatched)
rm(PMTmatched1)

intervals4 <- list()
for ( i in c(1:length(genes$start))){
  intervals4[[i]] <- cbind(genes$beforewindow[[i]], genes$start[[i]])
}


intervals5 <- list()


for (i in c(1:length(intervals4))){
  intervals5[[i]] <- cut(CHHvec, breaks = c(1, intervals4[[i]], max(CHHvec))) %>% na.omit()
}

for (i in c(1:length(intervals5))){
  levels(intervals5[[i]]) <- c("no1", "yes", "no2")
}


intervals6 <- lapply(intervals5, droplevels.factor, exclude=c("no1", "no2"))

intervals6 <- lapply(intervals6, na.omit)

for (i in c(1:length(intervals6))){
  genes[i,9] <- length(intervals6[[i]])
}


#help <- list()
#for (i in c(1:length(intervals1))){
#  help[[i]] <- c(1, intervals1[[i]], max(mCHHvec))
#}

#test = cbind(1:200)
#testhist <- hist(test, breaks = c(50,200), plot = FALSE)




#permout <- seq(1, 10)
#permout <- as.data.frame(permout)
#test <- seq(2, ncol(permutations))

#for ( i in seq(2, ncol(permutations)) ){
#  test[[i]] <- head(permutations[, i])
#}

#
# Code graveyard (nothing beyond here matters)

#permutationstest1 <- sfLapply(permutationstest, cbind(CHHvec, permutationstest[c(1:length(permutationstest))])

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

