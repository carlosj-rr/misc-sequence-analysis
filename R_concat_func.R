#!/usr/bin/R

#Function to perform concatenations in R. INPUT ALIGNMENTS MUST BE IN PHYLIP.

#Inputs:
#A. 'someGeneList', a table of alignment filenames that looks like this:

#       fileName
#1 EOG8B00NV.phy
#2 EOG8B00NW.phy
#3 EOG8B00P0.phy

#B. 'sppClTable', a table of species names and the clade they belong to. Must look like this:
# NOTE: Important - species names must be identical through all alignments.
# NOTE: If the cladeName is not known or not important, just fill the column with "Unknown".

#  cladeName sppName
#1   Clade1 Cholhoff
#2   Clade1 Dasynove
#3   Clade2 Echitelf
#4   Clade2 Loxoafri
#5   Clade2 Proccape
#6   Clade3 Ailumela

concatInPhy<-function(someGeneList,sppClTable,outName) {
 outName<-as.character(outName)
 firstAli<-toupper(read.dna(file=as.character(someGeneList$fileName[1]),format="interleave",as.character=TRUE))
 aliTaxa<-rownames(firstAli)

 taxaMissing<-as.character(sppClTable$sppName[!sppClTable$sppName %in% aliTaxa])

 if (length(taxaMissing)) {
   missingDataSeq<-rep("?",ncol(firstAli))
   newTaxaList<-c(aliTaxa,taxaMissing)

   for (i in 1:length(taxaMissing)) {
     firstAli<-rbind(firstAli,missingDataSeq)
   }
   rownames(firstAli)<-newTaxaList
 }

 newOrder<-sapply(sppClTable$sppName,reorganize,firstAli)
 counter<-1
 perc<-round(counter/nrow(someGeneList)*100,digits=0)
 print(paste("Adding ",as.character(someGeneList$fileName[1])," (",perc,"%)",sep=""),quote=FALSE)
 bigAli<-firstAli[newOrder,]

 for (i in someGeneList$fileName[2:nrow(someGeneList)])
 {
   counter<-counter+1
   aliFile<-as.character(i)
   newAli<-toupper(read.dna(file=aliFile,format="interleaved",as.character=TRUE))

   aliTaxa<-rownames(newAli)
   taxaMissing<-as.character(sppClTable$sppName[!sppClTable$sppName %in% aliTaxa])

   if (length(taxaMissing)) {
     missingDataSeq<-rep("?",ncol(newAli))
     newTaxaList<-c(aliTaxa,taxaMissing)

     for (i in 1:length(taxaMissing)) {
       newAli<-rbind(newAli,missingDataSeq)
     }
     rownames(newAli)<-newTaxaList
   }

   newOrder<-sapply(sppClTable$sppName,reorganize,newAli)
   newAli<-newAli[newOrder,]
   perc<-round(counter/nrow(someGeneList)*100,digits=0)
   print(paste("Adding ",aliFile," (",perc,"%)",sep=""),quote=FALSE)
   bigAli<-cbind(bigAli,newAli)

 }

  print(c("Writing concatenated file to ",outName),quote=FALSE)
  write.dna(bigAli,file=outName,format="interleaved")
}

reorganize<-function(sppClSpp,dataset) {
  return(grep(sppClSpp,rownames(dataset)))
}


# WATCH OUT: FUNCTIONS AFTER THIS LINE ARE NOT WORKING PROPERLY...
concatInFas<-function(someGeneList,sppClTable,outName) { #*****STILL NOT PROPERLY WORKING******* The functions for importing phy and fas are not the same when the "as.character=TRUE" argument is given, on phy the data gets imported as a character matrix (works well), and on fasta it gets imported as a list of character vectors, I must find a way to convert the list of character vectors into a char matrix.
 outName<-as.character(outName)
 firstAli<-read.dna(file=as.character(someGeneList$fileName[1]),format="fasta",as.character=TRUE)
 firstAli<-toupper(firstAli) #THIS command doesn't work on the character list
 aliTaxa<-rownames(firstAli)

 taxaMissing<-as.character(sppClTable$sppName[!sppClTable$sppName %in% aliTaxa])

 if (length(taxaMissing)) {
   missingDataSeq<-rep("?",ncol(firstAli))
   newTaxaList<-c(aliTaxa,taxaMissing)

   for (i in 1:length(taxaMissing)) {
     firstAli<-rbind(firstAli,missingDataSeq)
   }
   rownames(firstAli)<-newTaxaList
 }

 newOrder<-sapply(sppClTable$sppName,reorganize,firstAli)
 counter<-1
 perc<-round(counter/nrow(someGeneList)*100,digits=0)
 print(paste("Adding ",as.character(someGeneList$fileName[1])," (",perc,"%)",sep=""),quote=FALSE)
 bigAli<-firstAli[newOrder,]

 for (i in someGeneList$fileName[2:nrow(someGeneList)])
 {
   counter<-counter+1
   aliFile<-as.character(i)
   newAli<-read.dna(file=aliFile,format="fasta",as.character=TRUE)
   newAli<-toupper(newAli)

   aliTaxa<-rownames(newAli)
   taxaMissing<-as.character(sppClTable$sppName[!sppClTable$sppName %in% aliTaxa])

   if (length(taxaMissing)) {
     missingDataSeq<-rep("?",ncol(newAli))
     newTaxaList<-c(aliTaxa,taxaMissing)

     for (i in 1:length(taxaMissing)) {
       newAli<-rbind(newAli,missingDataSeq)
     }
     rownames(newAli)<-newTaxaList
   }

   newOrder<-sapply(sppClTable$sppName,reorganize,newAli)
   newAli<-newAli[newOrder,]
   perc<-round(counter/nrow(someGeneList)*100,digits=0)
   print(paste("Adding ",aliFile," (",perc,"%)",sep=""),quote=FALSE)
   bigAli<-cbind(bigAli,newAli)

 }

  print(c("Writing concatenated file to ",outName),quote=FALSE)
  write.dna(bigAli,file=outName,format="interleaved")
}

reorganize<-function(sppClSpp,dataset) {
  return(grep(sppClSpp,rownames(dataset)))
}


#UNTESTED-----SUGGESTED USE FOR MULTIPLE CONCATENATIONS:
#If many concatenations will be made, make a table of the filenames of the files that will later be the "someGeneList", like this:

#> head(masterList)
#         concListFileName
#1 ConcEdRev3-100_all.list
#2  ConcEdRev3-10_all.list
#3  ConcEdRev3-11_all.list
#4  ConcEdRev3-12_all.list
#5  ConcEdRev3-13_all.list
#6  ConcEdRev3-14_all.list

#And finally use the other function below: massiveConc; with lapply, like this:

#>lapply(masterList$concListFileName,massiveConc)


massiveConc<-function(concListFileName) {

 someGeneList<-read.csv(file=as.character(concListFileName),header=FALSE,col.names="Ali")
 outName<-paste(strsplit(as.character(concListFileName),".",fixed=TRUE)[[1]][1],".phy",sep="")
 concat(someGeneList,outName)

}
