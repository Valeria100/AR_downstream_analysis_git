# READ DATA --------------------------------------------------------------------------------------------------------------------
#Read the summary data set (explanatory) using the package gdata to read .xls files.
library(gdata)
summary <- read.xls("massietimecourse/Summary.xls", sheet=1, check.names=FALSE)
#As the file takes long to read it is saved in the working directory.
save(summary, file="summary.RData")
load("summary.RData")

#File has a double header. I contruct the correct column names.
cols <- c("Probes", rep("VCap",9), rep("LNCap",9), rep("VCap",2), rep("LNCap",2), rep("VCapLNCap",2))
sum_col <- gsub("[.]","",as.character(t(summary[1,])))
sum_col <- gsub(" ", "_", sum_col)
sum_col <- gsub("[?]", "", sum_col)

summary_colnames <- paste(cols, sum_col, sep=".")

#Assign the correct colnames to the summary file.
colnames(summary) <- summary_colnames

#Shift the rownames to later avoid indeces problems.
rownames(summary) <- c(0:(nrow(summary)-1))
summary <- summary[-1,]

#Create a file with only the probe type, id and corresponding hgnc symbol.
probes <- summary[,c(1,2,10,11,19)]

#For all gene I combine both the hgnc name and the illumina probe id. For later use.
# probes_VCap <- paste(probes[,3], probes[,2], sep=".") #"\n"
# probes_LNCap <- paste(probes[,5], probes[,4], sep=".")

#Read sample files and assign it to the corresponding file names
titles <- c("VCapAndrogen1","VCapAndrogen2","VCapControl1","VCapControl2",
            "LNCapAndrogen1","LNCapAndrogen2","LNCapControl1","LNCapControl2")

for(i in 1:length(titles)){
  assign(titles[i],read.table(paste("massietimecourse/",titles[i],".csv",sep=""),header=TRUE,check.names=FALSE,sep=","))
}

#As shown in summary file the last elements are LNCap Only and omitted in the VCap files.
# To even the rows number of the files I add NAs.
# Initially VCap files are 48802 long. After they are 56546 (same as LNCap files).
even_rows <- function(Data, Probes){
  matrix_na <- matrix(NA, length(which(Probes[,1]=="LNCap Only")), ncol(Data))
  colnames(matrix_na) <- colnames(Data)
  even_matrix <- rbind(Data, matrix_na)
  return(even_matrix)
}

VCapAndrogen1 <- even_rows(VCapAndrogen1, probes)
VCapAndrogen2 <- even_rows(VCapAndrogen2, probes)

VCapControl1 <- even_rows(VCapControl1, probes)
VCapControl2 <- even_rows(VCapControl2, probes)

#Convert time in the column titles into numeric and store them in vectors
time_titles <- c("time_VCapAndrogen1","time_VCapAndrogen2","time_VCapControl1","time_VCapControl2",
                 "time_LNCapAndrogen1","time_LNCapAndrogen2","time_LNCapControl1","time_LNCapControl2")

for(i in 1:length(time_titles)){
  assign(time_titles[i],as.numeric(gsub("Time","",colnames(get(titles[i])))))
}
